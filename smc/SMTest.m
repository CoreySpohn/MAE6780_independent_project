function SMTest

close all

%	Alphabetical List of Initial Conditions
	alpha =	    0;	% Angle of attack, deg	(relative to air mass)
	beta =		0;		% Sideslip angle, deg	(relative to air mass)
	cm =		0.25;	% Longitudinal center-of-mass location, % mac/100
	CONTROL = 	0;		% Feedback control ON (= 1) or OFF (= 0)
	dA =		0;		% Aileron angle, deg
	dAS =		0;		% Asymmetric spoiler angle, deg
	dE =		0;	% Elevator angle, deg
	dR =		0;		% Rudder angle, deg
	dF = 		0;		% Flap setting, deg
	dS = 		-1; %-1;%-1.948;	% Stabilator setting, deg
	dT = 		0.5;	% Throttle setting, % / 100
	GEAR = 		0;		% Landing gear DOWN (= 1) or UP (= 0)
	h =			5000;	% Altitude above Sea Level, m
	hdot =		0;		% Altitude rate, m/s
	LINEAR = 	0;		% Linear model flag (= 1 to calculate F and G)
	p =			0;		% Body-axis roll rate, deg/s
	phi =		0;		% Body roll angle wrt earth, deg
	psi =		0;		% Body yaw angle wrt earth, deg
	q	=		0;		% Body-axis pitch rate, deg/sec
	r =			0;		% Body-axis yaw rate, deg/s
	SIMUL =		1;		% Flight path flag (= 1 for nonlinear simulation)
	SPOIL =		0;		% Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
	tf =		30;		% Final time, sec
	ti = 		0;		% Initial time, sec
	theta =		alpha;	% Body pitch angle wrt earth, deg
	TRIM = 		1;		% Trim flag (= 1 to calculate trim)
	V =			200;	% True Air Speed, TAS, m/s	(relative to air mass)
	xe =		-500000;		% Initial longitudinal position, m
	ye = 		0;		% Initial lateral position, m
	ze = 		-h;		% Initial vertical position, m
		
%	Initial Conditions depending on prior initial conditions

	[airDens,airPres,temp,soundSpeed] = Atmos(h);

	gamma	=	57.29578 * atan(hdot / sqrt(V^2 - hdot^2));
						% Inertial vertical flight path angle, deg
	qbar	= 	0.5 * airDens * V^2;	
						% Dynamic pressure, N/m^2
	IAS		=	sqrt(2 * qbar / 1.225);
						% Indicated Air Speed, m/s
	Mach	= 	V / soundSpeed;	
						% Mach number
                        
%	State Vector and Control Initialization
	phir	=	phi * .01745329;
	thetar	=	theta * .01745329;
	psir	=	psi * .01745329;

	windb	=	WindField(-ze,phir,thetar,psir);
	alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;
    
    global xn0 un

	xn0	=	[V * cos(alphar) * cos(betar) - windb(1)
			V * sin(betar) - windb(2)
			V * sin(alphar) * cos(betar) - windb(3)
			xe
			ye
			ze
			p * .01745329
			q * .01745329
			r * .01745329
			phir
			thetar
			psir];
	
	un	=	[dE * .01745329
			dA * .01745329
			dR * .01745329
			dT
			dAS * .01745329
			dF * .01745329
			dS * .01745329];
						
%	Initial-Condition Perturbation (Test Inputs)
	delx	=	[0;0;0
				0;0;0
				0.00;0.00;0.00
				0.00;0.00;0];
				
	delu	=	[0;0;0;0;0;0;0];
        
%   Initial & final time
    t0 = 0;
    tf = 6000;

%   Linearized dynamics (no hurricane)
    global A B
    global hurricane
    hurricane = 0;
    thresh = 0.1 * ones(19, 1);
    xj = [xn0; un];
	xdotj = LinModel(t0,xj);
	[dFdX,~] = numjac('LinModel', t0, xj, xdotj, thresh, [], 0);
    A = dFdX(1:12,1:12);
	B = dFdX(1:12,[13:16 18:19]);
    
%   True initial state
    x0 = xn0 + delx;
    
%   Simulation steps
    N = 1000;
    T = linspace(t0, tf, N);
    
%   Simulate hurricane
    global hurr_para hurr_Z hurr_T
    hurr_para.linvel = 0; %m/s
    hurr_para.angvel = 0; % rad/s % This will be a very small number, but it will cause the hurricane path to arc
    hurr_para.Vmax = 252; %km/hr The units on this are important!
    hurr_para.Rmax = 47; % km The units on this are important!
    hurr_para.xmax = 20000;
    hurr_para.ymax = 20000;
    hurr_para.scalefactor = 100;
%   Generate Noise
    noise=make_noise(100,100);
    hurr_para.noise = noise;
    z_hurr = [ 0, 0, deg2rad(45)]; % Initial hurricane conditions
%   North position of center of mass WRT Earth, xe, m
%	East position of center of mass WRT Earth, ye, m
%   Angle of Hurricane Path with respect to the Earth Frame's East Axis
    global counter
    counter = 0;
    [hurr_T, hurr_Z] = ode45(@HurricaneEOM, T, z_hurr);

%   Integrate true motion (yes hurricane)
    hurricane = 1;
    [T, Xt] = ode23t(@SMEoM, T, x0);
    
%   Nominal states
    Xn = zeros(N, 12);
    for k = 1:N
        Xn(k,:) = NomState(T(k))';
    end
    
%   Fullscreen figure for plots
    %figure('units','normalized','outerposition',[0 0 1 1])
    
%   Plot planar motion
    subplot(2, 2, 1)
    plot(Xn(:,4), Xn(:,5), 'b-', Xt(:,4), Xt(:,5), 'r:', ...
        Xn(1,4), Xn(1,5), 'b*', Xt(1,4), Xt(1,5), 'r*')
    legend('Desired','True','Location','best') 
    xlabel('x (m)')
    ylabel('y (m)')
    axis equal
    
%   Plot altitude
    subplot(2, 2, 2)
    plot(T, -Xn(:,6), 'b-', T, -Xt(:,6), 'r:')
    legend('Desired','True','Location','best')
    xlabel('Time (s)')
    ylabel('Altitude (m)')
    
%   Plot velocity
    subplot(2, 2, 3)
    plot(T, Xn(:,1), 'b-', T, Xt(:,1), 'r:')
    legend('Desired','True','Location','best')
    xlabel('Time (s)')
    ylabel('Axial Velocity (m/s)')
    
%   Plot roll/pitch/yaw
    subplot(2, 2, 4)
    plot(T, rad2deg(Xt(:,10)), T, rad2deg(Xt(:,11)), T, rad2deg(Xt(:,12)))
    legend('Roll','Pitch','Yaw','Location','best')
    xlabel('Time (s)')
    ylabel('Angle (deg)')

%   Save plots
    PrepFigPresentation(gcf)
    saveas(gcf, 'sm_flight.png')

%--------------------------------------------------------------------------

    Vmax = hurr_para.Vmax;
    Rmax = hurr_para.Rmax*1000; % Convert from km to m
    
    [X,Y] = meshgrid(-500:100:500);
    for i = 1:length(X(1,:))
        for j = 1:length(Y(:,1))
            r = (X(1,i)^2+Y(j,1)^2)^.5;
            if r<Rmax
               Vr = Vmax*(r/Rmax)^(3/2);%Scalar multiplier of the vector field
            else
                Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
            end
            dx(j,i) = Vr*(Y(j,1));        %/(X.^2+Y.^2));
            dy(j,i) = Vr*(-X(1,i));       %/(X.^2+Y.^2));
        end
    end
    
    figure
    % contour(X,Y,dx,dy)
    hold on
    quiver(X,Y,dx,dy,'AutoScaleFactor',1)
    plot(Xn(:,4)/1000, Xn(:,5)/1000, 'b-', ...
         Xt(:,4)/1000, Xt(:,5)/1000, 'r:', 'LineWidth', 1.5)
    legend('Wind','Desired Trajectory', 'Actual Trajectory')
    hold off
    axis equal
    xlabel('x [km]')
    ylabel('y [km]')

    %PrepFigPresentation(gcf)
    saveas(gcf, 'sm_hurricane.png')
    
    % Plot RPY
    figure
    plot(T, rad2deg(Xt(:,10)), T, rad2deg(Xt(:,11)), T, rad2deg(Xt(:,12)))
    legend('Roll','Pitch','Yaw','Location','best')
    xlabel('Time (s)')
    ylabel('Angle (deg)')
    PrepFigPresentation(gcf)
    saveas(gcf, 'sm_rpy.png')
    
end

function xd = SMEoM(t, x)

    % Constants
    Kr = diag([100, 100, 100]);
    Kth = diag([1, 1, 0.001]);
    L = [0.001 * ones(3,1); 1 * ones(3,1)];
    
    % Nominal state
    xn = NomState(t);
    
    % Sliding-mode control
    global u un
    u = SMC(t, x, xn, un, Kr, Kth, L);
    
    % State derivative
    xd = EoM(t, x);

end

function xn = NomState(t)

    global xn0
    xn = xn0;
    H = DCM(xn(10), xn(11), xn(12))';
    xn(4:6) = xn(4:6) + t * H * xn(1:3);

%     global hurr_Z hurr_T
% 
%     z_hurr = interp1(hurr_T, hurr_Z, t, 'spline');
%     
%     xn = [zeros(3,1);...
%                   z_hurr(2);... % Hurricane Center Position, North
%                   z_hurr(1);... % Hurricane Center Position, East
%                   -3000;... % Desired Flight Height
%                   zeros(6,1)];

%     R = 10000;
%     V = 240;
%     z = -9000;
%     th0 = pi;
%     w = -V/R;
%     
%     th = w * t + th0;
%     
%     xn = [V; 0; 0; R*cos(th); R*sin(th); z; 0; 0; w; 0; 0; ...
%           th+sign(w)*(pi/2)];
    
    
end
