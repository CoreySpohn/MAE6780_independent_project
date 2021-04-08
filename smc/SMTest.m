function SMTest

    close all

    %	Alphabetical List of Initial Conditions
	alpha =		3.63;	% Angle of attack, deg	(relative to air mass)
	beta =		0;		% Sideslip angle, deg	(relative to air mass)
	cm =		0.25;	% Longitudinal center-of-mass location, % mac/100
	CONTROL = 	0;		% Feedback control ON (= 1) or OFF (= 0)
	dA =		0;		% Aileron angle, deg
	dAS =		0;		% Asymmetric spoiler angle, deg
	dE =		0;	% Elevator angle, deg
	dR =		0;		% Rudder angle, deg
	dF = 		0;		% Flap setting, deg
	dS = 		-1.948;	% Stabilator setting, deg
	dT = 		0.1919;	% Throttle setting, % / 100
	GEAR = 		0;		% Landing gear DOWN (= 1) or UP (= 0)
	h =			9150;	% Altitude above Sea Level, m
	hdot =		0;		% Altitude rate, m/s
	LINEAR = 	0;		% Linear model flag (= 1 to calculate F and G)
	p =			0;		% Body-axis roll rate, deg/s
	phi =		0;		% Body roll angle wrt earth, deg
	psi =		45;		% Body yaw angle wrt earth, deg
	q	=		0;		% Body-axis pitch rate, deg/sec
	r =			0;		% Body-axis yaw rate, deg/s
	SIMUL =		1;		% Flight path flag (= 1 for nonlinear simulation)
	SPOIL =		0;		% Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
	tf =		30;		% Final time, sec
	ti = 		0;		% Initial time, sec
	theta =		alpha;	% Body pitch angle wrt earth, deg
	TRIM = 		0;		% Trim flag (= 1 to calculate trim)
	V =			244;	% True Air Speed, TAS, m/s	(relative to air mass)
	xe =		0;		% Initial longitudinal position, m
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
						
%	Initial-Condition Perturbation (Test Inputs)
	delx	=	[10;10;1
				10000;20000;100
				0.01;0.01;0.01
				0.01;0.01;-0.1];
				
	delu	=	[0;0;0;0;0;0;0];

%	State Vector and Control Initialization
	phir	=	phi * .01745329;
	thetar	=	theta * .01745329;
	psir	=	psi * .01745329;

	windb	=	WindField(-ze,phir,thetar,psir);
	alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;

%   Nominal initial state
    global xn0
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
        
%   Nominal control
    global un
	un	=	[dE * .01745329
			dA * .01745329
			dR * .01745329
			dT
			dAS * .01745329
			dF * .01745329
			dS * .01745329];
        
%   Initial & final time
    t0 = 0;
    tf = 600;
        
%   Linearized dynamics
    global B
    thresh = 0.1 * ones(19, 1);
    xj = [xn0; un];
	xdotj = LinModel(t0,xj);
	[dFdX,~] = numjac('LinModel', t0, xj, xdotj, thresh, [], 0);
	B = dFdX(1:12,13:19);
    
%   True initial state
    x0 = xn0 + delx;
    
%   Simulation steps
    N = 1000;
    T = linspace(t0, tf, N);

%   Integrate true motion
    [T, X] = ode23s(@SMEoM, T, x0);
    
%   Nominal states
    Xn = zeros(N, 12);
    for k = 1:N
        Xn(k,:) = NomState(T(k))';
    end
    
%   Plot planar motion
    figure
    plot(Xn(:,4), Xn(:,5), 'b-', X(:,4), X(:,5), 'r:', ...
        Xn(1,4), Xn(1,5), 'b*', X(1,4), X(1,5), 'r*')
    legend('Desired','True') 
    xlabel('x (m)')
    ylabel('y (m)')
    axis equal
    PrepFigPresentation(gcf)
    saveas(gcf, 'position.png')
    
%   Plot attitude
    figure
    plot(T, -Xn(:,6), 'b-', T, -X(:,6), 'r:')
    legend('Desired','True')
    xlabel('Time (s)')
    ylabel('Altitude (m)')
    PrepFigPresentation(gcf)
    saveas(gcf, 'altitude.png')
    

end

function xd = SMEoM(t, x)

    % Constants
    tau = 100;
    L = 10;
    ls = 10;
    
    % Nominal state
    xn = NomState(t);
    
    % Sliding-mode control
    global u un
    u = SMC(t, x, xn, un, tau, L, ls);
    
    % State derivative
    xd = EoM(t, x);

end

function xn = NomState(t)

    global xn0
    xn = xn0;
    H = DCM(xn(10), xn(11), xn(12))';
    xn(4:6) = xn(4:6) + t * H * xn(1:3);

end