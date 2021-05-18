%	FLIGHT --  6-DOF Trim, Linear Model, and Flight Path Simulation 
%	August 4, 1999   
%	===============================================================
%	Copyright 1993-1999 by ROBERT F. STENGEL.  All rights reserved.

	close all; clearvars; clc
    set(0,'DefaultAxesFontSize',17);
	global GEAR CONTROL SPOIL u x V parhis %windb

%	This is the EXECUTIVE FILE.  It contains the Main Program, which:
%		Defines initial conditions
%		Calculates longitudinal trim condition
%		Calculates stability-and-control derivatives
%		Simulates flight path using nonlinear equations of motion
		
%	Functions used by FLIGHT:
%		AeroModel.m		Aerodynamic coefficients of the aircraft, thrust model,
%						and geometric and inertial properties
%		Atmos.m			Air density, sound speed
%		ControlSystem.m	Control law
%		DCM.m			Direction-cosine matrix
%		EoM.m			Equations of motion for integration
%		LinModel.m		Equations of motion for linear model definition
%		TrimCost.m		Cost function for trim solution
%		WindField.m		Wind velocity components

%	DEFINITION OF THE STATE VECTOR
%		x(1) = 		Body-axis x inertial velocity, ub, m/s
%		x(2) =		Body-axis y inertial velocity, vb, m/s
%		x(3) =		Body-axis z inertial velocity, wb, m/s
%		x(4) =		North position of center of mass WRT Earth, xe, m
%		x(5) =		East position of center of mass WRT Earth, ye, m
%		x(6) =		Negative of c.m. altitude WRT Earth, ze = -h, m
%		x(7) =		Body-axis roll rate, pr, rad/s
%		x(8) =		Body-axis pitch rate, qr, rad/s
%		x(9) =		Body-axis yaw rate, rr,rad/s
%		x(10) =		Roll angle of body WRT Earth, phir, rad
%		x(11) =		Pitch angle of body WRT Earth, thetar, rad
%		x(12) =		Yaw angle of body WRT Earth, psir, rad
	
%	DEFINITION OF THE CONTROL VECTOR
%		u(1) = 		Elevator, dEr, rad
%		u(2) = 		Aileron, dAr, rad
%		u(3) = 		Rudder, dRr, rad
%		u(4) = 		Throttle, dT, %
%		u(5) =		Asymmetric Spoiler, dASr, rad
%		u(6) =		Flap, dFr, rad
%		u(7) =		Stabilator, dSr, rad

%	BEGINNING of MAIN PROGRAM
%	=========================

	'FLIGHT'

%	Alphabetical List of Initial Conditions
	alpha =		0;	% Angle of attack, deg	(relative to air mass)
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
	h =			5000;	% Altitude above Sea Level, m
	hdot =		0;		% Altitude rate, m/s
	LINEAR = 	1;		% Linear model flag (= 1 to calculate F and G)
	p =			0;		% Body-axis roll rate, deg/s
	phi =		0;		% Body roll angle wrt earth, deg
	psi =		45;		% Body yaw angle wrt earth, deg
	q	=		0;		% Body-axis pitch rate, deg/sec
	r =			0;		% Body-axis yaw rate, deg/s
	SIMUL =		0;		% Flight path flag (= 1 for nonlinear simulation)
	SPOIL =		0;		% Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
	tfinal =		30;		% Final time, sec
	ti = 		0;		% Initial time, sec
	theta =		alpha;	% Body pitch angle wrt earth, deg
	TRIM = 		0;		% Trim flag (= 1 to calculate trim)
	V =			240;	% True Air Speed, TAS, m/s	(relative to air mass)
	xe =		0;		% Initial longitudinal position, m
	ye = 		-50000;		% Initial lateral position, m
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
	delx	=	[0;0;0
				0;0;0
				0;0;0
				0;0;0];
				
	delu	=	[0;0;0;0;0;0;0];

%	State Vector and Control Initialization
	phir	=	phi * .01745329;
	thetar	=	theta * .01745329;
	psir	=	psi * .01745329;
    alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;
    
%   Initialize State Before the First Loop
    x	=	[V * cosd(alpha) * cosd(beta)
        V * sind(beta)
        V * sind(alpha) * cosd(beta)
        xe
        ye
        ze
        p * .01745329
        q * .01745329
        r * .01745329
        phir
        thetar
        psir]
    
    u	=	[dE * .01745329
        dA * .01745329
        dR * .01745329
        dT
        dAS * .01745329
        dF * .01745329
        dS * .01745329]
    
    %	Trim Calculation (for Steady Level Flight at Initial V and h)
	if TRIM >= 1
		'TRIM'
		parhis	=	[];
		optPar				=	[-.05;.19;.06];
		options				=	foptions;
		[optPar,options]	=	fmins('TrimCost',optPar,options);
%		Optimizing dSr, dT, Theta, and Error Cost
		optPar
		J					=	options(8);
		iterations			=	options(10)
		index	=	[1:iterations];
		figure
		subplot(1,2,1)
		plot(index,parhis(1,:),index,parhis(2,:),index,parhis(3,:))
		xlabel('Iterations'), ylabel('Trim Parameters'), grid

		subplot(1,2,2)
		semilogy(index,parhis(4,:))
		xlabel('Iterations'), ylabel('Trim Cost'), grid
		u	=	[u(1)
				u(2)
				u(3)
				optPar(2)
				u(5)
				u(6)
				optPar(1)]
		format long			
		x	=	[V * cos(optPar(3))
				x(2)
				V * sin(optPar(3))
				x(4)
				x(5)
				x(6)
				x(7)
				x(8)
				x(9)
				x(10)
				optPar(3)
				x(12)]
		format short
    end

%   Start of the Iterative Loop
    Tstop = 500; % Total Simulation Time
    Ts = 1; % Sample Time
    p = 50; % Predictive Horizon
    c = 50; % Control Horizon
    tarray = (0:Ts:p*Ts); % Need to simulate hurricane through the predictive horizon each loop
    state = x'; % Set up variable to contain the state variable after each loop, put nominal in to start
    control = [0, 0, 0, 250, 0, 0, 0]; % Set up variable to contain the control variables after each loop
    z_total = []; % Set up variable to contain the hurricane track after each loop
    ref_total = [5000,x(5),x(6)]; % Start trajectory at aircraft starting position
    % We assume we have perfect knowledge of the states 
    C = eye(12);
    D = zeros(12,7);

%	Stability-and-Control Derivative Calculation
   	if LINEAR >= 1
		'LINEAR'
		thresh	=	[.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1];
		xj		=	[x;u];
		xdotj		=	LinModel(ti,xj);
		[dFdX,fac]	=	numjac('LinModel',ti,xj,xdotj,thresh,[],0);
		Fmodel		=	real(dFdX(1:12,1:12));
		Gmodel		=	real(dFdX(1:12,13:19));
		save Fmodel
		save Gmodel
        % Add three extra columns to Gmodel for the unmeasured disturbances
        disturbance_model = [eye(3);zeros(9,3)];
        Gmodel = [Gmodel, disturbance_model];
        D = zeros(12,10);
        CSTR = ss(Fmodel, Gmodel, C, D);
        CSTR = setmpcsignals(CSTR,'MV',[1 2 3 4 5 6 7],'UD',[8 9 10]);
        mpc_obj = mpc(CSTR, Ts,p,c);
        % Previous State of the aircraft
        %mpc_obj.Model.Nominal.Y = [1000, 1000, 20000, 0, 0, 0]; %Not Fully Observable State Space Model
        mpc_obj.Model.Nominal.Y = state(end,:); % Use last state, no disturbances included
        % Starting control input
        mpc_obj.Model.Nominal.U = [control(end,:) 0 0 0]; % Use last control
        % Set Weights for Q and R
        mpc_obj.Weights.MV = [.7 .7 .7 1 .7 .7 .7];
        mpc_obj.Weights.MVRate = [10 10 10 10 10 10 10];
%         mpc_obj.Weights.OV = [.7 0.001 0.001 .95 .95 .95 0.001 0.001 0.001 .1 .1 .1];
        mpc_obj.Weights.OutputVariables = [.1 10 0.1 10 10 1 0.0001 0.0001 0.0001 .1 .1 .1];
    %     mpc_obj.Weights.ECR = 100000;
    else
        mpc_obj = nlmpc(12,12,7);%,'MV',[1 2 3 4 5 6 7],'UD',[1 2 3]);
        mpc_obj.Ts = Ts;
        mpc_obj.PredictionHorizon = p;
        mpc_obj.ControlHorizon = c;
        mpc_obj.Model.StateFcn = @(x,u) EOM(t,x,u);
        % Set Weights for Q and R
        mpc_obj.Weights.ManipulatedVariables = [.9 .9 .9 1 .9 .9 .9];
    %     mpc_obj.Weights.MVRate = [0.1 0.1 0.1 0.1 0.1 0.1 0.1];
        mpc_obj.Weights.OutputVariables = [0.99 0.1 0.1 .5 1 1 0.0001 0.0001 0.0001 .1 .1 .1];
        mpc_obj.Weights.ECR = 1000;
    end
  
%     % Throttle max and min
%     % Control Surface Angles limited to +or- 30 deg = 0.5236
    for k = 1:7
        if k == 4
            mpc_obj.MV(k).Max = Inf;
            mpc_obj.MV(k).Min = 0;
            mpc_obj.MV(k).ScaleFactor = 100; 
%             mpc_obj.MV(k).MaxECR = 0;
%             mpc_obj.MV(k).MinECR = 0;
        else
            mpc_obj.MV(k).Max = deg2rad(60);%0.5236;
            mpc_obj.MV(k).Min = deg2rad(-60);%-0.5236;
%             mpc_obj.MV(k).MaxECR = 0;
%             mpc_obj.MV(k).MinECR = 0;
            % Set Span to Better Scaling of Hessian
            mpc_obj.MV(k).ScaleFactor = 1;

        end
    end
% 
%     % Output Minimums i.e. We cannot fly into the ground, We cannot fly
%     % faster than 400 mph, pitch, roll, and yaw are limited as well
    mpc_obj.OV(1).Max = 500; % maximum axial velocity, 400 mph = 175 m/s
    mpc_obj.OV(1).MaxECR = 1;
    hurr_para.maxVelAircraft = 175;
    mpc_obj.OV(1).Min = 0;
%     mpc_obj.OV(1).MaxECR = 0.1;
%     mpc_obj.OV(1).MinECR = 0;
%     mpc_obj.OV(1).ScaleFactor = 100;
% %     mpc_obj.OV(2).MaxECR = 10; % maximum velocity, 400 mph = 175 m/s
% %     mpc_obj.OV(2).MinECR = -10;
%     mpc_obj.OV(2).ScaleFactor = 10;
%     mpc_obj.OV(3).ScaleFactor = 10;
    mpc_obj.OV(4).ScaleFactor = 10;
    mpc_obj.OV(5).ScaleFactor = 100;
%     mpc_obj.OV(6).Max = 0; % meter, maximum operational range
%     mpc_obj.OV(6).Min = -16000;
    hurr_para.goalaltitude = -5000;
    mpc_obj.OV(6).ScaleFactor = 10;
%     mpc_obj.OV(6).MinECR = 0;
% %     mpc_obj.OV(7).ScaleFactor = 1;
% %     mpc_obj.OV(8).ScaleFactor = .1;
% %     mpc_obj.OV(9).ScaleFactor = 1;
%     mpc_obj.OV(10).ScaleFactor = 10;
% %     mpc_obj.OV(11).ScaleFactor = 1;
%     mpc_obj.OV(12).ScaleFactor = 10;
    mpc_obj.OV(10).Max = 3*pi/2; % Roll rad, large plane will not barrel roll
    mpc_obj.OV(10).Min = -3*pi/2;
    mpc_obj.OV(10).MaxECR = 0;
    mpc_obj.OV(10).MinECR = 0;
    mpc_obj.OV(11).Max = 3*pi/8; % Pitch rad, the plane cannot flip over the minor axis
    mpc_obj.OV(11).Min = -pi/2;
    mpc_obj.OV(11).MaxECR = 0;
    mpc_obj.OV(11).MinECR = 0;
    mpc_obj.OV(12).Max = 0.5236; % Yaw rad, the plane will not spin
    mpc_obj.OV(12).Min = -0.5236;
    mpc_obj.OV(12).MaxECR = 0;
    mpc_obj.OV(12).MinECR = 0;

    
%     mpc_refsignal = zeros(11, 7);
%     if i == 1
%         xmpc = mpcstate(mpc_obj);
%         xmpc.Plant = state(end,:);        
%     end    

%   Hurricane Initial Conditions and Parameters
    hurr_para.linvel = 4.9; %m/s
    hurr_para.angvel = 0; % rad/s % This will be a very small number, but it will cause the hurricane path to arc
    hurr_para.Vmax = 252; %km/hr The units on this are important!
    hurr_para.Rmax = 5;%47; % km The units on this are important!
    hurr_para.xmax = 20000;
    hurr_para.ymax = 20000;
    hurr_para.scalefactor = 10;
    science_traj = false;
    
%   Create a Noise WindField, only sliding mode will do this outside the
%   loop

    % Set up ODE 45 outside the iterative loop,
    hurr_prop = @(t,z) HurricaneEOM2(t,z,hurr_para);
    hurr_science = @(t,z) CircularReferenceTrajectoryEOM(t,z,hurr_para);
%     n = 100;
%     tf = 10; % This should equal the control horizon
%     tarray = linspace(0,tf,n);
    smallnumber=1e-10;
    opts=odeset('Abstol',smallnumber,'RelTol',smallnumber);

    z_hurr = [ 15000, 5000, deg2rad(45)]; % Initial hurricane conditions
%   North position of center of mass WRT Earth, xe, m
%	East position of center of mass WRT Earth, ye, m
%   Angle of Hurricane Path with respect to the Earth Frame's East Axis
     
%     mpc_options = mpcmoveopt();
%     mpc_options.MVMax = [0.5236 0.5236 0.5236 Inf 0.5236 0.5236 0.5236];
%     mpc_options.MVMin = [-0.5236 -0.5236 -0.5236 0 -0.5236 -0.5236 -0.5236];
%     mpc_options.OutputMax = [175 10 Inf Inf Inf 0 Inf Inf Inf 3*pi/2 pi/2 Inf];
%     mpc_options.OutputMin = [0 -10 -Inf -Inf -Inf -9000 -Inf -Inf -Inf -3*pi/2 -pi/2 -Inf];
    
%    Switch betweemn sim and mpcmove
    options = mpcsimopt();
    options.RefLookAhead = 'on';
    options.MDLookAhead = 'off';
    options.Constraints = 'on';
    options.OpenLoop = 'off';
    

%for i = 1:round(Tstop/Ts)+1
for i = 1:round(Tstop/c) % Number of Iterations if using Sim
    % Calculate current output before applying wind disturbances
%     y_mpcmove = C*x+D*u;
    
    % Calculate the Current Wind Disturbances
	windb	=	HurricaneWindField(x,z_hurr,hurr_para);

    mod1 = windb(1)*tf([1],[1 0]);
    mod2 = windb(2)*tf([1],[1 0]);
    mod3 = windb(3)*tf([1],[1 0]);
    mod_comb = eye(3)*[mod1;mod2;mod3];%[eye(3)*[mod1;mod2;mod3],zeros(3,9)];
    setindist(mpc_obj,'model',mod_comb);

%     % Add the Unmeasured Disturbance Going into Plant
% 	x	=	x - [windb(1);...
%                  windb(2);...
% 			     windb(3);...
%                  zeros(9,1)];
             
    mpc_obj.Model.Nominal.Y = state(end,:); % Use last state, no disturbances included
    mpc_obj.Model.Nominal.U = [control(end,:) 0 0 0]; % Use last control

%   Huricane Modeling
    % Run an ODE 45 Solver to calculate the hurricane trjectory over the
    % current control horizon  

    zarray = [ref_total(end,1)*ones(1,Ts*p);...
                linspace(ref_total(end,2),ref_total(end,2)+200*p*Ts,p*Ts);...
                ref_total(end,3)*ones(1,Ts*p)];
    index = c*Ts;
    z_hurr = z_hurr;
    z_total = [z_total;z_hurr]; % Store up through the control horizon
    ref_total = [ref_total;zarray(1:3,index)'];
    ref_signal = [240*ones(1,length(zarray(1,2:end)));...
                  zeros(2,length(zarray(1,2:end)));...
                  zarray(1,2:end);... % Hurricane Center Position, North
                  zarray(2,2:end);... % Hurricane Center Position, East
                  zarray(3,2:end);... % Desired Flight Height
                  zeros(6,length(zarray(1,2:end)))]';


    % Attempt with Simulation Function
    % sim(mpc_obj, 1, ref_signal, [], options)
    [y, t, u_t, xp, xc, output_options] = sim(mpc_obj, p+1, ref_signal, [], options);
    index = find(t==(c*Ts));
    u = u_t(index,:)';
    x = y(index,:)';
%     x = Fmodel*x+Gmodel*u;
    state = [state;y(2:index,:)];
    control = [control;u_t(2:index,:)];
    
    % Attempt with MPC Move with Explicit Defintion of Time
%     [u,Info] = mpcmove(mpc_obj, xmpc, y_mpcmove, ref_signal, [], mpc_options);
%     index = find(Info.Topt==c);
%     x = xmpc.Plant; % Info.Xopt(index,1:12)';
%     state = [state;Info.Xopt(1:index,1:12)];
%     control = [control;u'];
end

%     LQR_sys = ss(Fmodel-Gmodel*K1, Gmodel, C, D);
%     LQR = sim(LQR_sys,zeros(length(tarray),7),tarray,[0, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0]')

%% Quiver plot vals
east_vals = -500000:50000:500000;
north_vals = -500000:50000:500000;
[X, Y, U, V] = HurricaneQuiver(east_vals, north_vals, z_hurr, hurr_para);

%% Adapting the stuff from SMC stuff
t_control = linspace(0,Tstop,length(control(:,1)))';
ref_t = linspace(0, Tstop, length(ref_total));
ref_v = 300*ones(length(ref_total),1);
%   Plot planar motion
    subplot(2, 2, 1)
    plot(ref_total(:,2), ref_total(:,1), 'b-', state(:,5), state(:,4), 'r:', ...
        ref_total(1,2), ref_total(1,1), 'b*', state(1,5), state(1,4), 'r*', 'LineWidth',3)
    legend('Desired','True','Location','best') 
    xlabel('x (m)')
    ylabel('y (m)')
    axis equal

%   Plot altitude
    subplot(2, 2, 2)
    plot(ref_t, -ref_total(:,3), 'b-', t_control, -state(:,6), 'r:', 'LineWidth',3)
    legend('Desired','True','Location','best')
    xlabel('Time (s)')
    ylabel('Altitude (m)')

%   Plot velocity
    subplot(2, 2, 3)
    plot(ref_t, ref_v(:,1), 'b-', t_control, state(:,1), 'r:', 'LineWidth',3)
    legend('Desired','True','Location','best')
    xlabel('Time (s)')
    ylabel('Axial Velocity (m/s)')
    
%   Plot roll/pitch/yaw
    subplot(2, 2, 4)
    plot(t_control, state(:,10), t_control, state(:,11), t_control, state(:,12), 'LineWidth',3)
    legend('Roll','Pitch','Yaw','Location','best')
    xlabel('Time (s)')
    ylabel('Angle (deg)')
    
%   Save plots
    saveas(gcf, 'mpc_flight.png')
%%
    % East vs North
    figure()
    %plot(y(:,2), y(:,1)) %Not Fully Observable State Space Model
    plot(state(:,5), state(:,4),'r-.','Linewidth',5)
    hold on
    plot(ref_total(:,2),ref_total(:,1),'b','Linewidth',2);
    %plot(LQR(:,5), LQR(:,4)) % LQR if we can get it working
    plot(z_total(:,1), z_total(:,2),'k','Linewidth',2)
    quiver(X, Y, U, V)
    xticks([-500000:10000:500000])
    yticks([-500000:10000:500000])
    xlim([-500000 500000])
    ylim([-500000 500000])
    title('Aircraft Trajectory to Hurricane Center')
    legend('Aircraft Trajectory','Reference Trajectory','Hurricane Trajectory','Location','SouthEast')
    xlabel('East [m]')
    ylabel('North [m]')
%     axis equal

    figure()
    %plot3(y(:,2), y(:,1), y(:,3)) %Not Fully Observable State Space Model
    plot3(state(:,5), state(:,4), -1*state(:,6),'r-.','Linewidth',2)
    hold on
    plot3(ref_total(:,2),ref_total(:,1),-ref_total(:,3),'b','Linewidth',2);
    plot3(z_total(:,1), z_total(:,2),3000*ones(length(z_total(:,2)),1),'k','Linewidth',2)
    
    title('Aircraft Trajectory to Hurricane Center')
    legend('Aircraft Trajectory','Hurricane Trajectory','Location','SouthEast')
    xlabel('East [m]')
    ylabel('North [m]')
    zlabel('Altitude [m]')
    
    t_state = linspace(0,Tstop,length(state(:,1)))';
    figure()
    plot(t_state,(state(:,10)),'Linewidth',2)
    hold on
    plot(t_state,(state(:,11)),'Linewidth',2)
    plot(t_state,(state(:,12)),'Linewidth',2)
    title('Aircraft Roll Pitch Yaw')
    legend('Roll','Pitch','Yaw')
    xlabel('Time [sec]')
    ylabel('Angle [deg]')
    
    t_control = linspace(0,Tstop,length(control(:,1)))';
    figure()
    subplot(2,1,1)
    stairs(t_control,rad2deg(control(:,1)))
    hold on
    stairs(t_control,rad2deg(control(:,2)))
    stairs(t_control,rad2deg(control(:,3)))
    stairs(t_control,rad2deg(control(:,5)))
    stairs(t_control,rad2deg(control(:,6)))
    stairs(t_control,rad2deg(control(:,7)))
    title('Control Surface Angles')
    xlabel('Time [sec]')
    ylabel('Angle [deg]')
    legend('Elevator','Aileron','Rudder','Asymmetric Spoiler','Flap','Stabilator')
    
    subplot(2,1,2)
    stairs(t_control,control(:,4))
    title('Throttle')
    xlabel('Time [sec]')
    ylabel('Throttle [T]')
    
	%Flight Path History
	if SIMUL >= 1
		tspan	=	[ti tfinal];
		xo		=	x + delx
		u		=	u + delu;
		[t,x]	=	ode23('EoM',tspan,xo);

		figure
		subplot(2,2,1)
		plot(t,x(:,1))
		xlabel('Time, s'), ylabel('Axial Velocity, m/s'), grid
		subplot(2,2,2)
		plot(t,x(:,2))
		xlabel('Time, s'), ylabel('Side Velocity, m/s'), grid
		subplot(2,2,3)
		plot(t,x(:,3))
		xlabel('Time, s'), ylabel('Normal Velocity, m/s'), grid
		subplot(2,2,4)
		plot(t,x(:,4))
		xlabel('Time, s'), ylabel('North, m'), grid
	
		figure
		subplot(2,2,1)
		plot(t,x(:,5))
		xlabel('Time, s'), ylabel('East, m'), grid
		subplot(2,2,2)
		plot(t,-x(:,6))
		xlabel('Time, s'), ylabel('Altitude, m'), grid
		subplot(2,2,3)
		plot(t,x(:,7) * 57.29578)
		xlabel('Time, s'), ylabel('Roll Rate, deg/s'), grid
		subplot(2,2,4)
		plot(t,x(:,8) * 57.29578)
		xlabel('Time, s'), ylabel('Pitch Rate, deg/s'), grid

		figure
		subplot(2,2,1)
		plot(t,x(:,9) * 57.29578)
		xlabel('Time, s'), ylabel('Yaw Rate, deg/s'), grid
		subplot(2,2,2)
		plot(t,x(:,10) * 57.29578)
		xlabel('Time, s'), ylabel('Roll Angle, deg'), grid
		subplot(2,2,3)
		plot(t,x(:,11) * 57.29578)
		xlabel('Time, s'), ylabel('Pitch Angle, deg'), grid
		subplot(2,2,4)
		plot(t,x(:,12) * 57.29578)
		xlabel('Time, s'), ylabel('Yaw Angle, deg'), grid
		
		figure
		subplot(2,2,1)
		plot(t,x(:,1),t,x(:,2),t,x(:,3))
		xlabel('Time, s'), ylabel('Velocity, m/s'), grid
		subplot(2,2,2)
		plot(t,x(:,4),t,x(:,5),t,-x(:,6))
		xlabel('Time, s'), ylabel('Position, m'), grid
		subplot(2,2,3)
		plot(t,x(:,7) * 57.29578,t,x(:,8) * 57.29578,t,x(:,9) * 57.29578)
		xlabel('Time, s'), ylabel('Angular Rate, deg/s'), grid
		subplot(2,2,4)
		plot(t,x(:,10) * 57.29578,t,x(:,11) * 57.29578,t,x(:,12) * 57.29578)
		xlabel('Time, s'), ylabel('Angle, deg'), grid
		
		figure
		subplot(2,2,1)
		plot(t,x(:,1),t,x(:,3) * 10)
		xlabel('Time, s'), ylabel('u, 10w, m/s'), grid
		subplot(2,2,2)
		plot(t,-x(:,6))
		xlabel('Time, s'), ylabel('Altitude, m'), grid
		subplot(2,2,3)
		plot(t,x(:,8) * 57.29578)
		xlabel('Time, s'), ylabel('Pitch Rate, deg/s'), grid
		subplot(2,2,4)
		plot(t,x(:,11) * 57.29578)
		xlabel('Time, s'), ylabel('Pitch Angle, deg'), grid		
		
		figure
		subplot(2,2,1)
		plot(t,x(:,2))
		xlabel('Time, s'), ylabel('Side Velocity, m/s'), grid
		subplot(2,2,2)
		plot(t,x(:,5))
		xlabel('Time, s'), ylabel('Crossrange, m'), grid
		subplot(2,2,3)
		plot(t,x(:,7) * 57.29578,t,x(:,9) * 57.29578)
		xlabel('Time, s'), ylabel('Roll and Yaw Rates, deg/s'), grid
		subplot(2,2,4)
		plot(t,x(:,10) * 57.29578,t,x(:,12) * 57.29578)
		xlabel('Time, s'), ylabel('Roll and Yaw Angles, deg'), grid	
	
		figure
		subplot(1,2,1)
		plot(x(:,7) * 57.29578,x(:,8) * 57.29578)
		xlabel('Roll Rate, deg/s'), ylabel('Pitch Rate, deg/s'), grid
		subplot(1,2,2)
		plot(x(:,9) * 57.29578,x(:,8) * 57.29578)
		xlabel('Yaw Rate, deg/s'), ylabel('Pitch Rate, deg/s'), grid
		
    end