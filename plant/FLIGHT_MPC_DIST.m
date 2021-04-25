%	FLIGHT --  6-DOF Trim, Linear Model, and Flight Path Simulation 
%	August 4, 1999   
%	===============================================================
%	Copyright 1993-1999 by ROBERT F. STENGEL.  All rights reserved.

%	close all; clearvars; clc
    set(0,'DefaultAxesFontSize',17);
	global GEAR CONTROL SPOIL u x V parhis

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
	LINEAR = 	1;		% Linear model flag (= 1 to calculate F and G)
	p =			0;		% Body-axis roll rate, deg/s
	phi =		0;		% Body roll angle wrt earth, deg
	psi =		0;		% Body yaw angle wrt earth, deg
	q	=		0;		% Body-axis pitch rate, deg/sec
	r =			0;		% Body-axis yaw rate, deg/s
	SIMUL =		0;		% Flight path flag (= 1 for nonlinear simulation)
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
	delx	=	[0;0;0
				0;0;0
				0;0;0
				0;0;0];
				
	delu	=	[0;0;0;0;0;0;0];

%	State Vector and Control Initialization
	phir	=	phi * .01745329;
	thetar	=	theta * .01745329;
	psir	=	psi * .01745329;
    
%   Initialize State Before the First Loop
    x	=	[V * cos(alpha) * cos(beta)
        V * sin(beta)
        V * sin(alpha) * cos(beta)
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

%   Hurricane Initial Conditions and Parameters
    hurr_para.linvel = 4.9; %m/s
    hurr_para.angvel = 0; % rad/s % This will be a very small number, but it will cause the hurricane path to arc
    hurr_para.Vmax = 252; %km/hr The units on this are important!
    hurr_para.Rmax = 47; % km The units on this are important!

    % Set up ODE 45 outside the iterative loop,
    hurr_prop = @(t,z) HurricaneEOM(t,z,hurr_para);
    n = 100;
    tf = 10; % This should equal the control horizon
    tarray = linspace(0,tf,n);
    smallnumber=1e-10;
    opts=odeset('Abstol',smallnumber,'RelTol',smallnumber);

    z_hurr = [ -5000, 5000, deg2rad(45)]; % Initial hurricane conditions
%   North position of center of mass WRT Earth, xe, m
%	East position of center of mass WRT Earth, ye, m
%   Angle of Hurricane Path with respect to the Earth Frame's East Axis

%   Start of the Iterative Loop
    loop_steps = 500; %Number of Control Horizon Iterations
    state = []; % Set up variable to contain the state variable after each loop 
    control = []; % Set up variable to contain the control variables after each loop
    z_total = []; % Set up variable to contain the hurricane track after each loop
for i = 1:loop_steps
    % Calculate the Current Wind Disturbances
	windb	=	HurricaneWindField(x,z_hurr,hurr_para);
	alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;

	x	=	x - [windb(1);...
                 windb(2);...
			     windb(3);...
                 zeros(9,1)]
	
	u	=	u

%	Trim Calculation (for Steady Level Flight at Initial V and h)
	if TRIM >= 1
		'TRIM'
		parhis	=	[];
		optPar				=	[-.05;.19;.06];
		options				=	foptions;
		[optPar,options]	=	fmins('TrimCost',optPar,options);
%		Optimizing dSr, dT, Theta, and Error Cost
		optPar
		J					=	options(8)
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

%	Stability-and-Control Derivative Calculation
   	if LINEAR >= 1
		'LINEAR'
		thresh	=	[.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1;.1];
		xj		=	[x;u];
		xdotj		=	LinModel(ti,xj);
		[dFdX,fac]	=	numjac('LinModel',ti,xj,xdotj,thresh,[],0);
		Fmodel		=	dFdX(1:12,1:12);
		Gmodel		=	dFdX(1:12,13:19);
		save Fmodel
		save Gmodel
    end
%     C = [1 1 1 1 1 1 1 1 1 1 1 1];

    %C = [0 0 0 1 0 0 0 0 0 0 0 0;
%          0 0 0 0 1 0 0 0 0 0 0 0;
%          0 0 0 0 0 1 0 0 0 0 0 0;
%          0 0 0 0 0 0 1 0 0 0 0 0;
%          0 0 0 0 0 0 0 1 0 0 0 0;
%          0 0 0 0 0 0 0 0 1 0 0 0];
    %D = zeros(6, 7);
    C = eye(12);
    D = zeros(12,7);
    CSTR = ss(Fmodel, Gmodel, C, D);

    
%     CSTR.InputName = {'T_c','C_A_i'};
%     CSTR.OutputName = {'T','C_A'};
%     CSTR.StateName = {'C_A','T'};
%     CSTR.InputGroup.MV = 1; % Manipulated variable
%     CSTR.InputGroup.UD = 2; % Unmeasured distubances
%     CSTR.OutputGroup.MO = 1; % Measured outputs
%     CSTR.OutputGroup.UO = 2; % Unmeasured outputs
    Ts = 1;
    mpc_obj = mpc(CSTR, Ts);
    % specify prediction horizon
    mpc1.PredictionHorizon = 2000; % sec
    % specify control horizon
    control_step = 10; % sec
    mpc1.ControlHorizon = control_step;
    % Starting location of the aircraft
    %mpc_obj.Model.Nominal.Y = [1000, 1000, 20000, 0, 0, 0]; %Not Fully Observable State Space Model
    mpc_obj.Model.Nominal.Y = [0, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0];
    % Starting control input
    mpc_obj.Model.Nominal.U = [0, 0, 0, 250, 0, 0, 0];
    % Throttle max and min
    % mpc_obj.MV(4).Max = 400; % This constraint should actually be on the
    % axial velocity state 1
    % mpc_obj.MV(4).Min = 0;
    % Control Surface Angles limited to +or- 30 deg = 0.5236
    mpc_obj.MV(1).Max = 0.5236;
    mpc_obj.MV(1).Min = -0.5236;
    mpc_obj.MV(2).Max = 0.5236;
    mpc_obj.MV(2).Min = -0.5236;
    mpc_obj.MV(3).Max = 0.5236;
    mpc_obj.MV(3).Min = -0.5236;
    mpc_obj.MV(5).Max = 0.5236;
    mpc_obj.MV(5).Min = -0.5236;
    mpc_obj.MV(6).Max = 0.5236;
    mpc_obj.MV(6).Min = -0.5236;
    mpc_obj.MV(7).Max = 0.5236;
    mpc_obj.MV(7).Min = -0.5236;
    % Output Minimums i.e. We cannot fly into the ground, We cannot fly
    % faster than 400 mph, pitch, roll, and yaw are limited as well
    mpc_obj.OV(1).Max = 175; % maximum axial velocity, 400 mph = 175 m/s
    mpc_obj.OV(1).Min = 0;
    mpc_obj.OV(6).Max = 9000; % meter, maximum operational range
    mpc_obj.OV(6).Min = 1000;
    mpc_obj.OV(10).Max = 3*pi/2; % Roll rad, large plane will not barrel roll
    mpc_obj.OV(10).Min = -3*pi/2;
    mpc_obj.OV(11).Max = pi/2; % Pitch rad, the plane cannot flip over the minor axis
    mpc_obj.OV(11).Min = -pi/2;
    % Set Weights for Q and R
    mpc_obj.Weights.MV = [.8 .8 .8 .8 .8 .8 .8];
    mpc_obj.Weights.MVRate = [0.1 0.1 0.1 0.1 0.1 0.1 0.1];
    mpc_obj.Weights.OV = [0 0 0 .6 .6 .1 0 0 0 .5 .5 .5];
    mpc_obj.Weights.ECR = 100000;
    mpc_refsignal = zeros(11, 7);
    
    % LQR for comparison
%     Q_lqr = [0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 .6 0 0 0 0 0 0 0 0;
%           0 0 0 0 .6 0 0 0 0 0 0 0;
%           0 0 0 0 0 .1 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 .1 0 0;
%           0 0 0 0 0 0 0 0 0 0 .1 0;
%           0 0 0 0 0 0 0 0 0 0 0 .1];
%     R_lqr = .1*eye(7);
%     K1 = lqr(CSTR,Q_lqr,R_lqr);
    

%   Huricane Modeling
    % Run an ODE 45 Solver to calculate the hurricane trjectory over the
    % current control horizon
    z = ode45(hurr_prop,tarray,z_hurr,opts);
    zarray = deval(z,tarray);
    z_total = [z_total;zarray'];
    z_hurr = zarray(:,end)'; % Current Hurricane State

    ref_signal = [zeros(3,1);...
                  z_hurr(2);... % Hurricane Center Position, North
                  z_hurr(1);... % Hurricane Center Position, East
                  3000;... % Desired Flight Height
                  zeros(6,1)]';
    
    options = mpcsimopt();
    options.RefLookAhead = 'off';
    options.MDLookAhead = 'off';
    options.Constraints = 'on';
    options.OpenLoop = 'off';
    %sim(mpc_obj, 1, ref_signal, [], options)
    [y, t, u_t, xp, xc, output_options] = sim(mpc_obj, 1, ref_signal, [], options);
    x = y'
    state = [state;y];
    u = u_t'
    control = [control;u_t];
end
    
%     LQR_sys = ss(Fmodel-Gmodel*K1, Gmodel, C, D);
%     LQR = sim(LQR_sys,zeros(length(tarray),7),tarray,[0, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 0, 0]')
    
%%
    % East vs North
    figure()
    %plot(y(:,2), y(:,1)) %Not Fully Observable State Space Model
    plot(state(:,5), state(:,4),'r-.','Linewidth',2)
    hold on
    %plot(LQR(:,5), LQR(:,4)) % LQR if we can get it working
    plot(z_total(:,1), z_total(:,2),'b','Linewidth',2)
    title('Aircraft Trajectory to Hurricane Center')
    legend('Aircraft Trajectory','Hurricane Trajectory','Location','SouthEast')
    xlabel('East [m]')
    ylabel('North [m]')
     
    figure()
    %plot3(y(:,2), y(:,1), y(:,3)) %Not Fully Observable State Space Model
    plot3(state(:,5), state(:,4), state(:,6),'r-.','Linewidth',2)
    hold on
    plot3(z_total(:,1), z_total(:,2),3000*ones(length(z_total(:,2)),1),'b','Linewidth',2)
    
    title('Aircraft Trajectory to Hurricane Center')
    legend('Aircraft Trajectory','Hurricane Trajectory','Location','SouthEast')
    xlabel('East [m]')
    ylabel('North [m]')
    zlabel('Altitude [m]')
    
    figure()
    plot(t,(y(:,10)),'Linewidth',2)
    hold on
    plot(t,(y(:,11)),'Linewidth',2)
    plot(t,(y(:,12)),'Linewidth',2)
    title('Aircraft Roll Pitch Yaw')
    legend('Roll','Pitch','Yaw')
    xlabel('Time [sec]')
    ylabel('Angle [deg]')
    
    
    figure()
    subplot(2,1,1)
    line(t,rad2deg(control(:,1)))
    hold on
    line(t,rad2deg(control(:,2)))
    line(t,rad2deg(control(:,3)))
    line(t,rad2deg(control(:,5)))
    line(t,rad2deg(control(:,6)))
    line(t,rad2deg(control(:,7)))
    title('Control Surface Angles')
    xlabel('Time [sec]')
    ylabel('Angle [deg]')
    legend('Elevator','Aileron','Rudder','Asymmetric Spoiler','Flap','Stabilator')
    
    subplot(2,1,2)
    line(t,control(:,4))
    title('Throttle')
    xlabel('Time [sec]')
    ylabel('Throttle [T]')
    
	%Flight Path History
	if SIMUL >= 1
		tspan	=	[ti tf];
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
