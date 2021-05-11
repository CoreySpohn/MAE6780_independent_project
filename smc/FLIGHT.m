%	FLIGHT --  6-DOF Trim, Linear Model, and Flight Path Simulation 
%	August 4, 1999   
%	===============================================================
%	Copyright 1993-1999 by ROBERT F. STENGEL.  All rights reserved.

	clear
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

	windb	=	WindField(-ze,phir,thetar,psir);
	alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;

	x	=	[V * cos(alphar) * cos(betar) - windb(1)
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
	
	u	=	[dE * .01745329
			dA * .01745329
			dR * .01745329
			dT
			dAS * .01745329
			dF * .01745329
			dS * .01745329];

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
		Fmodel		=	dFdX(1:12,1:12)
		Gmodel		=	dFdX(1:12,13:19)
		save Fmodel
		save Gmodel
	end

%	Flight Path History
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