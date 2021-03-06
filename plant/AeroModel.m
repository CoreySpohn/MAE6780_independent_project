	function [CD,CL,CY,Cl,Cm,Cn,Thrust]	=	AeroModel(x,u,Mach,alphar,betar,V)
%	Aerodynamic Coefficients of the Aircraft, Thrust Model,
%	and Geometric and Inertial Properties

%	Business Jet -- Low-Angle-of-Attack, Mach-Dependent Model
%	Landing Gear: Up or Down (GEAR  = 1 or 0)
%	Flap Setting, u(6): 0? or 38? (0 or 0.66 rad)
%	Symmetric Spoiler: Deployed or Closed (SPOIL = 1 or 0)

	global m Ixx Iyy Izz Ixz S b c GEAR CONTROL SPOIL

%	Typical Mass and Inertial Properties
	m		=	4536;				% Mass, kg
	Ixx		=	35926.5;			% Roll Moment of Inertia, kg-m^2
	Iyy		=	33940.7;			% Pitch Moment of Inertia, kg-m^2
	Izz		=	67085.5;			% Yaw Moment of Inertia, kg-m^2
	Ixz		=	3418.17;			% Nose-high(low) Product of Inertia, kg-m^2
	
%	Geometric Properties
	c		=	2.14;				% Mean Aerodynamic Chord, m
	b		=	10.4;				% Wing Span, m
	S		=	21.5;				% Reference Area, m^2
	ARw		=	5.02;				% Wing Aspect Ratio
	taperw	=	0.507;				% Wing Taper Ratio
	sweepw	=	13 * .01745329;		% Wing 1/4-chord sweep angle, rad
	ARh		=	4;					% Horizontal Tail Aspect Ratio
	sweeph	=	25 * .01745329;		% Horiz Tail 1/4-chord sweep angle, rad
	ARv		=	0.64;				% Vertical Tail Aspect Ratio
	sweepv	=	40 * .01745329;		% Vert Tail 1/4-chord sweep angle, rad
	lvt		=	4.72;				% Vert Tail Length, m
	
%	Thrust Properties
	StaticThrust	=	26243.2;	% Static Thrust @ Sea Level, N
	
%	Current Thrust
	[airDens,airPres,temp,soundSpeed] = Atmos(-x(6));
	Thrust			=	u(4) * StaticThrust * (airDens / 1.225)^0.7;
									% Thrust at Altitude, N
	
%	Current Mach Effects, normalized to Test Condition B (Mach = 0.1734)
	PrFac			=	1 / (sqrt(1 - Mach^2) * 1.015);	
									% Prandtl Factor
	WingMach		=	1 / ((1 + sqrt(1 + ((ARw/(2 * cos(sweepw)))^2) ...
						* (1 - Mach^2 * cos(sweepw)))) * 0.268249);
									% Modified Helmbold equation
	HorizTailMach	=	1 / ((1 + sqrt(1 + ((ARh/(2 * cos(sweeph)))^2) ...
						* (1 - Mach^2 * cos(sweeph)))) * 0.294539);
									% Modified Helmbold equation
	VertTailMach	=	1 / ((1 + sqrt(1 + ((ARv/(2 * cos(sweepv)))^2) ...
						* (1 - Mach^2 * cos(sweepv)))) * 0.480338);
									% Modified Helmbold equation
	
%	Current Longitudinal Characteristics
%	====================================

%	Lift Coefficient
	CLo		=	0.1095;				% Zero-AoA Lift Coefficient (B)
	if GEAR >= 1
		CLo	=	CLo - 0.0192;		% Gear-down correction
	end
	if u(6) >= 0.65
		CLo	=	CLo + 0.5182;		% 38?-flap correction
	end
	if SPOIL >= 1
		CLo	=	CLo - 0.1897;		% 42?-Symmetric Spoiler correction
	end	
	
	CLar	=	5.6575;				% Lift Slope (B), per rad
	if u(6) >= 0.65
		CLar	=	CLar - 0.0947;
	end

	CLqr		=	4.231 * c / (2 * V);
									% Pitch-Rate Effect, per rad/s
	
	CLdSr	=	1.08;				% Stabilator Effect, per rad
	if u(6) >= 0.65
		CLdSr	=	CLdSr - 0.4802;	% 38?-flap correction
	end
	
	CLdEr	=	0.5774;				% Elevator Effect, per rad
	if u(6) >= 0.65
		CLdEr	=	CLdEr - 0.2665;	% 38?-flap correction
	end

	CL	=	CLo + (CLar*alphar + CLqr*x(8) + CLdSr*u(7) + CLdEr*u(1)) ...
			* WingMach;
								% Total Lift Coefficient, w/Mach Correction
	
%	Drag Coefficient
	CDo		=	0.0255;				% Parasite Drag Coefficient (B)
	if GEAR >= 1
		CDo	=	CDo + 0.0191;		% Gear-down correction
	end
	if u(6) >= 0.65
		CDo	=	CDo + 0.0836;		% 38?-flap correction
	end
	if SPOIL >= 1
		CDo	=	CDo + 0.0258;		% 42?-Symmetric Spoiler correction
	end	

	epsilon	=	0.0718;				% Induced Drag Factor
	if u(6) >= 0.65
		epsilon	=	0.079;			% 38?-flap correction
	end

	CD	=	CDo * PrFac + epsilon * CL^2;
							% Total Drag Coefficient, w/Mach Correction
	
%	Pitching Moment Coefficient
	Cmo		=	0;					% Zero-AoA Moment Coefficient (B)
	if GEAR >= 1
		Cmo	=	Cmo + 0.0255;		% Gear-down correction
	end
	if u(6) >= 0.65
		Cmo	=	Cmo - 0.058;		% 38?-flap correction
	end
	if SPOIL >= 1
		Cmo	=	Cmo - 0.0154;		% 42?-Symmetric Spoiler correction
	end	
	
	Cmar	=	-1.231;				% Static Stability (B), per rad
	if u(6) >= 0.65
		Cmar	=	Cmar + 0.0138;
	end

	Cmqr		=	 -18.8 * c / (2 * V);
							% Pitch-Rate + Alpha-Rate Effect, per rad/s
	
	CmdSr	=	-2.291;				% Stabilator Effect, per rad
	if u(6) >= 0.65
		CmdSr	=	CmdSr + 0.121;	% 38?-flap correction
	end
	
	CmdEr	=	-1.398;				% Elevator Effect, per rad
	if u(6) >= 0.65
		CmdEr	=	CmdEr + 0.149;	% 38?-flap correction
	end

	Cm	=	Cmo + (Cmar*alphar + Cmqr*x(8) + CmdSr*u(7) + CmdEr*u(1)) ...
			* HorizTailMach;
					% Total Pitching Moment Coefficient, w/Mach Correction
	
%	Current Lateral-Directional Characteristics
%	===========================================

%	Side-Force Coefficient
	CYBr	=	-0.7162;			% Side-Force Slope (B), per rad
	if u(6) >= 0.65
		CYBr	=	CYBr + 0.0826;
	end

	CYdAr	=	-0.00699;			% Aileron Effect, per rad
	
	CYdRr	=	0.1574;				% Rudder Effect, per rad
	if u(6) >= 0.65
		CYdRr	=	CYdRr - 0.0093;	% 38?-flap correction
	end
	
	CYdASr	=	0.0264;				% Asymmetric Spoiler Effect, per rad
	if u(6) >= 0.65
		CYdASr	=	CYdASr + 0.0766;	
									% 38?-flap correction
	end

	CY	=	(CYBr*betar + CYdRr*u(3)) * VertTailMach ... 
			+ (CYdAr*u(2) + CYdASr*u(5)) * WingMach;
						% Total Side-Force Coefficient, w/Mach Correction

%	Yawing Moment Coefficient
	CnBr	=	0.1194;				% Directional Stability (B), per rad
	if u(6) >= 0.65
		CnBr	=	CnBr - 0.0092;
	end

	Cnpr	=	CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V));				
									% Roll-Rate Effect, per rad/s
	
	Cnrr	=	(-2 * (lvt / b) * CnBr * VertTailMach - 0.1 * CL^2) ...
				* (b / (2 * V));				
									% Yaw-Rate Effect, per rad/s

	CndAr	=	0;					% Aileron Effect, per rad
	if u(6) >= 0.65
		CndAr	=	CndAr + 0.0028;
	end
	
	CndRr	=	-0.0713;			% Rudder Effect, per rad
	if u(6) >= 0.65
		CndRr	=	CndRr - 0.0185;	% 38?-flap correction
	end
	
	CndASr	=	-0.0088;			% Asymmetric Spoiler Effect, per rad
	if u(6) >= 0.65
		CndASr	=	CndASr - 0.0106;	
									% 38?-flap correction
	end

	Cn	=	(CnBr*betar + CndRr*u(3)) * VertTailMach ...
			+ Cnrr * x(9) + Cnpr * x(7) ...
			+ (CndAr*u(2) + CndASr*u(5)) * WingMach;
					% Total Yawing-Moment Coefficient, w/Mach Correction

%	Rolling Moment Coefficient
	ClBr	=	-0.0918;				% Dihedral Effect (B), per rad
	if u(6) >= 0.65
		ClBr	=	ClBr - 0.0092;
	end

	Clpr	=	-CLar * (1 + 3 * taperw)/(12 * (1 + taperw)) ...
				* (b / (2 * V));				
										% Roll-Rate Effect, per rad/s
	
	Clrr	=	(CL * (1 + 3 * taperw)/(12 * (1 + taperw)) ...
				* ((Mach * cos(sweepw))^2 - 2) / ((Mach * cos(sweepw))^2 - 1)) ...
				* (b / (2 * V));				
										% Yaw-Rate Effect, per rad/s

	CldAr	=	0.07686;				% Aileron Effect, per rad
	if u(6) >= 0.65
		CldAr	=	CldAr + 0.00589;
	end
	
	CldRr	=	0.01208;				% Rudder Effect, per rad
	if u(6) >= 0.65
		CldRr	=	CldRr + 0.01115;	% 38?-flap correction
	end
	
	CldASr	=	-0.01496;				% Asymmetric Spoiler Effect, per rad
	if u(6) >= 0.65
		CldASr	=	CldASr - 0.02376;	
									% 38?-flap correction
	end

	Cl	=	(ClBr*betar + CldRr*u(3)) * VertTailMach ... 
			+ Clrr * x(9) + Clpr * x(7) ...
			+ (CldAr*u(2) + CldASr*u(5)) * WingMach;
					% Total Rolling-Moment Coefficient, w/Mach Correction
