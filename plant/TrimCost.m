	function J = TrimCost(optPar)
%	Cost Function for Longitudinal Trim in Steady Level Flight
%	Copyright 1993-1999 by ROBERT F. STENGEL.  All rights reserved.

	global m Ixx Iyy Izz Ixz S b c GEAR CONTROL u x V parhis

	R	=	[1 0 0
			0 1 0
			0 0 1];

% Optimization Vector:
%	1 = Stabilator, rad
%	2 = Throttle, %
%	3 = Pitch Angle, rad
	optPar;
				
	u	=	[u(1)
			u(2)
			u(3)
			optPar(2)
			u(5)
			u(6)
			optPar(1)];
				
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
			x(12)];
	
	xdot	=	EoM(1,x);
	xCost	=	[xdot(1)
				xdot(3)
				xdot(8)];
	J		=	xCost' * R * xCost;
	pars	=	[optPar;J];
	parhis	=	[parhis pars];
	
	


