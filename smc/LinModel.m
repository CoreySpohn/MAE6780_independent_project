	function xdotj = LinModel(tj,xj)
%	Equations of Motion for Linear Model (Jacobian) Evaluation,
%	with dummy state elements added for controls

	global u
	
	x		=	xj(1:12);
	u		=	xj(13:19);

	xdot	=	EoM(tj,x);
	xdotj	=	[xdot;0;0;0;0;0;0;0];
