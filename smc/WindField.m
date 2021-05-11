	function windb = windField(height,phir,thetar,psir)
%	Wind Field Interpolation for 3-D wind as a Function of Altitude
%	Copyright 1993-1999 by ROBERT F. STENGEL.  All rights reserved.

	windh	=	[-10 0 100 200 500 1000 2000 4000 8000 16000];	% Wind-Height, m
	windx	=	[0 0 0 0 0 0 0 0 0 0];	% Northerly wind, m/s
	windy	=	[0 0 0 0 0 0 0 0 0 0];	% Easterly wind, m/s
	windz	=	[0 0 0 0 0 0 0 0 0 0];	% Vertical wind. m/s
	
	winde	=	[interp1(windh,windx,height)
				interp1(windh,windy,height)
				interp1(windh,windz,height)];	% Earth-relative frame
	HEB		=	DCM(phir,thetar,psir);
	windb	=	HEB * winde;					% Body-axis frame


