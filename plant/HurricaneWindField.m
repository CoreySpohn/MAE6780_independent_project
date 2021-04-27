function windb = HurricaneWindField(x,z,hurr_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hurricane Wind Field Calculation
%
%   Corey Spohn and Rachel Oliver
%   MAE6780 - Multivariable Controls
%
%   Unlike the normal wind field calculation for the flight model, this
%   needs to be updated at each time step and the inertial position of the
%   aircraft and the hurricane matters. The hurricane current conditions
%   are propagated through a simple euler's equation and the previous time
%   step. The inertial position of the hurricane is an input into this
%   function to create a windfield that is centered on the hurricane. The
%   inertial position of the aircraft with relation to the hurricane
%   determines the gradient of the windfield in the ECI frame.
%
%
%   Inputs:
%       x   State of the Aircraft.
%       z   State of the Hurricane.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack Hurricane Properties
Vmax = hurr_para.Vmax;
Rmax = hurr_para.Rmax*1000; % Convert from km to m

% Assign Radial Velocity
% Calculate the euclidean distance between the hurricane center and the
% aircraft position
r = ((x(5)-z(1))^2+(x(4)-z(2))^2)^.5;
theta = atan2((x(4)-z(2)),(x(5)-z(1)));

% Calculate the magnitude of the radial velocity, Vr, at the radial distance
if r<Rmax
    Vr = Vmax*(r/Rmax)^(3/2); %Scalar multiplier of the vector field
else
    Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
end

% Create Noise of +-5m/s for each portion of velocity
noise = -Vr/2 + (Vr/2+Vr/2)*rand(1,1);

% Use Trigonometry to Split Vr into x and y components
% dx = Vr*0.277778*cos(theta)+noise(1); %Convert from km/h to m/s
% dy = Vr*0.277778*sin(theta)+noise(2); %Convert from km/h to m/s
Vt = cross([(Vr*.277778+noise)*(x(5)-z(1))/abs((x(5)-z(1))),(Vr*.277778+noise)*(x(4)-z(2))/abs((x(4)-z(2))),0],[0,0,1]);

% Cross with <0,0,1> to get the tangential component
% Vt = cross([dx,dy,0],[0,0,1]);

% Calculate the magnitude of the vertical velocity
% This seems to be a function of distance from the eye wall and height in
% the hurricane - Initially set to be 15 m/s universally
noise = -5 + (5+5)*rand(1,1);
dz = -15+noise; %m/s

% [X,Y] = meshgrid(-240:20:240);
% for i = 1:length(X(1,:))
%     for j = 1:length(Y(:,1))
%         r = (X(1,i)^2+Y(j,1)^2)^.5;
%         if r<Rmax
%             Vr = Vmax*(r/Rmax)^(3/2);%Scalar multiplier of the vector field
%         else
%             Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
%         end
%         vect = Vr*[(X(1,i)-z(1)),(Y(j,1)-z(2)),0];
%         dxdy = cross(vect,[0,0,1]);
%         dx(j,i) = dxdy(1);        %/(X.^2+Y.^2));
%         dy(j,i) = dxdy(2);       %/(X.^2+Y.^2));
%     end
% end
% 
% %
% 
% figure
% contour(X,Y,dx,dy)
% hold on
% quiver(X,Y,dx,dy,'AutoScaleFactor',1)
% hold off
% axis equal
% title('Hurricane Wind Field')
% xlabel('East [km]')
% ylabel('North [km]')

% Assume the wind is uniform in height
    windh	=	[-10 0 100 200 500 1000 2000 4000 8000 16000];	% Wind-Height, m
	windx	=	[Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2)];	% Northerly wind, m/s
	windy	=	[Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1)];	% Easterly wind, m/s
	windz	=	[dz dz dz dz dz dz dz dz dz dz];	% Vertical wind. m/s
    height = abs(x(6));
	winde	=	[interp1(windh,windx,height)
				interp1(windh,windy,height)
				interp1(windh,windz,height)];	% Earth-relative frame
	HEB		=	DCM(x(10),x(11),x(12));
	windb	=	HEB * winde;					% Body-axis frame
