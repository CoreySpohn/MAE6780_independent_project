
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
Vmax = hurr_para.Vmax; % km/h
Rmax = hurr_para.Rmax*1000; % Convert from km to m
xmax = hurr_para.xmax;
ymax = hurr_para.ymax;
noise = hurr_para.noise; % Currently generated inside the TestSMC file
sf = hurr_para.scalefactor;

% Assign Radial Velocity
% Calculate the euclidean distance between the hurricane center and the
% aircraft position
r = ((x(5)-z(1))^2+(x(4)-z(2))^2)^.5;
theta = atan2((x(4)-z(2)),(x(5)-z(1)));

% Calculate the magnitude of the radial velocity, Vr, at the radial distance
if r<Rmax
    Vr = Vmax*(r/Rmax)^(3/2); %Scalar multiplier of the vector field
    % Eye Wall Altitude Scaling Factors
    scale = [0 .92 1.1 1.15 1.21 1.15 1.05 1 1 1];
else
    Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
    % Outer Vortex Altitude Scaling Factors
    scale = [0 .8 .94 1.01 1.08 1.06 1.05 1 1 1];
end

% generate noise that is large enough to account for the sampling location
noise_pt = get_noise(x(5),x(4),noise,xmax,ymax,sf);
if isnan(noise_pt)
    noise_pt = 0;
end
%noise_pt


% Use Trigonometry to Split Vr into x and y components
if x(5) == z(1) && x(4) == z(2)
    Vt = cross([0,0,0],[0,0,1]);
elseif x(5) == z(1)
    Vt = cross([0,((Vr+noise_pt)*.277778*(x(4)-z(2)))/abs((x(4)-z(2))),0],[0,0,1]);
elseif x(4) == z(2)
    Vt = cross([(Vr+noise_pt)*.277778*(x(5)-z(1))/abs((x(5)-z(1))),0,0],[0,0,1]);
else  
    Vt = cross([(Vr+noise_pt)*.277778*(x(5)-z(1))/abs((x(5)-z(1))),((Vr+noise_pt)*.277778*(x(4)-z(2)))/abs((x(4)-z(2))),0],[0,0,1]);
end
    
% Cross with <0,0,1> to get the tangential component
% Vt = cross([dx,dy,0],[0,0,1]);

% Calculate the magnitude of the vertical velocity
% Most hurricanes need very little wind shear to form so
dz = 0;

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
	windx	=	scale.*[Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2) Vt(2)];	% Northerly wind, m/s
	windy	=	scale.*[Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1) Vt(1)];	% Easterly wind, m/s
	windz	=	[dz dz dz dz dz dz dz dz dz dz];	% Vertical wind. m/s
    height = abs(x(6));
	winde	=	[interp1(windh,windx,height)
				interp1(windh,windy,height)
				interp1(windh,windz,height)];	% Earth-relative frame
	HEB		=	DCM(x(10),x(11),x(12));
	windb	=	HEB * winde;					% Body-axis frame
