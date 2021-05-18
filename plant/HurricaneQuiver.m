function [X, Y, U, V] = HurricaneQuiver(x_grid, y_grid, hurr_state, hurr_para)
% Unpack the hurricane's properties
Vmax = hurr_para.Vmax;
Rmax = hurr_para.Rmax*1000;
xmax = hurr_para.xmax;
ymax = hurr_para.ymax;
sf = hurr_para.scalefactor;

% Get the location of the hurricane's location
x_h = hurr_state(1);
y_h = hurr_state(2);

% Create the U and V vectors
X = zeros(length(x_grid), length(y_grid));
Y = zeros(length(x_grid), length(y_grid));
U = zeros(length(x_grid), length(y_grid));
V = zeros(length(x_grid), length(y_grid));

for i = 1:length(x_grid)
    x = x_grid(i);
    for j = 1:length(y_grid)
        % calculate distance from the center of the hurricane
        y = y_grid(j);
        r = sqrt((x - x_h)^2 + (y- y_h)^2);
        theta = atan2((y-y_h),(x-x_h));
        if r < Rmax
            Vr = Vmax*(r/Rmax)^(3/2); %Scalar multiplier of the vector field
            % Eye Wall Altitude Scaling Factors
            scale = [0 .92 1.1 1.15 1.21 1.15 1.05 1 1 1];
        else
            Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
            % Outer Vortex Altitude Scaling Factors
            scale = [0 .8 .94 1.01 1.08 1.06 1.05 1 1 1];
        end
        X(i, j) = x;
        Y(i, j) = y;
        U(i, j) = Vr*sin(theta);
        V(i, j) = -Vr*cos(theta);
    end
end
