%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hurricane Dynamics and Vector Field
%
%   Corey Spohn and Rachel Oliver
%   MAE6780 - Multivariable Controls
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars; clc
% Hurricane Properties
Vmax = 252; %km/h
Rmax = 47; %km

[X,Y] = meshgrid(-240:20:240);
for i = 1:length(X(1,:))
    for j = 1:length(Y(:,1))
        r = (X(1,i)^2+Y(j,1)^2)^.5;
        if r<Rmax
            Vr = Vmax*(r/Rmax)^(3/2);%Scalar multiplier of the vector field
        else
            Vr = Vmax*(2*Rmax*r)/(r^2+Rmax^2);
        end
        dx(j,i) = Vr*(Y(j,1));        %/(X.^2+Y.^2));
        dy(j,i) = Vr*(-X(1,i));       %/(X.^2+Y.^2));
    end
end

%

figure
% contour(X,Y,dx,dy)
% hold on
quiver(X,Y,dx,dy,'AutoScaleFactor',1)
hold off
axis equal
title('Hurricane Wind Field')
xlabel('East [km]')
ylabel('North [km]')