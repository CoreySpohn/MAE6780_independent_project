function z_hurr_dot = HurricaneSensingTrajectoryEOM(t,z_hurricane_prev,hurr_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hurricane Dynamics
%
%   Corey Spohn and Rachel Oliver
%   MAE6780 - Multivariable Controls
%
%   To calculate the desired trajectory of a sensing mission at the maximum
%   velocity radius, essentially sampling the hurricane eye. We will switch
%   to this reference trajectory once the aircraft reaches the Rmax from
%   the hurricane center.
%
%
%   Inputs:
%       z_hurricane_prev    Previous inertial frame state of the hurricane
%                           in east and north directions as well as the 
%                           direction of travel with relation to the fixed 
%                           inertial frame [East , North , theta]
%       hurr_para           Struture containing the hurricane parameters.
%                           The only one of use here is the inertial
%                           velocity of the hurricane. This is treated as a
%                           a constant in the hurricane fixed frame. There
%                           is an linear and angular velocity component.
%                           You can think of this as similar to a
%                           non-holonomic robot motion. Set the angular
%                           velocity to 0 for stright line motion.
%       t_step              This is the control horizon for the MPC
%                           controller in this scenario. Though it could be
%                           any time step that you wish to conduct euler's
%                           method.       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack Hurricane Parameters
v = hurr_para.linvel; v_th = hurr_para.angvel; Rmax = hurr_para.Rmax*1000;

% Unpack Previous Hurricane State
z = z_hurricane_prev(1:3)'; theta = z_hurricane_prev(3);

% Unpack Previous Reference State
% z_ref_pos = z_hurricane_prev(7:9)'; 
% z_ref_vel = z_hurricane_prev(4:6)'; 
z_ref_pos = z_hurricane_prev(4:6)'; 

% Create DCM to convert movement in hurricance fixed frame back to the
% inertial one
R_HI = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

if v_th ==0
    z_dot = R_HI*[v;0;v_th];
else
    % Use definition of an arc to determine the components of change 
    z_dot = R_HI*[2*(v/v_th)*sin(v_th/2)*sin((pi-v_th)/2);2*(v/v_th)*sin(v_th/2)*cos((pi-v_th)/2);v_th];
end

% The aircraft will now travel on an arc in either the positive or negative
% direction around the hurricane. So it needs to adjust with the hurricane
% movement and use the formula for arc length s = r*theta to find the new
% position of the aircraft with respect to the hurricane center
r = Rmax;
%s_dot = norm(z_ref_vel(1:2));
s_dot = hurr_para.maxVelAircraft;
% Angle Between Hurricane Center and Aircraft Position
theta = atan2((z_ref_pos(1)-z(2)),(z_ref_pos(2)-z(1)));
theta_dot = s_dot/Rmax;
theta_new = theta + theta_dot;
new_pos_y = sin(theta_new)*Rmax;
new_pos_x = cos(theta_new)*Rmax;
local_dot = [(z(1)+new_pos_x) - z_ref_pos(2);...
             (z(2)+new_pos_y) - z_ref_pos(1);...
             0];


% Update reference trajectory by translating with the hurricane center and moving in the negative gradient
z_ref_pos_dot = z_dot + local_dot;
    
z_hurr_dot = [z_dot; z_ref_pos_dot(2); z_ref_pos_dot(1); z_ref_pos_dot(3)];
end