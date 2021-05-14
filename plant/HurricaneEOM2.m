function z_hurr_dot = HurricaneEOM(t,z_hurricane_prev,hurr_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hurricane Dynamics
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
v = hurr_para.linvel; v_th = hurr_para.angvel; goal_alt = hurr_para.goalaltitude;

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
    z_dot = R_HI*[2*(v/v_th)*sin(v_th/2)*sin((pi-v_th)/2);2*(v/v_th)*sin(v_th/2)*cos((pi-v_th)/2);vth];
end

% The reference trajectory will be a straight line drawn from the aircraft
% current position to the hurricane center.
z_new = z' + z_dot;
%r_dot = norm(z_ref_vel(1:2)); % Velocity Magnitude of aircraft
r_dot = hurr_para.maxVelAircraft;
direction = [(z_new(1) - z_ref_pos(2))/norm([z_new(1),z_ref_pos(2)]); (z_new(2) - z_ref_pos(1))/norm([z_new(2),z_ref_pos(1)]); (goal_alt - z_ref_pos(3))/norm([goal_alt,z_ref_pos(3)])];
z_ref_pos_dot = r_dot*direction;

z_hurr_dot = [z_dot; z_ref_pos_dot(2); z_ref_pos_dot(1); z_ref_pos_dot(3)];
end