function z_hurr_dot = HurricaneSensingTrajectoryEOM(t,z_hurricane_prev,hurr_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hurricane Dynamics
%
%   Corey Spohn and Rachel Oliver
%   MAE6780 - Multivariable Controls
%
%   To calculate the desired trajectory of a sensing mission within one
%   vector field of the hurricane, the aircraft 
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
v = hurr_para.linvel; v_th = hurr_para.angvel;

% Unpack Previous Hurricane State
z = z_hurricane_prev(1:3)'; theta = z_hurricane_prev(3);

% Unpack Previous Reference State
z_ref = z_hurricane_prev(4:6)';

% Calculate Windfield at Previous Hurricane State and 
windb	=	HurricaneWindField3Point(z_ref,z,hurr_para);

% Create DCM to convert movement in hurricance fixed frame back to the
% inertial one
R_HI = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

if v_th ==0
    z_dot = R_HI*[v;0;v_th];
else
    % Use definition of an arc to determine the components of change 
    z_dot = R_HI*[2*(v/v_th)*sin(v_th/2)*sin((pi-v_th)/2);2*(v/v_th)*sin(v_th/2)*cos((pi-v_th)/2);v_th];
end

% Update reference trajectory by translating with the hurricane center and moving in the negative gradient
z_ref_dot = z_dot-windb;
    
z_hurr_dot = [z_dot;z_ref_dot];
end