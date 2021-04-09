function usm = SMC(t, x, xn, un, tau, L, ls)
% Sliding-mode control
%   t   - Time
%   x   - State
%   xn  - Nominal state
%   un  - Nominal control
%   tau - Time decay constant on sliding surface
%   L   - Sliding-mode control scale
%   ls  - Length scale

    % Linearized dynamics
    global A B
%     thresh = 0.1 * ones(19, 1);
%     xj = [x; un];
% 	xdotj = LinModel(t,xj);
% 	[dFdX,~] = numjac('LinModel', t, xj, xdotj, thresh, [], 0);
%     A = dFdX(1:12,1:12);
% 	B = dFdX(1:12,13:19);
    
    
    % Sliding-surface matrices
    S = KinMat(x, tau, ls);
    Sn = KinMat(xn, tau, ls);
    
    % Sliding function
    sig = S*x - Sn*xn;
    
    % Sliding-mode control
    v = lsqminnorm(S*B, -L*sigmoid(sig));
    %usm(1:6) = un(1:6) + v;
    %usm(7) = un(7);
    usm = un + v;
    
end

% Sliding-surface matrix
function S = KinMat(x, tau, ls)

    phi = x(10);
    theta = x(11);
    psi = x(12);

    H = DCM(phi, theta, psi);
    
    L = [1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
         0, cos(phi),            -sin(phi);
         0, sin(phi)*sec(theta),  cos(phi)*sec(theta)];

    S = zeros(6, 12);
    
    S(1:3, 1:3)  = tau * H' / ls;
    S(1:3, 4:6)  = eye(3) / ls;
    S(4:6, 7:9)  = tau * L;
    S(4:6,10:12) = eye(3);
    
end

% Sigmoid function
function s = sigmoid(z)

    s = (2/pi) * atan(100*z);

end