function usm = SMC(t, x, xn, un, tau, L, ls)
% Sliding-mode control
%   t   - Time
%   x   - State
%   xn  - Nominal state
%   un  - Nominal control
%   tau - Time decay constant on sliding surface
%   L   - Sliding-mode control scale
%   ls  - Length scale

    % Linearized model
    global B
    
    % Sliding-surface matrices
    S = KinMat(x, tau, ls);
    Sn = KinMat(xn, tau, ls);
    
    % Sliding function
    sig = S*x - Sn*xn;
    
    % Sliding-mode control
    usm = un + lsqminnorm(S*B, -L*sigmoid(sig));
    
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

    s = (2/pi) * atan(z*100);

end