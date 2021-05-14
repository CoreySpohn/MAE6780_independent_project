function usm = SMC(t, x, xn, un, Kr, Kth, L)
% Sliding-mode control
%   t   - Time
%   x   - State
%   xn  - Nominal state
%   un  - Nominal control
%   Kr, Kth - Constant matrices
%   L   - Sliding-mode control scale

    % Linearized dynamics
    global A B
%     thresh = 0.1 * ones(19, 1);
%     xj = [x; un];
% 	xdotj = LinModel(t,xj);
% 	[dFdX,~] = numjac('LinModel', t, xj, xdotj, thresh, [], 0);
%     A = dFdX(1:12,1:12);
% 	B = dFdX(1:12,13:19);
    
    
    % Sliding-surface matrices
    S = KinMat(x, Kr, Kth);
    Sn = KinMat(xn, Kr, Kth);
    
    % Sliding function
    sig = S*x - Sn*xn;
    
    % Sliding-mode control
    v = lsqminnorm(S*B, -L*sigmoid(sig));
    usm(1:4) = un(1:4) + v(1:4);
    usm(5) = un(5);
    usm(6:7) = un(6:7) + v(5:6);
    
end

% Sliding-surface matrix
function S = KinMat(x, Kr, Kth)

    phi = x(10);
    theta = x(11);
    psi = x(12);

    H = DCM(phi, theta, psi);
    
    L = [1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
         0, cos(phi),            -sin(phi);
         0, sin(phi)*sec(theta),  cos(phi)*sec(theta)];

    S = zeros(6, 12);
    
    S(1:3, 1:3)  = H';
    S(1:3, 4:6)  = Kr;
    S(4:6, 7:9)  = L;
    S(4:6,10:12) = Kth;
    
end

% Sigmoid function
function s = sigmoid(z)

    s = (2/pi) * atan(10*z);
    
end