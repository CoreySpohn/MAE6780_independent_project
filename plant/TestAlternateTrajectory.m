%   Hurricane Initial Conditions and Parameters
    hurr_para.linvel = 4.9; %m/s
    hurr_para.angvel = -.0001; % rad/s % This will be a very small number, but it will cause the hurricane path to arc
    hurr_para.Vmax = 252; %km/hr The units on this are important!
    hurr_para.Rmax = 47; % km The units on this are important!

    % Set up ODE 45 outside the iterative loop,
    hurr_prop = @(t,z) HurricaneSensingTrajectoryEOM(t,z,hurr_para);
%     n = 100;
%     tf = 10; % This should equal the control horizon
%     tarray = linspace(0,tf,n);
    smallnumber=1e-10;
    opts=odeset('Abstol',smallnumber,'RelTol',smallnumber);

    z_hurr = [ -5000, 5000, deg2rad(45),-2000,2000,3000]; % Initial hurricane conditions
    
    Tstop = 30000; % Total Simulation Time
    n = 500; % number of samples
    tarray = linspace(0,Tstop,n); 
    
    z = ode45(hurr_prop,tarray,z_hurr,opts);
    zarray = deval(z,tarray);
    
    hurr = zarray(1:2,:)';
    ref_traj = zarray(4:5,:)';
    
    figure()
    plot(hurr(:,1),hurr(:,2))
    hold on
    plot(ref_traj(:,1),ref_traj(:,2))
    title('Field Following Trajectory')
    legend('Hurricane Path','Reference Trajectory')
    xlabel('East [m]')
    ylabel('North [m]')
    