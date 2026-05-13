function [sys,x0,str,ts] = Error_Constrains(t,x,u,flag)
    switch flag
        case 0
            [sys,x0,str,ts] = mdlInitializeSizes;
        case 3
            sys = mdlOutputs(t,x,u);
        case {1,2,4,9}
            sys = [];
        otherwise
            error(['Unhandled flag = ',num2str(flag)]);
    end
end

function [sys,x0,str,ts] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;  % sig1, sig2, th1, th2, th1_dot, th2_dot
    sizes.NumInputs      = 5;  % E, phi_e
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    %% --- 常数定义 ---
    E     = u_in(1);
    phi_e = u_in(2);

    % 性能函数参数 (表2-2)
    Tf = 20;                % 预设收敛时间
    th10 = 10; th1_inf = 0.2;  % 位置误差边界
    th20 = 10; th2_inf = 0.05; % 角度误差边界
    k1 = 10; l1 = 300;
    k2 = 10; l2 = 300;

    %% --- 1. 性能函数 theta_i (式2-26) ---
    if t < Tf  
        exp_term = exp(1-Tf/(Tf-t));
        theta1 = (th10 - th1_inf) * exp_term + th1_inf;
        theta2 = (th20 - th2_inf) * exp_term + th2_inf;

        % theta 导数 (需要传递到虚拟控制律)
        theta1_dot = -((th10 - th1_inf) * Tf / ((Tf - t)^2 )) * exp_term;
        theta2_dot = -((th20 - th2_inf) * Tf / ((Tf - t)^2 )) * exp_term;
    else
        % t >= Tf: 保持在稳态值
        theta1 = th1_inf;
        theta2 = th2_inf;
        theta1_dot = 0;
        theta2_dot = 0;
    end

    %% --- 2. 误差转换 sigma_i (式2-29) ---
   
    sigma1 = (k1*E)/sqrt(E^2+l1 );  
    sigma2 = (k2*phi_e)/sqrt(phi_e^2+l2);

    sys = [sigma1; sigma2; theta1; theta2; theta1_dot; theta2_dot];
end
