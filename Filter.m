function [sys,x0,str,ts] = Filter(t,x,u,flag)
    switch flag
        case 0
            [sys,x0,str,ts] = mdlInitializeSizes;
        case 1
            sys = mdlDerivatives(t,x,u);
        case 3
            sys = mdlOutputs(t,x,u);
        case {2,4,9}
            sys = [];
        otherwise
            error(['Unhandled flag = ',num2str(flag)]);
    end
end

function [sys,x0,str,ts] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 2;  % alpha_bar_u, alpha_bar_r
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;  % dot_alpha_bar_u, dot_alpha_bar_r, alpha_bar_u, alpha_bar_r, u_e, r_e
    sizes.NumInputs      = 13; % 来自Transform(6) + VirtualCtrl(7)
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);

    % 初始状态设为期望轨迹初速度 (u_d=1.5, r_d=0.04)
    x0  = [1.5; 0.04];
    str = [];
    ts  = [0 0];
end

%% 核心：预定义时间滤波器导数 (式2-43)
function alpha_bar_dot = filter_dynamics(alpha, alpha_bar, Tp, gamma)
    e_f = alpha - alpha_bar;
    % sign(x)*|x|^p 结构 — 确保实数运算
    sig_minus = sign(e_f) * abs(e_f)^(1 - gamma);
    sig_plus  = sign(e_f) * abs(e_f)^(1 + gamma);
    alpha_bar_dot = (pi / (gamma * Tp)) * (sig_minus + sig_plus);
end

function sys = mdlDerivatives(t,x,u_in)
    alpha_u = u_in(7);
    alpha_r = u_in(8);
    a1 = u_in(9);
    a2 = u_in(10);
    c1 = u_in(11);
    c2 = u_in(12);
    phi_e = u_in(13);
    Tp = 20; gamma = 0.35;

    e_fu = alpha_u - x(1);
    e_fr = alpha_r - x(2);
    
    c1_safe = max(c1, 1e-6);
    c2_safe = max(c2, 1e-6);


    % 式(2-43) 完整滤波器：预定义时间项 + 耦合补偿(a/c) + 阻尼(-0.5*e_f)
    alpha_bar_u_dot = filter_dynamics(alpha_u, x(1), Tp, gamma)...
                      + a1*cos(phi_e)/c1_safe - 0.5*e_fu;
    alpha_bar_r_dot = filter_dynamics(alpha_r, x(2), Tp, gamma) ...
                      + a2/c2_safe - 0.5*e_fr;

    % 导数限幅：防止 c→0 时 a/c 补偿项过大导致数值爆炸
    dot_max = 1e4;
    alpha_bar_u_dot = max(min(alpha_bar_u_dot, dot_max), -dot_max);
    alpha_bar_r_dot = max(min(alpha_bar_r_dot, dot_max), -dot_max);

    sys = [alpha_bar_u_dot; alpha_bar_r_dot];
end

function sys = mdlOutputs(t,x,u_in)
    u_real = u_in(4);
    r_real = u_in(6);

    alpha_u = u_in(7);
    alpha_r = u_in(8);
    a1 = u_in(9);
    a2 = u_in(10);
    c1 = u_in(11);
    c2 = u_in(12);
    phi_e = u_in(13);
    Tp = 20; gamma = 0.35;

    % 零保护
    c1_safe = max(c1, 1e-6);
    c2_safe = max(c2, 1e-6);
    e_fu = alpha_u - x(1);
    e_fr = alpha_r - x(2);

    % 式(2-43) 完整滤波器输出（与导数函数保持一致）
    alpha_bar_u_dot = filter_dynamics(alpha_u, x(1), Tp, gamma) ...
                      + a1*cos(phi_e)/c1_safe - 0.5*e_fu;
    alpha_bar_r_dot = filter_dynamics(alpha_r, x(2), Tp, gamma) ...
                      + a2/c2_safe - 0.5*e_fr;

    % 导数限幅（与导数函数保持一致）
    dot_max = 1e4;
    alpha_bar_u_dot = max(min(alpha_bar_u_dot, dot_max), -dot_max);
    alpha_bar_r_dot = max(min(alpha_bar_r_dot, dot_max), -dot_max);

    % 速度误差 (式2-44)
    u_e = x(1) - u_real;
    r_e = x(2) - r_real;

    sys = [alpha_bar_u_dot; alpha_bar_r_dot; x(1); x(2); u_e; r_e];
end
