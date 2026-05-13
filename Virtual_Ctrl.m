function [sys,x0,str,ts] = Virtual_Ctrl(t,x,u,flag)
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
    sizes.NumOutputs     = 7;  % alpha_u, alpha_r, a1, a2, c1, c2, phi_e
    sizes.NumInputs      = 17; % 期望轨迹6 + 误差2 + v_bar + 约束6+phid，dphid
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    %% --- 1. 信号解包 ---
    % 期望轨迹信息 (u_in 1-6)
    x_d      = u_in(1);
    y_d      = u_in(2); 
    phi_d    = u_in(3);
    xd_dot   = u_in(4);
    yd_dot   = u_in(5);
    phid_dot = u_in(6);

    % 误差信息 (u_in 7-8)
    E     = u_in(7);
    phi_e = u_in(8);
    phid = u_in(9);
    dphid=u_in(10);
    % 辅助状态 (u_in 9)
    v_bar = u_in(11);

    % 误差约束信息 (u_in 10-15)
    sig1    = u_in(12);
    sig2    = u_in(13);
    th1     = u_in(14);
    th2     = u_in(15);
    th1_dot = u_in(16);
    th2_dot = u_in(17);
  

    %% --- 2. 控制参数 (表2-2) ---
    Tp = 20; gamma = 0.35;
    K1 = 0.1; K2 = 1;
    k1 = 10; k2 = 10; l1 = 300; l2 = 300;  % l1,l2 与 Error_Constrains 保持一致 (论文表2-2)

    rho_u = sig1 / th1;
    rho_r = sig2 / th2;

    %% --- 3. 计算辅助参数 a, c, d (式2-40) ---
    a1 = rho_u * k1 * l1 * th1;
    a2 = rho_r * k2 * l2 * th2;

    c1 = (1-rho_u^2)^2*th1^2*(E^2+l1)^1.5;
    c2 = (1-rho_r^2)^2*th2^2*(phi_e^2+l2)^1.5;

    d1 = rho_u*th1_dot*sig1*(E^2+l1)^1.5;
    d2 = rho_r*th2_dot*sig2*(phi_e^2+l2)^1.5;

    %% --- 保护未使用
     % a1_safe = sign(a1) * max(abs(a1), 1e-10);
     % a2_safe = sign(a2) * max(abs(a2), 1e-10);

    %% --- 6. Lyapunov 函数 V (式2-30, 分母加极小量防除零) ---
    Vu = 0.5*rho_u^2/(1-rho_u^2);
    Vr = 0.5*rho_r^2/(1-rho_r^2);

    %% --- 7. 虚拟控制律 alpha (式2-42) ---
    term_Vu = (pi / (gamma * Tp)) * (Vu^(1 - 0.5*gamma) + Vu^(1 + 0.5*gamma));
    term_Vr = (pi / (gamma * Tp)) * (Vr^(1 - 0.5*gamma) + Vr^(1 + 0.5*gamma));
  

    alpha_u = (xd_dot*cos(phid)+yd_dot*sin(phid)-v_bar*sin(phi_e) + ...
               (c1/a1)*term_Vu-(d1/a1)+K1*E)/cos(phi_e);

    alpha_r = dphid+c2/a2*term_Vr-d2/a2+K2*phi_e;

    sys = [alpha_u; alpha_r; a1; a2; c1; c2; phi_e];
end
