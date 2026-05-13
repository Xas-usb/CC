function [sys,x0,str,ts] = Aux_System(t,x,u,flag)
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
    sizes.NumContStates  = 2;  % Lambda_u, Lambda_r
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2;
    sizes.NumInputs      = 4;  % u_e, r_e, delta_tau_u, delta_tau_r
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [0.02; 0.02];  % 文献初值
    str = [];
    ts  = [0 0];
end

function sys = mdlDerivatives(t,x,u_in)
    u_e         = u_in(1);
    r_e         = u_in(2);
    delta_tau_u = u_in(3);
    delta_tau_r = u_in(4);

    Lambda_u = x(1);
    Lambda_r = x(2);

    % 控制参数
    Tp = 20; gamma = 0.35; mu = 0.01; ell_Lambda = 0.01;
    zeta_u = 1; zeta_r = 1;  

    m11 = 25.8; m22 = 33.8; m23 = 1.0948; m32 = 1.0948; m33 = 2.76;
    xi_u = 1 / m11;
    xi_r = m22 / (m22*m33 - m23*m32);

    % 预定义时间收敛项 (式2-48)
    term_PT_u = (pi/(gamma*Tp)) * ...
        (sign(Lambda_u)*abs(Lambda_u)^(1-gamma) + sign(Lambda_u)*abs(Lambda_u)^(1+gamma));
    term_PT_r = (pi/(gamma*Tp)) * ...
        (sign(Lambda_r)*abs(Lambda_r)^(1-gamma) + sign(Lambda_r)*abs(Lambda_r)^(1+gamma));

    dL_u = -ell_Lambda * u_e - term_PT_u;
    dL_r = -ell_Lambda * r_e - term_PT_r;

    % 饱和补偿项：仅当 |Lambda| >= mu 时激活 (式2-48下半部分)
    % 注意：分母 Lambda^2 加 1e-6 保护，防止 Lambda→0 时除法爆炸
    eps_div = 1e-6;  %

    if abs(Lambda_u) >= mu
        comp_u = (zeta_u * xi_u * u_e * delta_tau_u * Lambda_u) / (Lambda_u^2 + eps_div);
        dL_u = dL_u + comp_u;
    end

    if abs(Lambda_r) >= mu
        comp_r = (zeta_r * xi_r * r_e * delta_tau_r * Lambda_r) / (Lambda_r^2 + eps_div);
        dL_r = dL_r + comp_r;
    end

    % 导数饱和限幅：防止瞬态大信号导致 Lambda 积分爆炸
    dL_max = 1e4;
    dL_u = max(min(dL_u, dL_max), -dL_max);
    dL_r = max(min(dL_r, dL_max), -dL_max);

    sys = [dL_u; dL_r];
end

function sys = mdlOutputs(t,x,u_in)
    sys = x;  % [Lambda_u; Lambda_r]
end
