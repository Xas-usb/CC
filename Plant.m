function [sys,x0,str,ts] = Plant(t,x,u,flag)
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
    sizes.NumContStates  = 6;  % x, y, phi, u, v, r
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;  % 同上
    sizes.NumInputs      = 5;  % tau_fu, tau_fr, tau_Du, tau_Dv, tau_Dr
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [-5; -5; 0; 0; 0; 0];  % 初始位置 (x,y)=(-20,10), 航向角0, 速度0
    str = [];
    ts  = [0 0];
end

%% 核心：USV 动力学导数 (式2-1, 2-2)
function sys = mdlDerivatives(t,x,u_in)
    pos_x = x(1); pos_y = x(2);
    phi   = x(3); u     = x(4); v = x(5); r = x(6);

    tau_fu = u_in(1); tau_fr = u_in(2);
    tau_Du = u_in(3); tau_Dv = u_in(4); tau_Dr = u_in(5);

    % 模型参数 (表2-1)
    m11 = 25.8; m22 = 33.8; m23 = 1.0948;
    m32 = 1.0948; m33 = 2.76;

    % 水动力阻尼 (速度相关)
    d11 = 5.8664*u^2 + 1.3274*abs(u) + 0.7225;
    d22 = 0.805*abs(r) + 36.2823*abs(v) + 0.8612;
    d23 = abs(r) + 0.845*abs(v) - 0.1079;
    d32 = -0.13*abs(r) - 5.0437*abs(v) - 0.1052;
    d33 = 0.75*abs(r) - 0.08*abs(v) + 1.9;

    detM = m22*m33 - m23*m32;

    % 运动学 (式2-1)
    dx = u*cos(phi) - v*sin(phi);
    dy = u*sin(phi) + v*cos(phi);
    dphi = r;

    % 动力学 (式2-2)
    du_dot = (1/m11) * (m22*v*r + m23*r^2 - d11*u + tau_Du + tau_fu);

    dv_dot = -(m23/detM) * (m11*u*v - m22*u*v - m23*u*r - d32*v - d33*r + tau_Dr + tau_fr) ...
             - (m33/detM) * (m11*u*r + d22*v + d23*r - tau_Dv);

    dr_dot = (m22/detM) * (m11*u*v - m22*u*v - m23*u*r - d32*v - d33*r + tau_Dr + tau_fr) ...
             + (m32/detM) * (m11*u*r + d22*v + d23*r - tau_Dv);

    sys = [dx; dy; dphi; du_dot; dv_dot; dr_dot];
end


function sys = mdlOutputs(t,x,u_in)
    
    sys = x;  % [x, y, phi, u, v, r]
end
