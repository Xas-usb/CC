function [sys,x0,str,ts] = Transform(t,x,u,flag)
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
    sizes.NumOutputs     = 6;  % x_bar, y_bar, phi, u, v_bar, r
    sizes.NumInputs      = 6;  % x, y, phi, u, v, r (来自 Plant)
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    x_real   = u_in(1);
    y_real   = u_in(2);
    phi_real = u_in(3);
    u_real   = u_in(4);
    v_real   = u_in(5);
    r_real   = u_in(6);

    % 坐标变换常数 hbar = m23/m22 (式2-3)
    m22 = 33.8;
    m23 = 1.0948;
    hbar = m23 / m22;

    % 坐标变换 (式2-3)
    x_bar = x_real + hbar * cos(phi_real);
    y_bar = y_real + hbar * sin(phi_real);
    v_bar = v_real + hbar * r_real;

    sys = [x_bar; y_bar; phi_real; u_real; v_bar; r_real];
end
