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
    sizes.NumOutputs     = 8;  % x_bar, y_bar, phi, u, v_bar, r
    sizes.NumInputs      = 8;  % x, y, phi, u, v, r, dx, dy
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    x_real = u_in(1);
    y_real = u_in(2);
    phi = u_in(3);
    u0 = u_in(4);
    v_real = u_in(5);
    r = u_in(6);
    dx1 = u_in(7);
    dy= u_in(8);

    %% 坐标变换常数 hbar = m23/m22 (式2-3)
    m22 = 33.8;
    m23 = 1.0948;
    hbar = m23 / m22;

    %% 坐标变换 
    x_bar = x_real + hbar * cos(phi);
    y_bar = y_real + hbar * sin(phi);
    v_bar = v_real + hbar * r;



    sys = [x_bar; y_bar; phi; u0; v_bar; r; dx1;dy];
end
