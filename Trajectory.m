function [sys,x0,str,ts] = Trajectory(t,x,u,flag)
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
    sizes.NumContStates  = 3;  % xd, yd, phid
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 6;  % xd, yd, phid, xd_dot, yd_dot, phid_dot
    sizes.NumInputs      = 0;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [0; 0; 0];  % 期望轨迹起点 (0,0), 航向角 0
    str = [];
    ts  = [0 0];
end

function sys = mdlDerivatives(t,x,u_in)
    phi_d = x(3);

    % 期望速度：ud=1.5 m/s, rd=0.04 rad/s → 圆形轨迹
    u_d = 1.5;
    v_d = 0;
    r_d = 0.04;

    % 运动学 (式2-72)
    xd_dot = u_d * cos(phi_d) - v_d * sin(phi_d);
    yd_dot = u_d * sin(phi_d) + v_d * cos(phi_d);
    phid_dot = r_d;

    sys = [xd_dot; yd_dot; phid_dot];
end

function sys = mdlOutputs(t,x,u_in)
    xd = x(1); yd = x(2); phi_d = x(3);

    u_d = 1.5; v_d = 0; r_d = 0.04;

    xd_dot = u_d * cos(phi_d) - v_d * sin(phi_d);
    yd_dot = u_d * sin(phi_d) + v_d * cos(phi_d);
    phid_dot = r_d;

    sys = [xd; yd; phi_d; xd_dot; yd_dot; phid_dot];
end
