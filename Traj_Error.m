function [sys,x0,str,ts] = Traj_Error(t,x,u,flag)
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
    sizes.NumOutputs     = 2;  % E, phi_e
    sizes.NumInputs      = 12; % 期望轨迹6 + 实际状态6
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    % 输入解包：期望轨迹 (1-6)
    x_d    = u_in(1);
    y_d    = u_in(2);
    phi_d  = u_in(3);  
    xd_dot = u_in(4);  
    yd_dot = u_in(5);  
    phid_dot = u_in(6); 

    % 输入解包：实际状态 (7-12, 来自Transform)
    x_bar   = u_in(7);
    y_bar   = u_in(8);
    phi= u_in(9);  
    u_real  = u_in(10);
    v_bar   = u_in(11); 
    r_real  = u_in(12); 

    % 大地坐标系下的位置误差 (式2-17)
    xe_bar = x_d - x_bar;
    ye_bar = y_d - y_bar;

    % 船体坐标系下的位置误差 (式2-18)
    xe_tilde =  cos(phi)*xe_bar + sin(phi)*ye_bar;
    ye_tidle = -sin(phi)*xe_bar + cos(phi)*ye_bar;

    %% 综合位置误差 E (标量)
    E = sqrt(xe_bar^2 + ye_bar^2);

    %% 角度误差 phi_e — 统一使用 atan2，规避 x_e=0 分支缺陷
    phi_e = atan2(ye_tidle, xe_tilde);              % 范围 (-pi, pi]
    phi_e = mod(phi_e + pi, 2*pi) - pi;   % 规约到 [-pi, pi)，避免 x_e 过零跳变

    sys = [E; phi_e];
end
