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
    sizes.NumOutputs     = 5;  
    sizes.NumInputs      = 14; 
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
    dxd = u_in(4);  
    dyd = u_in(5);  
    phid_dot = u_in(6); 

    % 输入解包：实际状态 (7-12, 来自Transform)
    x_bar   = u_in(7);
    y_bar   = u_in(8);
    phi= u_in(9);  
    u_real  = u_in(10);
    v_bar   = u_in(11); 
    r  = u_in(12); 
    dx1 = u_in(13);
    dy = u_in(14);

    % 大地坐标系下的位置误差 (式2-17)
    xe_bar = x_d-x_bar;
    ye_bar = y_d-y_bar;

    % 船体坐标系下的位置误差 (式2-18)
    xe_tilde =  cos(phi)*xe_bar + sin(phi)*ye_bar;
    ye_tilde = -sin(phi)*xe_bar + cos(phi)*ye_bar;
 %% 坐标变换求导
    m22 = 33.8;
    m23 = 1.0948;
    hbar = m23 / m22;
    dphi = r;
    phi_e = atan2(ye_tilde,xe_tilde);
    dx_bar=dx1-hbar*sin(phi)*dphi;
    dy_bar=dy+hbar*cos(phi)*dphi;

    dxe_bar=dxd-dx_bar;
    dye_bar=dyd-dy_bar;

    dxe_tilde=-sin(phi)*dphi*xe_bar+cos(phi)*dxe_bar+cos(phi)*dphi*ye_bar+sin(phi)*dye_bar;
    dye_tilde=-cos(phi)*dphi*xe_bar-sin(phi)*dxe_bar-sin(phi)*dphi*ye_bar+cos(phi)*dye_bar;

    dphie=(dye_tilde*xe_tilde-ye_tilde*dxe_tilde)/(xe_tilde^2+ye_tilde^2);

    phid=phi+phi_e;
    dphid=dphi+dphie;


    %% 综合位置误差 E (标量)
    E = sqrt(xe_bar^2 + ye_bar^2);
     
    sys = [E; phi_e;phid;dphid;v_bar];
end
