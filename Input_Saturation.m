function [sys,x0,str,ts] = Input_Saturation(t,x,u,flag)
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
    sizes.NumOutputs     = 2;  % tau_su, tau_sr (饱和后)
    sizes.NumInputs      = 2;  % taucu, taucr (来自 Main_Controller)
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u_in)
    tau_u = u_in(1);
    tau_r = u_in(2);

    %% --- 执行器限幅参数 ---
    % 推进器：正向推力 0~500 N 
    tau_u_max =  500;
    tau_u_min = 0;   

    % 舵机：±50 N·m
    tau_r_max =  50;
    tau_r_min = -50;

    %% --- 纵向饱和 (式2-15) ---
    tau_bu = (tau_u_max + tau_u_min)/2 + ((tau_u_max - tau_u_min)/2) * sign(tau_u);
    tau_su = tau_bu * tanh(2*t) * tanh(tau_u / tau_bu);

    %% --- 转向饱和 (双向对称) ---
    tau_br = (tau_r_max+tau_r_min)/2+((tau_r_max-tau_r_min)/2)*sign(tau_r);  % 取正幅值作为缩放基准
    tau_sr = tau_br * tanh(2*t) * tanh(tau_r / tau_br);

    sys = [tau_su; tau_sr];
end
