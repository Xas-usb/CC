function [sys,x0,str,ts] = Disturb(t,x,u,flag)
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
    sizes.NumOutputs     = 3;  % tau_Du, tau_Dv, tau_Dr
    sizes.NumInputs      = 0;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [0 0];
end

function sys = mdlOutputs(t,x,u)
    
    tau_Du = -10 + 4*sin(0.5*t)*cos(0.5*t)-6*cos(t)*cos(0.5*t);
    tau_Dv =  5*sin(0.1*t);
    tau_Dr =  8*sin(1.1*t)*cos(0.3*t);

    sys = [tau_Du; tau_Dv; tau_Dr];
end
