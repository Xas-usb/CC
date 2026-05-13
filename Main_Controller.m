function [sys,x0,str,ts] = Main_Controller(t,x,u,flag)
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
    sizes.NumContStates  = 4;  
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2;   % tau_u, tau_r (未饱和)
    sizes.NumInputs      = 23;  % Filter(6) + Error约束(6);phi_e;E;phid;dphid;vbar;F1F2;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);

    x0  = [0.02; 0.02; 0; 0];
    str = [];
    ts  = [0 0];  % 连续时间
end

function sys = mdlDerivatives(t,x,u_in)
  
    dalphau_bar = u_in(1);
    dalphar_bar = u_in(2);
    alpha_bar_u = u_in(3);
    alpha_bar_r = u_in(4);
    u_e         = u_in(5);
    r_e         = u_in(6);
    sig1 = u_in(7);
    sig2 = u_in(8);
    th1  = u_in(9);
    th2  = u_in(10);

    E     = u_in(13); 
    phi_e = u_in(14);
    phid = u_in(15);
    dphid = u_in(16);
    v_bar = u_in(17);
    F1 = u_in(18);
    F2 = u_in(21);
   
    u1 = alpha_bar_u - u_e;
    r = alpha_bar_r - r_e;
   
    m11=25.8;
    m22=33.8;
    m23=1.0948;
    m32=1.0948;
    m33=2.76;
    hbar = m23 / m22;
    gamma = 0.35; Tp=20;
    l1 = 0.5; l2 = 0.5;
    Psi1 = 20; Psi2 = 20;

    d11=5.8664*u1^2+1.3274*abs(u1)+0.7225;
    d22=0.805*abs(r)+36.2823*abs(v_bar)+0.8612;
    d23=0.845*abs(v_bar)+3.45*abs(r)-0.1079;
    d32=-5.0437*abs(v_bar)-0.13*abs(r)-0.1052;
    d33=0.75*abs(r)-0.08*abs(v_bar)+1.9;


    sigamu=1/m11;
    sigamr=m22/(m22*m33-m23*m32);
    f1=sigamu*(m22*(v_bar-hbar*r)*r+m23*r^2-d11*u1);
    delta=m22*m33-m23*m32;
    f3=((m11*m22-m22^2)*u1*(v_bar-hbar*r)+(m11*m32-m23*m22)*u1*r- ...
    (d33*r+d32*(v_bar-hbar*r))*m22+(d23*r+d22*(v_bar-hbar*r))*m23)/delta;

    klu = 0.01; le = 300;
    ke = 10; kf = 10;
    klr = 0.01; lf = 300;
    Pi1=0.5; Pi2=0.5;
    K3=30; K4=30;

etaE=ke*E/sqrt(E^2+le);  
etaf=kf*phi_e/sqrt(phi_e^2+lf);

pu=etaE/th1;
pr=etaf/th2;

a1=pu*ke*le*th1;
a2=pr*kf*lf*th2;

c1=(1-pu^2)^2*th1^2*(E^2+le)*sqrt(E^2+le);
c2=(1-pr^2)^2*th2^2*(phi_e^2+lf)*sqrt(phi_e^2+lf);

PTSU=pi/(gamma*Tp)*(sign(u_e)*(abs(u_e))^(1-gamma)+sign(u_e)*(abs(u_e))^(1+gamma));
PTSR=pi/(gamma*Tp)*(sign(r_e)*(abs(r_e))^(1-gamma)+sign(r_e)*(abs(r_e))^(1+gamma));
if t>=0 && t<=50
    hfu=1;
    kfu=0;
else
    hfu=0.6+0.6*exp(-0.8*t);
    kfu=0.5+0.8*cos(t);

end

if t>=0 && t<=70
    hfr=1;
    kfr=0;
else
    hfr=0.4+0.6*exp(-0.2*t);
    kfr=0.6+0.8*cos(t);
end

taocu=(dalphau_bar-F1+x(3)*tanh(u_e/Pi1)+PTSU+(a1/c1)*cos(phi_e)+K3*u_e-klu*x(1))/(hfu*sigamu);
taocr=(dalphar_bar-F2+x(4)*tanh(r_e/Pi2)+PTSR+a2/c2+K4*r_e-klr*x(2))/(hfr*sigamr);

    tau_u_max =  500;
    tau_u_min = 0;   

    % 舵机：±50 N·m
    tau_r_max = 50;
    tau_r_min = -50;

    taomu=(tau_u_max+tau_u_min)/2+((tau_u_max-tau_u_min)/2)*sign(taocu);
    taomr=(tau_r_max+tau_r_min)/2+((tau_r_max-tau_r_min)/2)*sign(taocr);
    taou=tanh(2*t)*taomu*tanh(1*taocu/(1*taomu));
    taor=tanh(2*t)*taomr*tanh(1*taocr/(1*taomr));
%% 饱和误差

    deltataou=taou-taocu;
    deltataor=taor-taocr;

   

%% 扰动
    taowu=-10+4*sin(0.5*t)*cos(0.5*t)-6*cos(t)*cos(0.5*t);
    taowv=5*sin(0.1*t);
    taowr=8*cos(0.3*t)*sin(1.1*t);

    du=sigamu*taowu;
    dr=sigamr*(taowr-hbar*taowv);


PTSLU=-pi/(gamma*Tp)*sign(x(1))*abs(x(1))^(1-gamma)-pi/(gamma*Tp)*sign(x(1))*abs(x(1))^(1+gamma);
PTSLr=-pi/(gamma*Tp)*sign(x(2))*abs(x(2))^(1-gamma)-pi/(gamma*Tp)*sign(x(2))*abs(x(2))^(1+gamma);

if abs(x(1)) >= 0.01
   dlambdau=-klu*u_e+(hfu*sigamu*u_e*deltataou)/(abs(x(1))^2)+PTSLU;
else
   dlambdau=-klu*u_e+PTSLU;
end

if abs(x(2)) >= 0.01
   dlambdar=-klr*r_e+(hfr*sigamr*r_e*deltataor)/(abs(x(2))^2)+PTSLr;
else
   dlambdar=-klr*r_e+PTSLr; 
end

PTSUP1=(2-gamma)*pi/(gamma*Tp)*sign(x(3))*abs(x(3))^(1-gamma)+(2+gamma)*pi/(gamma*Tp)*sign(x(3))*abs(x(3))^(1+gamma);
PTSRP2=(2-gamma)*pi/(gamma*Tp)*sign(x(4))*abs(x(4))^(1-gamma)+(2+gamma)*pi/(gamma*Tp)*sign(x(4))*abs(x(4))^(1+gamma);

drho1=Psi1*(tanh(u_e/l1)*u_e-PTSUP1);
drho2=Psi2*(tanh(r_e/l2)*r_e-PTSRP2);


sys = [dlambdau;dlambdar;drho1;drho2];
end

%% ===== 控制律输出 =====
function sys = mdlOutputs(t,x,u_in)
    
    dalphau_bar = u_in(1);
    dalphar_bar = u_in(2);
    alpha_bar_u     = u_in(3);
    alpha_bar_r     = u_in(4);
    u_e             = u_in(5);
    r_e             = u_in(6);

    sig1 = u_in(7);
    sig2 = u_in(8);
    th1  = u_in(9);
    th2  = u_in(10);

    E     = u_in(13); 
    phi_e = u_in(14);
    phid = u_in(15);
    dphid = u_in(16);

    v_bar = u_in(17);
    F1 = u_in(18);
    F2 = u_in(19);
   
    u1 = alpha_bar_u - u_e;
    r = alpha_bar_r - r_e;

    Tp = 20; gamma = 0.35;
    K3 = 30; K4 = 30;
    le = 300; lf = 300;
    ke = 10; kf = 10;
%% 故障
if t>=0 && t<=50
    hfu=1;
    kfu=0;
else
    hfu=0.6+0.6*exp(-0.8*t);
    kfu=0.5+0.8*cos(t);

end

if t>=0 && t<=70
    hfr=1;
    kfr=0;
else
    hfr=0.4+0.6*exp(-0.2*t);
    kfr=0.6+0.8*cos(t);
end
%% 参数
    m11 = 25.8; m22 = 33.8; m23 = 1.0948; 
    m32 = 1.0948; m33 = 2.76;
    hbar = m23 / m22;
    xi_u = 1 / m11;
    xi_r = m22 / (m22*m33 - m23*m32);

d11=5.8664*u1^2+1.3274*abs(u1)+0.7225;
d22=0.805*abs(r)+36.2823*abs(v_bar)+0.8612;
d23=0.845*abs(v_bar)+3.45*abs(r)-0.1079;
d32=-5.0437*abs(v_bar)-0.13*abs(r)-0.1052;
d33=0.75*abs(r)-0.08*abs(v_bar)+1.9;

Pi1 = 0.5; 
Pi2 = 0.5;
klu=0.01;
klr=0.01;

sigamu=1/m11;
sigamr=m22/(m22*m33-m23*m32);
f1=sigamu*(m22*(v_bar-hbar*r)*r+m23*r^2-d11*u1);
delta=m22*m33-m23*m32;
f3=((m11*m22-m22^2)*u1*(v_bar-hbar*r)+(m11*m32-m23*m22)*u1*r- ...
    (d33*r+d32*(v_bar-hbar*r))*m22+(d23*r+d22*(v_bar-hbar*r))*m23)/delta;
    
etaE=ke*E/sqrt(E^2+le);  
etaf=kf*phi_e/sqrt(phi_e^2+lf);

pu=etaE/th1;
pr=etaf/th2;

a1=pu*ke*le*th1;
a2=pr*kf*lf*th2;

c1=(1-pu^2)^2*th1^2*(E^2+le)*sqrt(E^2+le);
c2=(1-pr^2)^2*th2^2*(phi_e^2+lf)*sqrt(phi_e^2+lf);

PTSU=pi/(gamma*Tp)*(sign(u_e)*(abs(u_e))^(1-gamma)+sign(u_e)*(abs(u_e))^(1+gamma));
PTSR=pi/(gamma*Tp)*(sign(r_e)*(abs(r_e))^(1-gamma)+sign(r_e)*(abs(r_e))^(1+gamma));


taocu=(dalphau_bar-F1+x(3)*tanh(u_e/Pi1)+PTSU+(a1/c1)*cos(phi_e)+K3*u_e-klu*x(1))/(hfu*sigamu);
taocr=(dalphar_bar-F2+x(4)*tanh(r_e/Pi2)+PTSR+a2/c2+K4*r_e-klr*x(2))/(hfr*sigamr);

 
tau_u_max =  500;
tau_u_min = 0;   

%% 输入饱和控制
tau_r_max = 50;
tau_r_min = -50;

taomu=(tau_u_max+tau_u_min)/2+((tau_u_max-tau_u_min)/2)*sign(taocu);
taomr=(tau_r_max+tau_r_min)/2+((tau_r_max-tau_r_min)/2)*sign(taocr);
taou=tanh(2*t)*taomu*tanh(1*taocu/(1*taomu));
taor=tanh(2*t)*taomr*tanh(1*taocr/(1*taomr));

taufu=hfu*taou+kfu;
taufr=hfr*taor+kfr;


taowu=-10+4*sin(0.5*t)*cos(0.5*t)-6*cos(t)*cos(0.5*t);
taowv=5*sin(0.1*t);
taowr=8*cos(0.3*t)*sin(1.1*t);
du=sigamu*taowu;
dr=sigamr*(taowr-hbar*taowv);

Ff1=f1+du+sigamu*kfu;
Ff2=f3+dr+sigamr*kfr;
    
    
sys = [taufu; taufr];
end
