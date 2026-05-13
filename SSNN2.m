function [sys,x0,str,ts] = SSNN2(t,x,u,flag)
switch flag
case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
% case 1
%     sys=mdlDerivatives(t,x,u);
case 3
    sys=mdlOutputs(t,x,u);
case {1,2,4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 3;
sizes.NumInputs      = 21;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [];
str = [];
ts=[-1 0];
% ts  = [0 0];

global c2 b2 w2 TK_PASS2 Ir2

c2=1*[ -4.2 -4 -3.8 -3.6 -3.4 -3.2 -3 -2.8 -2.6 -2.4 -2.2 -2 -1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 ;    
       -4.2 -4 -3.8 -3.6 -3.4 -3.2 -3 -2.8 -2.6 -2.4 -2.2 -2 -1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 ];

b2=5*ones(size(c2,1),size(c2,2));

w2=zeros(1,size(c2,2));      %权值初值
TK_PASS2=0;     %记录过去时刻
Ir2=ones(1,size(c2,2));   %分裂衰减速度系数

function sys=mdlOutputs(t,x,u)
%% 逼近fr
global c2 b2 w2 TK_PASS2 Ir2
deta_TK2=t-TK_PASS2;   %这个时刻减前一个时刻的时间
TK_PASS2=t;

r=u(19);
alphar=u(8);
re=u(6);%用于自适应调节权值=re

%% 神经输入
Vinput2=[r alphar]';   %NN输入
phi2=[];     %神经元激活
for j=1:size(c2,2)
    expsum2=0;
    for i=1:length(Vinput2)
        expsum2=expsum2+(Vinput2(i)-c2(i,j))^2/(b2(i,j)^2);   %数  神经元激活
    end
    phi2=[phi2 exp(-expsum2)];  %行   神经元激活
end
%% 模糊化

obs_f2 = w2 * phi2';  %权值*神经元=估计函数

%% 更新律
k_w2=35;

w_new2=-k_w2*(phi2*re+0.01*w2)* deta_TK2 + w2;   %更新权值

w2 = w_new2; %更新权值

%% 除去神经元
Eth2=0.5;
EPSW2=[];  %剔除标志列
SPSW2=[];  %保留标志列

for i=1:size(phi2,2)
if phi2(i)<=Eth2      %如果激活值小于设定值   激活效果太差  去除神经元
        Ir2(i)=Ir2(i)*0.5;     %给一个标志位衰减，不能立马除去神经元，不然可能一下除完了
else
    Ir2(i)=1;
end
end

for i=1:size(phi2,2)
    if Ir2(i)<=0.1      %衰减标志位小于某个值，可以清除这个神经元了
%     if phi(i)<=Eth
        EPSW2=[EPSW2 i];  %要剔除的神经元序号
    else
        SPSW2=[SPSW2 i];  %要保留的神经元序号
    end
end

if size(EPSW2,2)~=0    %检测到要剔除
    if size(SPSW2,2)==0   %如果没有保留神经元，最少留一个
        SPSW2=1;    %如果没有保留神经元，最少留一个
    end
    c2=c2(:,SPSW2);
    b2=b2(:,SPSW2);
    Ir2=Ir2(:,SPSW2);
    sum_Ew2=0;   %初始剔除总和值
    for i=1:size(EPSW2,2)
        sum_Ew2= sum_Ew2 + w2(:,EPSW2(i)) .* phi2(EPSW2(i));  %剔除总和值
    end
    w2(:,SPSW2(1)) = w2(:,SPSW2(1)) + sum_Ew2 ./ phi2(SPSW2(1));  %剔除总和值分布到第一个神经元
    w2=w2(:,SPSW2);
    phi2=phi2(:,SPSW2);
end
%%  分裂神经元
Gth2=0.9;
[phi_max2,phimax_idx2] = max(phi2);  %求最大NN激活值
if phi_max2<=Gth2   %如果小于某个值，最大激活值那个神经元分裂
    c2=[c2 (Vinput2+c2(:,phimax_idx2))/2];
%     b=[b ones(size(c,1),1)];
    b2=[b2 b2(:,phimax_idx2)];
    w2=[w2 zeros(1,1)];
    Ir2=[Ir2 1];
end


mw2=(w2*w2')^0.5;

m=size(c2,2);

sys(1)=obs_f2;
sys(2)=size(c2,2);
sys(3)=mw2;

















  
