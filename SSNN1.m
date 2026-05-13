function [sys,x0,str,ts] = SSNN1(t,x,u,flag)
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

global c1 b1 w1 TK_PASS1  Ir1 


c1=1*[ -4.2 -4 -3.8 -3.6 -3.4 -3.2 -3 -2.8 -2.6 -2.4 -2.2 -2 -1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 ;    
       -4.2 -4 -3.8 -3.6 -3.4 -3.2 -3 -2.8 -2.6 -2.4 -2.2 -2 -1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 ];

b1=5*ones(size(c1,1),size(c1,2));

w1=rand(1,size(c1,2));      %权值初值
TK_PASS1=0;     %记录过去时刻
Ir1=ones(1,size(c1,2));   %分裂衰减速度系数


function sys=mdlOutputs(t,x,u)
global c1 b1 w1 TK_PASS1 Ir1 

deta_TK1=t-TK_PASS1;   %这个时刻减前一个时刻的时间
TK_PASS1=t;

u1=u(17);
alphau=u(7);
ue=u(5);%用于自适应调节权值=ue

%% 神经输入
Vinput1=[u1' alphau']';   %NN输入

phi1=[];     %神经元激活
for j=1:size(c1,2)
    expsum1=0;
    for i=1:length(Vinput1)
        expsum1=expsum1+(Vinput1(i)-c1(i,j))^2/(1.0*b1(i,j)^2);   %数  神经元激活
    end
    phi1=[phi1 exp(-expsum1)];  %行 神经元激活
end
%% 模糊化

obs_f1 = w1 * phi1';  %权值*神经元=估计函数

%% 更新律
k_w1=20;%60  

w_new1=(-k_w1*(ue*phi1+0.03*w1)) * deta_TK1 + w1; 

w1 = w_new1; %更新权值
%% 除去神经元
Eth1=0.02; %0.2
EPSW1=[];  %剔除标志列
SPSW1=[];  %保留标志列

for i=1:size(phi1,2)
if phi1(i)<=Eth1      %如果激活值小于设定值   激活效果太差  去除神经元
        Ir1(i)=Ir1(i)*0.5;     %给一个标志位衰减，不能立马除去神经元，不然可能一下除完了
else
    Ir1(i)=1;
end
end

for i=1:size(phi1,2)
    if Ir1(i)<=0.1      %衰减标志位小于某个值，可以清除这个神经元了
%     if phi(i)<=Eth
        EPSW1=[EPSW1 i];  %要剔除的神经元序号
    else
        SPSW1=[SPSW1 i];  %要保留的神经元序号
    end
end

if size(EPSW1,2)~=0    %检测到要剔除
    if size(SPSW1,2)==0   %如果没有保留神经元，最少留一个
        SPSW1=1;    %如果没有保留神经元，最少留一个
    end
    c1=c1(:,SPSW1);
    b1=b1(:,SPSW1);
    Ir1=Ir1(:,SPSW1);
    sum_Ew1=0;   %初始剔除总和值
    for i=1:size(EPSW1,2)
        sum_Ew1= sum_Ew1 + w1(:,EPSW1(i)) .* phi1(EPSW1(i));  %剔除总和值
    end
    w1(:,SPSW1(1)) = w1(:,SPSW1(1)) + sum_Ew1 ./ phi1(SPSW1(1));  %剔除总和值分布到第一个神经元
    w1=w1(:,SPSW1);
    phi1=phi1(:,SPSW1);
end
%%  分裂神经元
Gth1=0.4;%0.9
[phi_max1,phimax_idx1] = max(phi1);  %求最大NN激活值
if phi_max1<Gth1   %如果小于某个值，最大激活值那个神经元分裂 （神经元激活值不够大，增加神经元）
    c1=[c1 (Vinput1+c1(:,phimax_idx1))/2];
    b1=[b1 b1(:,phimax_idx1)];
    w1=[w1 zeros(1,1)];
    Ir1=[Ir1 1];
end
%%
mw1=(w1*w1')^0.5;
 
sys(1)=obs_f1;
sys(2)=size(c1,2);
sys(3)=mw1;







  
