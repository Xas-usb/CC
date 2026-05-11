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
    % --- RBF 神经网络参数 ---
    N_rbf1 = 15;  % 纵向通道神经元数
    N_rbf2 = 15;  % 转首通道神经元数

    sizes = simsizes;
    sizes.NumContStates  = N_rbf1 + N_rbf2;  % PTSNN 自适应权重作为连续状态
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2;   % tau_u, tau_r (未饱和)
    sizes.NumInputs      = 16;  % Filter(6) + Error约束(6) + phi_e + E + Aux(2)
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);

    % NN 权重初始化为 0 (论文: Ŵ1(0)=Ŵ2(0)=0)
    x0  = zeros(N_rbf1 + N_rbf2, 1);
    str = [];
    ts  = [0 0];  % 连续时间
end

%% ===== RBF 核函数 =====
function Z = rbf_kernel(x, centers, width)
    % x:   标量输入
    % centers: N x 1
    % width:  标量 (所有神经元同宽)
    % Z:     N x 1, 高斯基函数输出
    Z = exp(-((x - centers).^2) / (width^2 + 1e-30));
end

%% ===== NN 状态导数 (自适应律 式2-57/58) =====
function sys = mdlDerivatives(t,x,u_in)
    % 权重状态解包
    N1 = 15; N2 = 15;
    W1 = x(1:N1);        % 纵向权重 Ŵ1
    W2 = x(N1+1:N1+N2);  % 转首权重 Ŵ2

    % 输入解包
    alpha_bar_u = u_in(3);
    alpha_bar_r = u_in(4);
    u_e         = u_in(5);
    r_e         = u_in(6);

    % 反算实际速度
    u_real = alpha_bar_u - u_e;
    r_real = alpha_bar_r - r_e;

    % --- 纵向 NN (Surge) ---
    % RBF centers: u ∈ [-2, 6] m/s, 15个均匀分布
    c1_vec = linspace(-2, 6, N1)';
    b1 = 1.0;  % 宽度参数
    Z1 = rbf_kernel(u_real, c1_vec, b1);

    Gamma1 = 20;   % 学习率 (论文表2-2)
    sigma1 = 0.03; % σ-修正因子
    W1_dot = Gamma1 * (-Z1 * u_e - sigma1 * W1);

    % --- 转首 NN (Yaw) ---
    % RBF centers: r ∈ [-0.5, 0.5] rad/s, 15个均匀分布
    c2_vec = linspace(-0.5, 0.5, N2)';
    b2 = 0.1;  % 宽度参数
    Z2 = rbf_kernel(r_real, c2_vec, b2);

    Gamma2 = 20;   % 学习率 (论文表2-2)
    sigma2 = 0.01; % σ-修正因子
    W2_dot = Gamma2 * (-Z2 * r_e - sigma2 * W2);

    % 权重导数限幅 (防止瞬态爆炸)
    dot_max = 1e4;
    W1_dot = max(min(W1_dot, dot_max), -dot_max);
    W2_dot = max(min(W2_dot, dot_max), -dot_max);

    sys = [W1_dot; W2_dot];
end

%% ===== 控制律输出 (式2-51, 含 PTSNN 补偿) =====
function sys = mdlOutputs(t,x,u_in)
    % --- 1. 信号解包 ---
    alpha_bar_u_dot = u_in(1);
    alpha_bar_r_dot = u_in(2);
    alpha_bar_u     = u_in(3);
    alpha_bar_r     = u_in(4);
    u_e             = u_in(5);
    r_e             = u_in(6);

    sig1 = u_in(7);
    sig2 = u_in(8);
    th1  = u_in(9);
    th2  = u_in(10);

    E     = u_in(13); %#ok<NASGU>
    phi_e = u_in(14);

    Lambda_u = u_in(15);
    Lambda_r = u_in(16);

    % --- 2. 控制参数 ---
    Tp = 20; gamma = 0.35;
    K3 = 30; K4 = 30;
    k1 = 10; k2 = 10; l1 = 300; l2 = 300;
    ell_Lambda = 0.01;

    zeta_u = 1; zeta_r = 1;

    m11 = 25.8; m22 = 33.8; m23 = 1.0948; m32 = 1.0948; m33 = 2.76;
    xi_u = 1 / m11;
    xi_r = m22 / (m22*m33 - m23*m32);

    % --- 3. 安全计算 rho ---
    th1_safe = max(th1, 1e-6);
    th2_safe = max(th2, 1e-6);
    rho_u = max(min(sig1 / th1_safe, 0.999), -0.999);
    rho_r = max(min(sig2 / th2_safe, 0.999), -0.999);

    % --- 4. 重算 a, c ---
    sig1_safe = max(min(sig1, k1*0.999), -k1*0.999);
    sig2_safe = max(min(sig2, k2*0.999), -k2*0.999);

    E_l1_term   = l1 / (1 - (sig1_safe/k1)^2);
    phi_l2_term = l2 / (1 - (sig2_safe/k2)^2);

    c1 = (1 - rho_u^2)^2 * th1^2 * (E_l1_term)^1.5;
    c2 = (1 - rho_r^2)^2 * th2^2 * (phi_l2_term)^1.5;
    c1_safe = max(c1, 1e-6);
    c2_safe = max(c2, 1e-6);

    a1 = rho_u * k1 * l1 * th1;
    a2 = rho_r * k2 * l2 * th2;

    % --- 5. PTSNN 输出 F̂1, F̂2 (式2-53) ---
    N1 = 15; N2 = 15;
    W1 = x(1:N1);
    W2 = x(N1+1:N1+N2);

    % 反算实际速度
    u_real = alpha_bar_u - u_e;
    r_real = alpha_bar_r - r_e;

    % RBF 基函数计算
    c1_vec = linspace(-2, 6, N1)';
    b1 = 1.0;
    Z1 = rbf_kernel(u_real, c1_vec, b1);

    c2_vec = linspace(-0.5, 0.5, N2)';
    b2 = 0.1;
    Z2 = rbf_kernel(r_real, c2_vec, b2);

    % 集总不确定性估计
    F1_hat = W1' * Z1;
    F2_hat = W2' * Z2;

    % --- 6. 控制律计算 (式2-51, 含 -F̂1/-F̂2) ---
    term_ue = (pi/(gamma*Tp)) * (sign(u_e)*abs(u_e)^(1-gamma) + sign(u_e)*abs(u_e)^(1+gamma));
    term_re = (pi/(gamma*Tp)) * (sign(r_e)*abs(r_e)^(1-gamma) + sign(r_e)*abs(r_e)^(1+gamma));

    tau_u = (alpha_bar_u_dot + term_ue + (a1/c1_safe)*cos(phi_e) + K3*u_e - ell_Lambda*Lambda_u - F1_hat) / (xi_u * zeta_u);
    tau_r = (alpha_bar_r_dot + term_re + (a2/c2_safe) + K4*r_e - ell_Lambda*Lambda_r - F2_hat) / (xi_r * zeta_r);

    sys = [tau_u; tau_r];
end
