function [CEI_vs_Cth, CEI_vs_t] = run_CEI_for_mode(tx_mode)

    % ---------- 1. 基本参数、网格、q(t) ----------
    dx = 5; dy = 5; dz = 0.5;
    x = -100:dx:100;  Nx = numel(x);
    y = -100:dy:100;  Ny = numel(y);
    z = 0:dz:5;       Nz = numel(z);

    dt      = 2;
    T_hours = 24;
    Nt      = T_hours*3600/dt;

    D = 1e-5;

    q = release_model(dt, Nt);   % 之前写的那个：总质量释放速率 [g/s]

    % ---------- 2. 阈值设置 ----------
    Cth_mol  = [1e4 1e6 1e8 1e10 1e12 1e14];
    MA = 152.149; NA = 6.022e23;
    conv_mol2g = MA/NA;
    Cth_mass = Cth_mol * conv_mol2g;   % [g/m^3]
    K = numel(Cth_mass);

    % ---------- 3. TX 布局 & 映射到网格 ----------
    N_tx_model = 100;  % 代表 TX 个数
    [x_tx, y_tx, z_tx] = tx_positions(tx_mode, N_tx_model);
    tx_cells = map_tx_to_grid(x_tx, y_tx, z_tx, x, y, z);

    % ---------- 4. CEI 累积变量 ----------
    C = zeros(Nx,Ny,Nz);          % 初始浓度
    V_cells = Nx*Ny*Nz;           % 网格总数（比体积差一个常数，对比例来说等价）

    sum_f = zeros(K,1);           % 24h 用：∑ f_n(Cth_k)
    sum_f_fixed = 0;              % 下图用：∑ f_n(Cth_fixed)
    CEI_vs_t = zeros(Nt,1);       % 存每个时间步的 CEI(t, Cth=1e14)

    idx_fixed = K;                % 这里最后一个就是 1e14

    % ---------- 5. 时间循环 ----------
    for n = 1:Nt
        t = n * dt;

        % 风场
        [vx,vy,vz] = wind_field(t, x,y,z);

        % 源项
        S = build_source(q(n), tx_cells, dx,dy,dz, Nx,Ny,Nz);

        % PDE 更新
        C = step_pde(C, vx,vy,vz, S, D, dx,dy,dz, dt);

        % ---- CEI 统计：f_n(Cth_k) ----
        for k = 1:K
            mask = (C >= Cth_mass(k));
            f_nk = sum(mask(:)) / V_cells;   % 区域比例
            sum_f(k) = sum_f(k) + f_nk;
            if k == idx_fixed
                sum_f_fixed = sum_f_fixed + f_nk;
                CEI_vs_t(n) = sum_f_fixed / n;   % 离散形式: (1/n)∑ f_k
            end
        end
    end

    % ---------- 6. 24h 的 CEI(Cth) ----------
    CEI_vs_Cth = sum_f / Nt;   % CEI(t_sim, Cth_k)

end