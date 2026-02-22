%% main_fig7_scenarios.m
% 在 4 种 wind scenario 下分别画 Fig.7（CEI vs Cth & CEI vs time）
% Wenhao 扩展版（支持风场创新）

clear; clc;

%% ----------------- 0. 全局参数 -----------------
dx = 5;  dy = 5;  dz = 0.5;
x  = -100:dx:100;   Nx = numel(x);
y  = -100:dy:100;   Ny = numel(y);
z  = 0:dz:5;        Nz = numel(z);

dt      = 2;        
T_hours = 24;
Nt      = T_hours*3600/dt;

D = 1e-5;

% ===== MeSA 释放模型 =====
q_total = release_model(dt, Nt);

% ===== CEI 阈值 =====
Cth_mol = [1e4 1e6 1e8 1e10 1e12 1e14];
MA = 152.149;  NA = 6.022e23;
conv_mol2g = MA / NA;
Cth_mass = Cth_mol * conv_mol2g;
K = numel(Cth_mass);
idx_fixed = K;

% ===== TX 布局 =====
tx_modes = {'central','uniform','corners','perimeter'};
num_modes = numel(tx_modes);

% ===== Wind Scenarios =====
scenario_list = 1:4;
scenario_names = {
    'Scenario 1: Baseline wind'
    'Scenario 2: Calm / sheltered'
    'Scenario 3: Strong wind'
    'Scenario 4: Gusty / noisy wind'
};


%% ============================================================
%    外层大循环：对每个 scenario 画一张 Fig.7
%% ============================================================

for s = 1:numel(scenario_list)

    scenario = scenario_list(s);
    fprintf("\n===== Running Fig.7 for Wind Scenario %d =====\n", scenario);

    % CEI 结果存储
    CEI_Cth_all = zeros(num_modes, K);
    CEI_t_all   = zeros(num_modes, Nt);

    %% -------------------------------------------------------
    %    遍历四种 TX 布局
    %% -------------------------------------------------------
    for m = 1:num_modes

        mode = tx_modes{m};
        fprintf("  TX mode = %s\n", mode);

        % -- TX 坐标 --
        N_tx_model = 100;
        [x_tx, y_tx, z_tx] = tx_positions(mode, N_tx_model);
        tx_cells = map_tx_to_grid(x_tx, y_tx, z_tx, x, y, z);

        % -- 初始化 --
        C = zeros(Nx,Ny,Nz);
        V_cells = Nx*Ny*Nz;
        sum_f     = zeros(K,1);
        sum_f_fix = 0;
        CEI_vs_t  = zeros(Nt,1);

        %% -------------------------------------
        %         时间推进主循环
        %% -------------------------------------
        for n = 1:Nt

            t = n * dt;

            % ---------------- Wind ----------------
            [vx,vy,vz] = wind_field(t, x, y, z, scenario);

            % ---------------- Source ---------------
            S = build_source(q_total(n), tx_cells, dx,dy,dz, Nx,Ny,Nz);

            % ---------------- PDE step -------------
            C = step_pde(C, vx,vy,vz, S, D, dx,dy,dz, dt);

            % ---------------- CEI 计算 ------------
            for k = 1:K
                mask = (C >= Cth_mass(k));
                f_nk = nnz(mask) / V_cells;
                sum_f(k) = sum_f(k) + f_nk;

                if k == idx_fixed
                    sum_f_fix   = sum_f_fix + f_nk;
                    CEI_vs_t(n) = sum_f_fix / n;
                end
            end
        end

        % ---- 仿真结束：算出 CEI ----
        CEI_Cth_all(m,:) = sum_f / Nt;
        CEI_t_all(m,:)   = CEI_vs_t;
    end


    %% ============================================================
    %                      绘图：Fig.7
    %% ============================================================

    figure;

    % 时间轴
    t_hours = (1:Nt)*dt/3600;

    % --- 上图：CEI vs Cth ---
    subplot(2,1,1);
    semilogx(Cth_mol, CEI_Cth_all(1,:), '-ob','LineWidth',1.2); hold on;
    semilogx(Cth_mol, CEI_Cth_all(2,:), '-or','LineWidth',1.2);
    semilogx(Cth_mol, CEI_Cth_all(3,:), '-og','LineWidth',1.2);
    semilogx(Cth_mol, CEI_Cth_all(4,:), '-x','LineWidth',1.2,'Color',[0.9 0.5 0]);

    grid on;
    xlabel('Threshold C_{th} (molecules/m^3)');
    ylabel('CEI(24h, C_{th})');
    legend(tx_modes,'Location','southwest');

    % --- 下图：CEI vs time ---
    subplot(2,1,2);
    plot(t_hours, CEI_t_all(1,:), 'b','LineWidth',1.2); hold on;
    plot(t_hours, CEI_t_all(2,:), 'r','LineWidth',1.2);
    plot(t_hours, CEI_t_all(3,:), 'g','LineWidth',1.2);
    plot(t_hours, CEI_t_all(4,:), 'Color',[0.9 0.5 0],'LineWidth',1.2);

    grid on;
    xlabel('Time (h)');
    ylabel('CEI(t, C_{th}=10^{14})');
    legend(tx_modes,'Location','southeast');

    sgtitle(sprintf("Fig.7 — %s", scenario_names{s}));

end