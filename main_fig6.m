%% ===========================
%  main_fig6_scenarios.m
%  在 4 种风场 scenario 下画 Fig.6（2x3），
%  每个 scenario 一张图
%% ===========================

clear; clc;

%% ---- Step 0：初始化网格、时间、参数 ----
dx = 5; dy = 5; dz = 0.5;
x = -100:dx:100;     Nx = numel(x);
y = -100:dy:100;     Ny = numel(y);
z = 0:dz:5;          Nz = numel(z);

dt      = 2;                   % s
T_hours = 24;
Nt      = T_hours*3600/dt;     % 43200 步

D = 1e-5;                      % m^2/s

% 需要输出的三个时刻
times_h   = [1, 11, 22];
idx_times = round(times_h*3600/dt);

% z 平均窗口 z∈[1.5, 2.5]
z0 = 2.0; dz_band = 0.5;
idx_z = find(z >= z0-dz_band & z <= z0+dz_band);

% 风场场景名字（自己在 wind_field.m 里定义的 1~4 场景）
scenario_list  = 1:4;
scenario_names = { ...
    'Scenario 1: Baseline open field', ...
    'Scenario 2: Sheltered / calm', ...
    'Scenario 3: Windy day', ...
    'Scenario 4: Gusty / unstable'};

%% ---- Step 1：加载 MeSA 释放模型 q(t) ----
q = release_model(dt, Nt);     % [g/s]

%% ---- Step 2：生成 TX 布局（central patch 和 four corners）----
N_tx_model = 100;   % 用 100 个代表点即可

[x_tx_c, y_tx_c, z_tx_c] = tx_positions('central', N_tx_model);
[x_tx_f, y_tx_f, z_tx_f] = tx_positions('corners', N_tx_model);

% 映射到网格点（注意 tx_cells 和风无关，所以外循环 scenario 不用重复算）
tx_cells_central = map_tx_to_grid(x_tx_c, y_tx_c, z_tx_c, x, y, z);
tx_cells_corners = map_tx_to_grid(x_tx_f, y_tx_f, z_tx_f, x, y, z);

%% ========= 外层循环：不同的风场 scenario =========
for s = 1:numel(scenario_list)
    scenario = scenario_list(s);
    fprintf('Running Fig.6 for wind scenario %d ...\n', scenario);

    % ------------------------
    % 为当前 scenario 准备存储
    % ------------------------
    Cxy_central = cell(1, numel(times_h));
    Cxy_corners = cell(1, numel(times_h));

    %% ---- Step 3：central patch 24h 仿真 ----
    C = zeros(Nx,Ny,Nz);

    for n = 1:Nt
        t = n * dt;

        % 关键修改：传入 scenario
        [vx,vy,vz] = wind_field(t, x, y, z, scenario);
        S          = build_source(q(n), tx_cells_central, dx,dy,dz, Nx,Ny,Nz);

        C = step_pde(C, vx,vy,vz, S, D, dx,dy,dz, dt);

        if any(n == idx_times)
            idx = find(n == idx_times);
            Cxy_central{idx} = mean(C(:,:,idx_z), 3);
        end
    end

    %% ---- Step 4：four corners 24h 仿真 ----
    C = zeros(Nx,Ny,Nz);

    for n = 1:Nt
        t = n * dt;

        % 同样传入 scenario
        [vx,vy,vz] = wind_field(t, x, y, z, scenario);
        S          = build_source(q(n), tx_cells_corners, dx,dy,dz, Nx,Ny,Nz);

        C = step_pde(C, vx,vy,vz, S, D, dx,dy,dz, dt);

        if any(n == idx_times)
            idx = find(n == idx_times);
            Cxy_corners{idx} = mean(C(:,:,idx_z), 3);
        end
    end

    %% ---- Step 5：画当前 scenario 的 Fig.6（2x3）----
    figure;
    for k = 1:3
        % 上面一行：central patch
        subplot(2,3,k);
        C_xy = Cxy_central{k};
        imagesc(x, y, C_xy'); set(gca,'YDir','normal');
        xlabel('x in m'); ylabel('y in m');
        title(sprintf('Central patch, t = %d h', times_h(k)));
        colorbar;

        % 下面一行：four corners
        subplot(2,3,k+3);
        C_xy = Cxy_corners{k};
        imagesc(x, y, C_xy'); set(gca,'YDir','normal');
        xlabel('x in m'); ylabel('y in m');
        title(sprintf('Four corners, t = %d h', times_h(k)));
        colorbar;
    end

    % 给整张图加一个总标题，说明是哪种风场
    sgtitle(scenario_names{s});
end