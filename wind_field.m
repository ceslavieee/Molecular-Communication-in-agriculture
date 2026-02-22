function [vx, vy, vz] = wind_field(t, x, y, z, scenario)
    if nargin < 5
        scenario = 1;
    end

    [X, Y, Z] = ndgrid(x, y, z);
    Tday = 24*3600;
    phase = 0;   % 可以以后用来移相

    % ==== 1. 根据 scenario 选择参数 ====
    switch scenario
        case 1  % Baseline: open field
            v_bar       = 0.5;
            alpha       = 0.5;
            sig2_speed  = 0.03;
            sig2_dir    = 1.5;
            theta0      = 0;

        case 2  % Sheltered / calm
            v_bar       = 0.2;
            alpha       = 0.3;
            sig2_speed  = 0.01;
            sig2_dir    = 0.3;
            theta0      = pi/4;

        case 3  % Windy
            v_bar       = 1.0;
            alpha       = 0.7;
            sig2_speed  = 0.02;
            sig2_dir    = 0.5;
            theta0      = 0;

        case 4  % Gusty / unstable
            v_bar       = 0.5;
            alpha       = 0.5;
            sig2_speed  = 0.08;
            sig2_dir    = 3.0;
            theta0      = 0;

        otherwise
            error('Unknown wind scenario');
    end

    % ==== 2. 日变化因子 ====
    d_t = sin(2*pi*(t+phase)/Tday);

    % ==== 3. 水平风速模长 ====
    sigma_speed = sqrt(sig2_speed);
    eta_speed   = sigma_speed * randn(size(X));
    vs = v_bar * (1 + alpha*d_t) .* (1 + eta_speed);

    % ==== 4. 风向 ====
    sigma_dir = sqrt(sig2_dir);
    eta_dir   = sigma_dir * randn(size(X));
    theta = 2*pi*(t+phase)/Tday + theta0 + eta_dir;

    % ==== 5. 分量 ====
    vx = vs .* cos(theta);
    vy = vs .* sin(theta);

    z_ref = 5;
    vz = 0.1 * vs .* log(1 + Z/z_ref);
end
