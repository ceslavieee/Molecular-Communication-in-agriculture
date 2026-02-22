function [q] = release_model(dt, Nt)
    % KP 参数（从论文图1拟合）
    k     = 0.4429;
    n_exp = 0.1789;
    M_inf = 200;          % g，总释放质量

    t_hours = (0:Nt)*dt/3600;       % 0~24h
    F = k * t_hours.^n_exp;
    F(F>1) = 1;
    M = M_inf * F;

    dM = diff(M);          % 每步释放的质量 (g)
    q  = dM / dt;          % g/s，长度 Nt
end
