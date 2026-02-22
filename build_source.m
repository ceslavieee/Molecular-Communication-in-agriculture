function S = build_source(q_total, tx_cells, dx, dy, dz, Nx, Ny, Nz)
    Vcell = dx * dy * dz;
    N_tx  = size(tx_cells,1);

    q_per_tx = q_total / N_tx;      % 每个代表点的质量速率

    S = zeros(Nx, Ny, Nz);
    for l = 1:N_tx
        ix = tx_cells(l,1);
        iy = tx_cells(l,2);
        iz = tx_cells(l,3);
        S(ix,iy,iz) = S(ix,iy,iz) + q_per_tx / Vcell;
    end
end
