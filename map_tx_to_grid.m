function tx_cells = map_tx_to_grid(x_tx, y_tx, z_tx, x, y, z)
    N = numel(x_tx);
    tx_cells = zeros(N,3);

    for n = 1:N
        [~, ix] = min(abs(x - x_tx(n)));
        [~, iy] = min(abs(y - y_tx(n)));
        [~, iz] = min(abs(z - z_tx(n)));
        tx_cells(n,:) = [ix, iy, iz];
    end

    % 去重（多个 microsphere 落在同一格只需要一个格点）
    tx_cells = unique(tx_cells, 'rows');
end
