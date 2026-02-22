function [x_tx, y_tx, z_tx] = tx_positions(mode, N)
    switch mode
        case 'central'     % 4m x 4m, centered at (0,0,0.5)
            x_tx = -2 + 4*rand(N,1);
            y_tx = -2 + 4*rand(N,1);
            z_tx = 0.5*ones(N,1);

        case 'uniform'     % 4m x 80m, centered at (0,0,0.5)
            x_tx = -40 + 80*rand(N,1);
            y_tx = -2  + 4*rand(N,1);
            z_tx = 0.5*ones(N,1);

        case 'corners'     % four 1m x 1m patches at (±39.5, ±1.5, 0.5)
            N4 = N/4;
            centers = [ 39.5,  1.5;
                        39.5, -1.5;
                       -39.5,  1.5;
                       -39.5, -1.5 ];
            x_tx = []; y_tx = []; z_tx = [];
            for c = 1:4
                cx = centers(c,1);
                cy = centers(c,2);
                x_tx = [x_tx; cx - 0.5 + rand(N4,1)*1];
                y_tx = [y_tx; cy - 0.5 + rand(N4,1)*1];
                z_tx = [z_tx; 0.5*ones(N4,1)];
            end

        case 'perimeter'   % 0.5m wide around the 4x80 stripe
            N4 = N/4;
            % top edge
            x1 = -40 + 80*rand(N4,1);
            y1 = 2 + 0.5*rand(N4,1);
            % bottom edge
            x2 = -40 + 80*rand(N4,1);
            y2 = -2 - 0.5*rand(N4,1);
            % left edge
            x3 = -40 - 0.5*rand(N4,1);
            y3 = -2 + 4*rand(N4,1);
            % right edge
            x4 = 40 + 0.5*rand(N4,1);
            y4 = -2 + 4*rand(N4,1);

            x_tx = [x1; x2; x3; x4];
            y_tx = [y1; y2; y3; y4];
            z_tx = 0.5*ones(N4*4,1);

        otherwise
            error('Unknown mode');
    end
end