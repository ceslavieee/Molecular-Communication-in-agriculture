function C_new = step_pde(C, vx, vy, vz, S, D, dx, dy, dz, dt)
    [Nx,Ny,Nz] = size(C);
    C_new = C;

    % ---- 扩散项 D∇²C (中心差分) ----
    Lap = zeros(size(C));
    Lap(2:end-1, 2:end-1, 2:end-1) = ...
        (C(3:end,   2:end-1, 2:end-1) - 2*C(2:end-1,2:end-1,2:end-1) + C(1:end-2,2:end-1,2:end-1))/dx^2 + ...
        (C(2:end-1, 3:end,   2:end-1) - 2*C(2:end-1,2:end-1,2:end-1) + C(2:end-1,1:end-2,2:end-1))/dy^2 + ...
        (C(2:end-1, 2:end-1, 3:end  ) - 2*C(2:end-1,2:end-1,2:end-1) + C(2:end-1,2:end-1,1:end-2))/dz^2;
    Lap = D * Lap;

    % ---- 对流项 -v·∇C (一阶迎风) ----
    ADV = zeros(size(C));
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                if vx(i,j,k) >= 0
                    dCdx = (C(i,j,k)-C(i-1,j,k))/dx;
                else
                    dCdx = (C(i+1,j,k)-C(i,j,k))/dx;
                end
                if vy(i,j,k) >= 0
                    dCdy = (C(i,j,k)-C(i,j-1,k))/dy;
                else
                    dCdy = (C(i,j+1,k)-C(i,j,k))/dy;
                end
                if vz(i,j,k) >= 0
                    dCdz = (C(i,j,k)-C(i,j,k-1))/dz;
                else
                    dCdz = (C(i,j,k+1)-C(i,j,k))/dz;
                end

                ADV(i,j,k) = -(vx(i,j,k)*dCdx + vy(i,j,k)*dCdy + vz(i,j,k)*dCdz);
            end
        end
    end

    % ---- 显式 Euler 更新 ----
    C_new = C + dt * (Lap + ADV + S);

    % ---- 边界条件 ----
    % ground z=0: absorbing
    C_new(:,:,1) = 0;

    % lateral/top: Neumann / outflow (zero-gradient)
    C_new(1,:,:)   = C_new(2,:,:);
    C_new(end,:,:) = C_new(end-1,:,:);
    C_new(:,1,:)   = C_new(:,2,:);
    C_new(:,end,:) = C_new(:,end-1,:);
    C_new(:,:,end) = C_new(:,:,end-1);
end
