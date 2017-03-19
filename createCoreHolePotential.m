function U = createCoreHolePotential(r,params)
    % r = position of the core hole r = (ix, iy)
    U = zeros(params.Nx, params.Ny, params.Nb, params.Nx, params.Ny, params.Nb);
    
    % transforming r to a lattice variable quantity
    r = round(r);
    r = [mod(r(1)-1,params.Nx)+1, mod(r(2)-1,params.Ny)+1];
    
    for ib = 1:params.Nb
        U(r(1), r(2), ib, r(1), r(2), ib) = params.U0;
    end
    
    U = tensor2mx(U,params);
end

% Tested and works!