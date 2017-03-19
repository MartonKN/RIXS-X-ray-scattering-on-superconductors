function mT = mx2tensor(m, params)
    mT = reshape(m, params.Nx, params.Ny, params.Nb, params.Nx, params.Ny, params.Nb);
end
