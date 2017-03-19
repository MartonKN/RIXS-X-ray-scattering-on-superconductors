function m = tensor2mx(mT, params)
    m = reshape(mT, params.Nx * params.Ny * params.Nb, params.Nx * params.Ny * params.Nb);
end

% Tested and works!