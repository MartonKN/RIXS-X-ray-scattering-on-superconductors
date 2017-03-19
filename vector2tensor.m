function vT = vector2tensor(v, params)
    vT = reshape(v, params.Nx, params.Ny, params.Nb);
end
