function v = tensor2vector(vT, params)
    v = reshape(vT, params.Nx * params.Ny * params.Nb,1);
end
