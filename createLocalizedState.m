function v = createLocalizedState(r, band, params)
    % Auxiliary function to create a single particle eigenstate, in the 
    % band 'band', and centered at site r.
    rn = [mod(round(r(1)-1), params.Nx)+1, mod(round(r(2)-1), params.Ny)+1];
    
    v = zeros(params.Nx,params.Ny,params.Nb);
    v(rn(1),rn(2),band) = 1;
    v = tensor2vector(v,params);
    
% Tested and works!
end