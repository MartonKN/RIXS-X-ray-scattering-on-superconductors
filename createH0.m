function H0 = createH0(params)
    H0 = zeros(params.Nx, params.Ny, params.Nb, params.Nx, params.Ny, params.Nb);
    
    % We choose each possible initial sites and scroll over the possible
    % hopping amplitudes from this site.
    [length_J,~] = size(params.J); 
    for ix1=1:params.Nx
        for iy1=1:params.Ny
            for i=1:length_J
                inew = spatialShift([ix1,iy1],params.J(i,1:2),params); 
                % the site [ix,iy] shifted by the hopping vector stored in 
                % the first two arguments of J(i,:).
                % J(i,:) = [Dx, Dy, band1, band2, Jvalue].
                H0(ix1,iy1,params.J(i,3), inew(1),inew(2),params.J(i,4)) = ...
                H0(ix1,iy1,params.J(i,3), inew(1),inew(2),params.J(i,4)) - ...
                params.J(i,5);
            end
        end
    end
    
    % Transforming H0 into matrix form
    H0 = tensor2mx(H0,params);
    
    % Determine the chemical potential
    E0 = sort(eig(H0));
    % (First give an estimate for the chemical potential)
    nOcc = 0.5*(1.0-params.doping) * params.Nx*params.Ny*params.Nb;
    mu0 = E0(floor(nOcc));
    % (Then determine the chemical potential accurately)
    fun = @(mu1) (filling(mu1,E0,params) - (1-params.doping));
    mu0 = fzero(fun,mu0);
    
    H0 = H0 - mu0 * eye(size(H0));
    
% Tested and works for a single band!
end


function n = filling(mu,E0,params)
    n = 2 * sum( 1./(1 + exp((E0-mu)/params.T)) ) / (params.Nx*params.Ny*params.Nb);
end