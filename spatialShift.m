function rshifted = spatialShift(r,d,params)
    % shift the vector r cyclically with the vector d on the lattice
    % 1,2,3,...,Nx and 1,2,3,...,Ny.
    rshifted = [mod(round(r(:,1)+d(:,1)-1),params.Nx)+1, mod(round(r(:,2)+d(:,2)-1),params.Ny)+1];
end

% Tested and works!