function V = spatialShiftTrf(d,params)
    % This unitary transform just performs a spatial shift in real space.
    % To get a core hole potential at site r you just do
    %   U0 = createCoreHolePotential([0,0], params);
    %   Vr = spatialShiftTrf(d,params);
    %   Ur = Vr*U0*Vr';
    % This gives the same result as Ur=createCoreHolePotential(r, params);
    x0=(1:params.Nx)'; 
    y0=(1:params.Ny)'; 
    x1 = spatialShift1D(x0,d(1),params.Nx);
    y1 = spatialShift1D(y0,d(2),params.Ny);
    
    V=zeros(params.Nx,params.Ny,params.Nb,params.Nx,params.Ny,params.Nb);
    for ib=1:params.Nb
        for ix=1:params.Nx
            for iy=1:params.Ny
                V(x1(ix),y1(iy),ib,ix,iy,ib) = 1;
            end
        end
    end
    
    V=tensor2mx(V,params);
    V=sparse(V);
end

function xshifted = spatialShift1D(x,dx,Nx)
    % shift x cyclically with dx on the lattice 1,2,3,...,Nx.
    xshifted = mod(round(x+dx-1),Nx)+1;
end


% Tested and works!