function SRIXS = calculateSRIXS_cluster(s,t,tau,rmn,params)
    % Values returned: SRIXS(t,tau,s,rmn,nu1,nu2,nu3,nu4,sf)
    % sf == 1 corresponds to the non-spin-flip intensity, whereas
    % sf == 2 to the spin-flip transition.
    %
    % Input parameters: 
    % - rmn: the positions of the core hole, e.g. [1,0; 2,1; 3,-1; 0,4]
    % - s: s values 
    % - t: t values
    % - tau: tau values
    % - params: physical parameters, as returned by the function
    %   loadParameters()
    % 
    % Band indices: nu == 1
    
    % Since there are huge arrays of unitary matrices that is the same for
    % t and tau, I do not want to calculate it twice. Thus, I just merge
    % the t and tau values into the array tAll, and use that. tAll will
    % then index these matrix arrays (eE0t, eH1t, eH1mnt) etc.
    % Then, the function tIndex will help helps me find the index of a 
    % given time value at question.
    tAll = union(t,tau); 
    
    [eE0t, eE0s, V0, eE1, V1, n0] = diagonalizeHamiltonians(s, tAll, params);
    
    % I precalculate V0' * expm(-1i * H1 * t) * V0, that will appear 
    % multiple times during the simulation
    V10 = V1'*V0;
    eH1t = zeros(length(eE0t(:,1)), length(eE0t(:,1)), length(tAll));
    for itAll=1:length(tAll)
        eH1t(:,:,itAll) = V10' * diag(eE1(:,itAll)) * V10;
    end
    
    
    % -----------------------------------------------------------
    % DEFINE SRIXS and IRIXS
    % -----------------------------------------------------------
    [lrmn, ~] = size(rmn);
    SRIXS = zeros(length(t), length(tau), length(s), lrmn,... 
                  params.Nb, params.Nb, params.Nb, params.Nb ,2);
    % SRIXS(t, tau, s, rmn, nu1, nu2, nu3, nu4, sf)
    % sf ==1 for non-spin-flip transitions (NSF)
    %    ==2 for spin-flip transitions (SF)
    % nu1, ..., nu4: band indices
    
    % -----------------------------------------------------------
    % DETERMINE SRIXS FOR ALL CORE HOLE POSITIONS rmn AND FOR 
    % TIMES t, tau AND s
    % -----------------------------------------------------------
    for imn = 1:lrmn
        % Unitary matrix for spatial shifts on the lattice r -> r + rmn.
        Vmn = spatialShiftTrf(rmn(imn,:),params); 
         
        % I precalculate V0' * expm(-i Hmn t) * V0, since it's used multiple
        % times. Here Hmn = Vmn * H1 * Vmn', with Hmn being the Hamiltonian 
        % H0, with the addition of a core hole at site rmn.
        V1mn0 = V1' * Vmn' * V0;
        eH1mnt = zeros(length(eE0t(:,1)), length(eE0t(:,1)), length(tAll));
        for itAll=1:length(tAll)
            eH1mnt(:,:,itAll) = V1mn0' * diag(eE1(:,itAll)) * V1mn0;
        end
        clear Vmn V1mn0;
        
        % To calculate SRIXS, we'll need to take matrix elements with the 
        % following vectors
        % vec0  = V0'*createLocalizedState([0,0], nu, params);
        % vecmn = V0'*createLocalizedState( rmn , nu, params);
        vec0  = zeros(params.Nx*params.Ny*params.Nb, params.Nb);
        vecmn = zeros(params.Nx*params.Ny*params.Nb, params.Nb);
        for nu=1:params.Nb
            vec0 (:,nu) = createLocalizedState( [0,0]  , nu, params);
            vecmn(:,nu) = createLocalizedState(rmn(imn,:), nu, params);
        end
        vec0  = V0' * vec0;
        vecmn = V0' * vecmn;
        
        
        % Iterate over all time arguments
        for it=1:length(t)
            for itau=1:length(tau)
                for is=1:length(s)
                    itAll   = tIndex(t(it),    tAll);
                    itauAll = tIndex(tau(itau),tAll);
                    % Since the matrix arrays are indexed by tAll, but my
                    % time variables are indexed by t and tau, I need to
                    % switch between these two languages.
                    
                    % Useful matrices
                    w0   = eE0t(:,itAll) .* conj(eE0t(:,itauAll)) .* eE0s(:,is);
                    wmn  = eH1t(:,:,itauAll) * diag(conj(eE0s(:,is))) * eH1mnt(:,:,itAll)';
                    f    = diag(1 - n0) + wmn * diag(w0.*n0);
                    detf = det(f);
                    
                    % Useful vectors
                    % fwmnvecmn = inv(f) * (wmn*vecmn);
                    % feH1tvec0 = inv(f) * (eH1t(:,:,itau) * vec0);
                    fwmnvecmn = wmn*vecmn;
                    feH1tvec0 = eH1t(:,:,itauAll) * vec0;
                    tmpMx = f\[fwmnvecmn,feH1tvec0];
                    fwmnvecmn = tmpMx(:,1:params.Nb);
                    feH1tvec0 = tmpMx(:,(params.Nb+1):(2*params.Nb));
                    
                    % Expectation values appearing in SRIXS.
                    c1 = zeros(params.Nb, params.Nb);
                    c2 = zeros(params.Nb, params.Nb);
                    c3 = zeros(params.Nb, params.Nb);
                    c4 = zeros(params.Nb, params.Nb);
                    
%                    if params.Nb > 1 
                        for nu1=1:params.Nb
                            for nu2=1:params.Nb
                                c1(nu1,nu2) = dot((1-n0) .* vec0(:,nu1), feH1tvec0(:,nu2));
                                c2(nu1,nu2) = dot((1-n0) .* (eH1t(:,:,itauAll) * (conj(eE0s(:,is)) .* vecmn(:,nu1))) , fwmnvecmn(:,nu2));
                                c3(nu1,nu2) = dot((1-n0) .* vec0(:,nu1), fwmnvecmn(:,nu2));
                                c4(nu1,nu2) = dot(conj(w0) .* n0 .* (eH1mnt(:,:,itAll)*vecmn(:,nu1)) , feH1tvec0(:,nu2));
                            end
                        end
            
                        % Determining SRIXS
                        for nu1=1:params.Nb
                            for nu2=1:params.Nb
                                for nu3=1:params.Nb
                                    for nu4=1:params.Nb
                                        % NSF case
                                        SRIXS(it,itau,is,imn,nu1,nu2,nu3,nu4,1) = ...
                                            detf^2 * (4*c1(nu1,nu2)*c2(nu3,nu4) + 2*c3(nu1,nu4)*c4(nu3,nu2));
                                        % SF case
                                        SRIXS(it,itau,is,imn,nu1,nu2,nu3,nu4,2) = ...
                                            detf^2 *  2*c3(nu1,nu4)*c4(nu3,nu2) ;
                                    end
                                end
                            end
                        end
%                     else
%                         c1 = dot((1-n0) .* vec0(:,1), feH1tvec0(:,1));
%                         c2 = dot((1-n0) .* (eH1t(:,:,itauAll) * (conj(eE0s(:,is)) .* vecmn(:,1))) , fwmnvecmn(:,1));
%                         c3 = dot((1-n0) .* vec0(:,1), fwmnvecmn(:,1));
%                         c4 = dot(conj(w0) .* n0 .* (eH1mnt(:,:,itAll)*vecmn(:,1)) , feH1tvec0(:,1));   
%                         % NSF case
%                         SRIXS(it,itau,is,imn,1,1,1,1,1) = ...
%                             detf^2 * (4*c1*c2 + 2*c3*c4);
%                         % SF case
%                         SRIXS(it,itau,is,imn,1,1,1,1,2) = ...
%                             detf^2 *  2*c3*c4;
%                     end
                end
            end
        end
    end
end


function it = tIndex(tValue,tAll)
    it = find(tAll == tValue);
end


function [eE0t, eE0s, V0, eE1, V1, n0] = ...
         diagonalizeHamiltonians(s, t, params)
    
    % Auxiliary function to calculate the following
    %   expm(-1i * t(n) * H0) = V0 * diag(eE0t(:,n)) * V0';
    %   expm(-1i * s(n) * H0) = V0 * diag(eE0s(:,n)) * V0';
    %   expm(-1i * t(n) * H1) = V1 * diag(eE1 (:,n)) * V1';
    %   inv(1 + expm(H0 / params.T)) = V0 * diag(n0) * V0';
    
    % Create Hamiltonians
    H0 = createH0(params);
    H1 = H0 + createCoreHolePotential([0,0], params);
    
    % Diagonalize the Hamiltonians
    % V0'*H0*V0 == diag(E0); V1'*H1*V1 == diag(E1);
    [V0, E0] = eig(H0); E0 = real(diag(E0));
    [V1, E1] = eig(H1); E1 = real(diag(E1));
    
    eE0s = zeros(length(E0),length(s));
    eE0t = zeros(length(E0),length(t));
    eE1  = zeros(length(E1),length(t)); 
    
    for is = 1:length(s)
        eE0s(:,is) = exp(-1i*E0*s(is));
    end
    for it = 1:length(t)
        eE0t(:,it) = exp(-1i*E0*t(it));
        eE1 (:,it) = exp(-1i*E1*t(it));
    end
    
    n0 = 1./(1 + exp(E0/params.T));
    
    % Tested and works!
end


function v = createLocalizedState(r, band, params)
    % Auxiliary function to create a single particle eigenstate, in the 
    % band 'band', and centered at site r.
    rn = [mod(round(r(1)-1), params.Nx)+1, mod(round(r(2)-1), params.Ny)+1];
    
    v = zeros(params.Nx,params.Ny,params.Nb);
    v(rn(1),rn(2),band) = 1;
    v = tensor2vector(v,params);
end