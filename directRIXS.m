function [IRIXS,SRIXS,neighborClasses,params] = directRIXS(params,fnameSRIXS,fnameIRIXS)
     % Currently this function only applies to a single band.
     % Values returned: IRIXS(Domega,omega,q,sf) and SRIXS(t,tau,s,rmn,sf), where
     % rmn is the index of the cell array neighborClasses.
     % sf == 1 corresponds to the non-spin-flip intensity, whereas
     % sf == 2 to the spin-flip transition.
     
     % -----------------------------------------------------------
     % LOAD BASIC PARAMETERS 
     % -----------------------------------------------------------
     omega  = params.omega;
     Domega = params.Domega;
     q = params.q;
     
     
     % -----------------------------------------------------------
     % LOAD THE PARAMETERS THAT ARE NEEDED IN EACH ITERATION
     % -----------------------------------------------------------
     [t, s, eE0t, eE0s, V0, eE1, V1, n0] = ...
                    diagonalizeHamiltonians(params);
     
     % I precalculate V0' * expm(-1i * H1 * t) * V0, that will appear 
     % multiple times during the simulation
     V10 = V1'*V0;
     eH1t = zeros(length(eE0t(:,1)), length(eE0t(:,1)), length(t));
     for it=1:length(t)
         eH1t  (:,:,it) = V10' * diag(eE1(:,it)) * V10;
     end
     

     % -----------------------------------------------------------
     % CREATE ALL NON-EQUIVALENT LATTICE SITE PAIRS AND THE 
     % ASSOCIATED QUANTITIES 
     % -----------------------------------------------------------
     neighborClasses = truncatedEquivalentNeighbors(params);
     
     
     % -----------------------------------------------------------
     % DEFINE SRIXS and IRIXS
     % -----------------------------------------------------------
     SRIXS = zeros(length(t), length(t), length(s),length(neighborClasses),2);
     % SRIXS(t,tau,s,rmn,sf)
     % Meaning of the variable sf: ==1 for non-spin-flip transitions (NSF)
     %                             ==2 for spin-flip transitions (SF)
     
     [length_q,~] = size(q);
     IRIXS = zeros(length(Domega), length(omega), length_q, 2);
     % IRIXS(Domega,omega,q,sf)
     
     
     % -----------------------------------------------------------
     % DETERMINE SRIXS FOR ALL INEQUIVALENT PAIRS rmn AND FOR 
     % TIMES t, tau AND s
     % -----------------------------------------------------------
     for imn = 1:length(neighborClasses)
         disp([num2str(imn), '/', num2str(length(neighborClasses))]);
         % Choose an element rmn from each class to represent it
         rmn = neighborClasses{imn}{1}; rmn = [rmn(1),rmn(2)];
         
         % Unitary matrix for spatial shifts on the lattice r -> r + rmn.
         Vmn = spatialShiftTrf(rmn,params); 
         
         % I precalculate V0' * expm(-i Hmn t) * V0, since it's used multiple
         % times. Here Hmn = Vmn * H1 * Vmn', with Hmn being the Hamiltonian 
         % H0, with the addition of a core hole at site rmn.
         V1mn0 = V1' * Vmn' * V0;
         eH1mnt = zeros(length(eE0t(:,1)), length(eE0t(:,1)), length(t));
         for it=1:length(t)
             eH1mnt(:,:,it) = V1mn0' * diag(eE1(:,it)) * V1mn0;
         end
         
         % To calculate SRIXS, we'll need to take matrix elements with the 
         % following vectors
         vec0  = V0' * createLocalizedState([0,0], 1, params);
         vecmn = V0' * createLocalizedState( rmn , 1, params);
         
         % Iterate over all time arguments
         for it=1:length(t)
             for itau=1:length(t)
                 for is=1:length(s)
                     % Useful matrices
                     w0   = eE0t(:,it) .* conj(eE0t(:,itau)) .* eE0s(:,is);
                     wmn  = eH1t(:,:,itau) * diag(conj(eE0s(:,is))) * eH1mnt(:,:,it)';
                     f    = diag(1 - n0) + wmn * diag(w0.*n0);
                     detf = det(f);
                     
                     % Useful vectors
                     fwmnvecmn = f \ (wmn*vecmn);
                     feH1tvec0 = f \ (eH1t(:,:,itau) * vec0);
                     
                     % Expectation values appearing in SRIXS.
                     c1 = dot((1-n0) .* vec0, feH1tvec0);
                     c2 = dot((1-n0) .* (eH1t(:,:,itau) * (conj(eE0s(:,is)) .* vecmn)) , fwmnvecmn);
                     c3 = dot((1-n0) .* vec0, fwmnvecmn);
                     c4 = dot(conj(w0) .* n0 .* (eH1mnt(:,:,it)*vecmn) , feH1tvec0);
                     
                     SRIXS(it,itau,is,imn,1) = detf^2 * (4*c1*c2 + 2*c3*c4); % NSF
                     SRIXS(it,itau,is,imn,2) = detf^2 *            2*c3*c4 ; % SF
                 end
             end
         end
     end

     save(fnameSRIXS,'t','s','neighborClasses','SRIXS');

     % Tested and works up to this point: it correctly reproduces SRIXS
     % for any values of t, tau, s and rmn.
     
     % -----------------------------------------------------------
     % TRANSFORM EVERYTHING INTO FOURIER SPACE (FREQ AND MOMENTUM)
     % -----------------------------------------------------------
     % We only need to perform integration over positive values of s, since
     % the integral over the negative s contour just gives the complex
     % conjugate of the one for positive values of s.
     % 
     % Note that here I won't perform a DFT here, which would restrict
     % me to consider only discrete values of omega and Domega, but I take
     % them as continuous variables. This amounts to assuming that the
     % signal will have decayed over times t>tFinal and s>sFinal. In the
     % former case it is trivially justified since the finite lifetime of
     % the core hole guarantees me an exp(-params.Gamma*(t + tau)) factor,
     % that cuts of my integrand at large distances. In case of the
     % integral over s is a more tricky thing, that needs to be checked at
     % the end.
     
     % FIRST, PERFORM THE TIME INTEGRAL
     % In order to correctly perform Romberg integration, the first
     % elements SRIXS(t,tau,s,rmn) in all directions need to be multiplied by
     % 0.5.
     SRIXStmp = SRIXS;
     SRIXStmp(1,:,:,:,:) = SRIXStmp(1,:,:,:,:) * 0.5;
     SRIXStmp(:,1,:,:,:) = SRIXStmp(:,1,:,:,:) * 0.5;
     SRIXStmp(:,:,1,:,:) = SRIXStmp(:,:,1,:,:) * 0.5;
     
     % Multiply SRIXS by the 
     % exp(-params.Gamma*(t+tau)) * exp(-(apparatusResolution*s)^2/2)
     % factor
     eGammat   = repmat(exp(-params.Gamma * t).', 1, length(t), length(s), length(neighborClasses), 2);
     eGammatau = repmat(exp(-params.Gamma * t)  , length(t), 1, length(s), length(neighborClasses), 2);
     eResolutions = shiftdim( repmat(exp(-(params.apparatusResolution * s).^2/2).', 1, length(neighborClasses), 2, length(t), length(t)), 3);
     SRIXStmp = SRIXStmp.*eGammat.*eGammatau.*eResolutions;
     
     IRIXS_real_space = zeros(length(Domega),length(omega),length(neighborClasses), 2);
     for iomega = 1:length(omega)
         % Create matrices that will multiply SRIXS by exp(1i*omega*(t-tau))
         eomegat   = repmat(exp(1i * omega(iomega) * t).', 1, length(t), length(s), length(neighborClasses), 2);
         eomegatau = repmat(exp(-1i* omega(iomega) * t)  , length(t), 1, length(s), length(neighborClasses), 2);
         eomegattau = eomegat.*eomegatau;
         for iDomega = 1:length(Domega)
             % Create the matrix which multiplies SRIXS by exp(-1i*Domega*s)
             eDomegas = shiftdim( repmat(exp(-1i*Domega(iDomega) * s).', 1, length(neighborClasses), 2, length(t), length(t)), 3);
             % Multiply SRIXS by the appropriate phase factors and
             % integrate (or rather sum up) over all three time directions.
             IRIXS_real_space(iDomega, iomega, :, :) = ...
                 reshape( ... 
                            sum(sum(sum( eomegattau.*eDomegas.*SRIXStmp ))), ...
                            [1,1,length(neighborClasses), 2]...
                        )...
                 * params.dt*params.dt*params.ds;
         end
     end
     
     
     
     % SECONDLY, PERFORM THE REAL SPACE INTEGRATION
     % IRIXS(Domega,omega,q) = sum(IRIXS_real_space(Domega,omega,rmn) * exp(1i*q*rmn)).
     for iq=1:length_q
         for imn = 1:length(neighborClasses)
             % Sum up all the phase factors from this class of neighbors
             emn = 0;
             for i2=1:length(neighborClasses{imn})
                 rmn = neighborClasses{imn}{i2};
                 emn = emn + exp(1i*q(iq,:)*rmn(1:2));
             end
             
             % Add the real space value of IRIXS from this class of neighbors
             % to the momentum space RIXS amplitude, multiplied by the
             % appropriate phase factor.
             IRIXS(:,:,iq,:) = IRIXS(:,:,iq,:) + IRIXS_real_space(:,:,imn,:) * emn;
         end
     end
     IRIXS = IRIXS / (params.Nx*params.Ny);
     
     
     % FINALLY, TAKE ONLY THE REAL PART OF IRIXS,
     % The imaginary part is removed by the contribution of the integral
     % over the negative s axis.
     IRIXS = 2.0 * real(IRIXS);
                                         
     save(fnameIRIXS,'Domega','omega','q','IRIXS');
end




function [t, s, eE0t, eE0s, V0, eE1, V1, n0] = ...
         diagonalizeHamiltonians(params)
    
    % Auxiliary function to calculate the following
    %   expm(-1i * t(n) * H0) = V0 * diag(eE0t(:,n)) * V0';
    %   expm(-1i * s(n) * H0) = V0 * diag(eE0s(:,n)) * V0';
    %   expm(-1i * t(n) * H1) = V1 * diag(eE1 (:,n)) * V1';
    %   inv(1 + expm(H0 / params.T)) = V0 * diag(n0) * V0';
    
    t = 0:params.dt:params.tFinal;
    s = 0:params.ds:params.sFinal;
    
    % Create Hamiltonians
    H0 = createH0(params);
    H1 = H0 + createCoreHolePotential([0,0], params);
    
    % Diagonalize the Hamiltonians
    % V0'*H0*V0 == diag(E0); V1'*H1*V1 == diag(E1);
    [V0, E0] = eig(H0); E0 = real(diag(E0));
    [V1, E1] = eig(H1); E1 = real(diag(E1));
    
    edE0t = exp(-1i*E0*params.dt);
    edE0s = exp(-1i*E0*params.ds);
    edE1  = exp(-1i*E1*params.dt);
    
    eE0s = ones(length(E0),length(s));
    eE0t = ones(length(E0),length(t));
    eE1  = ones(length(E1),length(t));
    
    
    for it=1:(length(t)-1)
        eE0t(:,it+1) = eE0t(:,it).*edE0t;
        eE1 (:,it+1) = eE1 (:,it).*edE1;
    end
    for is=1:(length(s)-1)
        eE0s(:,is+1) = eE0s(:,is).*edE0s;
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
    
% Tested and works!
end