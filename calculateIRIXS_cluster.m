function [IRIXS, omega, Domega, q, params] = ... 
         calculateIRIXS_cluster(filenameSRIXS,filenameIRIXS,params)
    % This function reads the SRIXS results from the file(s) filenameSRIXS
    % and determines their contribution to the total IRIXS signal.
    % NOTE: flienameSRIXS can cover many files at once, e.g. 
    % filenameSRIXS = '*.mat'.
    %
    % Input: 
    %    - filenameSRIXS: name of the file containing SRIXS, t, s, 
    %      neighborClasses and params (with the .mat file extension)
    %    - filenameIRIXS: file in which to save the IRIXS intensities
    %       (with the .mat file extension)
    % 
    % Output: IRIXS(Domega,omega,q,polarizations)
    
    omega  = params.omega;
    Domega = params.Domega;
    q  = params.q;
    lq = size(q,2);
    neighborClasses = truncatedEquivalentNeighbors(params);    
    IRIXS = zeros(length(Domega),length(omega),lq,2);
    
    % Loop through each filename and get the results of the given run
    SRIXSfiles = dir(filenameSRIXS);
    for id = 1:length(SRIXSfiles)
        fname = SRIXSfiles(id).name;
        try
            load(fname);
        catch
            error(['Error: could not load ',fname]);
        end
        [noTasksOfCore,~] = size(SRIXS_cluster_results);
        for iTask = 1:noTasksOfCore
            t   = SRIXS_cluster_results{iTask,1};
            tau = 0:params.dt:params.tFinal;
            s   = SRIXS_cluster_results{iTask,2};
            rmn = SRIXS_cluster_results{iTask,3};
            SRIXS = SRIXS_cluster_results{iTask,4};
            
            IRIXS = IRIXS + calculateIRIXS_contribution(t,tau,s,rmn,neighborClasses,SRIXS,params);
            disp('.');
        end
    end
    
    % Save
    save(filenameIRIXS, 'params', 'IRIXS', 'omega', 'Domega', 'q', '-v7.3');
end

function IRIXS = calculateIRIXS_contribution(t,tau,s,rmn,neighborClasses,SRIXS,params)
    % The RIXS intensity IRIXS is a many-component sum. This function
    % determines the contributions of the SRIXS amplitudes that had been
    % determined in the range 
    % SRIXS(t,tau,s,rmn,nu1,nu2,nu3,nu4,sf)
    % by distribute_SRIXS_tasks.m, with
    %   t   == tValues(itMin:itMax)
    %   tau == tauValues(itauMin:itauMax)
    %   s   == sValues(isMin:isMax)
    %   rmn == nieghborClasses{imn}{1}.' == [ix1,iy1; ix2,iy2; ...] row vectors
    
    omega = params.omega;
    Domega = params.Domega;
    q = params.q;
    length_q = size(q,2);    
    IRIXS = zeros(length(Domega), length(omega), length(q), 2);
    
    % 1/2 factors of romberg integration
    if abs(t(1)) <10^(-10)
        SRIXS(1,:,:,:,:,:,:,:,:) = 0.5 * SRIXS(1,:,:,:,:,:,:,:,:);
    end
    if abs(tau(1)) <10^(-10)
        SRIXS(:,1,:,:,:,:,:,:,:) = 0.5 * SRIXS(:,1,:,:,:,:,:,:,:);
    end
    if abs(s(1)) <10^(-10)
        SRIXS(:,:,1,:,:,:,:,:,:) = 0.5 * SRIXS(:,:,1,:,:,:,:,:,:);
    end
    
    % Perform frequency summation first
    sSRIXS = size(SRIXS);
    lrmn = size(rmn,1);
    SRIXSfreq = zeros(length(Domega),length(omega),lrmn,params.Nb,params.Nb,params.Nb,params.Nb,2);
    for iomega = 1:length(omega)
        % Performing the time integrals for t, tau and s, with a single
        % (omega, Domega) pair of frequencies
        eiot   = exp(( 1i*omega(iomega) - params.Gamma)*t);
        eiotau = exp((-1i*omega(iomega) - params.Gamma)*tau);
        
        eiot   = reshape(eiot,[length(t),1,1]);
        eiotau = reshape(eiotau,[1,length(tau),1]);
        
        eiotMx   = repmat(eiot  , [1, sSRIXS(2:length(sSRIXS))]);
        eiotauMx = repmat(eiotau, [sSRIXS(1), 1, sSRIXS(3:length(sSRIXS))]);
        for iDomega = 1:length(Domega)
            eios   = exp( -1i*Domega(iDomega)*s) .* exp(-(params.apparatusResolution*s).^2 / 2);
            eios   = reshape(eios,[1,1,length(s)]);
            eiosMx   = repmat(eios  , [sSRIXS(1:2), 1, sSRIXS(4:length(sSRIXS))]);
            
            SRIXSfreqTmp = sum(sum(sum(SRIXS .* eiotMx .* eiotauMx .* eiosMx, 3), 2), 1) * params.dt^2 * params.ds;
            SRIXSfreq(iDomega,iomega,:,:,:,:,:,:) = reshape(SRIXSfreqTmp, [1,1,sSRIXS(4:length(sSRIXS))]);
            % SRIXSfreq(iomega,Domega,rmn,nu1,nu2,nu3,nu4,sf) has been summed up
            % for t, tau and s time indices. Now, it only remains to sum up
            % the positions and the band indices
        end
    end
    clear eiotMx eiotauMx eiosMx;
    
    % Perform momentum space summation
    for iq = 1:length_q
        for i1 = 1:lrmn
            imn = find_rmnIndices(rmn(i1,:),neighborClasses,params);
            imn = imn(1);
            emn = 0;
            for i2 = 1:length(neighborClasses{imn})
                r = neighborClasses{imn}{i2};
                emn = emn + exp(1i*q(:,iq).'*r);
            end
            IRIXS(:,:,iq,:) = IRIXS(:,:,iq,:) + emn * reshape(SRIXSfreq(:,:,i1,1,1,1,1,:), size(IRIXS(:,:,iq,:)));
        end
    end
    
    % We did not perform summation over the negative values of s, since
    % those just give the complex conjugate of the contributions for 
    % positive s. 
    IRIXS = 2.0 * real(IRIXS) / (params.Nx*params.Ny);
    
    
    % FINALLY, ASSUMING THAT THE RESOLUTION OF OUR t TIME PARAMETER IS GOOD
    % ENOUGH, WE CAN TRUNCATE THE HIGHER ORDER HARMONICS
    % If you assume that the t resolution is so good, that a linear
    % approximation between these points gives back the time-dependence of
    % SRIXS well, then the actual Fourier transform of the time continous
    % SRIXS is given by
    % IRIXS * 2*(1-cos(omega*params.dt)./(omega*params.dt)).^2.
    %
    % omega_tmp = omega + (omega==0)*10^(-4);
    % regulator_vec = 2*(1-cos(omega_tmp*params.dt))./(omega_tmp*params.dt).^2;
    % regulator_mx = repmat(regulator_vec, [length(Domega), 1, length_q, 2]);
    % IRIXS = regulator_mx .* IRIXS;
    %
    % I will comment it out for now, however, to be sure about whether my
    % time resolution is good enough.
end

function imn = find_rmnIndices(rmn,neighborClasses,params)
    [lmn,lmn2] = size(rmn);
    if lmn2 ~=2
        rmn = rmn.';
        lmn = lmn2;
    end
    lnC = length(neighborClasses);
    imn = zeros(lmn,2);
    for i0=1:lmn
        found = 0;
        for i1=1:lnC
            for i2=1:length(neighborClasses{i1})
                v = rmn(i0,:) - neighborClasses{i1}{i2}.';
                if mod(v(1),params.Nx) == 0 && mod(v(2),params.Ny) == 0
                    imn(i0,:) = [i1,i2];
                    found = 1;
                    break;
                end
                if(found)
                    break;
                end
            end
        end
    end
end
