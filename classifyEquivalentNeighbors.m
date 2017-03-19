function latticeSiteClasses = classifyEquivalentNeighbors(params)
    % latticeSiteClasses is a class, containing classes of sites. Each site
    % in these classes can be connected by symmetry transformations of the
    % lattice with other sites in its class. The sites are given as vectors
    % [ix; iy], where ix, iy are displacements in {1,2,...,-1,0}
    
    latticeSites = cell(params.Nx*params.Ny,1);
    i = 1;
    for ix=1:params.Nx
        for iy=1:params.Ny
            latticeSites{i} = [ix;iy];
            i = i+1;
        end
    end
    
    latticeSiteClasses = cell(0);
    while ~isempty(latticeSites)
        v = latticeSites{1};  % Take a vector from the set of all lattice 
                              % site points.
        
        c = cell(0);    % Create a class (c) for that vector and all other 
                        % vectors that can be reached from that vector 
                        % using symmetry transformations.
        c{1} = v;
        latticeSites=removeCellElement(latticeSites,v);
        
        % Create all points that can be reached from class c using the
        % symmetry transformations of the lattice.
        % Perform a symmetry transformation on all elements in c. Do this
        % until you find no new points that can be moved from the cell
        % latticeSites to c.
        while true
            lc = length(c);
            for i1=1:length(c)
                % Perform all possible symmetry transformations on all
                % elements in c.
                newVectors = performSymmetryTransformations(c{i1}, params);
                [~,l] = size(newVectors);
                for i2=1:l
                    c = joinVectorToCell(c, newVectors(:,i2));
                    latticeSites = removeCellElement(latticeSites, newVectors(:,i2));
                end
            end
            if length(c) == lc
                % In case that there's no new sites that can be reached by
                % symmetry transforms, BREAK.
                break;
            end
        end
        
        % Add the newly formed class of lattice sites to
        % latticeSiteClasses.
        latticeSiteClasses{length(latticeSiteClasses) + 1} = c;
    end
    
    % Transform all vectors into displacements on the lattice
    % {1,2,...,params.Nx-1,params.Nx} -> {1,2,...,-1,0}
    for i1=1:length(latticeSiteClasses)
        for i2=1:length(latticeSiteClasses{i1})
            v = latticeSiteClasses{i1}{i2};
            v(1) = mod(v(1), params.Nx) - params.Nx*(params.Nx>v(1) & v(1)>params.Nx/2);
            v(2) = mod(v(2), params.Ny) - params.Ny*(params.Ny>v(2) & v(2)>params.Ny/2);
            latticeSiteClasses{i1}{i2} = v;
        end
    end
end

function mxOfNewVectors = performSymmetryTransformations(xy,params)
    % This function creates each symmetry transformation given in the class
    % params.symmetries on the vector xy = [dx, dy], and puts the
    % resulting vectors into the array mxOfNewVectors.
    [~,~,noSymms] = size(params.symmetriesLatt);
    
    mxOfNewVectors = zeros(2, noSymms);
    for i1 = 1:noSymms
        sigmaLatt = params.symmetriesLatt(:,:,i1);
        mxOfNewVectors(:,i1) = transformSite(xy, sigmaLatt, params);
    end
end

function xyNew = transformSite(xy, sigmaLatt, params)
    % Transforms a site xy = [dx, dy] according to the
    % transformation 'sigma'. sigmaLatt acts on the [dx, dy] lattice
    % vector.
    xyNew = cyclicVector(sigmaLatt * xy, params);
end

function vCyclic = cyclicVector(v, params)
    % If one of the components of vector v doesn't fit in the range
    % 1,2,...,params.Nx and 1,2,...params.Ny, respectively, then we reshift
    % it appropriately.
    vCyclic = [mod(round(v(1) - 1), params.Nx) + 1; ...
               mod(round(v(2) - 1), params.Ny) + 1];
end

function cnew=joinVectorToCell(c,v)
    % Cell c contains vectors. If v is one of these, then nothing happens,
    % otherwise it is added to c.
    ismember = 0;
    for i=1:length(c)
        if c{i} == v,
            ismember = 1;
            break;
        end
    end
    cnew = c;
    if ismember == 0
        cnew{length(c)+1} = v;
    end
end

function cnew=removeCellElement(c,v)
    ismember = 0;
    for i=1:length(c)
        if c{i} == v,
            ismember = 1;
            break;
        end
    end
    if ismember == 1
        cnew = {c{1:i-1}, c{(i+1):length(c)}};
    else
        cnew = c;
    end
end