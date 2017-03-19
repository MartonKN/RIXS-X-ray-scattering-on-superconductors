function truncatedClasses = truncatedEquivalentNeighbors(params)
    % Create equivalent classes of points, connected by symmetry of the
    % lattice as given by the function classifyEquivalentStates(), and
    % put those sites into the same class at the end of the cell array
    % truncatedClasses, whose distance from the origin
    % (0,0) == (params.Nx,params.Ny) is greater than the radius 
    % params.neighbor_cutoff.
    % 
    % When performing the RIXS calculation, the algorithm will choose the
    % first element of each class and perform the calculation only on that
    % site. This shall give the same result as for all elements of the
    % class except for the last class in 'truncatedClasses' where all sites
    % that are further from the origin than params.neighbor_cutoff are
    % dumped. Therefore, by choosing the first element of this class, we
    % indeed make an assumption that the RIXS contribution from these sites
    % will be roughly the same. 
    %
    % Importantly, truncatedEquivalentNeighbors() also sorts the classes in
    % their increasing distance from the origin. Therefore, choosing the
    % first element of the last class in 'truncatedClasses' means that we
    % represent all sites of this class with the one that is the closest to
    % the origin.
    
    truncatedClasses = cell(0);
    remainingSites = cell(0);
    latticeSiteClasses = classifyEquivalentNeighbors(params);
    latticeSiteClasses = sortSites(latticeSiteClasses);
    
    for i=1:length(latticeSiteClasses)
        v = latticeSiteClasses{i}{1}; 
        % Choose the first element of this class, they have the same length
        % anyway.
        if norm(v(1:2),2) <= params.neighbor_cutoff
            truncatedClasses{length(truncatedClasses)+1} = latticeSiteClasses{i};
        else
            for i1 = 1:length(latticeSiteClasses{i})
                remainingSites{length(remainingSites)+1} = latticeSiteClasses{i}{i1};
            end
        end
    end
    if ~isempty(remainingSites)
        truncatedClasses{length(truncatedClasses) + 1} = remainingSites;
    end
end

function orderedLatticeSiteClasses = sortSites(latticeSiteClasses)
    % Sorts the neighboring sites according to which has the smallest
    % distance from the origin
    orderedLatticeSiteClasses = {};
    tmpLatticeSiteClasses = latticeSiteClasses;
    
    while ~isempty(tmpLatticeSiteClasses)
        posSmallestElement = 0;
        valSmallestElement = 10^10;
        for i1 = 1:length(tmpLatticeSiteClasses)
            tmp = norm(tmpLatticeSiteClasses{i1}{1}(1:2), 2);
            if tmp < valSmallestElement
                posSmallestElement = i1;
                valSmallestElement = tmp;
            end
        end
        orderedLatticeSiteClasses{length(orderedLatticeSiteClasses) + 1} = tmpLatticeSiteClasses{posSmallestElement};
        tmpLatticeSiteClasses = { tmpLatticeSiteClasses{1 : (posSmallestElement-1)}, ...
                                  tmpLatticeSiteClasses{(posSmallestElement+1) : length(tmpLatticeSiteClasses)} };
    end
end