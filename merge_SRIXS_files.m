function [SRIXS,s,t,neighborClasses,params] = merge_SRIXS_files(filename, numberOfCores, params)
    s = 0:params.ds:params.sFinal;
    t = 0:params.dt:params.tFinal;
    neighborClasses = truncatedEquivalentNeighbors(params);
    
    SRIXS = zeros(length(t), length(t), length(s), length(neighborClasses),... 
                  params.Nb, params.Nb, params.Nb, params.Nb, 2);
    
    checkTable = zeros(length(t), length(s), length(neighborClasses));
    % We put a 1 into each entry if the SRIXS calculation at the 
    % corresponding s and rmn values was carried out. If at the end the
    % matrix is not all ones, there was a problem.
    
    for coreID = 1:numberOfCores
        try 
            coreIDStr = sprintf(['%0',num2str(floor(log10(numberOfCores))+1),'d'],coreID);
            load([filename,'_',coreIDStr,'.mat']);
            % These files contain the cell array 'SRIXS_cluster_results', 
            % containing the results of the SRIXS simulation on the given 
            % node in the form {s, rmn, SRIXS}.
        catch
            disp(['Error: no result from core #',num2str(coreID),'!!']);
        end
        [noTasksOfCore,~] = size(SRIXS_cluster_results);
        for i1 = 1:noTasksOfCore
            tValues= SRIXS_cluster_results{i1,1};
            sValue = SRIXS_cluster_results{i1,2};
            rmn    = SRIXS_cluster_results{i1,3};
            is  = find_s(sValue, s);
            imn = find_rmn(rmn, neighborClasses);
            itMin = find_s(tValues(1), t);
            itMax = find_s(tValues(length(tValues)), t);
            SRIXS(itMin:itMax,:,is,imn,:,:,:,:,:) = SRIXS_cluster_results{i1,4}(1:length(tValues),:,1,1,:,:,:,:,:);
            checkTable(itMin:itMax,is,imn) = ones(length(tValues),1,1);
        end
    end
    
    % Check if all entries of the matrix checkTable is 1
    % If not, find which entries were not found
    if ~isequal(checkTable, ones(length(t), length(s), length(neighborClasses)))
        missingElements = [];
        for it=1:length(t)
            for is=1:length(s)
                imn = find(checkTable(it,is,:)==0);
                if ~isempty(imn)
                    missingElements = [missingElements;[it,is,imn]];
                end
            end
            disp('Error: the following (t, s, rmn) pairs are missing:');
            disp(missingElements);
        end
    else
        % If everything was fine, we delete the cluster simulation files,
        % and put SRIXS and the other important variables into
        % filename.mat.
        save tmpFileSRIXS.mat SRIXS s t neighborClasses params -v7.3;
        movefile('tmpFileSRIXS.mat', [filename,'.mat'])
        for coreID = 1:numberOfCores
            coreIDStr = sprintf(['%0',num2str(floor(log10(numberOfCores))+1),'d'],coreID);
            delete([filename,'_',coreIDStr,'.mat']);
        end
    end
end

function index = find_s(sVaule, s)
    % Find the index of sValue in the vector s.
    index = find(sVaule == s);
end

function index = find_rmn(rmn, neighborClasses)
    % Find which class of neighborClasses is rmn from
    index = 0;
    for i1=1:length(neighborClasses)
        for i2=1:length(neighborClasses{i1})
            if isequal(rmn, neighborClasses{i1}{i2}.')
                index = i1;
                break;
            end
        end
        if index
            break;
        end
    end
end