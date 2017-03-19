function distribute_SRIXS_tasks(filenameStub, numberOfCores, coreID, params)
    % SRIXS_cluster_results{i1,1} = t
    % SRIXS_cluster_results{i1,2} = s
    % SRIXS_cluster_results{i1,3} = rmn
    % SRIXS_cluster_results{i1,4} = SRIXS

    coreIDStr = sprintf(['%0',num2str(floor(log10(numberOfCores))+1),'d'],coreID);
    filenameSRIXS = [filenameStub,'_',coreIDStr,'.mat'];
    filenameIRIXS = [filenameStub,'_IRIXS_',coreIDStr,'.mat'];
    filenameXAS   = [filenameStub,'_XAS_',coreIDStr,'.mat'];

    % If the files already exist, we do nothing. Otherwise we run the
    % simulation.
    allFiles = dir();
    SRIXS_file_exists = 0;
    IRIXS_file_exists = 0;
    XAS_file_exists = 0;
    
    for id = 1:length(allFiles)
        fname = allFiles(id).name;
        if numel(strfind(fname, filenameSRIXS))>0
            SRIXS_file_exists = 1;
            try
               	[~] = load(fname);
            catch
                SRIXS_file_exists = 0;
            end
        end
        if numel(strfind(fname, filenameIRIXS))>0
            IRIXS_file_exists = 1;
            try
               	[~] = load(fname);
            catch
                IRIXS_file_exists = 0;
            end
        end
        if numel(strfind(fname, filenameXAS))>0
            XAS_file_exists = 1;
            try
               	[~] = load(fname);
            catch
                XAS_file_exists = 0;
            end
        end        
    end
    
    % Run SRIXS part
    if ~SRIXS_file_exists
        tau = 0:params.dt:params.tFinal;

        % tasksOfCore{taskID, 1} = tRange
        % tasksOfCore{taskID, 2} = s
        % tasksOfCore{taskID, 1} = rmn
        tasksOfCore = getTasksOfCore(params, numberOfCores, coreID);
        [noTasksOfCore,~] = size(tasksOfCore);

        SRIXS_cluster_results = cell(noTasksOfCore,4);
        for i1 = 1:noTasksOfCore
            t = tasksOfCore{i1,1};
            s = tasksOfCore{i1,2};
            rmn = tasksOfCore{i1,3};
            SRIXS = calculateSRIXS_cluster(s,t,tau,rmn,params);
            SRIXS_cluster_results{i1,1} = t;
            SRIXS_cluster_results{i1,2} = s;
            SRIXS_cluster_results{i1,3} = rmn;
            SRIXS_cluster_results{i1,4} = SRIXS;
            disp('.');
        end
        save(filenameSRIXS, 'SRIXS_cluster_results');
        
        % Making sure the files were saved appropriately
        for i1=1:100
            succesfulSRIXSSave = 1;
            try
               	load filenameSRIXS;
            catch
                succesfulSRIXSSave = 0;
            end
            if succesfulSRIXSSave
                break;
            end
        end
    end
        
    % Determine IRIXS contribution of the current SRIXS parameters
    if ~IRIXS_file_exists || (~SRIXS_file_exists && succesfulSRIXSSave)
        calculateIRIXS_cluster(filenameSRIXS,filenameIRIXS,params);
    end
    
    % Determine XAS contribution of the current SRIXS parameters
    if ~XAS_file_exists   || (~SRIXS_file_exists && succesfulSRIXSSave)
        calculateXAS_cluster(filenameSRIXS,filenameXAS,params);
    end
end

function tasksOfCore = getTasksOfCore(params, numberOfCores, coreID)
    % tasksOfCore{taskID, 1} = tRange
    % tasksOfCore{taskID, 2} = s
    % tasksOfCore{taskID, 1} = rmn
    
    sValues = 0:params.ds:params.sFinal;
    tValues = 0:params.dt:params.tFinal;
    neighborClasses = truncatedEquivalentNeighbors(params);
    
    lmn = length(neighborClasses);
    ls = length(sValues);
    lt = length(tValues);
    
    if numberOfCores > (lt * lmn * ls)
        error(['Too many jobs (max number is ', num2str(lt * lmn * ls),')']);
    end
    
    % Let's determine how many pieces we are going to cut our tValues to
    % (We need to keep as big chunks of t values together as possible
    % for low computational load)
    % Then, we put the different chunks of tValues into the cell tRanges.
    n = ceil(numberOfCores/(lmn*ls));
    m = floor(lt / n);
    tRange = cell(n,1);
    for i1=1:(n-1)
        tRange{i1} = tValues(((i1-1)*m + 1):(i1*m));
    end
    tRange{n} = tValues(((n-1)*m + 1):lt);
    
    % Number of all computational units and their number per core
    nCU = lmn*ls*n;
    compUnitsPerNode = ceil(nCU / numberOfCores);

    
    % Create an array of the (t, s, rmn) pairs.
    compUnits = zeros(nCU,3);
    iCU = 1;
    for imn=1:length(neighborClasses)
        for is=1:length(sValues)
            for it=1:n
                compUnits(iCU,:) = [it,is,imn];
                iCU = iCU+1;
            end
        end
    end
    
    taskIDs = find(mod((1:nCU) - 1, numberOfCores) + 1 == coreID);
    tasksOfCore = cell(length(taskIDs),3);
    itask = 1;
    for iCU=taskIDs
        tasksOfCore{itask,1} = tRange{compUnits(iCU,1)};
        tasksOfCore{itask,2} = sValues(compUnits(iCU,2));
        tasksOfCore{itask,3} = neighborClasses{compUnits(iCU,3)}{1}.';
        itask = itask+1;
    end
end
