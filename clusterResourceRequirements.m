function [idealNumberOfCores, noTasksPerCore, ...
          memoryRequirementPerCore, timeRequirementPerCore] = ...
          clusterResourceRequirements(estimatedNumberOfCores, params)
    
    [idealNumberOfCores, noTasksPerCore] = checkNumberOfCores(estimatedNumberOfCores, params);
    disp(['Ideal number of cores: ', num2str(idealNumberOfCores)]);
    disp(['Number of SRIXS runs per core: ', num2str(noTasksPerCore)]);
    
    % Estimate memory requirements
    memoryRequirementPerCore = estimatedMemoryRequirement(params, noTasksPerCore);
    disp(['Estimated memory: ', num2str(memoryRequirementPerCore), ' MB']);
    
    % Estimate memory requirements
    timeRequirementPerCore = estimatedTimeRequirement(params, noTasksPerCore);
    hrs = floor(timeRequirementPerCore / 3600);
    mins = floor((timeRequirementPerCore - 3600*hrs) / 60);
    secs = mod(timeRequirementPerCore, 60);
    disp(['Estimated time: ', num2str(hrs), 'h ', ...
                              num2str(mins), 'm ', ...
                              num2str(secs), 's ',...
                              '(', num2str(timeRequirementPerCore), ' s)']);
end

function mem = estimatedMemoryRequirement(params, noTasksOfCore)
    % Returns the estimated memory requirements in MB
    Nd = params.Nx*params.Ny*params.Nb*2;
    t = 0:params.dt:params.tFinal;
    
    % Number of complex numbers stored roughly
    mem = 2*Nd^2*length(t) + 6*Nd^2 + 11*Nd + noTasksOfCore*length(t)^2*params.Nb^2*2;
    % ~ Nd: eE0t, eE0s, eE1, n0, w0, vec0, vecmn, tmpMx, fwmnvecmn, feH1tvec0
    % ~ Nd^2: V0, V1, V10, V1mn0, f, wmn
    % ~ Nd^2 * length(t): eH1t, eH1mnt
    % ~ length(t)^2 * params.Nb^2 * 2: SRIXS
    
    % Each complex number is stored on 16 bytes
    mem = mem * 16;
    
    % We calculate in kBytes
    mem = mem / 1024^2;
    
    % Add 80%, plus 200 MB, just to be sure.
    mem = round(mem * 1.80 + 200);
end

function timeNeeded = estimatedTimeRequirement(params, noTasksOfCore)
    % Very crude estimate of the time requirements in seconds on a single 
    % core machine, based on the simulations run on my laptop.
    Nd = params.Nx*params.Ny*params.Nb*2;
    t = 0:params.dt:params.tFinal;
    
    timeNeeded = (30 + 2.5 * length(t)^2) * (0.3 + (Nd/500)^3) * noTasksOfCore;
    
    % Add 50% plus 5 minutes just to be sure
    timeNeeded = round(timeNeeded + 300, -1);
end


function [idealNumberOfCores,noTasksPerCore] = checkNumberOfCores(estimatedNumberOfCores, params)
    % Checks if the number of cores requested is far from ideal, i.e. too
    % many cores do nothing at all.
    s = 0:params.ds:params.sFinal;
    t = 0:params.dt:params.tFinal;
    neighborClasses = truncatedEquivalentNeighbors(params);
    
    noCompTasks = length(s)*length(t)*length(neighborClasses); 
    % Number of all computational tasks to be carried out
    
    idealNumberOfCores = ceil(noCompTasks / round(noCompTasks / estimatedNumberOfCores));
    noTasksPerCore = ceil(noCompTasks/idealNumberOfCores);
end