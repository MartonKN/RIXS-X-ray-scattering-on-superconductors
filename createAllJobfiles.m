noCores = 98;
timeReq = 2*24*3600; % in sec
memReq  = 5000; % in MB
jobNameStub = 'YBCO';

loadParametersFileNames = dir('loadParameters*.m');

for jobNo=1:length(loadParametersFileNames)
    jobName = [jobNameStub,sprintf(['%0',num2str(floor(log10(length(loadParametersFileNames)))+1),'d'],jobNo)];
    createJobFile(jobName, loadParametersFileNames(jobNo).name, noCores, timeReq, memReq);
end