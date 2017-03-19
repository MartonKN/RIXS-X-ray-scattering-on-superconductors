function createJobFile(filename, paramsFilename, noCores, timeRequested, memoryRequested)
    jobfilename = [filename,'_jobfile.slurm'];
    fid = fopen(jobfilename, 'w');
    
    timeMins = ceil(timeRequested/60); % seconds -> minutes
    days = floor(timeMins / (60*24));
    hrs  = floor((timeMins - days*60*24)/60);
    mins = timeMins - days*60*24 - hrs*60;
    
    fprintf(fid, '#!/bin/bash \n\n');
    
    fprintf(fid, ['#SBATCH --job-name=', filename, '        \t\t\t\t\t\t # name of the job\n']);
    fprintf(fid, '#SBATCH --array=1-%d                                \t # number of jobs needed\n',noCores);
    fprintf(fid, '#SBATCH -n 1                                        \t # number of cores per job\n');
    fprintf(fid, '#SBATCH -N 1                                        \t # the n jobs need to be on N nodes\n');
    fprintf(fid, '#SBATCH -t %d-%d:%d                               \t\t # time requirement (D-HH:MM)\n',days,hrs,mins);
    fprintf(fid, '#SBATCH --mem-per-cpu=%d                            \t # required memory (MB)\n', memoryRequested);
    fprintf(fid, '#SBATCH -p serial_requeue                           \t # where to submit jobs (serial_requeue is the fastest)\n');
    fprintf(fid, ['#SBATCH --error=',filename,'_%%a.err                 \t # stderr: %%a = ID of job in the array, %%A = job array ID\n']);
    fprintf(fid, ['#SBATCH --output=',filename,'_%%a.out                \t # stdout\n']);
    fprintf(fid, '# #SBATCH --mail-type=END                             \t # email notifications (BEGIN, END, FAIL, ALL)\n');
    fprintf(fid, '# #SBATCH --mail-user=kanasz.nagy.marton@gmail.com    \t # your email\n\n');

    fprintf(fid, 'hostname                                            \t # print hostname into the output file\n\n');
    
    fprintf(fid, 'module load matlab/R2014a-fasrc01                   \t # load matlab\n');
    fprintf(fid, ['matlab -nosplash -nodesktop -r "distribute_SRIXS_tasks(''', filename,''', ', num2str(noCores),', ${SLURM_ARRAY_TASK_ID}, ', paramsFilename,')"\n']);
    fprintf(fid, '                                                    \t # start program\n\n\n');
    
    fprintf(fid, '# Run this job \n');
    fprintf(fid, ['# sbatch ', jobfilename, ' \n']);
    
    fclose(fid);
end