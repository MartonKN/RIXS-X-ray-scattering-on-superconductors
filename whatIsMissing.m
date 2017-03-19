function missingIDs = whatIsMissing(filename, noCores, firstID, lastID)
    allFiles = dir();
    IDs = firstID:1:lastID;
    % This array contains all job IDs that should have been run
    % If we do find such a file in the directory, then we set the
    % corresponding ID to zero.
    
    strIDs = cell(1,length(IDs)); % ID strings
    for i1 = 1:length(IDs)
        strIDs{i1} = sprintf(['%0',num2str(floor(log10(noCores))+1),'d'], IDs(i1));
    end
    
    for id = 1:length(allFiles)
        fname = allFiles(id).name;
        if numel(strfind(fname, filename))>0
            for i2 = 1:length(IDs)
                if numel(strfind(fname, strIDs{i2}))>0
                    IDs(i2) = 0;
                end
            end
        end
    end
    
    
    missingIDs = IDs(find(IDs>0));
end