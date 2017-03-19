function findCorruptedFiles(filenames)
    allFileNames = dir(filenames);
    for id = 1:length(allFileNames)
        fname = allFileNames(id).name;
        try
            load(fname);
        catch
            disp(['Error: could not load ',fname]);
            delete(fname);
        end
    end
end

