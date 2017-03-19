function [XAS, omega, params] = ...
         sum_XAS_contributions(XAS_input_filenames, XAS_output_filename, params)
    % After calculateXAS_cluster() has determined the XAS contributions
    % of several outputs of the cluster, we just need to sum them up.
    % XAS_filenames shall refer to many files at once, like '*XAS*.mat'
    %
    % Input:
    %   - sf == 1 corresponds to the non-spin-flip intensity, whereas
    %     sf == 2 to the spin-flip transition.
    %
    % Output: XAS(omega,sf)
    
    XAS_total = zeros(length(params.omega), 2);
    
    allFileNames = dir(XAS_input_filenames);
    for id = 1:length(allFileNames)
        fname = allFileNames(id).name;
        try
            load(fname);
        catch
            error(['Error: could not load ',fname]);
        end
        XAS_total = XAS_total + XAS;
    end
    
    XAS = XAS_total;
    clear IRIXS_total;
    
    save(XAS_output_filename, 'params', 'XAS', 'omega', '-v7.3');

end
