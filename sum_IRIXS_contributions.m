function [IRIXS, omega, Domega, q, params] = ...
         sum_IRIXS_contributions(IRIXS_input_filenames, IRIXS_output_filename, params)
    % After calculateIRIXS_cluster() has determined the IRIXS contributions
    % of several outputs of the cluster, we just need to sum them up.
    % IRIXS_filenames shall refer to many files at once, like '*IRIXS*.mat'
    %
    % Input:
    %   - sf == 1 corresponds to the non-spin-flip intensity, whereas
    %     sf == 2 to the spin-flip transition.
    %
    % Output: IRIXS(Domega,omega,q,sf)
    
    [~,length_q] = size(params.q);
    IRIXS_total = zeros(length(params.Domega), length(params.omega), length_q, 2);
    
    allFileNames = dir(IRIXS_input_filenames);
    for id = 1:length(allFileNames)
        fname = allFileNames(id).name;
        try
            load(fname);
        catch
            error(['Error: could not load ',fname]);
        end
        IRIXS_total = IRIXS_total + IRIXS;
    end
    
    IRIXS = IRIXS_total;
    clear IRIXS_total;
    
    save(IRIXS_output_filename, 'params', 'IRIXS', 'omega', 'Domega', 'q', '-v7.3');

end
