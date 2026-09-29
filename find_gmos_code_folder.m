function codeDir = find_gmos_code_folder(thisPipelineDir)
%FIND_GMOS_CODE_FOLDER Try to locate the user's GMOS Code folder.

    candidates = { ...
        fullfile(thisPipelineDir, 'Code'), ...
        fullfile(thisPipelineDir, '..', 'Code'), ...
        fullfile(thisPipelineDir, 'Code_GMOS', 'Code'), ...
        fullfile(thisPipelineDir, '..', 'Code_GMOS', 'Code'), ...
        fullfile(thisPipelineDir, '..', '..', 'Code'), ...
        thisPipelineDir};

    codeDir = '';
    for k = 1:numel(candidates)
        c = candidates{k};
        if isfolder(c) && (isfile(fullfile(c, 'pick_po_from_JPL.m')) || ...
                isfile(fullfile(c, 'L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv')))
            codeDir = c;
            return;
        end
    end

    codeDir = thisPipelineDir;
    warning('Could not automatically find a GMOS Code folder. Add your Code folder to the MATLAB path manually if needed.');
end
