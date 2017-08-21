% Author: Halle R. Dimsdale-Zucker
initialize_ABCDCon

analysis_name = 'univariate_sanityCheck';
contrast_lbl = 'allRHitsxFHits_Miss_scon0005_N23';

fpath = fullfile(analMRIDir,analysis_name, contrast_lbl,'coordinates_p001_k88_no_hdr.csv');

if exist(fpath,'file')
    coord = csvread(fpath);
    [oneline, cellarray]=cuixuFindStructure(coord);
    
    % save out the labels
    cell2csv_modeoptions(fullfile(dropbox_dir,'writeups', 'figures','univariate_allRHitsxFHits_Misses_p001_k88_xjview_lbls.csv'), cellarray)
else
    fprintf('Results file %s does not exist.', fpath)
end
