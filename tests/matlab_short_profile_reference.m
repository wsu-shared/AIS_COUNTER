function matlab_short_profile_reference(basepath, clicks, outfile)
% MATLAB_SHORT_PROFILE_REFERENCE  Where the original stops instead of measuring.
%
% ``ais_auto.m``'s sliding mean runs ``lv_ss(i) = mean([lv_smooth(1:i) lv_smooth(i:i+d)])``
% for ``i <= d``, so it reads index ``2d``. On a profile shorter than that MATLAB raises
% "Index exceeds the number of array elements" and the script produces no length at all --
% it is not a bad measurement, it is the absence of one. ``aiscounter`` clamps the window
% instead (there is nothing to be faithful to) and must flag the row rather than report the
% clamped number as a measurement; this fixture is what pins that boundary.
%
% *clicks* is an N-by-2 matrix of 1-based (column, row) clicks. Each one is run through
% ``matlab_rethreshold_reference`` and recorded as either its result or the error it raised,
% so the outcome is ground truth either way.
%
% Regenerate with:
%   /Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('tests'); \
%       matlab_short_profile_reference('<base>', [184 521; 556 1023], 'out.json')"

results = cell(1, size(clicks,1));
tmp = [tempname '.json'];

for k = 1:size(clicks,1)
    entry = struct('click_col', clicks(k,1), 'click_row', clicks(k,2));
    try
        matlab_rethreshold_reference(basepath, 0, clicks(k,1), clicks(k,2), tmp);
        raw = fileread(tmp);
        entry.ok      = true;
        entry.error   = '';
        entry.result  = jsondecode(raw);
        entry.n_pix   = entry.result.n_pix;
        entry.lngth   = entry.result.lngth;
    catch err
        % The whole point of the fixture: a click the original cannot answer.
        entry.ok    = false;
        entry.error = err.message;
    end
    results{k} = entry;
end

if exist(tmp, 'file'); delete(tmp); end

fid = fopen(outfile, 'wt');
fprintf(fid, '%s', jsonencode(results));
fclose(fid);
fprintf('WROTE %s  (%d click(s))\n', outfile, numel(results));
end
