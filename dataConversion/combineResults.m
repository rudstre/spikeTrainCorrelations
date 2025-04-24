function res = combineResults(resPath)
% combineResults: Combines multiple result files located in a directory into
% aggregated arrays of correlation measures and their corresponding lags.

% If no input is given, prompt the user to select a directory
if nargin == 0
    resPath = uigetdir;
end

% Get the list of files in the specified directory
files = dir(resPath);

% Loop through all files starting from index 3 because:
% files(1) = '.' and files(2) = '..' are directory references, not data files
for i = 1:length(files)
    [~,~,ext] = fileparts(files(i).name);

    if ~strcmp(ext,'.mat')
        continue
    end

    fprintf('Loading file %d of %d...\n', i, length(files))
    % Load the 'results' structure from the .mat file
    data = load(fullfile(resPath, files(i).name)).results;

    % Extract z-score arrays for negative and positive lags
    zn = data.zscoreNeg;
    zp = data.zscorePos;
    ln = data.lagNeg;
    lp = data.lagPos;

    % If this is the first file, initialize the cumulative arrays
    if ~exist('zneg','var')
        [zneg, zpos, lneg, lpos] = deal(NaN(size(zn)));
    end

    znNaN = all(isnan(cat(4,zneg, zn)),4);
    zpNaN = all(isnan(cat(4,zpos, zp)),4);
    lpNaN = all(isnan(cat(4,lpos, lp)),4);
    lnNaN = all(isnan(cat(4,lneg, ln)),4);

    % Add current file's data to the cumulative arrays
    if isempty(zn)
        continue
    end
    zneg = nansum(cat(4,zneg,zn),4); zneg(znNaN) = NaN;
    zpos = nansum(cat(4,zpos,zp),4); zpos(zpNaN) = NaN;
    lneg = nansum(cat(4,lneg,ln),4); lneg(lnNaN) = NaN;
    lpos = nansum(cat(4,lpos,lp),4); lpos(lpNaN) = NaN;

end

% Adding reverse pairs to matrix
[x, y, z] = ind2sub(size(zneg), find(isnan(zneg)));
for i = 1:length(x)
    xi = x(i);
    yi = y(i);
    zi = z(i);

    if any(~isnan([zneg(xi, yi, zi), zpos(xi, yi, zi)]))
        continue
    end

    % Adjust symmetry: if zneg at (xi, yi) is zero, try to mirror from (yi, xi).
    % Similarly handle zpos, lneg, and lpos to maintain consistent symmetrical structure.
    % Note: lneg and lpos get their sign flipped for the reversed pair.
    zpos(xi, yi, zi) = zneg(yi, xi, zi);
    zneg(xi, yi, zi) = zpos(yi, xi, zi);
    lpos(xi, yi, zi) = -lneg(yi, xi, zi);
    lneg(xi, yi, zi) = -lpos(yi, xi, zi);
end

% Store the final combined results into the output structure
res.zneg = zneg;
res.zpos = zpos;
res.lneg = lneg;
res.lpos = lpos;

end
