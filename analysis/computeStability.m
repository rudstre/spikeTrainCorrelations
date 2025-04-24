function [observed_medDiff, nullDist, p_value] = computeStability(data)
[nx, ny, nt] = size(data);

% Flatten the (nx, ny) pairs into a single dimension for convenience
pairs_data = convertPairsFromSqFm3d(data);  % Now size is [nx*ny, 3]

for i = 2:size(pairs_data,2)
    currentData = pairs_data(:, (i-1 : i));

    % Compute observed variance for each pair (row), then median across all pairs
    observed_diffs = abs(diff(currentData, [], 2)); % variance along the time dimension
    observed_medDiff = nanmedian(observed_diffs);

    % Set the number of permutations
    nPerm = 10000;
    nullDist = nan(nPerm, 1);

    for p = 1:nPerm
        clc
        fprintf('Null %d of %d...\n',p,nPerm);
        randshift = round(length(observed_diffs) * rand());
        shuffled_data = currentData;
        shuffled_data(:,2) = circshift(shuffled_data(:,2),randshift);

        % Compute variance for the shuffled data
        shuffled_diffs = abs(diff(shuffled_data, [], 2));

        % Compute the median variance for this permutation
        nullDist(i,p) = nanmedian(shuffled_diffs);
    end
    hold on
    barWithError([2*(i-1), 2*(i-1) + 1],...
        [mean(nullDist(i,:)),observed_medDiff],...
        [nanstd(nullDist(i,:)), ...
           nanstd(observed_diffs) / sqrt( sum(isnan(observed_diffs)) - 1 )], ...
        'grouped')
end

title('Variation in z-score compared to shuffled distribution for Animal T362')
xlabel('Median variance in z-score across time')
ylabel('Probability')
tset