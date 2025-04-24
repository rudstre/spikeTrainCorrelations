function groups = generatePairGroups(ncells, nworkers)
    % Generate all possible pairs from the cells
    pairs = nchoosek(1:ncells, 2);
    
    % Total number of pairs
    total_pairs = size(pairs, 1);
    
    % Initialize groups for workers
    groups = cell(1, nworkers);  % Each worker gets a cell array of pairs
    
    % Distribute the pairs among workers
    pairs_per_worker = floor(total_pairs / nworkers);  % Base number of pairs each worker gets
    remainder = mod(total_pairs, nworkers);  % Extra pairs to distribute
    
    pair_index = 1;  % Start index for assigning pairs
    
    % Assign pairs to workers
    for i = 1:nworkers
        % Calculate number of pairs for this worker (handling remainder)
        num_pairs = pairs_per_worker + (i <= remainder);
        
        % Assign pairs to this worker
        if pair_index + num_pairs - 1 <= total_pairs
            groups{i} = pairs(pair_index:pair_index + num_pairs - 1, :);
            pair_index = pair_index + num_pairs;  % Update pair index
        else
            % If we run out of pairs, assign empty group
            groups{i} = [];
        end
    end
    
    % Display results
    for i = 1:nworkers
        if ~isempty(groups{i})
            fprintf('Worker %d is assigned %d pairs\n', i, size(groups{i}, 1));
        else
            fprintf('Worker %d is assigned 0 pairs\n', i);
        end
    end
end
