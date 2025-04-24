function spikeTrainCorrelation(optPath, workerId, totalWorkers)
% Spike train cross-correlation with LRU cache for overlapping units
% Args:
%   optPath: Path to the options .mat file
%   workerId: Worker ID for parallel processing
%   totalWorkers: Total number of workers

% Initialize global variables
global startTime currentWorkerId
startTime = tic;
currentWorkerId = workerId;

% Load configuration options
options = load(optPath).opt;
pairsToProcess = options.pairs;
binSize = 4; % number of weeks per time bin
centralWindowSize_ms = 5; % range around which to detect peaks
baselineMaxLag_ms = 50;

% Assign pairs to the current worker
pairGroups = generatePairGroups(max(pairsToProcess(:)), totalWorkers);
workerPairs = pairGroups{workerId};
options.pairs = workerPairs;

% Exit early if no pairs assigned
if isempty(workerPairs)
    logMessage('No pairs to process. Skipping...\n');
    saveResults(options, workerId, [], [], [], [], 'No pairs to process.');
    return;
end

% Open the spikes file and initialize the cache
logMessage('Opening spikes matfile: %s\n', options.spike_path);
spikeCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
lruQueue = []; % LRU queue for cache eviction
cacheLimit = 2;

logMessage('Processing %d pairs out of %d total pairs.\n', size(workerPairs, 1), size(options.pairs, 1));

% Process each pair assigned to this worker
for pairIdx = 1:size(workerPairs, 1)
    logMessage('Progress: Pair %d of %d\n', pairIdx, size(workerPairs, 1));
    neuronPair = workerPairs(pairIdx, :);

    % Load or retrieve spike trains for the neuron pair
    [spikesNeuron1, lruQueue] = loadFromCache(spikeCache, lruQueue, neuronPair(1), options.spike_path, cacheLimit);
    [spikesNeuron2, lruQueue] = loadFromCache(spikeCache, lruQueue, neuronPair(2), options.spike_path, cacheLimit);

    if pairIdx == 1
        numSessions = size(spikesNeuron1.spikes,2);
        numTimeGroups = ceil(numSessions / binSize);
        numNeurons = max(pairsToProcess(:));
        zscorePosMax = NaN(numNeurons, numNeurons, numTimeGroups);
        zscoreNegMax = NaN(numNeurons, numNeurons, numTimeGroups);
        lagPosMax = NaN(numNeurons, numNeurons, numTimeGroups);
        lagNegMax = NaN(numNeurons, numNeurons, numTimeGroups);
        binStarts = (1:binSize:numSessions)';
        sessionBins = [binStarts,binStarts + binSize - 1];
        sessionBins(sessionBins > numSessions) = numSessions;
        baselineMaxLag = baselineMaxLag_ms/milliseconds(spikesNeuron2.tbin);
        centralWindowSize = centralWindowSize_ms/milliseconds(spikesNeuron2.tbin);
    end

    % Process in time bins
    for timeGroupIdx = 1:numTimeGroups
        % Define the time segment
        firstSession = sessionBins(timeGroupIdx,1);
        lastSession = sessionBins(timeGroupIdx,2);

        % Extract spikes for the current segment
        segmentSpikesNeuron1 = [spikesNeuron1.spikes{firstSession:lastSession}];
        segmentSpikesNeuron2 = [spikesNeuron2.spikes{firstSession:lastSession}];

        % Compute cross-correlation
        [crossCorr, lagValues_samp] = xcorr(segmentSpikesNeuron2, segmentSpikesNeuron1, baselineMaxLag);
        lagValues = lagValues_samp * milliseconds(spikesNeuron2.tbin);

        % Define baseline lag ranges
        baselinePosLags = 10:baselineMaxLag;
        baselineNegLags = -baselinePosLags(end:-1:1);
        baselineIdxs = ismember(lagValues, [baselineNegLags, baselinePosLags]);

        % Calculate baseline statistics
        baselineValues = crossCorr(baselineIdxs);
        % Under the Poisson model, the mean of the baseline is our lambda
        lambda = mean(baselineValues);
        if lambda*centralWindowSize < 2
            continue
        end

        % If max is at 0 and super high, probably actually the same unit
        [maxCorr,lagMax] = max(crossCorr);
        if lagValues_samp(lagMax) == 0 && maxCorr > 2*lambda
            continue  
        end

        % Positive lags
        zscorePosMax(neuronPair(1), neuronPair(2), timeGroupIdx) = ...
            computeExtrema(crossCorr, lagValues, [1 centralWindowSize], lambda, 2);

        % Negative lags
        zscoreNegMax(neuronPair(1), neuronPair(2), timeGroupIdx) = ...
            computeExtrema(crossCorr, lagValues, [-centralWindowSize -1], lambda, 2);
    end
end

% Save results
saveResults(options, workerId, zscorePosMax, zscoreNegMax, lagPosMax, lagNegMax);
logMessage('Results saved successfully.\n');
end


function zFinal = computeExtrema(ccf, lags, lagRange, lambda, ~)
% computeExtremaWindow
%
%   [zFinal, lagFinal] = computeExtremaWindow(ccf, lags, lagRange, lambda, ~)
%
%   Sums all counts within lagRange once and computes a single z‐score
%   against Poisson(expectedSum).  No bias corrections.

if sum(ccf(iswithin(lags,-1,1))) < 1
    % collision correction
    onePos = abs(lagRange) == 1;
    sgn = lagRange(onePos);
    lagRange(onePos) = sgn + 2 * sgn;
end

% restrict to specified lag range
idx = iswithin(lags, lagRange(:));
ccf_sub  = ccf(idx);

% total counts in window and expected sum
totalCount   = sum(ccf_sub);
expectedSum  = numel(ccf_sub) * lambda;

% helper for asymptotic z‐approximation
zApprox = @(x,lam) sqrt(-2 * (-lam - x*log(x/lam) + x - 0.5*log(2*pi*x)));

% compute raw z
if totalCount > expectedSum
    % excitation (upper‐tail)
    p = 1 - poisscdf(totalCount-1, expectedSum);
    if p > 0
        zFinal = norminv(1 - p);
    else
        zFinal = zApprox(totalCount, expectedSum);
    end

elseif totalCount < expectedSum
    % inhibition (lower‐tail)
    p = poisscdf(totalCount, expectedSum);
    if p > 0
        zFinal = norminv(p);
    else
        zFinal = -zApprox(max(totalCount,1), expectedSum);
    end

else
    zFinal = 0;
end

end

function [spikes, queue] = loadFromCache(cache, queue, unitId, spikePath, maxCacheSize)
% Retrieve spikes from cache or load from file if not cached
try
    if isKey(cache, unitId)
        spikes = cache(unitId);
        queue(queue == unitId) = [];
        queue = [queue, unitId]; % Update LRU queue
    else
        spikes = load(fullfile(spikePath, sprintf('spikes_%d.mat', unitId))).sp;
        cache(unitId) = spikes;
        queue = [queue, unitId];
        if numel(queue) > maxCacheSize
            remove(cache, queue(1)); % Evict least recently used item
            queue(1) = [];
        end
    end
catch ME
    error('Error loading spikes for unit %d: %s', unitId, ME.message);
end
end

function saveResults(options, workerId, zscorePos, zscoreNeg, lagPos, lagNeg, message)
% Save results to a .mat file
results.options = options;
results.zscorePos = zscorePos;
results.zscoreNeg = zscoreNeg;
results.lagPos = lagPos;
results.lagNeg = lagNeg;
if exist('message', 'var')
    results.message = message;
end
if ~isfolder(options.save_path)
    mkdir(options.save_path);
end
save(fullfile(options.save_path, sprintf('results_%s_%d.mat', options.name, workerId)), 'results', '-v7.3');
end

function logMessage(message, varargin)
% Log worker-specific messages with timestamps
global startTime currentWorkerId;
fprintf('[%.2fs] [Worker %d] %s', toc(startTime), currentWorkerId, sprintf(message, varargin{:}));
end
