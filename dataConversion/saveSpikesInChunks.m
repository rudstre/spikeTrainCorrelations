function saveSpikesInChunks(timeResolution, spikeTrainPath)
% SAVESPIKESINCHUNKS Processes spike trains into binned spikes and saves results session-by-session.
%
% Usage:
%   saveSpikesInChunks()
%   saveSpikesInChunks(timeResolution, spikeTrainPath)
%
% Inputs:
%   timeResolution - (optional) Time resolution for binning spike times (default = 1 ms).
%   spikeTrainPath - (optional) Path to .mat file containing 'SessionSplitSpikeTrains'.

% Handle defaults and file prompts
if nargin < 1 || isempty(timeResolution)
    timeResolution = milliseconds(1);
end
if nargin < 2 || isempty(spikeTrainPath)
    spikeTrainPath = fileSelector('Select path to spikes');
end

% Load spike data
S = load(spikeTrainPath, 'SessionSplitSpikeTrains');
if ~isfield(S, 'SessionSplitSpikeTrains')
    error('The selected file does not contain ''SessionSplitSpikeTrains''.');
end
spikeTrain = S.SessionSplitSpikeTrains;
[numNeurons, numSessions] = size(spikeTrain);

% Select output directory
outputDir = uigetdir([], 'Select output directory for saved spike files');
if isequal(outputDir, 0)
    error('No output directory selected. Aborting.');
end

% Step 1: Process each session and save per-neuron data
for s = 1:numSessions
    spikes = spikeTimesToBins(spikeTrain(:, s), timeResolution);

    for n = 1:numNeurons
        clc
        fprintf('Computing spike time array for session %d of %d...\nNeuron %d of %d...\n', ...
                s, numSessions, n, numNeurons);

        fp = fullfile(outputDir, sprintf('spikes_%d.mat', n));
        varname = sprintf('spikes%d_%d', n, s);
        sessionData.(varname) = spikes{n}; 
        
        if isfile(fp)
            save(fp, '-struct', 'sessionData', '-append');
        else
            save(fp, '-struct', 'sessionData');
        end

        clear sessionData  % Clear struct to avoid appending multiple vars
    end
end

% Step 2: Concatenate across sessions for each neuron
for n = 1:numNeurons
    clc
    fprintf('Concatenating spike time array for neuron %d of %d...\n', n, numNeurons);

    fp = fullfile(outputDir, sprintf('spikes_%d.mat', n));
    vars = load(fp);
    sp.tbin = timeResolution;
    sp.spikes = cell(1, numSessions);

    for s = 1:numSessions
        varname = sprintf('spikes%d_%d', n, s);
        sp.spikes{s} = vars.(varname);
    end

    save(fp, 'sp');  % Overwrite with clean structure
    clear vars sp
end
end