function saveSpikesInChunks(timeResolution, spikeTrainPath)
% SAVESPIKESINCHUNKS Loads a compiled spike train, processes it into binned spikes,
% and saves the results in separate .mat files.
%
% Usage:
%   saveSpikesInChunks()                % Prompts user to select a .mat file and a directory to save outputs.
%   saveSpikesInChunks(timeResolution, spikeTrainPath)
%                                       % Uses the provided time resolution and .mat file path.
%
% Inputs:
%   timeResolution - (optional) Time resolution for binning spike times. Defaults to 1 ms.
%   spikeTrainPath - (optional) Path to the .mat file containing 'SessionSplitSpikeTrains'.
%
% Outputs:
%   spikes - Processed spike times (binned), returned as a cell array.
%
% The function performs the following steps:
% 1. Loads the 'SessionSplitSpikeTrains' variable from the specified .mat file.
% 2. Converts spike times into binned spikes using 'spikeTimesToBins'.
% 3. Prompts the user to select an output directory to save .mat files.
% 4. Saves each neuron's spikes as a separate .mat file, appending session data incrementally.

% If no input is given, prompt the user to select the file containing 'SessionSplitSpikeTrains'
if nargin < 2
    spikeTrainPath = fileSelector('Select path to spikes'); % User-defined file selection function
end

if nargin < 1
    timeResolution = milliseconds(1); % Default time resolution to 1 ms
end

% Load 'SessionSplitSpikeTrains' from the specified file
S = load(spikeTrainPath, 'SessionSplitSpikeTrains');
if ~isfield(S, 'SessionSplitSpikeTrains')
    error('The selected file does not contain ''SessionSplitSpikeTrains''.');
end

spikeTrain = S.SessionSplitSpikeTrains;
numSessions = size(spikeTrain, 2); % Total number of sessions
numNeurons = size(spikeTrain, 1); % Total number of neurons

% Prompt user to select the output directory
outputDir = uigetdir([], 'Select output directory for saved spike files');
if isequal(outputDir, 0)
    error('No output directory selected. Aborting.');
end

% Step 1: Process spike times for each session
for s = 1:numSessions
    % Convert spike times to binned spikes for this session
    spikes = spikeTimesToBins(spikeTrain(:, s), timeResolution);

    for n = 1:numNeurons
        clc
        fprintf('Computing spike time array for session %d of %d...\nNeuron %d of %d...\n', s, numSessions, n, numNeurons);

        % Generate file path and variable name for this neuron
        fp = fullfile(outputDir, sprintf('spikes_%d.mat', n));
        varname = sprintf('spikes%d_%d', n, s); % Create a unique variable name per session and neuron

        % Save the current session's spikes for this neuron
        eval(sprintf('%s = spikes(n);', varname)); % Dynamically assign spikes to the variable
        if isfile(fp)
            save(fp, varname,'-append'); % Append to the neuron's file
        else
            save(fp, varname);
        end
        clear(varname); % Clear the variable to free memory
    end
end

% Step 2: Concatenate spikes across sessions for each neuron
for n = 1:numNeurons
    % Load neuron-specific file
    fp = fullfile(outputDir, sprintf('spikes_%d.mat', n));
    vars = load(fp); % Load all session data for this neuron

    % Initialize the structure with time resolution
    sp.tbin = timeResolution;

    for s = 1:numSessions
        clc
        fprintf('Concatenating spike time array for neuron %d of %d...\n', n, numNeurons);

        % Dynamically load session data and assign it to the structure
        varname = sprintf('spikes%d_%d', n, s);
        eval(sprintf('sp.spikes(s) = vars.%s;', varname));
    end

    % Save the concatenated structure back to the file
    save(fp, 'sp');
    clear('vars'); % Clear temporary variables to free memory
end
end
