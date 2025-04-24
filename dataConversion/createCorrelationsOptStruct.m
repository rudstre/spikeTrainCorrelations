function opt = createCorrelationsOptStruct(pairs,varargin)
% Define default values
defaultBinning = 4;
defaultCentralWindow = 5;
defaultThrSpikes = 0;
defaultSpikePath = '';
defaultSavePath = '';
defaultMaxLag = 50;

% Set up the input parser
p = inputParser;
addRequired(p, 'pairs');
addParameter(p, 'binning', defaultBinning, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(p, 'central_window', defaultCentralWindow, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(p, 'thr_spikes', defaultThrSpikes, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(p, 'max_lag', defaultMaxLag, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(p, 'spike_path', defaultSpikePath, @(x) ischar(x) || isstring(x));
addParameter(p, 'save_path', defaultSavePath, @(x) ischar(x) || isstring(x));

% Parse inputs
parse(p, pairs, varargin{:});

% Create the options structure
opt = struct();
opt.pairs = pairs;
opt.binning = p.Results.binning;
opt.central_window = p.Results.central_window;
opt.thr_spikes = p.Results.thr_spikes;
opt.spike_path = p.Results.spike_path;
opt.max_lag = p.Results.max_lag;
if isempty(opt.spike_path), opt.spike_path = uigetdir(); end
opt.save_path = p.Results.save_path;
if isempty(opt.save_path), opt.save_path = uigetdir(); end
opt.name = input('Enter name: ', 's');

% Save the structure to the specified save path
opt_path = fullfile(uigetdir, sprintf('opt_%s.mat', opt.name));
save(opt_path,'opt');

% Print a detailed summary to the user
fprintf('\n=========== Options Summary ===========\n');
fprintf('Name: %s\n', opt.name);
fprintf('Spike Path (Cluster): %s\n', opt.spike_path);
fprintf('Save Path: %s\n', opt.save_path);
fprintf('Binning: %d\n', opt.binning);
fprintf('Central Window: %d\n', opt.central_window);
fprintf('Max Lag: %d\n', opt.max_lag);
fprintf('Threshold Spikes: %d\n', opt.thr_spikes);
fprintf('Options saved to: %s\n', opt_path);
fprintf('=======================================\n\n');
end