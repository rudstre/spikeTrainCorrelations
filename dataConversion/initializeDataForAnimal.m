function opt = initializeDataForAnimal(animalDataPath)
if nargin < 1
    animalDataPath = fileSelector('Select path to animal .mat workspace');
end

load(animalDataPath,'SessionSplitSpikeTrains');
pairs = nchoosek(1:length(SessionSplitSpikeTrains),1);

spike_path = input('Enter path where spikes will be on cluster: ', 's');
save_path = input('Enter path where results will be saved on cluster: ', 's');
opt = createCorrelationsOptStruct(pairs,'spike_path',spike_path,'save_path',save_path);

