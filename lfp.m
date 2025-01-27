
clear; close all

params = struct();
% Create a dialog box with dropdown
options = {'Force plate', 'Beam', 'Gait'};
[params.exptType(1), selection(1)] = listdlg('PromptString', 'Select an option:','SelectionMode', 'single','ListString', options,'ListSize', [200, 75]);

options = {'knockdown','sema3','dcnDbs','corticalDbs'};
[params.exptType(2), selection(2)] = listdlg('PromptString', 'Select an option:','SelectionMode', 'single','ListString', options,'ListSize', [200, 75]);

clear options

%% 

% Extract recording channels
tempParams = inputdlg({'Number of recording channels: '; 'Ch #1: '; 'Ch #2: '; 'Peak Hz range: '; 'NumEpochs'},'Parameters',1,{'2';'3';'8';'15 25';'1'});
params.chs(1) = str2double(tempParams{2});
params.chs(2) = str2double(tempParams{3});
params.chs = params.chs(1:str2double(tempParams{1}));
params.extractFreq = str2num(tempParams{4});
params.numEpochs = str2double(tempParams{5});

% Extract data from NS2 file
openNSx;
try
    params.filename = NS2.MetaTags.Filename;
    fpTemp = [NS2.ElectrodesInfo.ElectrodeID];
    rawData{1} = double(NS2.Data');
    params.sf = NS2.MetaTags.SamplingFreq;
catch
    params.filename = NS3.MetaTags.Filename;

    fpTemp = [NS3.ElectrodesInfo.ElectrodeID];
    rawData{1} = double(NS3.Data');
    params.sf = NS3.MetaTags.SamplingFreq;
end

for i = 1:length(params.chs)
    tempCh = find(fpTemp==params.chs(i));
    signal{1}(:,i) = rawData{1}(:,tempCh);
end

% Create parameters for computing PSD
params.timeWndw = 20;
params.timeframe = params.timeWndw;
params.cohSec = 1;
params.timeshift = 0;

clear NS2 tempParams tempCh exptType

if params.exptType(1)==1
    params.chs(end+1) = find(fpTemp==141);
    signal{2} = rawData{1}(:,params.chs(end));
end

[psdData,params] = computePsdV2(signal,params);

clear fpTemp

save([params.filename '.mat'])
