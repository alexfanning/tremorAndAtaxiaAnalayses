%   PSD analysis of behavior
%
%   Written by Alex Fanning on 4/19/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

% Input parameters
params = struct();

params.chRef1 = 0;

params.BW = 0; 
params.timecut = 0; 
params.timeshift = 0; 
params.trimStartTime = 0; 
params.trimStopTime = 390; 

params.sf = 1000; 
params.timeframe = 20; 
params.cohSec = 1;

params.plotLimits = [0,50,0,200];

% Define the frequency range for normalizing data
params.sumWndw = [35,45];

% Extract parameters from user
% prompt = {'Force plate 1 (141 or 144) and ref: ','Force plate rig 2 (143): ','Extract wndw: ',...
%     '1st epoch times (e.g. 1 200): ','2nd epoch times (e.g. 201 400): ','3rd epoch times (e.g. 401 600): ','Individual tremor events: '};
% dlgtitle = 'Parameters';
% defaults = {'141', '141','15 25','1 600' '601 1200','0','0'};
% tempParams = inputdlg(prompt,dlgtitle,1,defaults);

tempParams = inputdlg({'Plate number (141 or 143): '; 'Extraction window (Hz): '; 'Number of epochs: '; 'Group (1 = Control, 2 = Expt group A, 3 = Expt group B): ';...
    'Subject number: '; 'Individual tremor events (0 = no, 1 = yes): '},'Parameters',1,{'141'; '15 25'; '3'; '1';'1';'0'});

params.chSig = str2double(tempParams{1});
% params.chSig2 = str2double(tempParams{2});
params.extractFreqRange = str2num(tempParams{2});
params.numEpochs = str2double(tempParams{3});
params.group = str2double(tempParams{4});
params.sxNum = str2double(tempParams{5});
params.indTremor = str2num(tempParams{6});

params.extractWndw(1,:) = str2num(tempParams{4});
params.extractWndw(2,:) = str2num(tempParams{5});
params.extractWndw(3,:) = str2num(tempParams{6});

tempInput = inputdlg({'Epoch 1 time:'},...
    'Parameters',1,{'141';'3'});

params.taskNumber = str2double(tempInput{1});
params.group = str2double(tempInput{2});
params.sxNum = str2double(tempInput{3});
params.chSig = str2double(tempParams{1});
params.chSig2 = str2double(tempParams{2});
params.extractFreqRange = str2num(tempParams{3});
params.extractWndw(1,:) = str2num(tempParams{4});
params.extractWndw(2,:) = str2num(tempParams{5});
params.extractWndw(3,:) = str2num(tempParams{6});
params.indTremor = str2num(tempParams{7});
clear prompt dlgtitle defaults tempParams

if params.chSig == params.chSig2
    params.numSxs = 1;
else
    params.numSxs = 2;
end

% params.wndwNum = nnz(params.extractWndw) / 2;

%% Load NS2 data

openNSx;
params.sourceFile = [NS2.MetaTags.Filename NS2.MetaTags.FileExt];
data{1} = double(NS2.Data');
NS2 = rmfield(NS2,'Data');
params.filename = [params.sourceFile(1:end-4) '.mat'];

%%%%%%%% Re-assign channelID to legacy variables %%%%%%%%%%
params.channelSig1 = find(NS2.MetaTags.ChannelID==params.chSig);
params.channelSig2 = find(NS2.MetaTags.ChannelID==params.chSig2);

if params.chRef1 == 0
    params.channelRef1 = 0; 
else
    params.channelRef1 = find(NS4.MetaTags.ChannelID==params.chRef1);
end

if params.chRef2 == 0
    params.channelRef2 = 0; 
else
    params.channelRef2 = find(NS2.MetaTags.ChannelID==params.chRef2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.title{1} = [sprintf('%d %s',params.chSig,params.sourceFile(1:end-4))];
params.title{2} = [sprintf('%d %s',params.chSig2,params.sourceFile(1:end-4))];

if params.timecut ~= 0
    data{1} = data{1}(params.trimStartTime*params.sf+1:params.trimStopTime*params.sf,:);
end

if params.channelRef1==0
    data{1,2} = data{1}(:,params.channelSig1);
else
    data{1,2} = data{1}(:,params.channelSig1) - data(:,params.channelRef1);
end

if params.channelRef2==0
    data{1,3} = data{1}(:,params.channelSig2);
else
    data{1,3} = data{1}(:,params.channelSig2) - data(:,params.channelRef2);
end    

%% Determine coherency between LFP and force plate data

params.N = params.sf*params.timeframe; 
params.freqStep = params.sf/params.N; 
params.n1s = params.sf;          
params.nCoh = params.n1s*params.cohSec; 
params.freqCohRes = params.sf/params.nCoh;
params.endNum = fix(length(data{1,1})/params.n1s) - params.timeframe-ceil(params.timeshift);
params.shiftedFrame = params.n1s*params.timeshift;

for i = 1:params.numSxs
    for s = 1:params.endNum
        
        y1 = data{1,i+1}((1+(s-1)*params.n1s):(params.N+(s-1)*params.n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
        
        [PxxY1{i}, params.freqPsY1{i}] = pwelch(y1,hanning(params.nCoh),1/2*params.nCoh,params.nCoh,params.sf); % PSD by Wilch method and Hanning window
        data{2,i}(:,s) = PxxY1{i};

        params.freqLabel(:,s) = params.freqPsY1{i};
        params.timeLabel(1:(params.nCoh/2+1),s) = s-1;
    
    end
end

clear s i PxxY1 y1

%% Extract tremor episodes

tremor = cell(params.numSxs); move = cell(params.numSxs);
epochIdxs = cell(params.numSxs); rawData = cell(params.numSxs,params.numEpochs);
if params.indTremor == 1

    for i = 1:params.numSxs
        for t = 1:params.numEpochs
            rawData{i,t} = data{1,i+1}(params.extractWndw(t,1)*params.sf:params.extractWndw(t,2)*params.sf);
            [tremor{i},move{i},epochIdxs{i},data] = extractDualFPtremor(rawData{i,t},data,params,tremor{i},move{i},epochIdxs{i},t);
            for k = 1:length(tremor{i}{t,5})
                tremAmp(i,t) = mean(tremor{i}{t,5}{k}(:));
                if k ~= length(tremor{i}{t,5})
                    iei{i}(k,t) = epochIdxs{i}{t,1}(k+1,1) - epochIdxs{i}{t,1}(k,2);
                end
            end
            tremDur(i,t) = mean(tremor{i}{t,10}(:));
            tremFreq(i,t) = nanmean(tremor{i}{t,9}(:));
            tremOccur(i,t) = length(tremor{i}{1,t}) / sum([length(move{i}{t}),length(tremor{i}{1,t})]);

        end
    end
end

% histogram(iei{1}(:,1),100)

%% Plot raw data

figure('Name','Raw data')
ax1 = subplot(2,1,1); hold on
plot(data{1,2}(:,1))
ylabel('\muV')
title(ax1,'Force plate 1')
set(gca,'FontSize',16)

ax2 = subplot(2,1,2); hold on
plot(data{1,3}(:,1))
ylabel('\muV')
title(ax2,'Force plate 2')
xlabel('Timepoints')
set(gca,'FontSize',16)
linkaxes([ax1 ax2],'x')

clear ax1 ax2

%% Plot PSDs

for i = 1:params.numSxs
    figure(i+1);
    gca3 = pcolor(params.timeLabel,params.freqLabel,data{2,i});
    set(gca3, 'LineStyle','none');
    colorbar('location','eastoutside');
    caxis([0 1000]);
    title(['PSD: ',params.title{i}]); 
    xlabel('Time (seconds)'); 
    ylabel('Frequency (Hz)'); 
    axis([-inf,inf,0,50]);
    set(gca,'FontSize',16)
    saveas(i+1,[params.title{i} '_PSD.fig']);
end

clear gca3

%% Compute mean PSD and normalize data

% Frequency band cutoffs
bands = [1 3; 4 7; 8 12; 13 30];

% Define the Endtime by subtract the timeframe (for internal consistency of previous version)
params.extractWndw(:,3) = params.extractWndw(:,2) - params.timeframe;
mPSD = cell(1); nMeanPsd = cell(1);

% Mean PSD and normalization
for t = 1:params.numSxs
    for i = 1:params.numEpochs
        mPSD{t}(:,i) = mean(data{2,t}(:,params.extractWndw(i,1):params.extractWndw(i,3)),2);
        params.sumPSD(t,i) = sum(mPSD{t}(params.sumWndw(1):params.sumWndw(2),i));
        nMeanPsd{t}(:,i) = mPSD{t}(:,i)/params.sumPSD(t,i);

        [maxPsd(t,i),maxPsdFreq(t,i)] = max(mPSD{t}(params.extractFreqRange(1):params.extractFreqRange(2),i));
        [nMaxPsd(t,i),nMaxPsdFreq(t,i)] = max(nMeanPsd{t}(params.extractFreqRange(1):params.extractFreqRange(2),i));

        maxFreq(t,i) = (params.extractFreqRange(1) + maxPsdFreq(t,i) - 1) * params.freqCohRes;
        nMaxFreq(t,i) = (params.extractFreqRange(1) + nMaxPsdFreq(t,i) - 1) * params.freqCohRes;

        nPsdAucAll{t,i} = cumtrapz(nMeanPsd{t}(1:35,i));
        for k = 2:length(nPsdAucAll{t,i})
            nPsdAucTotal{t}(i,k) = nPsdAucAll{t,i}(k) - nPsdAucAll{t,i}(k-1);
        end
        nPsdAuc(t,i) = trapz(nMeanPsd{t}(params.extractFreqRange(1):params.extractFreqRange(2),i));
        nPsdAucPct(t,i) = (nPsdAuc(t,i) / nPsdAucAll{t,i}(end)) * 100;

        for m = 1:4
            nPsdBands{t}(i,m) = sum(nMeanPsd{t}(bands(m,1):bands(m,2),i));
            nPsdBandAuc{t}(i,m) = sum(nPsdAucTotal{t}(i,bands(m,1):bands(m,2)));
            nPsdBandsAucPct{t}(i,m) = (nPsdBandAuc{t}(i,m) / nPsdAucAll{t,i}(end)) * 100;
        end

        pk40hz(t,i) = max(mPSD{t}(35:45,i));
    end
end

%% Plot normalized data

figure(4);
params.legData = {[num2str(params.extractWndw(1,1)) '-' num2str(params.extractWndw(1,2))],[num2str(params.extractWndw(2,1)) '-' num2str(params.extractWndw(2,2))],...
            [num2str(params.extractWndw(3,1)) '-' num2str(params.extractWndw(3,2))]};
params.col = {'k','r','b'};
for t = 1:params.numSxs
    subplot(2,1,t); hold on
    for i = 1:params.numEpochs
        plot(params.freqPsY1{t},mPSD{t}(:,i),'Color',params.col{i},'LineWidth',1);
    end
    title(['mean PSD: ' params.title{t}]);
    xlabel('Frequency (Hz)'); 
    ylabel('Mean PSD (dB/Hz)');
    axis([params.plotLimits(1:3),max(maxPsd(t,:))*1.2]);
    if t == 1
        legend(params.legData{1:params.numEpochs})
    end
    set(gca,'FontSize',16)
end

figure(5);
for t = 1:params.numSxs
    subplot(2,1,t); hold on
    for i = 1:params.numEpochs
        plot(params.freqPsY1{t},nMeanPsd{t}(:,i),'Color',params.col{i},'LineWidth',1);
    end
    title(['Normalized PSD: ' params.title{t}]);
    xlabel('Frequency (Hz)'); 
    ylabel('Norm PSD (dB/Hz)');
    xlim([0 50])
    %axis([params.plotLimits(1:3),1*1.2]);
    if t == 1
        legend(params.legData{1:params.numEpochs})
    end
    set(gca,'FontSize',16)
end

saveas(4,[params.title{1} '_normPSD.fig']);

%% Write data to Excel

params.cellFormat{1} = ['A', num2str(1)]; params.cellFormat{2} = ['F', num2str(1)];
params.cellFormat{3} = ['B', num2str(2)]; params.cellFormat{4} = ['G', num2str(2)];
params.cellFormat{5} = ['B', num2str(8)]; params.cellFormat{6} = ['G', num2str(8)];
params.cellFormat{7} = ['B', num2str(10)]; params.cellFormat{8} = ['G', num2str(10)];
params.cellFormat{9} = ['B', num2str(14)]; params.cellFormat{10} = ['G',num2str(14)];
sumStats{1} = {'normMaxPSDandFreq','normPSD','maxPSDandFreq','meanPSD'};
sumStats{2} = {'Force plate 1','Force plate 2'};
sumStats{3} = {'Peak nPSD';'Peak nPSD freq.';'Amp. (uVs)';'Dur. (ms)';'Auc'; 'Auc Pct';'Occurrence';'nPsdBands';'nPsdBandsAucPct'};
sumStats{4} = {'Peak nPSD';'Peak nPSD freq.';'Auc'; 'Auc Pct';'nPsdBands';'nPsdBandsAucPct'};
params.file = [params.filename(1:end-4) '.xlsx'];
for i = 1:length(sumStats{1})
    for t = 1:params.numSxs
        writecell(sumStats{2}(t),params.file,'Sheet',sumStats{1}{i},'Range',params.cellFormat{t})
        if params.indTremor == 1
            writematrix([nMaxPsd(t,:); nMaxFreq(t,:); tremAmp(t,:); tremDur(t,:); nPsdAuc(t,:); nPsdAucPct(t,:); tremOccur(t,:)],params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+2})
            writematrix([nPsdBands{t}; nPsdBandsAucPct{t}],params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+6})
            writematrix(pk40hz,params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+8})
            writecell(sumStats{3},params.file,'Sheet',sumStats{1}{1},'Range','A2:A10')
        else
            writematrix([nMaxPsd(t,:); nMaxFreq(t,:); nPsdAuc(t,:); nPsdAucPct(t,:)],params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+2})
            writematrix([nPsdBands{t}; nPsdBandsAucPct{t}],params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+4})
            writematrix(pk40hz,params.file,'Sheet',sumStats{1}{1},'Range',params.cellFormat{t+8})
            writecell(sumStats{4},params.file,'Sheet',sumStats{1}{1},'Range','A2:A7')
        end
        writematrix(nMeanPsd{t},params.file,'Sheet',sumStats{1}{2},'Range',params.cellFormat{t+2})
        writematrix([maxPsd(t,:); maxFreq(t,:)],params.file,'Sheet',sumStats{1}{3},'Range',params.cellFormat{t+2})
        writematrix(mPSD{t},params.file,'Sheet',sumStats{1}{4},'Range',params.cellFormat{t+2})
    end
end
% 
% params.xlCellId = {['B' num2str(params.sxNum)],['H' num2str(params.sxNum)]; ['B' num2str(params.sxNum+24)],['H' num2str(params.sxNum+24)];...
%     ['B' num2str(params.sxNum+49)],['C' num2str(params.sxNum+49)]; ['F' num2str(params.sxNum+49)],['G' num2str(params.sxNum+49)]};
% params.xlCellId2 = {['A' num2str(1)];['A' num2str(2)];['A' num2str(25)];['A' num2str(50)];['E' num2str(50)];['G' num2str(1)]};
% params.xlTabNames = {'Ch 37','Ch 42','Ch 47','Ch 35','Ch 49','Ch 44'};
% params.xlTabNames2 = {'Ch X','Ch Y','Ch Z'};
% params.dataNames = {'Ctrl';'nPsdBands';'nPsdBandsAucPct';'nMaxPsd';'nMaxPsdFreq';'Ataxia'};
% for k = 1:2
%     for m = 1:3
%         writematrix(nPsdBands{k}(m,:),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames{k*m},'Range',params.xlCellId{1,params.group})
%         writematrix(nPsdBandsAucPct{k}(m,:),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames{k*m},'Range',params.xlCellId{2,params.group})
%         writematrix(nMaxPsd{k}(m),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames{k*m},'Range',params.xlCellId{3,params.group})
%         writematrix(nMaxPsdFreq{k}(m),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames{k*m},'Range',params.xlCellId{4,params.group})
%         for d = 1:6
%             writecell(cellstr(params.dataNames{d}),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames{k*m},'Range',params.xlCellId2{d})
%         end
%     end
% end
% for k = 1:3
%     writematrix(data2export.nPsdBands{j}(k,:),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{1,params.group})
%     writematrix(data2export.nPsdBandsAucPct{j}(k,:),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{2,params.group})
%     writematrix(data2export.nMaxPsd{j}(k),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{3,params.group})
%     writematrix(data2export.nMaxPsdFreq{j}(k),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId{4,params.group})
%     for d = 1:6
%         writecell(cellstr(params.dataNames{d}),[params.sheetNames{j} '_summaryStats.xlsx'],'Sheet',params.xlTabNames2{k},'Range',params.xlCellId2{d})
%     end
% end


clear i t k maxPsdFreq nMaxPsdFreq
save(params.filename)
