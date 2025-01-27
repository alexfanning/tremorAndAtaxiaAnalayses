function tremor = extractTremor(data,parameters,file) 

%   Calculates change in fluorescence divided by isobestic fluorescence
%
%   Uses 20th percentile value calculated across a sliding window
%   on the raw fluorescence
%
%   f = raw fluorescence signal
%   window = width of the sliding Window
%   f0 = baseline fluorescence
%
%   Alex Fanning, September 2022
% *************************************************************************

prompt = {'Window size (ms): ', 'Bandpass lower limit: ', 'Bandpass upper limit: '};
dlgtitle = 'Sliding window';
default = {'60', '10', '30'};
paramsTemp = inputdlg(prompt,dlgtitle,1,default);
parameters(1).wndwSize = str2double(paramsTemp{1});
parameters(1).bpRange = [str2double(paramsTemp{2}) str2double(paramsTemp{3})];

% Bandpass filter the data for a flat baseline
bpVec = bandpass(data,parameters(1).bpRange,1000);

% Loop through each sample to calculate f0
% [f0,g0] = slideWndw(parameters,data);
[f0bp,g0bp] = slideWndw(parameters,bpVec);

%% Calculate difference between bounds and thresholds
% Bounds difference and tremor start/stop threshold
bndDiff = NaN(1,length(f0bp));
bndDiff = g0bp - f0bp;

% baseline threshold
bndDiffThr = prctile(bndDiff,35);

% Find threshold for movement event detection based on cumulative sum
movThr = prctile(bndDiff,80);

%% Check bandpass filtering in relation to raw data

figure('Name','Raw data and bounds')
ax1 = subplot(2,1,1); hold on
plot(data)
ylabel('\muVs')
set(gca,'FontSize',14)

ax2 = subplot(2,1,2); hold on
plot(bpVec)
linkaxes([ax1 ax2], 'x')
xlabel('Time (ms)')
ylabel('\muVs')
set(gca,'FontSize',14)

%% Plot raw data, bandpass vector, and difference of bounds

figure('Name','Bandpass data and bounds')
ax = subplot(3,1,1); hold on
plot(data)
ylabel('\muVs')
set(gca,'FontSize',14)

ax1 = subplot(3,1,2); hold on
plot(bpVec)
plot(f0bp)
plot(g0bp)
ylabel('\muVs')
set(gca,'FontSize',14)

ax2 = subplot(3,1,3); hold on
plot(bndDiff)
yline(bndDiffThr,'LineWidth',2,'Color','r')
yline(movThr,'LineWidth',2)
linkaxes([ax ax1 ax2], 'x')
xlabel('Time (ms)')
ylabel('Difference')
legend('','Tremor endpts threshold','Movement threshold')
set(gca,'FontSize',14)
%set(ax1,'ylim',[-50 50])
%set(ax2,'ylim',[0 

%% Plot histogram of difference of bounds with threshold indicated

figure('Name','Movement threshold'); hold on
hist(bndDiff,1000,'k')
xline(movThr,'r','LineWidth',1)
xlabel('Difference value')
ylabel('Number of datapoints')
set(gca,'FontSize',16)
title('Difference between bounds')

%% Extract potential tremor epoch start and end times

% Find points above movement threshold
tempVec = zeros(1,length(bndDiff));
for i = 1:length(bndDiff)
    if bndDiff(i) > movThr
        tempVec(i) = 1;
    end
end

% Find start and end times for tremor epochs
a = 1;
for j = 1:length(tempVec)
    if tempVec(j) == 1 && tempVec(j-1) == 0
        tempTimepts(a,1) = j;
        
        % Find start timept
        for k = 1:j
            if bndDiff(k) < bndDiffThr
                startTemp = k;
            end
        end

        if exist('startTemp')
            segStart(a) = startTemp;

        % Find stop timept
            for ii = j:length(bndDiff)
                if bndDiff(ii) < bndDiffThr
                    segStop(a) = ii;
                    break
                end
            end
        end

        a = a + 1;
    end
end

segStart = unique(segStart);
segStop = unique(segStop);

if segStart(1) == 0
    segStart = segStart(2:end);
    segStop = segStop(2:end);
end

% Create cell array structure to hold tremor epochs
for m = 1:length(segStop)
    allEpochs{m} = data(segStart(m):segStop(m));
end

%% Determine frequency of each potential tremor event

proceed = 'n';
while proceed == 'n'
    % Parameters for finding peaks
    prompt = {'minPeakProm: ','minPeakInt: ','maxPeakWidth: ','minPeakWidth: '};
    dlgtitle = 'FindPeaks params';
    default = {'3500','25','75','5'};
    tempParams = inputdlg(prompt,dlgtitle,1,default);
    parameters.minPeakProm = str2double(tempParams{1});
    parameters.minPeakDist =  str2double(tempParams{2});
    parameters.maxPeakWidth =  str2double(tempParams{3});
    parameters.minPeakWidth =  str2double(tempParams{4});
    
    % Find peaks of each tremor event
    for i = 1:length(allEpochs)
        if length(allEpochs{i}) >= 100
            [pks,idxs,w,p] = findpeaks(allEpochs{i},1,'MinPeakProminence',parameters.minPeakProm,'MinPeakDistance',parameters.minPeakDist,'Annotate','extents','MaxPeakWidth',parameters.maxPeakWidth,'MinPeakWidth',parameters.minPeakWidth);
            peaks{i} = pks; locs{i} = idxs; widths{i} = w; otherPeaks{i} = p;
        end
        % For trough detection: Flip the sign of each tremorEpoch and rerun
        % peak detection to get locs and peaks
    end
    
    % Plot individual tremor events
    counter = 1;
    for t = 1:2
        figure()
        for i = counter:counter + 20 -1
            nexttile; hold on
            plot(allEpochs{i},'LineWidth',2,'Color','k')
            peakIdx = locs{i}; peakAmp = peaks{i};
            for t = 1:length(locs{i})
                scatter(peakIdx(t),peakAmp(t),100,'LineWidth',2)
            end
        end
        counter = 21;
    end

    prompt = {'Keep peaks? (y/n): '};
    dlgtitle = 'Peak analysis';
    proceed = inputdlg(prompt,dlgtitle);
    proceed = proceed{1};
end

%%%%%%% INVERT tremorEpochs DATA AND TAKE PEAKS TO CAPTURE NEGATIVE PEAKS
%%%%%%% FOR AMPLITUDE MEASUREMENTS

%% Delinate movement from tremor events

a = 1;
for j = 1:length(peaks)
    if length(peaks{j}) > 2
        newEpochs{a} = allEpochs{j};
        newPeaks{a} = peaks{j};
        newLocs{a} = locs{j};
        newAbsPks{a} = otherPeaks{j};
        idxHold(a) = j;

        % Find time difference between peaks
        for t = 1:length(newLocs{a}) - 1
            tremPeakDiff{a}(t) = newLocs{a}(t+1) - newLocs{a}(t);
        end
        a = a + 1;
    end
end

% Calculate frequency of each tremor episode
tremFreq = cellfun(@(x) 1/(mean(x)/1000),tremPeakDiff,'UniformOutput',false);

b = 1;
for i = 1:length(tremFreq)
    if tremFreq{i} > 25 || tremFreq{i} < 15
        tremFreq{i} = {};
    else
        tremFreqActual(b) = tremFreq{i};
        b = b + 1;
    end
end

tremFreqIdxs = find(~cellfun(@isempty,tremFreq));
tremIdxs = idxHold(tremFreqIdxs);
move = allEpochs;
move(tremIdxs) = [];

tremor{1} = newEpochs(tremFreqIdxs);
tremor{2} = newPeaks(tremFreqIdxs);
tremor{3} = newLocs(tremFreqIdxs);

tremor{4} = newAbsPks{1};
for i = 2:length(newAbsPks)
    tremor{4} = cat(1,tremor{4},newAbsPks{i});
end

tremor{5} = tremor{2}{1};
for i = 2:length(tremor{2})
    tremor{5} = cat(1,tremor{5},tremor{2}{i});
    tremor{6}(i) = length(tremor{2}{i}); 
end

if length(tremor{3}) >= 20
    counter = 1;
    for t = 1:2
        while counter < length(tremor{1})
            figure()
            for i = counter:counter + 20 - 1
                nexttile; hold on
                plot(tremor{1}{i},'LineWidth',2,'Color','k')
                peakIdx = tremor{3}{i}; peakAmp = tremor{2}{i};
                for b = 1:length(tremor{3}{i})
                    scatter(peakIdx(b),peakAmp(b),100,'LineWidth',2)
                end
            end
            counter = 21;
        end
    end
end

%%  Tremor epoch summary statistics

figure(); hold on
histogram(tremor{4},30)
xlabel('Peak amplitude (\muVs)')
ylabel('# of events')
set(gca,'FontSize',16)

figure(); hold on
histogram(tremor{6})
xlabel('# of peaks')
ylabel('# of epochs')
set(gca,'FontSize',16)

figure(); hold on
histogram(tremFreqActual,15)
xlabel('Tremor frequency')
ylabel('# of epochs')
set(gca,'FontSize',16)

save([file, '_tremorEpochs.mat'])

%% may need to align events based on first peak using circshift

%% may need to take subset of non-tremor movements (e.g. 1-2 hz oscillations with multiple peaks)

%% plot summary statistics for all tremor epochs (i.e. amplitude, duration)
