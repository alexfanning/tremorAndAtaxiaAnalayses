function [tremor,move,epochIdxs,data] = extractDualFPtremor(signal,data,parameters,tremor,move,epochIdxs,type) 

%   Extracts movement and tremor epochs.
% 
%   Upper and lower bounds of the
%   bandpass filtered data are calculated based on percentiles and the
%   difference of the bounds vector is used to calculate movement and
%   baseline thresholds.
%
%   Alex Fanning, September 2022
% *************************************************************************

prompt = {'Window size (ms): ', 'diff lower limit: ', 'diff upper limit: ','Bandpass lower limit: ', 'Bandpass upper limit: ','Baseline thresh: ', 'Movement thresh: '};
dlgtitle = 'Sliding window';
default = {'60','25','75','10','30','40','90'};
paramsTemp = inputdlg(prompt,dlgtitle,1,default);
parameters.wndwSize = str2double(paramsTemp{1});
parameters.diffRange = [str2double(paramsTemp{2}) str2double(paramsTemp{3})];
parameters.bpRange = [str2double(paramsTemp{4}) str2double(paramsTemp{5})];
parameters.baseThr = str2double(paramsTemp{6});
parameters.moveThr = str2double(paramsTemp{7});

% Bandpass filter the data for a flat baseline
signal = detrend(signal);
bpVec = bandpass(signal,parameters.bpRange,parameters.sf);
bp40hzVec = bandpass(signal,[35 45],parameters.sf);

% Loop through each sample to calculate f0
[f0bp,g0bp] = slideWndw(parameters,bpVec);
data{3,1} = bpVec;

%% Difference of bounds and tremor start/stop threshold

goodThresh = 'n';
while goodThresh == 'n'
    bndDiff = NaN(1,length(f0bp));
    bndDiff = g0bp - f0bp;
    bndDiff = smooth(bndDiff,100);

    % baseline threshold
    bndDiffThr = prctile(bndDiff,parameters(1).baseThr);

    % Threshold for movement event detection
    movThr = prctile(bndDiff,parameters(1).moveThr);

    % Plot raw data, bandpass vector, and difference of bounds

    figure('Name','Bounds and difference vectors')
    ax = subplot(3,1,1); hold on
    plot(signal)
    ylabel('\muVs')
    ylim([-1000 1000])
    title('Raw')
    set(gca,'FontSize',16)

    ax1 = subplot(3,1,2); hold on
    plot(bpVec)
    plot(f0bp)
    plot(g0bp)
    ylabel('\muVs')
    title('Bandpass vector and bounds')
    ylim([-250 250])
    legend('','Upper bnd','Lower bnd')
    set(gca,'FontSize',16)

    ax2 = subplot(3,1,3); hold on
    plot(bndDiff)
    yline(bndDiffThr,'LineWidth',2,'Color','r')
    yline(movThr,'LineWidth',2)
    linkaxes([ax ax1 ax2], 'x')
    xlabel('Timepoints')
    ylabel('Diff')
    title('Diff of bounds')
    legend('','Tremor endpts thresh','Movement thresh')
    ylim([0 500])
    set(gca,'FontSize',16)

    prompt = 'Good thresholds? (y/n): ';
    goodThresh = input(prompt,'s');
    if goodThresh == 'y'
        continue
    else
        prompt = {'Window size (ms): ', 'diff lower limit: ', 'diff upper limit: ','Bandpass lower limit: ', 'Bandpass upper limit: ','Baseline thresh: ', 'Movement thresh: '};
        dlgtitle = 'Sliding window';
        default = {'60','25','75','10','30','40','98'};
        paramsTemp = inputdlg(prompt,dlgtitle,1,default);
        parameters(1).wndwSize = str2double(paramsTemp{1});
        parameters(1).diffRange = [str2double(paramsTemp{2}) str2double(paramsTemp{3})];
        parameters(1).bpRange = [str2double(paramsTemp{4}) str2double(paramsTemp{5})];
        parameters(1).baseThr = str2double(paramsTemp{6});
        parameters(1).moveThr = str2double(paramsTemp{7});
    end
end

close all

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
for j = 2:length(tempVec)
    if tempVec(j) == 1 && tempVec(j-1) == 0
        tempTimepts(a,1) = j;
        
        % Find start timept
        for k = 1:j
            if bndDiff(k) < bndDiffThr
                startTemp = k;
            end
        end
        if exist('startTemp','var')
            segStart(a) = startTemp;
        else
            continue
        end

        % Find stop timept
        for ii = j:length(bndDiff)
            if bndDiff(ii) < (bndDiffThr)
                segStop(a) = ii;
                break
            end
        end

        a = a + 1;
    end
end

epochEndPts(:,1) = unique(segStart)';
if length(unique(segStart)) ~= length(unique(segStop))
    epochEndPts(end) = [];
end
epochEndPts(:,2) = unique(segStop)';

% tempEndPts(:,1) = unique(segStart)';
% if length(unique(segStart)) ~= length(unique(segStop))
%     tempEndPts(end) = [];
% end
% tempEndPts(:,2) = unique(segStop)';

% a = 1;
% for i = 1:length(tempEndPts)-1
%     if tempEndPts(i+1,1) - tempEndPts(i,2) < 150
%         epochEndPts(a,1) = tempEndPts(i,1);
%         epochEndPts(a,2) = tempEndPts(i+1,2);
%         a = a + 1;
%     else
%         epochEndPts(a,1) = tempEndPts(i,1);
%         epochEndPts(a,2) = tempEndPts(i,2);
%         a = a + 1;
%     end
% end

% Create cell array structure to hold tremor epochs
allEpochs = cell(1);
for m = 1:length(epochEndPts)
    allEpochs{m} = bpVec(epochEndPts(m,1):epochEndPts(m,2));
    allEpochs40hz{m} = bp40hzVec(epochEndPts(m,1):epochEndPts(m,2));
end

%% Determine frequency of each potential tremor event

proceed = 'n';
while proceed == 'n'
    % Parameters for finding peaks
    prompt = {'minPeakProm: ','minPeakInt: ','maxPeakWidth: ','minPeakWidth: '};
    dlgtitle = 'FindPeaks params';
    default = {'75','25','100','5'};
    tempParams = inputdlg(prompt,dlgtitle,1,default);
    parameters(1).minPeakProm = str2double(tempParams{1});
    parameters(1).minPeakDist =  str2double(tempParams{2});
    parameters(1).maxPeakWidth =  str2double(tempParams{3});
    parameters(1).minPeakWidth =  str2double(tempParams{4});
    
    % Find peaks of each tremor event
    for i = 1:length(allEpochs)
        if length(allEpochs{i}) >= 25
            [pks,idxs,w,p] = findpeaks(allEpochs{i},1,'MinPeakProminence',parameters(1).minPeakProm,'MinPeakDistance',parameters(1).minPeakDist,'Annotate','extents','MaxPeakWidth',parameters(1).maxPeakWidth,'MinPeakWidth',parameters(1).minPeakWidth);
            peaks{i} = pks; locs{i} = idxs; widths{i} = w; otherPeaks{i} = p;

            [pks40,idxs40,w40,p40] = findpeaks(allEpochs40hz{i},1,'MinPeakProminence',5);
            peaks40{i} = pks40; locs40{i} = idxs40;
            otherPeaks40{i} = p40;
        end
        % For trough detection: Flip the sign of each tremorEpoch and rerun
        % peak detection to get locs and peaks
    end
    
    % Plot individual tremor events
    tempFigNum = [14,15,16];
    endNum = length(allEpochs);
    if length(allEpochs) <= 40 && length(allEpochs) > 20
        tempFigs = 2; a = 2;
        counter = [1, 21, 20, length(allEpochs)];
    elseif length(allEpochs) <= 20
        tempFigs = 1; a = 1;
        counter = [1, length(allEpochs)];
    else
        tempFigs = 3; a = 3;
        counter = [1, 21, 41, 20, 40, 60];
    end

    for t = 1:tempFigs
        fig = figure(tempFigNum(t));
        for i = counter(t):counter(t+a)
            nexttile; hold on
            plot(allEpochs{i},'LineWidth',2,'Color','k')
            peakIdx = locs{i}; peakAmp = peaks{i};
            for j = 1:length(locs{i})
                scatter(peakIdx(j),peakAmp(j),100,'LineWidth',2)
            end
        end
    end

%     prompt = {'Keep peaks? (y/n): '};
%     dlgtitle = 'Peak analysis';
%     proceed = inputdlg(prompt,dlgtitle);
    prompt = 'Keep peaks? (y/n): ';
    proceed = input(prompt,'s');
    close(figure(14),figure(15),figure(16))

end

%%%%%%% INVERT tremorEpochs DATA AND TAKE PEAKS TO CAPTURE NEGATIVE PEAKS
%%%%%%% FOR AMPLITUDE MEASUREMENTS

%% Delinate movement from tremor events

% Find time difference between peaks for each epoch
a = 1; b = 1;
for j = 1:length(peaks)

    if length(peaks{j}) > 2

        for t = 1:length(locs{j}) - 1
            tremPeakDiff{j}(t) = locs{j}(t+1) - locs{j}(t);
        end

    end

end

% Calculate frequency of each epoch
tremFreq = cellfun(@(x) 1/(mean(x)/parameters(1).sf),tremPeakDiff,'UniformOutput',false);

% Remove epochs that do not match the tremor frequency range
b = 1; c = 1;
for i = 1:length(tremFreq)
    if tremFreq{i} > 23 || tremFreq{i} < 17
        tremFreq{i} = {};
    else
        tremFreqActual(b) = tremFreq{i};
        b = b + 1;
    end
end

for i = 1:length(tremFreq)
    if isempty(tremFreq{i})
    elseif isnan(tremFreq{i})
        tremFreq{i} = {};
    end
end

% Grab tremor and movement indexes
tremIdxs = find(~cellfun(@isempty,tremFreq));
allIdxs = 1:length(allEpochs);
moveIdxs = setdiff(allIdxs,tremIdxs);

epochIdxs{type,1} = epochEndPts(tremIdxs,:);
epochIdxs{type,2} = epochEndPts(moveIdxs,:);
tremor{type,10} = (epochIdxs{type,1}(:,2) - epochIdxs{type,1}(:,1));

% Separate out movements from tremor epochs
moveTemp = allEpochs;
moveTemp(tremIdxs) = [];
move{type} = moveTemp;
% move{type,2} = peaks(moveIdxs);

tremor{type,1} = allEpochs(tremIdxs);
tremor{type,2} = peaks(tremIdxs);
tremor{type,3} = locs(tremIdxs);
tremor{type,4} = widths(tremIdxs);
tremor{type,5} = otherPeaks(tremIdxs);
tremor{type,11} = allEpochs40hz(tremIdxs);
tremor{type,12} = peaks40(tremIdxs);
tremor{type,13} = otherPeaks40(tremIdxs);
tremor{type,14} = locs40(tremIdxs);

numPlots = length(tremor{type,1});
if numPlots <= 20
    numPlots = 1;
    numSegs = length(tremor{type,1});
elseif numPlots <= 40
    numPlots = 2;
    numSegs = [20,length(tremor{type,1}) - 20];
else
    numPlots = 3;
    numSegs = [20,20,length(tremor{type,1}) - 40];
end

counter = 1;
for t = 1:numPlots

    figure('Name','Tremor epochs')
    for i = counter:counter + numSegs(t) - 1
        nexttile; hold on
        plot(tremor{type,1}{i},'LineWidth',2,'Color','k')
        peakIdx = tremor{type,3}{i}; peakAmp = tremor{type,2}{i};
        for b = 1:length(tremor{type,3}{i})
            scatter(peakIdx(b),peakAmp(b),100,'LineWidth',2)
        end
    end
    xlabel('Timepoints')
    ylabel('\muVs')
    title('Tremor epochs')

    counter = counter + 20;
end

close all

%% Shift each epoch to be aligned on first peak

minLoc = 1000;
for i = 1:length(tremor{type,1})
    for k = 1:length(tremor{type,3}{i})
        tempMin(i) = min(tremor{type,3}{i});
        if tempMin(i) < minLoc
            minLoc = tempMin(i);
            minLocIdx = i;
        end
    end
    tremorDur(i) = length(tremor{type,2}{i}); 
end
tremor{type,6} = tremorDur;

for i = 1:length(tremor{type,1})
    if exist('minLocIdx')
        if i ~= minLocIdx
            shiftAmt = minLoc - tremor{type,3}{i}(1);
            tempShift = circshift(tremor{type,1}{i},shiftAmt);
            tremor{type,7}{i} = tempShift;
            tremor{type,8}(i) = shiftAmt;
        end
    end
end

tremor{type,9} = tremFreqActual;

if ~isempty(tremor{type,7})
    figure('Name','Peak-shifted tremor epochs'); hold on
    for i = 1:length(tremor{type,7})
        plot(tremor{type,7}{i})
    end
    xlabel('Timepts')
    ylabel('\muVs')
    title('Peak-aligned tremor epochs')
    set(gca,'FontSize',16)
end

close all

%% may need to take subset of non-tremor movements (e.g. 1-2 hz oscillations with multiple peaks)
