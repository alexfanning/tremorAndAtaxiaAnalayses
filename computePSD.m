%
%   Computes PSD
%
%   Written by Alex Fanning, 5/2/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataOut,parameters,count] = computePSD(dataIn,parameters,dataOut,timeWndw,sfNum,sigNum,expt,count)

parameters(1).timeframe = timeWndw;
parameters(1).cohSec = 1;
parameters(1).timeshift = 0;

n = parameters(sfNum).sf * parameters(1).timeframe;
freqStep = parameters(sfNum).sf / n; 
n1s = parameters(sfNum).sf;          
nCoh = n1s * parameters(1).cohSec; 
freqCohRes = parameters(sfNum).sf / nCoh;
shiftedframe = n1s * parameters(1).timeshift;

if sigNum ~= 5
    signal = dataIn{1,sigNum}(:,3);
    endNumber = fix(length(dataIn{sigNum}) / n1s) - parameters(1).timeframe - ceil(parameters(1).timeshift);
else
    signal = dataIn;
    endNumber = fix(length(dataIn) / n1s) - parameters(1).timeframe - ceil(parameters(1).timeshift);
end

%% Compute PSD

freqLabel = NaN(nCoh/2+1,endNumber);
timeLabel = NaN(nCoh/2+1,endNumber);
psd = NaN(nCoh/2+1,endNumber);

for s = 1:endNumber
    
    y1 = signal((1 + (s-1) * n1s):(n+(s-1)*n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point

    [PxxY, freqPsY] = pwelch(y1,hanning(nCoh),1/2*nCoh,nCoh,parameters(sfNum).sf); % PSD by Wilch method and Hanning window
    psd(:,s) = PxxY;

    freqLabel(:,s) = freqPsY;
    timeLabel(1:(nCoh/2+1),s) = s-1;

end

%% Normalize PSD

if expt == 1
    prompt = {'Lower bound: ','Upper bound: '};
    dlgtitle = 'Extraction window';
    default = {'15','25'};
    tempParams = inputdlg(prompt,dlgtitle,1,default);
    parameters.extractWndw = [str2double(tempParams{1}) str2double(tempParams{2})];
elseif expt == 2
    if count == 1
        prompt = {'Lower bound: ','Upper bound: '};
        dlgtitle = 'Extraction window';
        default = {'8','13'};
        tempParams = inputdlg(prompt,dlgtitle,1,default);
    else
        tempParams = {'8','13'};
    end
    parameters.extractWndw = [str2double(tempParams{1}) str2double(tempParams{2})];
elseif expt == 3
%     if count == 1
%         prompt = {'Lower bound: ','Upper bound: '};
%         dlgtitle = 'Extraction window';
%         default = {'1','9'};
%         tempParams = inputdlg(prompt,dlgtitle,1,default);
%     else
%         tempParams = {'1','9'};
%     end
    tempParams = {parameters.bpRange(1),parameters.bpRange(2)};
    parameters.extractWndw = [tempParams{1} tempParams{2}];
end
if expt == 3
    sumWndw = [1 35];
else
    sumWndw = [37 43];
end

if expt == 1 || expt == 2
    dataOut.meanPSD{count(1)}(:,count(2)) = mean(psd,2);
    sumPSD = sum(dataOut.meanPSD{count(1)}(sumWndw(1):sumWndw(2),count(2)));
    dataOut.normPSD{count(1)}(:,count(2)) = dataOut.meanPSD{count(1)}(:,count(2)) / sumPSD;
elseif expt == 3
    dataOut.mPSD{count(1)}(:,count(2)) = nanmean(psd,2);
    sumPSD = sum(dataOut.mPSD{count(1)}(sumWndw(1):sumWndw(2),count(2)));
    dataOut.nMeanPsd{count(1)}(:,count(2)) = dataOut.mPSD{count(1)}(:,count(2)) / sumPSD;

    chnkTime = 1;
    for i = 1:parameters.numEpochs
        try
            chnkTime(i+1) = floor(parameters.chunkSize{count(1)}(i) + dataOut.idxs{count(1),count(2)}(1)/parameters.sf);
            dataOut.mPsdChnk{count(1),count(2)}(:,i) = mean(psd(:,chnkTime(i)+1:chnkTime(i+1)),2);
        catch
            chnkTime(i+1) = endNumber;
            dataOut.mPsdChnk{count(1),count(2)}(:,i) = mean(psd(:,chnkTime(i)+1:chnkTime(i+1)),2);
        end

        sumPsdChnk = sum(dataOut.mPsdChnk{count(1),count(2)}(sumWndw(1):sumWndw(2),i));
        dataOut.nPsdChnk{count(1),count(2)}(:,i) = dataOut.mPsdChnk{count(1),count(2)}(:,i) / sumPsdChnk;
    end
elseif expt == 4
    meanPSD = mean(psd,2);
    sumPSD = sum(meanPSD(sumWndw(1):sumWndw(2)));
    normPSD = meanPSD / sumPSD;
end

% CREATE EPOCHS OF MEAN, SUM, AND NORM
% 

if expt == 3

    % Frequency band cutoffs
    bands = [1 3; 4 7; 8 12; 13 30];

    % Find maximum value
    [dataOut.maxPSD{count(1)}(count(2)),dataOut.maxPsdFreq{:,count(1)}(count(2))] = max(dataOut.mPSD{count(1)}(parameters.extractWndw(1):parameters.extractWndw(2),count(2)));
    dataOut.maxPsdFreq{:,count(1)}(count(2)) = dataOut.maxPsdFreq{:,count(1)}(count(2)) + parameters.extractWndw(1) - 1;
    [dataOut.nMaxPsd{count(1)}(count(2)), dataOut.nMaxPsdFreq{count(1)}(count(2))] = max(dataOut.nMeanPsd{count(1)}(parameters.extractWndw(1):parameters.extractWndw(2),count(2)));
    dataOut.nMaxPsdFreq{count(1)}(count(2)) = dataOut.nMaxPsdFreq{count(1)}(count(2)) + parameters.extractWndw(1) - 1;

    freqMaxPSD = (parameters.extractWndw(1) + dataOut.maxPsdFreq{:,count(1)}(count(2)) - 1);
    % nmaxPSDbaseFreq = dataOut.mPSD(parameters.extractWndw(1)+ dataOut.maxPSDidx{count(1)}(count(2),(1)) - 1) / sumPSD;

    % Find max PSD and frequency and AUC for each frequency band
    dataOut.nPsdAucAll{count(1)}(count(2),:) = cumtrapz(dataOut.nMeanPsd{count(1)}(1:35,count(2)));
    for i = 2:length(dataOut.nPsdAucAll{count(1)}(count(2),:))
        dataOut.nPsdAucTotal{count(1)}(count(2),i-1) = dataOut.nPsdAucAll{count(1)}(count(2),i) - dataOut.nPsdAucAll{count(1)}(count(2),i-1);
    end
    dataOut.nPsdAuc{count(1)}(:,count(2)) = trapz(dataOut.nMeanPsd{count(1)}(parameters.extractWndw(1):parameters.extractWndw(2),count(2)));
    dataOut.nPsdAucPct{:,count(1)}(count(2)) = (dataOut.nPsdAuc{count(1)}(:,count(2)) / dataOut.nPsdAucAll{count(1)}(count(2),end)) * 100;

    for t = 1:4
        dataOut.nPsdBands{count(1)}(count(2),t) = sum(dataOut.nMeanPsd{count(1)}(bands(t,1):bands(t,2),count(2)));
        dataOut.nPsdBandAuc{count(1)}(count(2),t) = sum(dataOut.nPsdAucTotal{count(1)}(count(2),bands(t,1):bands(t,2)));
        dataOut.nPsdBandsAucPct{count(1)}(count(2),t) = (dataOut.nPsdBandAuc{count(1)}(count(2),t) / dataOut.nPsdAucAll{count(1)}(count(2),end)) * 100;
    end

    % Compute maximum PSD and frequency for each epoch
    for i = 1:parameters.numEpochs
        [dataOut.nPsdChnkMax{count(1)}(count(2),i), dataOut.nPsdChnkMaxFreq{count(1)}(count(2),i)] = max(dataOut.nPsdChnk{count(1),count(2)}(parameters.extractWndw(1):parameters.extractWndw(2),i));
        dataOut.nPsdChnkMaxFreq{count(1)}(count(2),i) = dataOut.nPsdChnkMaxFreq{count(1)}(count(2),i) + parameters.extractWndw(1) - 1;
    end

elseif expt == 1 || expt == 2
    [dataOut.maxPSD{count(1)}(:,count(2)),dataOut.maxPSDidx{count(1)}(:,count(2))] = max(dataOut.meanPSD{count(1)}(parameters.extractWndw(1):parameters.extractWndw(2),count(2)));
    dataOut.maxPSDidx{count(1)}(:,count(2)) = dataOut.maxPSDidx{count(1)}(:,count(2)) + parameters(1).extractWndw(1) - 1;
    [dataOut.normPSDmax{count(1)}(:,count(2)), dataOut.normPSDmaxIdx{count(1)}(:,count(2))] = max(dataOut.normPSD{count(1)}(parameters.extractWndw(1):parameters.extractWndw(2),count(2)));
    dataOut.normPSDmaxIdx{count(1)}(:,count(2)) = dataOut.normPSDmaxIdx{count(1)}(:,count(2)) + parameters.extractWndw(1) - 1;

    dataOut.freqMaxPSD{count(1)}(:,count(2)) = (parameters.extractWndw(2) + dataOut.maxPSDidx{count(1)}(:,count(2)) - 2) * freqCohRes;
    dataOut.nmaxPSDbaseFreq{count(1)}(:,count(2)) = dataOut.meanPSD{count(1)}(parameters.extractWndw(2) + dataOut.maxPSDidx{count(1)}(count(2)) - 1) / sumPSD;
end

%% Plot PSDs

if expt == 1
    figure();
    gca3 = pcolor(timeLabel,freqLabel,psd);
    set(gca3, 'LineStyle','none');
    colorbar('location','eastoutside');
    clim([0 1000]);
    title('PSD');
    xlabel('Time (s)'); 
    ylabel('Frequency (Hz)'); 
    axis([-inf,inf,0,50]);
    set(gca,'FontSize',16);
    saveas(gcf,[parameters(1).recName 'timeFreqPsd.tif']);
elseif expt == 4
    subNames = {'Cb (left)','Cb (right)','M1 (left)','M1 (right)','S1 (left)','S1 (right)'};
    figure(6)
    subplot(3,2,count(1))
    gca3 = pcolor(timeLabel,freqLabel,psd);
    set(gca3, 'LineStyle','none');
    colorbar('location','eastoutside');
    clim([0 1000]);
    title('PSD');
    if count(1) == 3
        ylabel('Frequency (Hz)');
    elseif count(1) == 5 || count(1) == 6
        xlabel('Time (s)');
    end
    title(subNames{count(1)})
    axis([-inf,inf,0,40]);
    if count(1) == 2
        sgtitle('PSD')
    end
    set(gca,'FontSize',16);
    if count(1) == 6
        saveas(gcf,[parameters(1).recName '_timeFreqPsd.tif']);
    end
elseif expt == 3
    titleNames = {'Action','Posture','Spiral','Tapping','Rest'};
    figure();
    gca3 = pcolor(timeLabel,freqLabel,psd);
    set(gca3, 'LineStyle','none');
    h = colorbar('location','eastoutside');
    psdLims = prctile(psd,99,'all');
    clim([0 psdLims]);
    xlabel('Time (s)'); 
    ylabel('Frequency (Hz)');
    ylabel(h,'Intensity (dB/Hz)','FontSize',16)
    axis([-inf,inf,0,50]);
    title([titleNames{count(1)} ': Mean PSD']);
    set(gca,'FontSize',16);
    saveas(gcf,[parameters(1).fileName '_timeFreqPsd.fig']);
end

plotLmts = [0,40,0,0.5];
if expt == 1
    figure(7);
    subplot(2,1,1); hold on
    h = plot(dataOut.meanPSD{count(1)}{count(2),count(3)}(count(4),:),'LineWidth',1);
    title('Mean PSD');
    ylabel('PSD (dB/Hz)');
    axis([plotLmts(1:3), dataOut.maxPSD{count(1)}{count(2),count(3)}(count(4),:)*2]);
    set(h(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)
    set(gca,'FontSize',16)
    
    subplot(2,1,2); hold on
    h = plot(dataOut.normPSD{count(1)}{count(2),count(3)}(count(4),:),'LineWidth',1);
    sgtitle('Mean PSD (normalized)');
    xlabel('Frequency (Hz)'); 
    ylabel('Norm PSD (dB/Hz)');
    axis([plotLmts(1:3), dataOut.normPSDmax{count(1)}{count(2),count(3)}(count(4),:)*2]);
    set(h(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)
    set(gca,'FontSize',16)
elseif expt == 2
    subNames = {'Cb (left)','Cb (right)','M1 (left)','M1 (right)','S1 (left)','S1 (right)'};
    col = cool(3);
    figure(7);
    subplot(2,3,count(1)); hold on
    h = plot(dataOut.meanPSD{count(1)}(:,count(2)),'LineWidth',1);
    if count(1) == 2
        sgtitle('Mean PSD');
    end
    title(subNames{count(1)})
    ylabel('PSD (dB/Hz)');
    axis([plotLmts(1:3), max(max(dataOut.maxPSD{count(1)}))*1.2]);
    set(h,'Color',col(count(2),:)); % Set the line 1 by RGB color [0 0 1] (blue)
    set(gca,'FontSize',16)
    if count(1) == 6
        legend('Pre','Stim','Post')
    end

    figure(8);
    subplot(2,3,count(1)); hold on
    h = plot(dataOut.normPSD{count(1)}(:,count(2)),'LineWidth',1,'Color',col(count(2),:));
    if count(1) ==  2
        sgtitle('Mean PSD (normalized)');
    end
    title(subNames{count(1)})
    xlabel('Frequency (Hz)');
    ylabel('Norm PSD (dB/Hz)');
    axis([plotLmts(1:3), dataOut.normPSDmax{count(1)}(:,count(2))*3]);
    set(h,'Color',col(count(2),:)); % Set the line 1 by RGB color [0 0 1] (blue)
    if count(1) == 6
        legend('Pre','Stim','Post')
    end
    set(gca,'FontSize',16)
elseif expt == 3
    plotLmts = [0,30,0,0.5];

    figure();
    subplot(2,1,1); hold on
    title([titleNames{count(1)} ': Mean PSD']);
    h = plot(dataOut.mPSD{1,count(1)}(:,count(2)),'LineWidth',1);
    ylabel('PSD (dB/Hz)');
    axis([plotLmts(1:3), dataOut.maxPSD{count(1)}(:,count(2))*2]);
    set(h(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)
    set(gca,'FontSize',16)

    subplot(2,1,2); hold on
    h = plot(dataOut.nMeanPsd{count(1)}(:,count(2)),'LineWidth',1);
    title('Mean PSD (normalized)');
    xlabel('Frequency (Hz)');
    ylabel('Norm PSD (dB/Hz)');
    axis([plotLmts(1:3), dataOut.nMaxPsd{1,count(1)}(:,count(2))*2]);
    set(h(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)
    set(gca,'FontSize',16)
    
end
saveas(gcf,[parameters(1).fileName '_meanPsd.fig'])

