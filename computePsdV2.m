%
%
%
%   Written by Alex Fanning on 1/2/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataOut,parameters] = computePsdV2(signal,parameters)

psd1 = cell(1);
n = parameters.sf * parameters.timeframe;
freqStep = parameters.sf / n;
n1s = parameters.sf;
nCoh = n1s * parameters.cohSec;
freqCohRes = parameters.sf / nCoh;
shiftedframe = n1s * parameters.timeshift;

endNum = fix(length(signal{1})/n1s) - parameters.timeframe-ceil(parameters.timeshift);

freqLabel = NaN(nCoh/2+1,endNum);
timeLabel = NaN(nCoh/2+1,endNum);
sumWndw = [5 55];

% Frequency band cutoffs
bands = [1 3; 4 7; 8 12; 13 30];

if parameters.exptType(1)==1
    psd2{1} = NaN(nCoh/2+1,endNum);
    iterNum = length(parameters.chs)-1;
    parameters.numChs = length(parameters.chs);
else
    iterNum = length(parameters.chs);
    parameters.numChs = length(parameters.chs);
end

for i = 1:iterNum
    psd1{i} = NaN(nCoh/2+1,endNum);
    for s = 1:endNum
    
        y1 = signal{1}((1+(s-1)*n1s):(n+(s-1)*n1s),i);  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
        [PxxY1, freqPsY1] = pwelch(y1,hanning(nCoh),1/2*nCoh,nCoh,parameters.sf); % PSD by Wilch method and Hanning window
        psd1{i}(:,s) = PxxY1;

        if parameters.exptType(1)==1
            y2 = signal{2}((1+(s-1)*n1s+shiftedframe):(n+(s-1)*n1s+shiftedframe));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
            [PxxY2, freqPsY2] = pwelch(y2,hanning(nCoh),1/2*nCoh,nCoh,parameters.sf); % PSD by Wilch method and Hanning window
            psd2{1}(:,s) = PxxY2;

            [Cy1y2,FreqCoh] = mscohere(y1,y2,hanning(nCoh),1/2*nCoh,nCoh,parameters.sf); % (y1,y2 coherece, with hamming window, take NCoh points, half points overlap, number of fft:NCoh, sampling frequency:fs
            coh{i}(:,s) = Cy1y2;
        end
    
        freqLabel(:,s) = freqPsY1;
        timeLabel(1:(nCoh/2+1),s) = s-1;
    
    end

    figure(i)
    gca3 = pcolor(timeLabel,freqLabel,psd1{i});
    set(gca3, 'LineStyle','none');
    colorbar('location','eastoutside');
    yMax = max(prctile(psd1{i},98));
    clim([0 yMax]);
    title(['PSD ch ' num2str(parameters.chs(i))]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    axis([-inf,inf,0,50]);
    set(gca,'FontSize',16);
    saveas(i,[parameters.filename '_ch' num2str(parameters.chs(i)) '_timeFreqPsd.tif']);

    if parameters.exptType(1)==1
        figure(iterNum+1)
        gca3 = pcolor(timeLabel,freqLabel,psd2{1});
        set(gca3, 'LineStyle','none');
        colorbar('location','eastoutside');
        yMax = max(prctile(psd2{1},98));
        clim([0 yMax]);
        title('PSD FP');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        axis([-inf,inf,0,50]);
        set(gca,'FontSize',16);
        saveas(iterNum+1,[parameters.filename '_fp_timeFreqPsd.tif']);
    
        figure(i+iterNum+1)
        gca3 = pcolor(timeLabel,freqLabel,coh{i});
        set(gca3, 'LineStyle','none');
        colorbar('location','eastoutside');
        yMax = max(prctile(coh{i},98));
        clim([0 yMax]);
        title('Coherence PSD');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        axis([-inf,inf,0,50]);
        set(gca,'FontSize',16);
        saveas(i+iterNum+1,[parameters.filename  '_ch' num2str(parameters.chs(i)) '_coh_timeFreqPsd.tif']);
    end

    for s = 1:parameters.numEpochs

        mPsd{1,i}(:,s) = mean(psd1{i},2);
        sumPsd(1) = sum(mPsd{1,i}(sumWndw(1):sumWndw(2),s));
        normPsd{1,i}(:,s) = mPsd{1,i}(:,s) / sumPsd(1);
        [maxPsd{1,i}(s),maxPsdFreq{1,i}(s)] = max(mPsd{1,i}(parameters.extractFreq(1):parameters.extractFreq(2),s));
        [nMaxPsd{1,i}(s),nMaxPsdFreq{1,i}(s)] = max(normPsd{1,i}(parameters.extractFreq(1):parameters.extractFreq(2),s));
        maxFreq{1,i}(s) = (parameters.extractFreq(1) + maxPsdFreq{1,i}(s) - 1) * freqCohRes;
        nMaxFreq{1,i}(s) = (parameters.extractFreq(1) + nMaxPsdFreq{1,i}(s) - 1) * freqCohRes;
        pk40hz{1,i}(s) = max(mPsd{1,i}(35:45,s));

        nPsdAucAll{1,i}(:,s) = cumtrapz(normPsd{1,i}(1:35,s));
        for k = 2:length(nPsdAucAll{1,i}(:,s))
            nPsdAucTotal{1,i}(s,k-1) = nPsdAucAll{1,i}(k,s) - nPsdAucAll{1,i}(k-1,s);
        end
        nPsdAuc{1,i}(s) = trapz(normPsd{1,i}(parameters.extractFreq(1):parameters.extractFreq(2),s));
        nPsdAucPct{1,i}(s) = (nPsdAuc{1,i}(s) / nPsdAucAll{1,i}(end,s)) * 100;

        for m = 1:4
            nPsdBands{1,i}(s,m) = sum(normPsd{1,i}(bands(m,1):bands(m,2),s));
            nPsdBandAuc{1,i}(s,m) = sum(nPsdAucTotal{1,i}(s,bands(m,1):bands(m,2)));
            nPsdBandsAucPct{1,i}(s,m) = (nPsdBandAuc{1,i}(s,m) / nPsdAucAll{1,i}(end,s)) * 100;
        end

        if parameters.exptType(1)==1
            mPsd{2,1}(:,s) = mean(psd2{1},2);
            sumPsd(2) = sum(mPsd{2,1}(sumWndw(1):sumWndw(2),s));
            normPsd{2,1}(:,s) = mPsd{2,1}(:,s) / sumPsd(2);
            [maxPsd{2,1}(s),maxPsdFreq{2,1}(s)] = max(mPsd{2,1}(parameters.extractFreq(1):parameters.extractFreq(2),s));
            [nMaxPsd{2,1}(s),nMaxPsdFreq{2,i}(s)] = max(normPsd{2,1}(parameters.extractFreq(1):parameters.extractFreq(2),s));
            maxFreq{2,1}(s) = (parameters.extractFreq(1) + maxPsdFreq{2,1}(s) - 1) * freqCohRes;
            nMaxFreq{2,1}(s) = (parameters.extractFreq(1) + nMaxPsdFreq{2,1}(s) - 1) * freqCohRes;
            pk40hz{2,1}(s) = max(mPsd{2,1}(35:45,s));
    
            nPsdAucAll{2,1}(:,s) = cumtrapz(normPsd{2,1}(1:35,s));
            for k = 2:length(nPsdAucAll{2,1}(:,s))
                nPsdAucTotal{2,1}(s,k-1) = nPsdAucAll{2,1}(k,s) - nPsdAucAll{2,1}(k-1,s);
            end
            nPsdAuc{2,1}(s) = trapz(normPsd{2,1}(parameters.extractFreq(1):parameters.extractFreq(2),s));
            nPsdAucPct{2,1}(s) = (nPsdAuc{2,1}(s) / nPsdAucAll{2,1}(end,s)) * 100;
            
            for m = 1:4
                nPsdBands{2,1}(s,m) = sum(normPsd{2,1}(bands(m,1):bands(m,2),s));
                nPsdBandAuc{2,1}(s,m) = sum(nPsdAucTotal{2,1}(s,bands(m,1):bands(m,2)));
                nPsdBandsAucPct{2,1}(s,m) = (nPsdBandAuc{2,1}(s,m) / nPsdAucAll{2,1}(end,s)) * 100;
            end

            mCoh{i}(:,s) = mean(coh{i},2);
            [maxCoh{i}(s),maxCohFreq{i}(s)] = max(mCoh{i}(parameters.extractFreq(1):parameters.extractFreq(2),s));
            cohFreq{i}(s) = (parameters.extractFreq(1) + maxCohFreq{i}(s) - 1) * freqCohRes;
        end

    end

end

%% Plot PSDs

figure();
plotLmts = [0,40,0,0.5];
colr = {[copper(parameters.numEpochs)];[summer(parameters.numEpochs)];[autumn(parameters.numEpochs)]};
for i = 1:iterNum
    count(1) = 1; count(2) = parameters.numEpochs+1; count(3) = parameters.numEpochs*2+1;
    for s = 1:parameters.numEpochs
        subplot(parameters.numChs,parameters.numEpochs,count(1)); hold on
        plot(mPsd{1,i}(:,s),'LineWidth',1,'Color',colr{i}(s,:));
        count(1) = count(1) + 1;
    end
    title('Mean PSD');
    ylabel('dB/Hz');
    axis([plotLmts(1:3), maxPsd{1,i}*2]);
    set(gca,'FontSize',12)

    for s = 1:parameters.numEpochs
        subplot(parameters.numChs,parameters.numEpochs,count(2)); hold on
        plot(normPsd{1,i}(:,s),'LineWidth',1,'Color',colr{i}(s,:));
        count(2) = count(2) + 1;
    end
    if parameters.exptType(1)==1
        plot(normPsd{2,1}(:,s),'LineWidth',1,'Color',colr{3}(s,:));
    end
    title('Normalized PSD');
    xlabel('Frequency (Hz)');
    ylabel('Norm PSD');
    axis([plotLmts(1:3), nMaxPsd{1,i}*2]);
    set(gca,'FontSize',12)
    if parameters.exptType(1)==1
        legend(['Ch ' num2str(parameters.chs(1))],['Ch ' num2str(parameters.chs(2))],'FP','Location','best','FontSize',10)
    else
        legend(['Ch ' num2str(parameters.chs(1))],['Ch ' num2str(parameters.chs(2))],'Location','best','FontSize',10)
    end

    if parameters.exptType(1)==1
        for s = 1:parameters.numEpochs
            subplot(parameters.numChs,parameters.numEpochs,count(3)); hold on
            plot(mCoh{i}(:,s),'LineWidth',1,'Color',colr{i}(s,:));
            count(3) = count(3) + 1;
        end
        title('Coherence PSD');
        xlabel('Frequency (Hz)');
        ylabel('PSD');
        axis([plotLmts(1:3), maxCoh{i}*2]);
        set(gca,'FontSize',12)
        legend(['Ch ' num2str(parameters.chs(1)) ' v fp' ],['Ch ' num2str(parameters.chs(2)) ' v fp' ],'Location','best','FontSize',10)
    end
end

saveas(gcf,[parameters.filename '_normPsd.tif']);

%% Export data
if parameters.exptType(1)==1
    dataOut = {mPsd; normPsd; mCoh; maxPsd; maxFreq; nMaxPsd; nMaxFreq; nPsdAuc; nPsdAucPct; pk40hz; maxCoh; cohFreq; nPsdBands; nPsdBandAuc; nPsdBandsAucPct};
else
    dataOut = {mPsd; normPsd; maxPsd; maxFreq; nMaxPsd; nMaxFreq; nPsdAuc; nPsdAucPct; pk40hz; nPsdBands; nPsdBandAuc; nPsdBandsAucPct};
end