
%clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters
CallExcel = 1;

ExtractFreqRange = [15,25]; % Define the frequency range to extract the data
ExtractWindows = [0,100;0,500]; % Define the extract windows (can be multiple)

nSig1PlotLimit = [0,100,0,0.1];
nSig2PlotLimit = [0,100,0,0.1];

% Define the frequency range for normalizing data
SumWindow = [5,55];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load excel parameters
if CallExcel==1
    A = xlsread('Input parameters3.xlsx','G2:O2');
    ExtractFreqRange=A(8:9); % Define the frequency range to extract the data
    ExtractWindows =[A(1:2);A(3:4);A(5:6)]; % Define the extract windows (can be multiple)
end

save temp.mat;
a_ExtractPara=load('temp.mat');
delete temp.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temp = dir('*.ns2.mat');
% load(temp.name);
%load Coherence_PSD_Profiles&Results.mat PSDY1 PSDY2 Coh FreqCoh a_InputParameters;
Title1=[sprintf('%d-%d vs %d-%d ',a_InputParameters.ChSig1,a_InputParameters.ChRef1,a_InputParameters.ChSig2,a_InputParameters.ChRef2),a_InputParameters.SourceFile];
FreqCohRes=1/a_InputParameters.CohSec; % Defined Frequency resolution

%%%%%%%%%%%%%%%%%%%%% Consistency check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FreqCohRes~=FreqCoh(2)
    errordlg('frequency resolution mismatch');
    return
end

ExtractWindows(:,3) = ExtractWindows(:,2)-a_InputParameters.timeframe; % Define the Endtime by subtract the timeframe (for internal consistency of previous version)
[WindowNumber,b] = size(ExtractWindows); % Define the extraction size
clear b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean PSD and normalization
SumWindowIndex=SumWindow/FreqCohRes+1; % calculating the normalizing value from the SumWindow)

if WindowNumber>=2
    for i=1:3
        meanPSDY1(:,i)=mean(PSDY1(:,ExtractWindows(i,1)+1:ExtractWindows(i,3)),2);
        meanPSDY2(:,i)=mean(PSDY2(:,ExtractWindows(i,1)+1:ExtractWindows(i,3)),2);
        meanCoh(:,i)=mean(Coh(:,ExtractWindows(i,1)+1:ExtractWindows(i,3)),2);

        sumPSDY1(i)=sum(meanPSDY1(SumWindowIndex(1):SumWindowIndex(2),i));
        sumPSDY2(i)=sum(meanPSDY2(SumWindowIndex(1):SumWindowIndex(2),i));

        nmeanPSDY1(:,i)=meanPSDY1(:,i)/sumPSDY1(i);
        nmeanPSDY2(:,i)=meanPSDY2(:,i)/sumPSDY2(i);
    end
else

    meanPSDY1(:,1)=mean(PSDY1(:,ExtractWindows(1,1)+1:ExtractWindows(1,3)),2);
    meanPSDY2(:,1)=mean(PSDY2(:,ExtractWindows(1,1)+1:ExtractWindows(1,3)),2);
    meanCoh(:,1)=mean(Coh(:,ExtractWindows(1,1)+1:ExtractWindows(1,3)),2);
    
    sumPSDY1=sum(meanPSDY1(SumWindowIndex(1):SumWindowIndex(2),1));
    sumPSDY2=sum(meanPSDY2(SumWindowIndex(1):SumWindowIndex(2),1));
    %sumCoh=sum(meanCoh(SumWindowIndex(1):SumWindowIndex(2),1));

    nmeanPSDY1=meanPSDY1/sumPSDY1;
    nmeanPSDY2=meanPSDY2/sumPSDY2;
end

%% Whole record average (Alex Fanning 10/2/22)
recAvgPSD1 = mean(PSDY1,2);
recAvgPSD2 = mean(PSDY2,2);
sumRecAvgPSD1 = sum(recAvgPSD1(SumWindowIndex(1):SumWindowIndex(2),1));
sumRecAvgPSD2 = sum(recAvgPSD2(SumWindowIndex(1):SumWindowIndex(2),1));
normRecAvgPSD1 = recAvgPSD1/sumRecAvgPSD1;
normRecAvgPSD2 = recAvgPSD2/sumRecAvgPSD2;

[maxPSDephys,indexMaxPSDephys] = max(recAvgPSD1(ExtractFreqRange(1):ExtractFreqRange(2),:));
normMaxPSDephys = maxPSDephys/sumRecAvgPSD1;
[maxPSDbeh,indexMaxPSDbeh] = max(recAvgPSD2(ExtractFreqRange(1):ExtractFreqRange(2),:));
normMaxPSDbeh = maxPSDbeh/sumRecAvgPSD2;
freqMaxEphys = (ExtractFreqRange(1)+indexMaxPSDephys-2) * FreqCohRes;
freqMaxRecBeh = (ExtractFreqRange(1)+indexMaxPSDbeh-2) * FreqCohRes;

%% 
%%%%%%%%%%%%%%%%%% Find max and freqeuncy %%%%%%%%%%%%%%%%%%%%%%%%
ExtractFreqIndex = ExtractFreqRange/FreqCohRes+1;
[maxPSDY1,indexPSDY1] = max(meanPSDY1(ExtractFreqIndex(1):ExtractFreqIndex(2),:));
nmaxPSDY1 = maxPSDY1/sumPSDY1;
[maxPSDY2,indexPSDY2] = max(meanPSDY2(ExtractFreqIndex(1):ExtractFreqIndex(2),:));
nmaxPSDY2 = maxPSDY2/sumPSDY2;
[maxCoh,indexCoh] = max(meanCoh(ExtractFreqIndex(1):ExtractFreqIndex(2),:));
FreqmaxPSDY1 = (ExtractFreqIndex(1)+indexPSDY1-2)*FreqCohRes;
FreqmaxPSDY2 = (ExtractFreqIndex(1)+indexPSDY2-2)*FreqCohRes;
FreqmaxCoh = (ExtractFreqIndex(1)+indexCoh-2)*FreqCohRes;

%%%%%%%%%%%%%%%%%% Find max according to baseline freqeuncy %%%%%%%%%%%%%%%%%%%%%%%%
nmaxPSDY1BaseFreq=meanPSDY1(ExtractFreqIndex(1)+indexPSDY1(1)-1,:)/sumPSDY1;
nmaxPSDY2BaseFreq=meanPSDY2(ExtractFreqIndex(1)+indexPSDY2(1)-1,:)/sumPSDY2;
maxCohBaseFreq=meanCoh(ExtractFreqIndex(1)+indexCoh(1)-1,:);


a_ExtractSummary=cell(14,WindowNumber+1);
a_ExtractSummary(:,1)={'Window','nPSDY1','nPSDY2','Coh','Y1Frq','Y2Freq','CohFreq','','nPSDY1_bfreq','nPSDY2_bfreq','Coh_bfreq','Y1Frq','Y2Freq','CohFreq'};
for i=WindowNumber:-1:1
    a_ExtractSummary(:,i+1)={sprintf('%d-%d',ExtractWindows(i,1), ExtractWindows(i,2)),nmaxPSDY1(i), nmaxPSDY2(i), maxCoh(i),FreqmaxPSDY1(i),FreqmaxPSDY2(i),FreqmaxCoh(i),'',nmaxPSDY1BaseFreq(i), nmaxPSDY2BaseFreq(i), maxCohBaseFreq(i),FreqmaxPSDY1(1),FreqmaxPSDY2(1),FreqmaxCoh(1)};
end

%clear Coh PSDY1 PSDY2 i;
save(params.filename)
%save Mean_PSD_new2.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data
figure(1);
subplot(3,1,1);
h11=plot(FreqCoh ,meanPSDY1);
title(['PSD (Sig1) ' Title1]);
xlabel('Frequency (Hz)'); 
ylabel('PSD (dB/Hz)');
xlim([0 30])
%axis([nSig1PlotLimit(1:3),max(maxPSDY1)*1.2]);
legend('1-200 s', '201-400 s','401-600 s')
%legend(a_ExtractSummary(1,2:WindowNumber+1));
set(h11(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)

subplot(3,1,2);
plot(FreqCoh ,meanPSDY2);
title('PSD (Sig2)'); 
xlabel('Frequency (Hz)'); 
ylabel('PSD (dB/Hz)');
axis([nSig2PlotLimit(1:3),max(maxPSDY2)*1.2]);

subplot(3,1,3);
plot(FreqCoh ,meanCoh);
title('PSD (Coh)'); 
xlabel('Frequency (Hz)'); 
ylabel('Coh');
xlim([nSig2PlotLimit(1),nSig2PlotLimit(2)]);
saveas(1,'PSD2.fig');
saveas(1,'PSD2.jpg');


figure(2);
subplot(3,1,1); hold on
h21=plot(FreqCoh ,nmeanPSDY1);
title(['Norm. PSD (Sig1) ' Title1]);
xlabel('Frequency (Hz)'); 
ylabel('Norm. PSD (dB/Hz)');
set(h21(1),'Color',[0,0,1]); % Set the line 1 by RGB color [0 0 1] (blue)
%plot(FreqCoh,normRecAvgPSD2,'LineWidth',2)
%axis(nSig1PlotLimit);
xlim([0 30])
legend(a_ExtractSummary(1,2:WindowNumber+1));

subplot(3,1,2);
plot(FreqCoh ,nmeanPSDY2);
title('Norm. PSD (Sig2)'); 
xlabel('Frequency (Hz)'); 
ylabel('Norm. PSD (dB/Hz)');
axis(nSig2PlotLimit);

subplot(3,1,3);
plot(FreqCoh ,meanCoh);
title('Norm. PSD (Coh)'); 
xlabel('Frequency (Hz)'); 
ylabel('Coh');
xlim([nSig2PlotLimit(1),nSig2PlotLimit(2)]);

saveas(2,'nPSD2.fig');
saveas(2,'nPSD2.jpg');


% xlswrite('PSDresults.xls',a_ExtractSummary,'Sheet1','A4');
% xlswrite('PSDresults.xls',{'FileName',a_InputParameters.SourceFile},'Sheet1','A1');
% xlswrite('PSDresults.xls',{'Title',Title1},'Sheet1','A2');
% xlswrite('PSDresults.xls',meanPSDY1,'meanPSDY1','A2');
% xlswrite('PSDresults.xls',meanPSDY2,'meanPSDY2','A2');
% xlswrite('PSDresults.xls',meanCoh,'meanCoh','A2');


