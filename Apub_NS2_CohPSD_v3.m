clear; close all;

%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%
params = struct();

params.chSig1 = 142;
params.chRef1 = 0;
params.chSig2 = 144;
params.chRef2 = 0;

params.BW = 0; 
params.timecut = 0; 
params.timeshift = 0; 
params.trimStartTime = 0; 
params.trimStopTime = 390; 

params.sf = 1000; 
params.timeframe = 20; 
params.cohSec = 1;

%% Extract parameters from excel sheet

params.xlData = xlsread('Input parameters3.xlsx','B2:E2');
prompt = {'Force plate channel # (141 or 144): '};
dlgtitle = 'Parameters';
defaults = {'141'};
tempParams = inputdlg(prompt,dlgtitle,1,defaults);
params.chSig1 = str2double(tempParams{2});
params.chSig2 = str2double(tempParams{1});
params.chRef1 = params.xlData(2) ;
params.chRef2 = params.xlData(4);
clear prompt dlgtitle defaults tempParams

%%%%%%%%  Load NS2 data %%%%%%%%%%%%%%%%%%%%%%

openNSx;
params.sourceFile = [NS2.MetaTags.Filename NS2.MetaTags.FileExt];
data = double(NS2.Data');
NS2 = rmfield(NS2,'Data');

%%%%%%%% Re-assign channelID to legacy variables %%%%%%%%%%
params.channelSig1 = find(NS2.MetaTags.ChannelID==params.chSig1);
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
Title1 = [sprintf('%d-%d vs %d-%d ',params.chSig1,params.chRef1,params.chSig2,params.chRef2), params.sourceFile];

if params.timecut ~= 0
    data = data(params.trimStartTime*params.sf+1:params.trimStopTime*params.sf,:);
end

if params.channelRef1==0
    Signal1 = data(:,params.channelSig1);
else
    Signal1 = data(:,params.channelSig1) - data(:,params.channelRef1);
end

if params.channelRef2==0
    Signal2 = data(:,params.channelSig2);
else
    Signal2 = data(:,params.channelSig2)-data(:,params.channelRef2);
end    

%% Determine coherency between LFP and force plate data

N = params.sf*params.timeframe; 
freqStep = params.sf/N; 
N1s = params.sf;          
NCoh = N1s*params.cohSec; 
FreqCohRes = params.sf/NCoh;
EndNumber = fix(length(data)/N1s) - params.timeframe-ceil(params.timeshift);
Shiftedframe = N1s*params.timeshift;

for s = 1:EndNumber
    
    y1 = Signal1((1+(s-1)*N1s):(N+(s-1)*N1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
    y2 = Signal2((1+(s-1)*N1s+Shiftedframe):(N+(s-1)*N1s+Shiftedframe));
    
    [Cy1y2,FreqCoh] = mscohere(y1,y2,hanning(NCoh),1/2*NCoh,NCoh,params.sf); % (y1,y2 coherece, with hamming window, take NCoh points, half points overlap, number of fft:NCoh, sampling frequency:params.sf
    Coh(:,s) = Cy1y2; 
    FreqLable(:,s) = FreqCoh;
    TimeLable(1:(NCoh/2+1),s) = s-1;
    
    [PxxY1, freqPsY1] = pwelch(y1,hanning(NCoh),1/2*NCoh,NCoh,params.sf); % PSD by Wilch method and Hanning window
    [PxxY2, freqPsY2] = pwelch(y2,hanning(NCoh),1/2*NCoh,NCoh,params.sf);
    PSDY1(:,s) = PxxY1;
    PSDY2(:,s) = PxxY2;

end

filename = [a_InputParameters.params.sourceFile '.mat'];

%% Extract tremor episodes

% rawData = data(:,1);
% params = struct();
% [tremEpochs] = extractTremor(rawData,params,params.sourceFile);

%% 

%save(filename)
%save Coherence_PSD_Profiles&Results.mat Coh PSDY1 PSDY2 FreqCoh a_InputParameters NS2 TimeLable FreqLable data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Drawing graphs of Time-Frequency Plot of Coherence
% Drawing graphs of Time-Frequency Plot of Power Spectrum Density

%load Coherence_PSD_Profiles&Results.mat TimeLable FreqLable Coh;
%Coherence: full frequency
figure(1);
gca=pcolor(TimeLable,FreqLable,Coh);
set(gca, 'LineStyle','none');
colorbar('location','eastoutside');
title(['Cohere ',Title1]); 
xlabel('Time (seconds)'); 
ylabel('Frequency (Hz)'); 
saveas(1,'Coherence_0to250Hz.fig');

figure(2);
gca=pcolor(TimeLable,FreqLable,Coh);
set(gca, 'LineStyle','none');
colorbar('location','eastoutside');
title(['Cohere ',Title1]); 
xlabel('Time (seconds)'); 
ylabel('Frequency (Hz)'); 
axis([-inf,inf,0,100]);
saveas(2,'Coherence_0to100Hz.fig');

figure(3);
gca=pcolor(TimeLable,FreqLable,Coh);
set(gca, 'LineStyle','none');
colorbar('location','eastoutside');
title(['Cohere ',Title1]); 
xlabel('Time (seconds)'); 
ylabel('Frequency (Hz)'); 
axis([-inf,inf,0,50]);
saveas(3,'Coherence_0to50Hz.fig');

figure(6);
gca3=pcolor(TimeLable,FreqLable,PSDY1);
set(gca3, 'LineStyle','none');
colorbar('location','eastoutside');
caxis([0 1000]);
title(['PSD Sig1 ',Title1]); 
xlabel('Time (seconds)'); 
ylabel('Frequency (Hz)'); 
axis([-inf,inf,0,100]);
saveas(6,'PSD_Singal1_0to100Hz.fig');

figure(8);
gca3=pcolor(TimeLable,FreqLable,PSDY2);
set(gca3, 'LineStyle','none');
colorbar('location','eastoutside');
caxis([0 1000]);
title(['PSD Sig2 ',Title1]); 
xlabel('Time (seconds)'); 
ylabel('Frequency (Hz)');
axis([-inf,inf,0,100]);
saveas(8,'PSD_Signal2_0to100Hz.fig');
close all;

figure(9);
subplot(3,1,1);
gca3=pcolor(TimeLable,FreqLable,PSDY1);
set(gca3, 'LineStyle','none');
colorbar('location','eastoutside');
caxis([0 1000]);
title([Title1]); 
ylabel('PSDY1(Hz)'); 
axis([-inf,inf,0,50]);

subplot(3,1,2);
gca3=pcolor(TimeLable,FreqLable,PSDY2);
set(gca3, 'LineStyle','none');
colorbar('location','eastoutside');
caxis([0 1000]);
ylabel('PSDY2(Hz)');
axis([-inf,inf,0,50]);

subplot(3,1,3);
gca=pcolor(TimeLable,FreqLable,Coh);
set(gca, 'LineStyle','none');
colorbar('location','eastoutside');
xlabel('Time (seconds)'); 
ylabel('Coh(Hz)'); 
axis([-inf,inf,0,50]);
saveas(9,'TFPLot_PSDY1Y2+Coh.fig');
saveas(9,'TFPLot_PSDY1Y2+Coh.jpg');

%% Plot raw LFP and behavioral data
% 
% figure('Name','Raw LFP and tremor')
% ax1 = subplot(3,1,1); hold on
% plot(data(:,1))
% ylabel('\muV','FontSize',18)
% title(ax1,'Thalamic LFP','FontSize',18)
% 
% ax2 = subplot(3,1,2); hold on
% plot(data(:,2))
% ylabel('\muV','FontSize',18)
% title(ax2,'Cerebellar LFP','FontSize',18)
% 
% ax3 = subplot(3,1,3); hold on
% plot(data(:,4))
% xlabel('Time (ms)','FontSize',18)
% ylabel('mV','FontSize',18)
% title(ax3,'Force plate','FontSize',18)
% linkaxes([ax1 ax2 ax3],'x')
% saveas(figure,'Raw.fig')
