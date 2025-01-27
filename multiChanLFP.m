%   Multi-channel LFP recording analysis
%
%   Written by Alex Fanning, November 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all

% Set channel numbers to analyze
params = {}; data2export = struct();
data2export.maxCoh = cell(1); data2export.maxCohIdx = cell(1); data2export.coh = cell(1);
params.prompt = {'Cb Ch: ','Cb Ch:  ','M1 ch: ','M1 ch:','S1 ch: ','S1 ch: '};
params.dlgtitle = 'Recording channels';
params.default = {'30','29','1','2','3','4'};
params.temp = inputdlg(params.prompt,params.dlgtitle,1,params.default);
params.chSig = str2double(params.temp);

[list,~] = getFolders(1);
recTable = table2array(readtable(list(2).name));
test = detrend(recTable(:,3));
vel(:,1) = cumtrapz(recTable(5:end,4),1);

params.sf = recTable(2,2);
params.numBlocks = 3;
params.recName = list(1).name;
data = cell(1);

params.blkLngth = floor((length(recTable)-3) / params.numBlocks);
params.epochLngth{1} = 1:params.blkLngth;
params.epochLngth{2} = params.blkLngth+1:params.blkLngth*2;
params.epochLngth{3} = params.blkLngth*2+1:params.blkLngth*3;
for i = 1:length(params.chSig)
    params.chIdx = find(params.chSig(i) == recTable(1,:));
    data{1,i} = recTable(4:end,params.chIdx);
    for m = 1:params.numBlocks
        data{2,i}(:,m) = data{1,i}(params.epochLngth{m});
        [data2export,params,params.count] = computePSD(data{2,i}(:,m),params,data2export,20,1,5,2,[i m]);
    end
    [~,params] = computePSD(data{1,i},params,data2export,20,1,5,4,[i m]);
end

for i = 1:length(params.chSig)
    for m = 1:params.numBlocks
        [data2export] = coherence(data,data2export,params,20,1,1);
        [data2export] = coherence(data,data2export,params,20,1,2);
    end
end

subNames{1} = {'Cb (left) vs. M1 (right)' 'Cb (left) vs. S1 (right)','Cb (left) vs. Cb (right)','M1 (right) vs. S1(right)'};
subNames{2} = {'Cb (right) vs. M1 (left)' 'Cb (right) vs. S1 (left)','Cb (right) vs. Cb (left)','M1 (left) vs. S1(left)'};
col = cool(3);
params.count = 1;
figure('Name','Coherence')
for m = 1:size(data2export.mCoh,1)
    for t = 1:size(data2export.mCoh,2)
        h(t) = subplot(size(data2export.mCoh,1),size(data2export.mCoh,2),params.count); hold on
        for i = 1:params.numBlocks
            params.yMax = max(max(data2export.mCoh{m,t}));
            plot(data2export.mCoh{m,t}(:,i),'LineWidth',1,'Color',col(i,:))
        end
        if t == 2
            ylabel('Coherence')
        elseif t == 3
            xlabel('Freq')
        end
        xlim([0 30])
        ylim([0 params.yMax])
        set(gca,'FontSize',16)
        title(subNames{m}{t})
        params.count = params.count + 1;
    end
    linkaxes([h(1) h(2) h(3) h(4)],'x')
end

summary = {'maxPSD';'maxPSDfreq';'normPSDmax';'nPSDmaxFreq';'maxCoh';'maxCohFreq'};
summaryData = {data2export.maxPSD; data2export.maxPSDidx; data2export.normPSDmax;data2export.normPSDmaxIdx;data2export.maxCoh;data2export.maxCohIdx};

names = {'CbL v M1R','CbL v S1R','CbL v CbR','CbR v M1L','CbR v S1L','CbR v CbL'};
% avgPSD = cell2mat(data2export.meanPSD);
% normPSD = cell2mat(data2export.normPSD);
coher = cat(2,data2export.coh{1},data2export.coh{2});
coher = cell2mat(coher);
maxCoh = cat(2,data2export.maxCoh(1,:),data2export.maxCoh(2,:));
maxCohIdx = cat(2,data2export.maxCohIdx(1,:),data2export.maxCohIdx(2,:));

writecell(summary,[params(1).recName '.xlsx'],'Sheet','summary','Range','A2');
writecell(names,[params(1).recName '.xlsx'],'Sheet','summary','Range','B1');
writecell(summaryData,[params(1).recName '.xlsx'],'Sheet','summary','Range','B2')

writecell(names,[params(1).recName '.xlsx'],'Sheet','meanPSD','Range','A1');
writematrix(data2export.meanPSD,[params(1).recName '.xlsx'],'Sheet','meanPSD','Range','A2')

writecell(names,[params(1).recName '.xlsx'],'Sheet','normPSD','Range','A1');
writematrix(data2export.normPSD,[params(1).recName '.xlsx'],'Sheet','normPSD','Range','A2')

writecell(names,[params(1).recName '.xlsx'],'Sheet','coherence','Range','A1');
writematrix(coher,[params(1).recName '.xlsx'],'Sheet','coherence','Range','A2')

save([params(1).recName '.mat'])
