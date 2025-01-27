
% Data Analysis For Rota Rod Experiments 
%
% Alex Fanning 5/26/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;  

filename = uigetfile('*.xlsx');
params = inputdlg({'Num trials: '; 'Num sessions: '},'Parameters',1,{'3';'6'});
numTrials = str2double(params{1});
numSess = str2double(params{2});

% sheetNames = {'Group 1 Experimental','Group 2 Control'};
sheetNames = sheetnames(filename);

numSheets = numel(sheetNames);
excelData = cell(1, numSheets);
data = cell(1, numSheets);
dataMean = cell(1, numSheets);
colAvgs = cell(1, numSheets);
dataStdErr = cell(1, numSheets);

%% Average trials across each session and across population

legendHold = 0;
%iterate through each sheet
for sheetIdx = 1:numSheets

    if strcmpi(sheetNames{sheetIdx} ,'Stats') || strcmpi(sheetNames{sheetIdx},'array4stats')
        legendHold = 1;
        continue
    else

        %Load data from each Excel sheet
        excelData{sheetIdx} = xlsread(filename,sheetIdx);
        numSxs(sheetIdx) = size(excelData{sheetIdx},2);
        %excelData{sheetIdx}(isnan(excelData{sheetIdx})) = [];
        excelData{sheetIdx} = reshape(excelData{sheetIdx},numTrials*numSess+1, numSxs(sheetIdx));
        data{sheetIdx} = excelData{sheetIdx}(2:end, :);
    
        %for loop iterating through every three rows
        count = 1;
        for i = 1:numTrials:numTrials*numSess
            dataMean{sheetIdx}(:,count) = mean(data{sheetIdx}(i:i+numTrials-1, :),1, 'omitnan');
            dataConcat{sheetIdx,count} = data{sheetIdx}(i:i+numTrials-1,:);
            colAvgs{sheetIdx}(count) = mean(dataConcat{sheetIdx,count},'all','omitnan');
            dataStdErr{sheetIdx}(count) = std(dataConcat{sheetIdx,count},0,'all','omitnan') / sqrt(numel(dataConcat{sheetIdx,count}));
            count = count + 1;
        end
    
        figure(sheetIdx); hold on;
            clrs = [0.16 0.57 0.84; 0.93 0.69 0.13; 0.77 0.06 0.06];
            plot(dataMean{sheetIdx}', 'LineWidth', 1, 'Marker', 'o','MarkerSize',10,'MarkerFaceColor',clrs(sheetIdx,:),'Color', clrs(sheetIdx,:));
            plot(colAvgs{sheetIdx}, 'LineWidth', 2, 'Marker', 'd','MarkerSize',10,'MarkerFaceColor',clrs(sheetIdx,:), 'Color', 'k');
            errorbar(colAvgs{sheetIdx}, dataStdErr{sheetIdx}, 'LineWidth', 2, 'Marker', 'o', 'Color', 'k');
            xlim([0.5 numSess+0.5])
            xlabel('Session');
            ylabel('Time on rotarod (seconds)');
            xticks(1:numSess)
            sessionLabels = cellstr(string(1:size(dataMean{sheetIdx}, 1)));
            if numSess > 5
                sessionLabels(6) = {'Retention'};
            end
            xticklabels(sessionLabels);
            set(gca,'FontSize',14)
            title('Rotarod performance', 'FontSize', 16);
        
            legend(sheetNames{sheetIdx},'FontSize', 12, 'Location', 'best');
            set(gcf, 'Position', [100, 100, 800, 500]);
    
    %     ax = gca;
    %     ax.Box = 'off';
    %     ax.LineWidth = 1.5;
    %     ax.FontSize = 10;
    %     ax.FontName = 'Times New Roman';
    %     ax.TickDir = 'out';
    %     ax.TickLength = [0.02, 0.02];
    %     ax.XAxis.MinorTick = 'off';
    %     ax.YAxis.MinorTick = 'off';
    end

end

figure(5); hold on;
for i = 1:numSheets
    if strcmpi(sheetNames{i},'Stats') || strcmpi(sheetNames{i},'array4stats')
        continue
    else
        plot(colAvgs{i}, 'LineWidth', 2, 'Marker', 'o', 'Color', clrs(i,:));
        errorbar(colAvgs{i}, dataStdErr{i}, 'LineWidth', 2, 'Marker', 'o','MarkerSize',10,'MarkerFaceColor',clrs(i,:),'Color', clrs(i,:));
    end
end
xlim([0.5 numSess+0.5])
xlabel('Session');
ylabel('Time on rotarod (seconds)');
xticks(1:numSess)
sessionLabels = cellstr(string(1:numSess));
if numSess > 5
    sessionLabels(6) = {'Retention'};
end
xticklabels(sessionLabels);
set(gca,'FontSize',14)
title('Rotarod performance', 'FontSize', 16);
if (sheetIdx == 3 && legendHold == 1) || (sheetIdx == 2 && legendHold == 0)
    legend(sheetNames{1},'',sheetNames{2},'FontSize', 12, 'Location', 'best');
elseif (sheetIdx == 3 && legendHold ~= 1) || (sheetIdx == 4 && legendHold == 1)
    legend(sheetNames{1},'',sheetNames{2},'',sheetNames{3},'FontSize', 12, 'Location', 'best');
end
set(gcf, 'Position', [100, 100, 800, 500]);
saveas(5,filename(1:end-5),'fig')

%%  Statistics
    
% Create repeated measures ANOVA model
for m = 1:size(dataConcat,2)
    count = 1;
    for i = 1:size(dataConcat,1)
        tempData = cat(1,dataConcat{i,m}(:));
        dataStat(count:count+numel(dataConcat{i,m})-1,m) = tempData;
        count = count + numel(dataConcat{i,m});
    end
end
group = repmat('A',numel(dataConcat{1}),1);
group(numel(dataConcat{1})+1:numel(dataConcat{1})+numel(dataConcat{2,1})) = repmat('B',numel(dataConcat{2,1}),1);
if (sheetIdx == 3 && legendHold == 0) || (sheetIdx == 4 && legendHold == 1) || (sheetIdx == 5 && legendHold == 1)
    group(numel(dataConcat{1})+numel(dataConcat{2,1})+1:numel(dataConcat{1})+numel(dataConcat{2,1})+numel(dataConcat{3,1})) = repmat('C',numel(dataConcat{3,1}),1);
end
time = (1:numSess)';
if numSess == 3
    t = table(group,dataStat(:,1),dataStat(:,2),dataStat(:,3),'VariableNames',{'group','a1','a2','a3'});
    rm = fitrm(t,'a1-a3 ~ group','WithinDesign',time);
elseif numSess == 5
    t = table(group,dataStat(:,1),dataStat(:,2),dataStat(:,3),dataStat(:,4),dataStat(:,5),'VariableNames',{'group','a1','a2','a3','a4','a5'});
    rm = fitrm(t,'a1-a5 ~ group','WithinDesign',time);
elseif numSess == 6
    t = table(group,dataStat(:,1),dataStat(:,2),dataStat(:,3),dataStat(:,4),dataStat(:,5),dataStat(:,6),'VariableNames',{'group','a1','a2','a3','a4','a5','a6'});
    rm = fitrm(t,'a1-a6 ~ group','WithinDesign',time);
end
ranovatbl = ranova(rm)
writetable(ranovatbl,filename,'Sheet','Stats','WriteRowNames',true)

trialsPerGroup = numSxs*numTrials;
tempOutput = dataStat';
for i = 1:numSess
    outputHolder(i,:) = [tempOutput(i,1:trialsPerGroup(1)) repmat(NaN,1,1) tempOutput(i,trialsPerGroup(1)+1:end)];
    if length(time) == 3
        output4stats(i,:) = [outputHolder(i,1:trialsPerGroup(1) + trialsPerGroup(2)+1) repmat(NaN,1,1) outputHolder(i,trialsPerGroup(1) + trialsPerGroup(2)+2:end)];
    end
end
writematrix(output4stats,filename,'Sheet','array4stats')

save([filename(1:end-5) '.mat'])


    % Add asterisks to indicate significance
% for sheetIdx = 1:numSheets
%     if p_values(sheetIdx) < 0.05
%         sig_indices = find(p_values < 0.05);
%         x_coords = 1:numSessions;
%         y_coords = colAvgs{sheetIdx};
%         
%         if ismember(sheetIdx, sig_indices)
%             significance_level = '***';
%         else
%             significance_level = '';
%         end
%         
%         text(x_coords, y_coords, significance_level, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
%     end
% end
