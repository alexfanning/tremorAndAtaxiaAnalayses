%
%
%
%
%

clear; close all

cameraType = 1;

if cameraType == 1
    sf = 60;
     
    % Extract data from dlc excel workbooks
    filename = uigetfile('*.xlsx');
    sheets = sheetnames(filename);
    if strcmpi(sheets{end},'Key metrics')
        sheets = sheets(1:end-1);
    end

    orderGiven = readcell(filename);
    orderGiven = orderGiven(1:4,7);
    order(1) = strmatch('C',orderGiven); % LF
    order(2) = strmatch('A',orderGiven); % RF
    order(3) = strmatch('D',orderGiven); % LH
    order(4) = strmatch('B',orderGiven); % RH

    cali = inputdlg('Calibration value: ','Scale movement');
    cali = str2double(cali{1});

    stanceWdthFL = NaN(20,10);
    stanceWdthHL = NaN(20,10);
    fpHpRight = NaN(20,10);
    fpHpLeft = NaN(20,10);
    fpTimeDiff = NaN(20,10);
    hpTimeDiff = NaN(20,10);
    LFRHtimeDiff = NaN(20,10);
    RFLHtimeDiff = NaN(20,10);
    strideDur = cell(1);
    cadence = cell(1);
    strideLngth = cell(1);
    normSpeed = cell(1);

    for j = 1:length(sheets)
        dataTable{j} = readmatrix(filename,"Sheet",j);
        num = 1;
        for t = order
            iter = 1;
            for i = t:4:size(dataTable{j},1)
                paw{num,j}(iter,:) = dataTable{j}(i,1:6);
                iter = iter + 1;
            end
            num = num + 1;
        end
    
        %
        for i = 1:4
            strideDur{i}(:,j) = NaN(20,1);
            cadence{i}(:,j) = NaN(20,1);
            strideLngth{i}(:,j) = NaN(20,1);
            normSpeed{i}(:,j) = NaN(20,1);

            pos{i,j} = paw{i,j}(:,2:3);
            
            % Compute velocity
            time = paw{i,j}(:,1) * ((1000/sf) * 0.001);
            % position (cms)
            position1 = (pos{i,j}(:,1) / cali);
            position2 = (pos{i,j}(:,2) / cali);

            dt = diff(time);
            dPos1 = diff(position1);
            dPos2 = diff(position2);
            v1 = dPos1 ./ dt;
            v2 = dPos2 ./ dt;
            
            % Midpoint times
            t_mid = time(1:end-1) + dt / 2;
            
            % Interpolation of velocity
            vel{i,j}(:,1) = abs(interp1(t_mid, v1, time, 'linear', 'extrap'));
            vel{i,j}(:,2) = abs(interp1(t_mid, v2, time, 'linear', 'extrap'));

            for ii = 2:size(pos{i,j},1)
                % Stride duration (ms)
                strideDur{i}(ii-1,j) = (paw{i,j}(ii,1) - paw{i,j}(ii-1,1)) * (1000/sf);
                % Cadence (m/s-1)
                cadence{i}(ii-1,j) = 1/(strideDur{i}(ii-1,j) / 1000);
                % Stride length (cm)
                strideLngth{i}(ii-1,j) = abs((pos{i,j}(ii) - pos{i,j}(ii-1))) / cali;
                % normalized speed
                normSpeed{i}(ii-1,j) = (strideLngth{i}(ii-1,j) * 0.01) / (strideDur{i}(ii-1,j) * 0.001);
            end

            minSize(i,j) = size(pos{i,j},1);
                    
                    % Y excursion aligned to swing phase
                    
                    % Forepaw planting relative to hindpaw planting (take forepaw last planting position and hindpaw first planting
                    % position difference -- ab = sqrt((a-b)*(a-b) +(a-b)*(a-b))
                    % first two abs is the x coordinate, 2nd two are y coordinate
                    % bd is the same as ab
                    % ab-bd is  ab/bd
        
                    % X and y excursion of each limb
                    %strideTrajectory{i,ii-1} = pos{i,j}(stridePosIdxs{i}(ii-1,2):stridePosIdxs{i}(ii,2));
                    %pos{i,j}(stridePosIdxs{i}(ii-1,2))
        end
    
%         minLng = cellfun(@(x) size(x,1),pos(:,j),'UniformOutput',false);
%         minLng = cell2mat(minLng);
        minLng(j) = min(minSize(:,j));
    
        % stance width (cm)
        for ii = 1:minLng(j)
            stanceWdthFL(ii,j) = abs((pos{1,j}(ii,2) - pos{2,j}(ii,2))) / cali;
            stanceWdthHL(ii,j) = abs((pos{3,j}(ii,2) - pos{4,j}(ii,2))) / cali;
        end
    
        % fp-hp
        if pos{2,j}(1,1) > pos{4,j}(1,1) && pos{2,j}(1,1) < pos{2,j}(2,1)
            [val,idx] = min(abs(pos{2,j}(1,1)-pos{4,j}(:,1)));
            direction = 1;
        elseif pos{2,j}(1,1) > pos{4,j}(1,1) && pos{2,j}(1,1) > pos{2,j}(2,1)
            [val,idx] = min(abs(pos{4,j}(1,1)-pos{2,j}(:,1)));
            direction = 2;
        elseif pos{4,j}(1,1) > pos{2,j}(1,1) && pos{4,j}(1,1) < pos{4,j}(2,1)
            [val,idx] = min(abs(pos{4,j}(1,1)-pos{2,j}(:,1)));
            direction = 3;
        else
            [val,idx] = min(abs(pos{2,j}(1,1)-pos{4,j}(:,1)));
            direction = 4;
        end

        if idx(1) >= 3
            idx(1) = idx(1)-1;
        end
    
        if pos{1,j}(1,1) > pos{3,j}(1,1) && pos{1,j}(1,1) < pos{1,j}(2,1)
            [val(2),idx(2)] = min(abs(pos{1,j}(1,1)-pos{3,j}(:,1)));
            direction(2) = 1;
        elseif pos{1,j}(1,1) > pos{3,j}(1,1) && pos{1,j}(1,1) > pos{1,j}(2,1)
            [val(2),idx(2)] = min(abs(pos{3,j}(1,1)-pos{1,j}(:,1)));
            direction(2) = 2;
        elseif pos{3,j}(1,1) > pos{1,j}(1,1) && pos{3,j}(1,1) < pos{3,j}(2,1)
            [val(2),idx(2)] = min(abs(pos{3,j}(1,1)-pos{1,j}(:,1)));
            direction(2) = 3;
        else
            [val(2),idx(2)] = min(abs(pos{1,j}(1,1)-pos{3,j}(:,1)));
            direction(2) = 4;
        end

        if idx(2) >= 3
            idx(2) = idx(2)-1;
        end
    
        for ii = 1:minLng(j)-1
            % fp-hp (cm)
            if direction(1) == 1 || direction(1) == 4
                fpHpRight(ii,j) = abs((pos{2,j}(ii,1) - pos{4,j}(ii+idx(1)-1,1))) / cali;
            elseif direction(1) == 2 || direction(1) == 3
                fpHpRight(ii,j) = abs((pos{2,j}(ii+idx(1)-1,1) - pos{4,j}(ii,1))) / cali;
            end

            if direction(2) == 1 || direction(2) == 4
                fpHpLeft(ii,j) = abs((pos{1,j}(ii,1) - pos{3,j}(ii+idx(2)-1,1))) / cali;
            elseif direction(2) == 2 || direction(2) == 3
                fpHpLeft(ii,j) = abs((pos{1,j}(ii+idx(2)-1,1) - pos{3,j}(ii,1))) / cali;
            end
        end 
    
        for ii = 1:minLng(j)
    
            % timing difference between forepaws planting (ms)
            fpTimeDiff(ii,j) = abs((paw{1,j}(ii,1) - paw{2,j}(ii,1)) * (1000/sf));
            % timing difference between hindpaws planting (ms)
            hpTimeDiff(ii,j) = abs((paw{3,j}(ii,1) - paw{4,j}(ii,1)) * (1000/sf));
    
            % timing difference between LF-RH and RF-LH
            LFRHtimeDiff(ii,j) = abs((paw{1,j}(ii,1) - paw{4,j}(ii,1)) * (1000/sf));
            RFLHtimeDiff(ii,j) = abs((paw{2,j}(ii,1) - paw{3,j}(ii,1)) * (1000/sf));
        end

    end
    
    for i = 1:4
        strideDurAvg(i) = mean(strideDur{i},'all','omitnan');
        strideDurCv(i) = std(strideDur{i},0,'all','omitnan') / strideDurAvg(i);
        cadenceAvg(i) = mean(cadence{i},'all','omitnan');
        cadenceCv(i) = std(cadence{i},0,'all','omitnan') / cadenceAvg(i);
        strideLngthAvg(i) = mean(strideLngth{i},'all','omitnan');
        strideLngthCv(i) = std(strideLngth{i},0,'all','omitnan') / strideLngthAvg(i);
        normSpeedAvg(i) = mean(normSpeed{i},'all','omitnan');
        normSpeedCv(i) = std(normSpeed{i},0,'all','omitnan') / normSpeedAvg(i);

        stanceWdthFLavg = mean(stanceWdthFL,'all','omitnan');
        stanceWdthFLcv = std(stanceWdthFL,0,'all','omitnan') / stanceWdthFLavg;
        stanceWdthHLavg = mean(stanceWdthHL,'all','omitnan');
        stanceWdthHLcv = std(stanceWdthHL,0,'all','omitnan') / stanceWdthHLavg;

        fpHpRightAvg = mean(fpHpRight,'all','omitnan');
        fpHpRightCv = std(fpHpRight,0,'all','omitnan') / fpHpRightAvg;
        fpHpLeftAvg = mean(fpHpLeft,'all','omitnan');
        fpHpLeftCv = std(fpHpLeft,0,'all','omitnan') / fpHpLeftAvg;

        fpTimeDiffAvg = mean(fpTimeDiff,'all','omitnan');
        fpTimeDiffCv = std(fpTimeDiff,0,'all','omitnan') / fpTimeDiffAvg;
        hpTimeDiffAvg = mean(hpTimeDiff,'all','omitnan');
        hpTimeDiffCv = std(hpTimeDiff,0,'all','omitnan') / hpTimeDiffAvg;

        LFRHtimeDiffAvg = mean(LFRHtimeDiff,'all','omitnan');
        LFRHtimeDiffCv = std(LFRHtimeDiff,0,'all','omitnan') / LFRHtimeDiffAvg;
        RFLHtimeDiffAvg = mean(RFLHtimeDiff,'all','omitnan');
        RFLHtimeDiffCv = std(RFLHtimeDiff,0,'all','omitnan') / RFLHtimeDiffAvg;

        velCompile{i,1} = NaN(20,20);
        velCompile{i,2} = NaN(20,20);
        for j = 1:length(sheets)
            velCompile{i,1}(1:size(vel{i,j},1),j) = vel{i,j}(:,1);
            velCompile{i,2}(1:size(vel{i,j},1),j) = vel{i,j}(:,2);
        end
        velAvg(i,1) = mean(velCompile{i,1},'all','omitnan');
        velAvg(i,2) = mean(velCompile{i,2},'all','omitnan');
        velCv(i,1) = std(velCompile{i,1},0,'all','omitnan') / velAvg(i,1);
        velCv(i,2) = std(velCompile{i,2},0,'all','omitnan') / velAvg(i,1);
    end

elseif cameraType == 2

    sf = 200;
    N = 19;
    formatSpec = '%s';
    fileId = fopen(filename);
    rowHeadings{1} = textscan(fileId,formatSpec,N,'delimiter',',');
    rowHeadings{2} = textscan(fileId,formatSpec,N,'delimiter',',');

    % Extract data from dlc excel workbooks
    filename = uigetfile('*.csv');
    dataTable{1} = readmatrix(filename);

    % Open calibration image and convert 1 cm to pixels
    img = imread(uigetfile('*.tif'));
    imshow(img)
    lineHolder = drawline;
    pos = lineHolder.Position;
    cali = diff(pos);
    cali = cali(1);
    clear lineHolder pos fileId
    close

    %
    count = 2; nextCell = 1;
    for i = 1:3:size(dataTable{1},2)-1
        pos{nextCell} = dataTable{1}(:,count:count+1);
        speed(nextCell) = (((pos{nextCell}(end,1) - pos{nextCell}(1,1)) / cali) *0.01) / (((1000/sf)/1000)*size(dataTable{1},1));
        for k = 1:2
            vel{nextCell}(:,k) = diff(pos{nextCell}(:,k));
            acc{nextCell}(:,k) = diff(vel{nextCell}(:,k));
        end
        strideHolder{nextCell} = find(vel{nextCell}(:,1) < 5);
        endpts = [];
        start = 1;
        ending = 1;
        for ii = 2:length(strideHolder{nextCell})
         
            if strideHolder{nextCell}(ii) == strideHolder{nextCell}(ii-1)+1
                ending = ii;
            else
                endpts = [endpts; start, ending];
                start = ii;
                ending = ii;
            end
        
        end
        consecTimepts{nextCell} = endpts;
        strideIdxs{nextCell} = strideHolder{nextCell}(consecTimepts{nextCell});
        if size(strideIdxs{nextCell},2) > 1
            stridePosIdxs{nextCell}(:,1) = strideIdxs{nextCell}(:,1);
            stridePosIdxs{nextCell}(:,2) = strideIdxs{nextCell}(:,2) + 1;
            for ii = 2:size(stridePosIdxs{nextCell},1)
                % Stride duration (ms)
                strideDur{nextCell}(ii-1) = (stridePosIdxs{nextCell}(ii,1) - stridePosIdxs{nextCell}(ii-1,1)) * (1000/sf);
                % Cadence (m/s-1)
                cadence{nextCell}(ii-1) = 1/(strideDur{nextCell}(ii-1) / 1000);
                % Stride length (cm)
                strideLngth{nextCell}(ii-1) = (pos{nextCell}(stridePosIdxs{nextCell}(ii,1)) - pos{nextCell}(stridePosIdxs{nextCell}(ii-1,1))) / cali;
                
                % Y excursion aligned to swing phase
                
                % Forepaw planting relative to hindpaw planting (take forepaw last planting position and hindpaw first planting
                % position difference -- ab = sqrt((a-b)*(a-b) +(a-b)*(a-b))
                % first two abs is the x coordinate, 2nd two are y coordinate
                % bd is the same as ab
                % ab-bd is  ab/bd
    
                % X and y excursion of each limb
                strideTrajectory{nextCell,ii-1} = pos{nextCell}(stridePosIdxs{nextCell}(ii-1,2):stridePosIdxs{nextCell}(ii,2));
                %pos{nextCell}(stridePosIdxs{nextCell}(ii-1,2))
            end
        end
        count = count + 3; nextCell = nextCell + 1;
    end

    % ADD STANCE DURATION
    
    stanceWidth(1) = mean((pos{1,j}(:,2) - pos{2,j}(:,2)) / cali);
    stanceWidth(2) = mean((pos{3,j}(:,2) - pos{4,j}(:,2)) / cali);

    % Plot locomotion
    plot(vel{1})
    scatter(dataTable(:,2),dataTable(:,3))
end


%% Export data
dataLabels = {'strideLength (cm)';'strideLengthCV';'strideDuration (ms)';'strideDurationCV';'cadence (m/s-1)';'cadenceCv'; 'speed (cm/s)'; 'speedCv'; 'normSpeed'; 'normSpeedCv (m/s)'};
data2export = {strideLngthAvg; strideLngthCv; strideDurAvg; strideDurCv; cadenceAvg; cadenceCv; velAvg(:,1)'; velCv(:,1)'; normSpeedAvg; normSpeedCv};

dataLabels2 = {'stanceWidthFLs (cm)';'stanceWidthFLsCV';'stanceWidthHLs (cm)';'stanceWidthHLsCV';...
    'fpHpRight (cm)';'fpHpRightCv';'fpHpLeft (cm)';'fpHpLeftCv';'RFLHtimeDiff (ms)';'RFLHtimeDiffCV';...
    'LFRHtimeDiff (ms)';'LFRHtimeDiffCv';'fpTimeDiff (ms)';'fpTimeDiffCV';'hpTimeDiff(ms)';'hpTimeDiffCV'};
data2export2 = {stanceWdthFLavg; stanceWdthFLcv; stanceWdthHLavg; stanceWdthHLcv; fpHpRightAvg; fpHpRightCv; fpHpLeftAvg; fpHpLeftCv;...
    RFLHtimeDiffAvg; RFLHtimeDiffCv; LFRHtimeDiffAvg; LFRHtimeDiffCv; fpTimeDiffAvg; fpTimeDiffCv; hpTimeDiffAvg; hpTimeDiffCv};

dataLabels3 = {'LF','RF','RH','LH'};

writecell(dataLabels2,filename,'Sheet','Key metrics','Range','A1')
writecell(data2export2,filename,'Sheet','Key metrics','Range','B1')

writecell(dataLabels,filename,'Sheet','Key metrics','Range','A20')
writecell(data2export,filename,'Sheet','Key metrics','Range','B20')

writecell(dataLabels3,filename,'Sheet','Key metrics','Range','B19')

%% Plot data

clear i idx ii iter j t

save([filename(1:end-5) '.mat'])
