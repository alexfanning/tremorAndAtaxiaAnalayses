function [dataOut] = coherence(dataIn,dataOut,parameters,timeFrame,sfNum,type)

n = parameters(sfNum).sf * timeFrame; 
n1s = parameters(sfNum).sf;          
nCoh = n1s * parameters(1).cohSec; 
shiftedframe = n1s * parameters(1).timeshift;
endNum = fix(length(dataIn{2,1}) / n1s) - timeFrame - ceil(parameters(1).timeshift);

for i = 1:4
    for m = 1:3

        if type == 1
            vec = {dataIn{2,1}(:,m) dataIn{2,4}(:,m) dataIn{2,6}(:,m) dataIn{2,2}(:,m)};
        elseif type == 2
            vec = {dataIn{2,2}(:,m) dataIn{2,3}(:,m) dataIn{2,5}(:,m) dataIn{2,1}(:,m)};
        end
    
        for s = 1:endNum
            if i <= 3
                y1 = vec{1}((1+(s-1)*n1s):(n+(s-1)*n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
                y2 = vec{i+1}((1+(s-1)*n1s+shiftedframe):(n+(s-1)*n1s+shiftedframe));
            elseif i == 4
                y1 = vec{2}((1+(s-1)*n1s):(n+(s-1)*n1s));  % Each iteration grabs 20s of data and the starting point is shifted 1 s to the right of the prior starting point
                y2 = vec{3}((1+(s-1)*n1s+shiftedframe):(n+(s-1)*n1s+shiftedframe));
            end
            
            [Cy1y2,FreqCoh] = mscohere(y1,y2,hanning(nCoh),1/2*nCoh,nCoh,parameters(sfNum).sf); % (y1,y2 coherece, with hamming window, take NCoh points, half points overlap, number of fft:NCoh, sampling frequency:fs
            coh(:,s) = Cy1y2;
        
        end
        if type == 1
            dataOut.mCoh{1,i}(:,m) = mean(coh,2);
        
            [dataOut.maxCoh{1}(i,m),dataOut.maxCohIdx{1}(i,m)] = max(dataOut.mCoh{1,i}(parameters(1).extractWndw(1):parameters(1).extractWndw(2),m));
            dataOut.maxCohIdx{1}(i,m) = dataOut.maxCohIdx{1}(i,m) + parameters.extractWndw(1) - 1;
        elseif type == 2
            dataOut.mCoh{2,i}(:,m) = mean(coh,2);
        
            [dataOut.maxCoh{2}(i,m),dataOut.maxCohIdx{2}(i,m)] = max(dataOut.mCoh{2,i}(parameters(1).extractWndw(1):parameters(1).extractWndw(2),m));
            dataOut.maxCohIdx{2}(i,m) = dataOut.maxCohIdx{2}(i,m) + parameters.extractWndw(1) - 1;
        end
    end
end
