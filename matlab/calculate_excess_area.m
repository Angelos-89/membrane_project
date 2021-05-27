%RETURN HERE
PATH = "/home/angelos-89/runs/long_runs/";
PINS = ["unpinned/","pinning/8%/","pinning/16%/","pinning/32%/"];
SRUNS = ["run0/","run1/"];
FILENAMES = [];
taus = [0.01,0.02,0.04,0.08,0.1,0.2,0.4,0.8,1.0,2.0,4.0,8.0,10.0,40.0,100.0,1000.0]';
CERTAIN_TAU_INDECES = [0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60];

TOWRITE = [];

for pin=1:length(PINS)

    MEANS_EXCESS_AREA=[];
    SDE_EXCESS_AREA=[];
    
    for jj = 1:length(CERTAIN_TAU_INDECES)

        CERTAIN_TAU_INDEX = CERTAIN_TAU_INDECES(jj);
        NUMS = [CERTAIN_TAU_INDEX,CERTAIN_TAU_INDEX+1 ...
                CERTAIN_TAU_INDEX+2,CERTAIN_TAU_INDEX+3];
        NUMS = string(NUMS);

        excessAreas = [];

        for run = 1:length(NUMS)

            FILENAMES = [];
            for i=1:length(SRUNS)
                FILENAMES = [FILENAMES, PATH+PINS(pin)+SRUNS(i)+"timeseries_"+NUMS(run)+".txt"]; 
            end

            rawArea = [];
            for i=1:length(SRUNS)
                areaData = importdata(FILENAMES(i)).data(:,3);
                rawArea = [rawArea; areaData];
            end
            eqArea = rawArea(2e3:end);
            maxInd = 2^13;
            eqArea = eqArea(1:maxInd);
            lenArea = length(eqArea);
            meanOfArea = mean(eqArea);

            rawPrjArea = [];
            for i=1:length(SRUNS)
                rawPrjAreaData = importdata(FILENAMES(i)).data(:,4);
                rawPrjArea = [rawPrjArea; rawPrjAreaData];
            end
            eqPrjArea = rawPrjArea(2e3:end);
            eqPrjArea = eqPrjArea(1:maxInd);
            lenPrjArea = length(eqPrjArea);
            meanOfPrjArea = mean(eqPrjArea);

            excessAreas = [excessAreas (meanOfArea - meanOfPrjArea)/meanOfArea];
        end
        format long
        MEANS_EXCESS_AREA = [MEANS_EXCESS_AREA; mean(excessAreas)];
        varExcessArea = var(excessAreas);
        SDE_EXCESS_AREA = [SDE_EXCESS_AREA; sqrt(varExcessArea)/sqrt(3)];
    end
    TOWRITE = [TOWRITE MEANS_EXCESS_AREA SDE_EXCESS_AREA];
    
end
%write to file
dataToWrite = [taus TOWRITE];
writematrix(dataToWrite,"excess-area-data.txt",'Delimiter','tab');
