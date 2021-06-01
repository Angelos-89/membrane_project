PATH = "/home/angelos-89/runs/long_runs/";
PINS = ["unpinned/","pinning/8%/","pinning/16%/","pinning/32%/"];
SRUNS = ["run0/","run1/","run2/","run3/","run4/","run5/"];
FILENAMES = [];
TAUS = [0.01,0.02,0.04,0.08,0.1,0.2,0.4,0.8,1.0,2.0,4.0,8.0,10.0,40.0,100.0,1000.0]';
CERTAIN_TAU_INDECES = [0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60];

TOWRITE = [];
start = 3*1e3;

for pin=1:length(PINS)

    MEANS_EXCESS_AREA=[];
    MIN_EXCESS_AREA = [];
    MAX_EXCESS_AREA = [];
    
    for jj = 1:length(CERTAIN_TAU_INDECES)

        CERTAIN_TAU_INDEX = CERTAIN_TAU_INDECES(jj);
        NUMS = [CERTAIN_TAU_INDEX, CERTAIN_TAU_INDEX+1, CERTAIN_TAU_INDEX+2, CERTAIN_TAU_INDEX+3];
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
            eqArea = rawArea(start:end);
            meanOfArea = mean(eqArea);

            rawPrjArea = [];
            for i=1:length(SRUNS)
                rawPrjAreaData = importdata(FILENAMES(i)).data(:,4);
                rawPrjArea = [rawPrjArea; rawPrjAreaData];
            end
            eqPrjArea = rawPrjArea(start:end);
            meanOfPrjArea = mean(eqPrjArea);

            excessAreas = [excessAreas , (meanOfArea - meanOfPrjArea)/meanOfArea];
        end
        format long
        MEANS_EXCESS_AREA = [MEANS_EXCESS_AREA; mean(excessAreas)];
        MIN_EXCESS_AREA = [MIN_EXCESS_AREA; min(excessAreas)];
        MAX_EXCESS_AREA = [MAX_EXCESS_AREA; max(excessAreas)];        
    end
    TOWRITE = [TOWRITE MEANS_EXCESS_AREA, MIN_EXCESS_AREA, MAX_EXCESS_AREA];
    
end
%write to file
dataToWrite = [TAUS TOWRITE];
writematrix(dataToWrite,"excess-area-data.txt",'Delimiter','tab');

%ffs fix this nonsense
infid = fopen('excess-area-data.txt', 'rt');
outfid = fopen('excess-area.data', 'wt');
fprintf( outfid, "%%taus \t Excess area(unpinned) \t Excess area min \t Excess area max \t Excess area(8) \t Excess area min \t Excess areamax \t Excess area(16) \t Excess area min \t Excess area max \t Excess area(32) \t Excess area min \t Excess area max \n");  %the text you are adding at the beginning
while true
    thisline = fgetl(infid);    %read line from input file
    if ~ischar(thisline); break; end    %reached end of file
    fprintf( outfid, '%s\n', thisline );   %write the line to the output file
end
fclose(infid);
fclose(outfid);
delete excess-area-data.txt;