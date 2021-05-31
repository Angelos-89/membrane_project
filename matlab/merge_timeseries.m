%% 
%Change the PATH to merge the series for all cases (pinned and unpinned membrane)
PATH = "/home/angelos-89/runs/long_runs/unpinned/"; %<---------change if needed
SRUNS = ["run0/","run1/","run2/","run3/"]; %<---------change if needed
FILENAMES = [];
NUM = 0;
%%
%Create one timeseries from the consecutive runs, drop the equilibration part of the first run.
for i=1:64
    areaSeries = [];
    prjAreaSeries = [];
    energySeries = [];
    for run = 1 : length(SRUNS)
       FILENAME = PATH + SRUNS(run) + "timeseries_" + string(i-1) +".txt"; %<---------change if needed
       areaFromIndividualRun = importdata(FILENAME).data(:,3);
       areaSeries = [areaSeries ; areaFromIndividualRun]; 

       prjAreaFromIndividualRun = importdata(FILENAME).data(:,4);
       prjAreaSeries = [prjAreaSeries ; prjAreaFromIndividualRun];

       energyFromIndividualRun = importdata(FILENAME).data(:,9);
       energySeries = [energySeries ; energyFromIndividualRun];
    end
    start = 5*1e3; %<---------change if needed
    areaSeries = areaSeries(start:end);
    prjAreaSeries = prjAreaSeries(start:end);
    energySeries = energySeries(start:end);
    excessAreaSeries = (areaSeries-prjAreaSeries)./areaSeries;
    %Write to .txt file
    dataToBeWritten = [areaSeries prjAreaSeries excessAreaSeries energySeries];
    SAVE_FILENAME = "merged_series_" + string(i-1)+ ".txt"; %<---------change if needed
    fileID = fopen(SAVE_FILENAME,'w');
    fprintf(fileID,'%s \t %s \t %s \t %s\n','Area','Projected Area','Excess Area','Energy');
    fprintf(fileID,'%f \t %f \t %f \t %f\n', dataToBeWritten);
    fclose(fileID);
    DATA = importdata(SAVE_FILENAME);
end