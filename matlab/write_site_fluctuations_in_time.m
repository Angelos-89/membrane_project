%{ 
    This script takes an extendible HDF5 dataset and
    reads the selected hyperslabs. Each hyperslab corresponds
    to the snapshot of the height field at a specific time.
    This slab is reshaped to a 1D array and then appended to a
    matrix that has as many rows, as time points (snapshots).
    Each column corresponds to the "time-evolution" of a specific
    site/node. The purpose is to find the autocorrelation function
    of every time series/column, and then average them to obtain the
    averaged over many points autocorrelation function of the height.
    Since this operation takes time, the final matrix is written to txt
    for later use. Reading this txt takes much less time.
 %}
%%-----------------These change every time-------------%%
FILENAME_TO_READ = "snapshots_0.h5";
FILENAME_TO_WRITE = 'timeSrData-tau-4-pin-8.txt';
totNodes = 6400;  % find this from PARAMS.txt
firstSlab = 2e4;  % find this by data inspection.
lastSlab = 3e4; % read this from stats.txt
%------------------------------------------------------%%
extend = lastSlab-firstSlab;
timeSrData = [];
for slab = firstSlab : lastSlab
    hyperSlab = Select_Hyperslab(FILENAME_TO_READ, slab);
    reshapedSlab = reshape(hyperSlab, 1, totNodes);    
    timeSrData = [timeSrData; reshapedSlab];
end
writematrix(timeSrData,FILENAME_TO_WRITE,'Delimiter','tab');
