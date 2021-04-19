% This script reads an HDF5 file that contains 2 datasets. The first one 
% is extendible and contains many height field snapshots. The second
% dataset contains metadata, such as the number of ghost points of the
% height fields. Then the script removes the ghost zones and selects the
% individual snapshots.

data = h5read("height_snapshots.h5","/data")';
metadata = h5read("height_snapshots.h5","/metadata")';

samples = metadata.Value(1);
nrows   = metadata.Value(7);
ncols   = metadata.Value(8);
nghost  = metadata.Value(9);

% Remove ghost zones.

data(: , 1:nghost)         = [];
data(: , end-nghost+1:end) = [];
data(1:nghost , :)         = [];
data(end-nghost+1:end , :) = [];

START = 0;
END   = 0;
SKIP  = 2*nghost-1;

for i=1:(samples-1)
   
    START = i*nrows + 1;
    END = START + SKIP;
    data(START : END , :) = [];
   
end

% Display/Select the height field snapshots.

START_ = 1;
END_   = nrows;

for i = 2:samples+1
    disp(data( START_ : END_ , :));
    START_ = START_ + nrows;
    END_   = START_ + nrows - 1;
end
