%{ 
    This function takes as input two arguments.
    The first argument is an HDF5 file containing 
    a root group.
    Inside the root group there are two datasets.
    The first contains many instances of the class 
    RectMesh.
    The second contains metadata such as the number 
    of ghost points.
    The second argument is the block/hyperslab number
    that you want to read and store in data. 
 %}

function hyperslab = Select_Hyperslab(H5data_FILENAME , N)

    % Open file and root group.
    file_id  = H5F.open(H5data_FILENAME);
    group_id = H5G.open(file_id,'/');
    
    % If there are no metadata written, default behaviour.
    if ~H5L.exists(group_id,'/metadata','H5P_DEFAULT')
        error("Dataset containing metadata is not found.");
    end
    
    % Read metadata.
    mdata  = h5read(H5data_FILENAME,"/metadata")'; 
    blocks = mdata.Value(1);
    nrows  = mdata.Value(7);
    ncols  = mdata.Value(8);
    nghost = mdata.Value(9);
    block_length = nrows + 2*nghost;
    block_width  = ncols + 2*nghost;
    
    % Bound check.
    if N > blocks
        error("The block you have chosen exceeds the number of stored blocks.");
    end
    
    property_list = 'H5P_DEFAULT';
    dataset_id = H5D.open(file_id,'/data'); %check here
    
    % Check existence of dataset containing many snapshots/blocks.
    if ~H5L.exists(group_id,'/data','H5P_DEFAULT')
        error("Dataset containing data is not found.");
    end
    
    dims = [block_length , block_width]; 
    memory_space_id = H5S.create_simple(2 , dims , []);
    file_space_id  = H5D.get_space(dataset_id);
    
    if N==1
        START_ = 1;
    else 
        START_ = (N-1)*block_length + 1;
    end
    
    % Select hyperslab and read it.
    offset =[START_-1 , 0];
    block = [block_length , block_width];
    H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',offset,[],[],block);
    hyperslab = H5D.read(dataset_id,'H5ML_DEFAULT',memory_space_id,file_space_id,property_list)';
   
    % Close resources.
    H5D.close(dataset_id);
    H5F.close(file_id);
end
