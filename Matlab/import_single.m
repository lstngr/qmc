function data = import_single(filename,nCols)
% Imports simulation "out" files. The function looks for nCols columns of 
% data and returns a cell of size 1*nCols storing the values found in the 
% columns of the out-file.
% Doesn't recognize the name of the columns, the user is responsible for
% ordering his observables in a convenient way.
formatSpec = repmat('%f', [1 nCols]);
fileID = fopen(filename,'r');
data = textscan(fileID,formatSpec,'CommentStyle','#','Delimiter',';');

end