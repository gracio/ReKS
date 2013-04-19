function [dataStruct] =  ReKS_main(dataStruct,varargin)
% run ReKS on a datastructure such as
%     dataStruct.data
%     dataStruct.headers
%     dataStruct.genes
% to use, run
%     dataStruct = ReKS_main(dataStruct);
%     dataStruct = ReKS_main(dataStruct,35);
%     dataStruct = ReKS_main(dataStruct,35,2);

n = 35; % default for minimum cluster size
dim = 1; % by default, variables to be clustered are ROWS (dimension 1 of the data matrix)

if nargin > 1 % if minim cluster size is supplied
    n = varargin{1};
    if nargin > 2 % if dimension is supplied
        dim = varargin{2};
    end
end


% if user wishes to use dimension 2, transpose the data matrix
if dim == 2
    dataStruct.data = dataStruct.data';
    fprintf('use transposed dimension\n')
end

% perform the data preparation step to obtain sparse symetrical affinity
% matrix. Default options searches for 15 closest neighbors on the
% covariance matrix, does not impute data, and uses 1-pearson correlation as the distance measure
[dataStruct.sparseSymA, dataStruct.sparseA, dataStruct.dist, dataStruct.data] = prepDataForSPC(dataStruct.data,15,0,'corr');

[dataStruct.ReKS] = runClustering(dataStruct,[],'ReKS',n, []); 

% finally, transpose the data matrix back to return to user
if dim==2
    dataStruct.data = dataStruct.data';
end
