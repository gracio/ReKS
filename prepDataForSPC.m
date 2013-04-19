function [S_sparse_sym, S_sparse, distances, data] = prepDataForSPC(data, topNneighbors,KNNimpute,distance)

S_sparse_sym=[];S_sparse=[];distances=[];

% filter genes if needed
if KNNimpute == 1
tic
% fprintf(['total # genes ' num2str(size(data,1)) ' \n']);
% 
% maskVal = genelowvalfilter(data);
% fprintf(['low value filter # genes kept ' num2str(sum(maskVal)) ' \n']);
% maskVar = genevarfilter(data);
% fprintf(['low variance filter # genes kept ' num2str(sum(maskVar)) ' \n']);

% impute data
data = knnimpute(data);
toc
end

tic
% find top neighbors
fprintf('calculating knn..\n')
[topNeighborsIdx,antiC] = knnsearch(data,data,'K',topNneighbors,'Distance',distance);
toc
fprintf('knnsearch done \n')


% calculate affinity matrix
C = 1-antiC;
sig = 1;
A = exp(( -1*( sin( acos(C)/2 )).^2 )/(sig^2) );  

% index into square form
tic
S_sparse = sparse(size(data,1),size(data,1));
for i=1:size(topNeighborsIdx,1)
    S_sparse(i,topNeighborsIdx(i,:)) = A(i,:);
end
toc
fprintf('affinity done \n')

fprintf('make affinity matrix symmetric \n')
% make symmetric
S_sparse_sym = S_sparse' .*  ( ( (S_sparse>0)+(S_sparse'>0) ) == 1) + S_sparse;

tic
fprintf('calculating distance matrix')
distances = pdist(data,'correlation');
distances = squareform(distances);
toc