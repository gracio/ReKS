function [dataStruct] =  ReKS_main(dataStruct,n,dim)
if dim == 2
    dataStruct.data = dataStruct.data';
    fprintf('use transposed dimension\n')
end
[dataStruct.sparseSymA, dataStruct.sparseA, dataStruct.dist, dataStruct.data] = prepDataForSPC(dataStruct.data,15,0,'corr');
[dataStruct.ReKS] = runClustering(dataStruct,[],'ReKS',n, []); 

if dim==2
    dataStruct.data = dataStruct.data';
end
%clusteMembership = treeStruct.groupMembership;