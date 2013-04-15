function [result] = runClustering(processedData,A2, clusteringMethod,minClusterSize, K)
% processed data is a structure with fields
% 'data','genes','headers','sparseSymA', 'sparseA'(?), 'A'(?), 'dist'.
% clustering Method = 'ReKS', 'TReKS',
% 'Hierarchical','HierarchicalSpectral','Kmeans','KmeansSpectral'. K is needed for
% hierarchical clustering, HierarchicalSpectral,  Kmeans and KmeansSpectral

tic
set(0,'RecursionLimit',1000)

fprintf('initializing...\n')

switch clusteringMethod
    case {'ReKS','TReKS'} % uses sparseSymA and data
        
        % initialize
        treeStruct.groupMembers = tree(1:length(processedData.sparseSymA)); % need recursion
        treeStruct.numDescendants = tree(length(processedData.sparseSymA)); % calculate along the way
        treeStruct.S = tree; % calculate along the way
        treeStruct.U = tree; % calculate along the way
        treeStruct.discComp = tree; % calculate along the way
        treeStruct.nodeID = tree(1);
        
        treeStruct.depthFromTop = tree(0); % calculate along the way, = number of parents, 0 for root node
        treeStruct.depthFromBottom = tree(0); % calculate along the way, = number of descendant generations, 0 for leaf nodes
        
        treeStruct.distFromParent = tree(0); % can calculate aftermath
        treeStruct.centroid = tree(mean(processedData.data,1)); % can calculate aftermath
        treeStruct.distFromRoot = tree(0);
        
        currentNodeID = 1;
        currentGroupGeneID = 1:length(processedData.sparseSymA);
        depthFromTop = 0;
        
        
        if strcmp(clusteringMethod,'ReKS') % do ReKS , BR Elapsed time is 441.962956 seconds.
            fprintf('performing recursive spectral clustering ReKS...\n')
            %[treeStruct,depthFromTop] = recursiveSpectralClustering(processedData.sparseSymA,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data);
            [treeStruct,depthFromTop] = ReKS(processedData.sparseSymA,A2,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data,minClusterSize-2);
            %[treeStruct,depthFromTop] = ReKS_lite(processedData.sparseSymA,A2,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data,minClusterSize-2);
            
        else % do T-ReKS, BR Elapsed time is 392.584826 seconds
            fprintf('performing transfer recursive spectral clustering, T-ReKS...\n')
            %[treeStruct,depthFromTop] = recursiveSpectralClustering_weighted(processedData.sparseSymA,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data);
            [treeStruct,depthFromTop] = TReKS(processedData.sparseSymA,A2,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data,minClusterSize);
            %[treeStruct,depthFromTop] = TReKS_lite(processedData.sparseSymA,A2,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,processedData.data,minClusterSize);
        end
        
        result.treeStruct = treeStruct;
        
        fprintf(['cutting by minimum cluster size ..\n'])
        % cut when cluster sizes gets smaller than minClustSize
        
        [temp] = cutByMinClust(treeStruct,tree(),minClusterSize,[],1,1);
        temp = sortrows(temp);
        result.groupMembership = temp(:,2);
        
        
    case 'Hierarchical' % uses data, K BR%Elapsed time is 90.037161 seconds.
        fprintf(['performing heirarchical clustering with maximum of ' num2str(K) ' clusters....\n'])
        result.groupMembership = clusterdata(processedData.data,'maxclust', K);
        
    case 'HierarchicalSpectral' % uses sparseSymA, K, %Elapsed time is 21.485911 seconds.
        
        fprintf('performing heirarchical clustering on top 3 eigenvectors....\n')
        [hierarchical_top3.discComp,hierarchical_top3.U,hierarchical_top3.S] = clusterKmeans(processedData.sparseSymA,[]);
        top3 = hierarchical_top3.U(:,2:4);
        hierarchical_top3.normed_top3 = top3./repmat((sum(top3.^2,2)).^0.5,1,3);
        result.groupMembership = clusterdata(hierarchical_top3.normed_top3,'maxclust',K);
        
    case 'Kmeans' % warning: time consuming
        fprintf('performing vanilla kmeans....\n') % needs data, K
        result.groupMembership = kmeans(processedData.data,K);
        
    case 'KmeansSpectral' % uses sparseSYmA, K, Elapsed time is 395.561184 seconds.
        fprintf('performing spectral kmeans....\n')
        [discComp,U,S] = clusterKmeans(processedData.sparseSymA,[],K);
        result.groupMembership = zeros(size(discComp,1),1);
        for i=1:size(discComp,2)
            result.groupMembership = result.groupMembership + discComp(:,i)*i;
        end
end

[result.groupID,result.groupSizes] = countOccurences(result.groupMembership);
result.k = length(result.groupID);

for i=1:length(result.groupID)
    result.groupMembers{i} = find(result.groupMembership == result.groupID(i));
end


toc