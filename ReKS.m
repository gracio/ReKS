function [treeStruct,depthFromTop] = ReKS(A,A2,treeStruct,currentNodeID,currentGroupGeneID,depthFromTop,data,minClustSize)
% first recursion: 14 sec for the svds part for top 20 eigenvectors, 18s including kmeans
depthFromTop = depthFromTop + 1;

% stopping criterion
if length(A) <= minClustSize
    
    % stop when there are less than two genes
    depthFromTop = depthFromTop - 1; % depth from top is decreased as return to the parent node
    return
    
else
    
    % perform spectral clustering
    [discComp,U,S] = clusterKmeans(A,A2);
    fprintf('did a clusterKmeans\n')
    depthFromTop
    % taking care of weird bug...some genes don't belong to any centroid
    check = sum(discComp,1);
    toDelete = find(check == 0);
    if ~isempty(toDelete)
        discComp(:,toDelete) = [];
    end
    
    % put the eigen values and eigen vectors to current node
    [treeStruct.U] = treeStruct.U.set(currentNodeID,U);
    [treeStruct.S] = treeStruct.S.set(currentNodeID,S);
    [treeStruct.discComp] = treeStruct.discComp.set(currentNodeID,discComp);
    
    % recurse on each cluster
    for i = 1 : size(discComp,2)
        
        groupMembers = find(discComp(:,i));
        
        % add members to this subnode
        [treeStruct.groupMembers, subNodeID] = treeStruct.groupMembers.addnode(currentNodeID,currentGroupGeneID(groupMembers));
        % add number of descendants to this subnode
        [treeStruct.numDescendants] = treeStruct.numDescendants.addnode(currentNodeID,length(groupMembers));
        % add nnodeID to this subnode
        [treeStruct.nodeID] = treeStruct.nodeID.addnode(currentNodeID,subNodeID);
        % add empty eigen values and eigen vector space holdesr to this subnode
        [treeStruct.U] = treeStruct.U.addnode(currentNodeID,[]);
        [treeStruct.S] = treeStruct.S.addnode(currentNodeID,[]);
        [treeStruct.discComp] = treeStruct.discComp.addnode(currentNodeID,[]);
        % add depth from top to this subnode
        [treeStruct.depthFromTop] = treeStruct.depthFromTop.addnode(currentNodeID,depthFromTop);
        % add depth from bottom to this subnode
        [treeStruct.depthFromBottom] = treeStruct.depthFromBottom.addnode(currentNodeID,0); % assign default depth from bottom = 0
        
        if ~isempty(data)
            % add new centroid to this subnode
            temp = mean(data(currentGroupGeneID(groupMembers),:),1);
            [treeStruct.centroid] = treeStruct.centroid.addnode(currentNodeID,temp);
            % add distance to parent node to this subnode
            temp = (sum((treeStruct.centroid.get(currentNodeID) - temp ).^2 ))^0.5;
            [treeStruct.distFromParent] = treeStruct.distFromParent.addnode(currentNodeID,temp );
            % add distance to root node to this subnode
            [treeStruct.distFromRoot] = treeStruct.distFromRoot.addnode(currentNodeID,treeStruct.distFromRoot.get(currentNodeID) + temp );
        end
        
        % carve out the submatrix to decompose
        Asub = A(groupMembers,groupMembers);
        if ~isempty(A2)
            A2temp = A2(groupMembers,groupMembers);
        else
            A2temp = [];
        end
        [treeStruct,depthFromTop] = ReKS(Asub,A2temp,treeStruct,subNodeID,currentGroupGeneID(groupMembers),depthFromTop,data,minClustSize);
        
        % depth from bottom is modified if a return occurs with a deeper
        % leaf
        if (treeStruct.depthFromBottom.get(subNodeID) + 1) > treeStruct.depthFromBottom.get(currentNodeID)
            [treeStruct.depthFromBottom] = treeStruct.depthFromBottom.set(currentNodeID,(treeStruct.depthFromBottom.get(subNodeID) + 1));
        end
        
    end
    depthFromTop = depthFromTop - 1; % depth from top is decreased as return to the parent node
    
end

