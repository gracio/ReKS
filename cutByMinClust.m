function [groupMembership, shortTree] = cutByMinClust(fullTree, shortTree, minClustSize,groupMembership,fullTree_currentNodeID,shortTree_currentNodeID)

%copy node to short tree
[shortTree] = shortTree.set(shortTree_currentNodeID,fullTree.numDescendants.get(fullTree_currentNodeID));

kids = fullTree.numDescendants.getchildren(fullTree_currentNodeID);
if isempty(kids)
    %return membership
    temp = [fullTree.groupMembers.get(fullTree_currentNodeID)' repmat(fullTree_currentNodeID,length(fullTree.groupMembers.get(fullTree_currentNodeID)),1) ];
    groupMembership = [groupMembership; temp];
    %copy node to short tree
    
else
    % look at the sizes of the descendant nodes
    for i=1:length(kids)
        kidsSize(i) = fullTree.numDescendants.get(kids(i));
    end
    isBigClust = (kidsSize >= minClustSize);
    
    if sum(isBigClust) == 0 % all clusters are small, output all the members here
        temp = [fullTree.groupMembers.get(fullTree_currentNodeID)' repmat(fullTree_currentNodeID,length(fullTree.groupMembers.get(fullTree_currentNodeID)),1) ];
        groupMembership = [groupMembership; temp];
        % copy node to short tree
        
    else % cluster is big enough to recurse down
        
        for i=1:length(kids)
            if isBigClust(i) % else recurse on childs with big enough cluster
                % add empty child nodes and obtain their ID
               [shortTree shortTree_kids(i)] = shortTree.addnode(shortTree_currentNodeID,0);
               [groupMembership, shortTree] = cutByMinClust(fullTree, shortTree, minClustSize,groupMembership,kids(i),shortTree_kids(i));
               
            else % and output the child nodes with cluster size too small
                temp = [fullTree.groupMembers.get(kids(i))' repmat(kids(i),length(fullTree.groupMembers.get(kids(i))),1) ];
                groupMembership = [groupMembership; temp];
                 [shortTree] = shortTree.addnode(shortTree_currentNodeID,fullTree.numDescendants.get(kids(i)));
                % copy node to short tree
            end
        end
    end
    
end