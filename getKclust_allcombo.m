function [ groupMembership, groupMembers,nodeMembers,firstChunk, bundled]  = getKclust_allcombo(treeStruct, k)
% get k clusters from a ReKS tree, using breadth first search nodes and the
% first combination encountered from the left that allows for k clusters.

% definite members are the ones that don't have children that split off
% early, and should be included as single clusters
% candidateMembers are the ones that should be grouped into fewer than
% their numbers, of which bundled are the child nodes that should be grouped together into one
% single cluster

% first chunck consists of all the guys that are single clusters (including
% definiteMembers, otherMembers for cases when k is obtained by just cutting at a level, and first combination of candidateSets
% there are 5 scenarios (k<min#leafs, k>max#leafs, k==1, k is done by cutting tree at a level, k is in between



BSI =  treeStruct.numDescendants.breadthfirstiterator;
BSIlevels = [treeStruct.depthFromTop.Node{treeStruct.numDescendants.breadthfirstiterator}];

% get total number of leafs at each tree level
for i=1:max(BSIlevels)
    numLeafs(i) = sum(BSIlevels == i);
end

% get total number of kids at each tree level
for i=1:length(BSI)
    numKids(i) = length(treeStruct.numDescendants.getchildren(BSI(i)));
end

% need to do this because some nodes end early
for i=0:max(BSIlevels)
    cumNoKids(i+1) = sum(BSIlevels==i & numKids==0);
end
cumNoKids = cumsum(cumNoKids);
cumNoKids(end) = [];
numKatlevel = numLeafs + cumNoKids; % number of k one would obtain at each level

% if k is 1 -> trivial case
if k==1
    fprintf('k==1\n');
    firstChunk = 1;
    bundled = [];
    
    % if k is smaller than number of leafs in the tree
elseif k < min(numKatlevel)
    fprintf('k is smaller than number of leafs in the tree\n')
    
    nodeIDs = treeStruct.nodeID.getchildren(1);
    [uqClusterMem,uqClusters] = orderNelementsIntoKgroups(min(numKatlevel),k);
    groupMembership = [];
    
    for m=1:size(uqClusterMem,1) % for all the combinations
        for n=1:max(uqClusterMem(m,:)) % for all the k labels in that combination
            group=find(uqClusterMem(m,:)== n); % the items in this combo that is of label n
            groupMembers{m}{n} = [];
            nodeMembers{m}{n} = [];
            
           for i=1:length(group)
                nodeMembers{m}{n} = [nodeMembers{m}{n} nodeIDs(group(i))];
                groupMembers{m}{n}  = [groupMembers{m}{n} treeStruct.groupMembers.Node{nodeIDs(group(i))}];
                groupMembership(m,treeStruct.groupMembers.Node{nodeIDs(group(i))}) = n;
            end
        end
    end
    
  return
    
elseif k > max(numKatlevel) % if k is larger than total number of leaf nodes
    fprintf(' k is larger than total number of leaf nodes\n');
    k = max(numKatlevel);
    
    level2focus = find(numKatlevel == k);
    definiteMembers = BSI(find((BSIlevels < level2focus) & (numKids ==0)));
    otherMembers = BSI(find(BSIlevels == level2focus));
    
    %
    firstChunk = [definiteMembers otherMembers];
    bundled = [];
    
    % if k happens to be obtained by just cutting across the tree
elseif ismember(k, numKatlevel)
    fprintf(['k is obtained by cutting at k==' num2str(find(numKatlevel == k)) '\n'])
    % output group structure right here, no more work to do.
    
    level2focus = find(numKatlevel == k);
    definiteMembers = BSI(find((BSIlevels < level2focus) & (numKids ==0)));
    otherMembers = BSI(find(BSIlevels == level2focus));
    
    %
    firstChunk = [definiteMembers otherMembers];
    bundled = [];
    
else % k is a number between two numbers in cumsum(numLeafs)
    fprintf(['k is a number between ' num2str(max(find(numKatlevel<k))+1)  '\n']);
    
    level2focus = max(find(numKatlevel<k))+1;
    definiteMembers = BSI(find((BSIlevels < level2focus) & (numKids ==0)));
    candidateMembers = BSI(find(BSIlevels == level2focus));
    
    % parents of the candidate members
    for i=1:length(candidateMembers)
        parents(i) = treeStruct.numDescendants.getparent(candidateMembers(i));
    end
    
    
    uqParents = unique(parents);
    for i=1:length(uqParents)
        for j=1:max(find(parents==uqParents(i)))
            parentsIncluded = uqParents(i+1:end);
            candidateSets{i,j} = [candidateMembers(1:j) parentsIncluded];
            bundled_iterate{i,j} = setdiff(setdiff(candidateMembers,candidateMembers(1:j)),candidateMembers(ismember(parents, parentsIncluded)));
            
        end
    end
    
    % produce all combinations that gives this k number
    [ii,jj] = find((cellfun(@length,candidateSets) + ~cellfun(@isempty,bundled_iterate) + length(definiteMembers)) == k);
    
    
    % get all combinations
    for m = 1:length(ii)
    firstChunk{m} = [definiteMembers candidateSets{ii(m),jj(m)}];
    bundled{m} = bundled_iterate{ii(m), jj(m)};
    end
    
end

% for the case when there is combinations
if iscell(firstChunk)
    
    for m=1:length(firstChunk{m})
        
    i=0;
    
    if ~isempty(firstChunk{m})
        for i=1:length(firstChunk{m})
            groupMembers{m}{i} = treeStruct.groupMembers.Node{firstChunk{m}(i)};
            groupMembership{m}(treeStruct.groupMembers.Node{firstChunk{m}(i)}) = firstChunk{m}(i);
        end
    end
    
    i=i+1;
    
    if ~isempty(bundled{m})
        groupMembers{i}=[];
        for j=1:length(bundled{m})
            groupMembers{m}{i} = [groupMembers{i} treeStruct.groupMembers.Node{bundled{m}(j)}];
            
        end
        if isempty(firstChunk{m})
            groupMembership{m}(groupMembers{i}) = 1;
        else
            groupMembership{m}(groupMembers{i}) = max(firstChunk{m})+1;
        end
        
    end
    
    end
    
else % if k<first tier, k>all, or k== a particular level
    
    i=0;
    
    if ~isempty(firstChunk)
        for i=1:length(firstChunk)
            groupMembers{i} = treeStruct.groupMembers.Node{firstChunk(i)};
            groupMembership(treeStruct.groupMembers.Node{firstChunk(i)}) = firstChunk(i);
            nodeMembers{i} = firstChunk(i);
        end
    end
    
    i=i+1;
    
    if ~isempty(bundled)
        groupMembers{i}=[];
        for j=1:length(bundled)
            groupMembers{i} = [groupMembers{i} treeStruct.groupMembers.Node{bundled(j)}];
            nodeMembers{i} = [nodeMembers{i} bundled(j)];
        end
        %?
        if isempty(firstChunk)
            groupMembership(groupMembers{i}) = 1;
        else
            groupMembership(groupMembers{i}) = max(firstChunk)+1;
        end
        
    end
    
end

groupMembers = {groupMembers};
nodeMembers = {nodeMembers};

end