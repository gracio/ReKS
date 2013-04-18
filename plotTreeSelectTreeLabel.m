function [] = plotTreeSelectTreeLabel(treeStruct, node2display, varargin)
% plot the tree given in treeStruct, displaying only labels of
% node2display. optional paramters is "contents" where one specifies
% exactly what is being displayed at each node.

treeCopy = treeStruct;
allNodes = treeCopy.breadthfirstiterator;
allNodes = setdiff(allNodes,node2display);

for i=1:length(allNodes)
    treeCopy = treeCopy.set(allNodes(i), []); % don't want to display everything else
end

% if user specifies what contents to be displayed, change to that
if ~isempty(varargin)
    for i=1:length(node2display)
        treeCopy = treeCopy.set(node2display(i),varargin{1}{i});
    end
end

figure
treeCopy.plot2