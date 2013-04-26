function [] = plotNodeAffinity(node2plot,treeStruct)


dummy = treeStruct.discComp.Node{node2plot};
order = [];
orderLab = {};
for i=1:size(dummy,2)
    order = [order; find(dummy(:,i))];
    orderLab = [orderLab; repmat({num2str(i)},sum(dummy(:,i)),1)];
end

figure;
imagesc(treeStruct.A.Node{node2plot}(order,order));

d=treeStruct.groupMembers.Node{node2plot}(order);

for i=1:length(d)
ddy{i} = [num2str(d(i)) ' (' orderLab{i} ')'];
ddx{i} = num2str(d(i));
end
set(gca,'XTick',1:length(d))
set(gca,'XTickLabel',ddx)
set(gca,'YTick',1:length(d))
set(gca,'YTickLabel',ddy)
title(['Affinity Map of Node ' num2str(node2plot)])