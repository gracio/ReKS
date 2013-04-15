function obj = chop(obj, node)
%% CHOP  Remove the target node and all subnodes from the given tree.

    iterator = obj.depthfirstiterator(node);
    obj.Parent(iterator) = [];
    obj.Node(iterator) = [];

end