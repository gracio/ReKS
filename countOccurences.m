function [uniqueItems,occ] = countOccurences(vec)
% vec is a vector with some number of unique items with different
% occurences 

tp = sort(vec);
[uniqueItems,Ifirst,J] = unique(tp,'first');
[uniqueItems,Ilast,J] = unique(tp,'last');
occ = (Ilast-Ifirst + 1);