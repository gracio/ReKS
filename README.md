Huang GT, Cunningham KI, Benos PV, Chennubhotla CS. Spectral clustering strategies for heterogeneous disease expression data. 
Pac Symp Biocomput. 2013:212-23

The main file takes a data structure with fields 'data','genes' and 'headers' as input. An example input file is provided (BR_Normal.mat). Future release will allow a tab delimited file as input.

Set path to  the TreeStructure class included, courtesy of http://www.mathworks.com/matlabcentral/fileexchange/35623-tree-data-structure-as-a-matlab-class. BSD license provided.

To run ReKS using the example .mat file, type:

>>load BR_Normal.mat;

>>clusterMembership = ReKS_main(BR_Tumor);

Future release will include visualization options and data export options.  


