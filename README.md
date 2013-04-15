Huang GT, Cunningham KI, Benos PV, Chennubhotla CS. Spectral clustering strategies for heterogeneous disease expression data. Pac Symp Biocomput. 2013:212-23.\n
The main file takes a data structure with fields 'data','genes' and 'headers' as input. \n
You will need to set path to include the TreeStructure class.\n
To run using the BR_Normal example \n

>> load BR_Normal.mat;\n
>> clusterMembership = ReKS_main(BR_Tumor);\n
