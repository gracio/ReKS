Huang GT, Charalambos A, Benos P. mirConnX: Condition-specific mRNA-microRNA network integrator. http://mirconnx.csb.pitt.edu/ Nucleic Acids Res. 2011 May 10; 1-8\n
The main file takes a data structure with fields 'data','genes' and 'headers' as input. \n
You will need to set path to include the TreeStructure class.\n
To run using the BR_Normal example \n

>> load BR_Normal.mat;\n
>> clusterMembership = ReKS_main(BR_Tumor);\n
