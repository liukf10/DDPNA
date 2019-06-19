## DDPNA Version 0.1.0 
-   The first version of DDPNA package.

## DDPNA Version 0.1.1 
-   fix some small error
-   add wgcnatest function 

## DDPNA Version 0.2.1 
-   fix ID_match function and allow user run it without blast+ software.
-   in Sample_ID_data.rda file, add more data to use in example.
-   add example in MaxquantDataconvert.Rd and ID_match.Rd
-   add dataStatInf function which can summrize the statatic information of data.
-   add DEP_Mod_HeatMap function which can get the DEP enrich fold in Module and plot a HeatMap
-   add rename_dupnewID function which can rename the duplicated newID.
-   add DEPsets function which get two or more IDsets interesection set and complementary set and give the colors.
-   add DEP_Mod_net_plot function which remove hubs which not in the IDsets and replot the PFG network

## DDPNA Version 0.2.1
-   fix a bug in Data_impute function, previous function will error when doclust.
