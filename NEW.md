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

## DDPNA Version 0.2.2
-   fix a bug in Data_impute function, previous function will error when doclust.
-   fix a bug in dataStatInf function when ID have "-". 190716
-   fix some bugs in Data_impute function. 190716
-   fix the order of the module in Module_Enrich function. 191105
- allow switch x and y axis in DEP_Mod_HeatMap function. 191118 

## DDPNA Version 0.2.3
-   fix the format of data "Sample_ID_data".200212

## DDPNA Version 0.2.4
-   fix some bugs in getmoduleHub and DEP_Mod_net_plot function.200327
-   fix example in DEP_Mod_net_plot.Rd which use T/F instead of TRUE/FALSE. 200327 

## DDPNA Version 0.2.5
-   fix bug in ID_match.20200602
-   add correlation P value adjust in getmoduleHub function. 20200610
-   fix a donntest error in IDmatch function. 20200626
-   fix example error in MaxQdataconvert function. 20200626
-   fix a problem of long system runtime in example of changedID and FCSenrichplot function. 20200626

## DDPNA Version 0.2.6
-   remove creat Plot Directories when output image.20200703

## DDPNA Version 0.2.7
-   change ftp website in ID_match.Rd and MaxQdataconvert.Rd. 20210422

## DDPNA Version 0.2.8
-   fix ID_match function occur error when contains metacharacters in GN. 20210528

## DDPNA Version 0.3.0
-   fix Data_impute function to save plot file in plot dir and add bicor method in outlier detect. 20210615,20210628,20210728
