## Resubmission
This is a resubmission. In this version I have:
*In in MaxQdataconvert.Rd, ID_match.Rd example, use tempdir() instead of getwd().
* Now we use ggfortify package to replace ggbiplot in modcomp function and delete the Vincent Q. Vu who have write ggbiplot package in github in Authors field. The previous version, we copy the ggbiplot function by reason of  CRAN not support install package from GitHub.
* Remove or replaced all cat() in code and remain the print() in code because it only for plotting.
* Instead of \dontrun{}, use \donttest{} in wgcnatest.Rd, DEP_Mod_net_plot.Rd, DEP_Mod_HeatMap.Rd, because these example take long time to run.
* Instead of \dontrun{}, use \donttest{} in MaxQdataconvert.Rd, ID_match.Rd, because the function can run without extra software but it will take long time to run in example.

Previous changed information
* delete “rm(list = ls())” in all example
* Add new function: dataStatInf,  DEP_Mod_HeatMap, rename_dupnewID, DEPsets, DEP_Mod_net_plot in this version.
* fix ID_match function and allow run it without extra software.
* add (old.wd <- getwd(); on.exit(setwd(old.wd))) inside the ID_match function to reset the default working directory. 
* Add single quotes in package names, software names and API names in title and description field.
* Add reference in DESCRIPTION file.
* Deleted extra directory and files.
* fix some words in DESCRIPTION and the check result  "et al" is abbreviation, which are not mis-spelled words.
* fix the grammatical mistake
* Add a new function wgcnatest into the package.

## Test environments
* local Win x86_64-w64-mingw32  R 3.5.1
* x86_64-w64-mingw32 R3.6.0 (win-builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC(on Rhub)
* Fedora Linux, R-devel, clang, gfortran(on Rhub)

## R CMD check results
There were no ERRORs or WARNINGs

There was 1 NOTE:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kefu Liu <liukefu19@163.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  al (15:25)
  et (15:22)
  

