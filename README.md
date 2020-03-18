# This code will recreate all plots and all analyses for the paper "The optimal size of an ex situ conservation population: a comparison among 11 taxa in 5 genera"
This project was made possible through the support of 
the Institute of Museum and Library Services [Grants MA-05-12-0336-12, MA-30-14-0123-14, MA-30-18-0273-18, and MG-30-16-0085-16] 
and the National Science Foundation (DEB 1050340, DBI 1203242 and DBI 1561346)
Fieldwork was supported by 
the Plant Exploration Fund, the Association of Zoological Horticulture, SOS—Save Our Species (Grant 2012A-035), 
and the Mohamed bin Zayed Species Conservation Fund (Projects 0925331, 12254271, and 162512606). 
Code primarily written by S. Hoban, with small contributions or comments from E Spence and E Schumacher
Note: some of this code was originally developed for the Flowers et al 2018 white ash genetic study

In this folder you will find data and files for recreating results in the paper "The optimal size of an ex situ conservation population: a comparison among 11 taxa in 5 genera".  Note to reviewers.  Only five species are included in this zip file.  The other files will be submitted to Dryad upon paper acceptance.  These will be requested to be on embargo for one year.  Nonetheless, all code should run for these five species with very minimal adjustment

The data files are of the format "Gs_total.gen" and "Gs_wild.gen" where G is the genus initial and s is the taxa initial.  For example "Qb_total.gen" for Quercus boyntonii.  All files are in genpop format (thus are delimited by POP).  "total" is all wild populations (each one as a separate pop) plus a population representing all ex situ accessions pooled together.  "wild" is all wild populations only.  Each data file is in a folder for that species (necessary because the sampling code will output a file that is named the same for all species so it has to be in its own folder to avoid being wiped each time as it loops over species.

Run the sampling code and ex situ vs. in situ comparison by running the Spp_sampling_for_pub.R file
Once this is complete, run the plotting code by running Spp_plotting.R
