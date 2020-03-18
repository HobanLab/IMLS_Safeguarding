# This code will recreate all plots and all analyses for the paper "The optimal size of an ex situ 
# conservation population: a comparison among 11 taxa in 5 genera"
# This project was made possible through the support of 
# the Institute of Museum and Library Services [Grants MA-05-12-0336-12, MA-30-14-0123-14, MA-30-18-0273-18, and MG-30-16-0085-16] 
# and the National Science Foundation (DEB 1050340, DBI 1203242 and DBI 1561346)
# Fieldwork was supported by 
# the Plant Exploration Fund, the Association of Zoological Horticulture, SOS—Save Our Species (Grant 2012A-035), 
# and the Mohamed bin Zayed Species Conservation Fund (Projects 0925331, 12254271, and 162512606). 
# Code primarily written by S. Hoban, with small contributions or comments from E Spence and E Schumacher
# Note: some of this code was originally developed for the Flowers et al 2018 white ash genetic study

######################################################################

# There are three sections in this file: Prep work, Part 1, and Part 2.  They are explained below
# To run this you will need a folder for each species and the files "Sp_total.gen" and "Sp_wild.gen"
# Within each folder, where Sp stands for species name, for example Qboyntonii or Mpyramidata
# Thus you have 11 folders (for 11 species) and withing each folder two gen data files
#
# SUMMARY PREP
# loads in libraries, defines species list (which is also the list of folders the data are found in)
# 
# SUMMARY PART 1: Code for testing different collection sizes (number of trees)
#
# There are several loops.  
# The outer loop is over the two options of what kind of rare alleles to consider 
# drop zero is to not drop any allles; drop 2 is to drop alleles present in two or less copies
# 	The next loop is over the number of species (e.g. eleven species
# 	Within this loop it loads in the data file (all in situ populations) and categorizes all alleles existing in the datasets
# 	e.g. as "global", "very common" etc.  The file does NOT include ex situ individuals, only in situ (wild)
#		 The next loop is the actual sampling, which is actually another nested loop
#		 There are a certain number of reps to include stochasticity (some thousands, usually)
#		 For each rep, there is a loop over all possible numbers of trees, randomly, up through the total population size
# 			For every sampling effort, for every rep, for every species, it records the 
# 			number of alleles captured in that set of trees
# 			This will later be converted to a proportion of genetic diversity and will be used for plotting all results,
# 			e.g how much genetic diversity is gained by adding more samples (separate .Rfile)
# There will be two output files- a "summ_results_tree.R".  and a "summ_results_tree_dr_0.R" file
# Each of these will be a 4 dimensional array with dimensions 
# 	(number of trees in the wild) x (number of allele categories+2) x (1) x (number of reps)
# Note there are 9 allele categories in this code even though only 5 will be used in this paper-
# the remaining four are for local and regional alleles that COULD be examined in future work
# Note also the first two columns in the array are blank; an example is below
# The columns are (blank) (blank) (global alleles captured) (very common) (common) (low freq) (rare)
# [2,]   NA   NA   36   24   32   12    0 	(ignore rest)
# [3,]   NA   NA   47   23   33   24    0 	(ignore rest)
#
# SUMMARY PART 2: Code for calculating current allele capture

# This has a similar loop structure over n_to_drop and species as noted in Part 1
# However it is reading in files that are all wild populations plus all ex situ individuals merged into one population
# And similar to above it categorizes alleles into categories
# And it is simply counting the number of alleles total and then the number captured in the ex situ population.  
# It records results in ex_situ_vs_in_situ.csv a single spreadsheet type file
#
# The code can be run using parallel processing or not
# It is currently set up to run parallel processing on 24 processors
# if NOT using this, just comment out the foreach and uncomment out the for loop code adjacent
#
####################################################################################################################


#####################
# PREP WORK
#####################

library(adegenet)
library(parallel);	library(doParallel) #will load foreach
library(abind)

#source("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Fa_sample_funcs.R")
source("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Fa_sample_funcs.R")
#source("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/Fa_sample_funcs.R")

species_names<-c("Hhannerae","Hwaimeae","Masheii","Mpyramidata","Pekmanii","Psargentii","Qboyntonii","Qgeorgiana","Qoglethorpensis","Zdecumbens","Zlucayana") 
region_makeup_list<-list(list(1:2,3:4),list(1:2,3),list(1,2),list(1:2,3:4),list(1,2),list(1:2,3),list(1:2,3,4,5:9),list(1:2,3,4,5:10),list(1:2,3,4,5:8),list(1,2),list(1:2,3))
alleles_existing_by_sp<-matrix(nrow=length(species_names),ncol=9)
	
	
##########################################################
#	PART 1- SAMPLING IN SITU POPS FOR NUMBER OF TREES
###########################################################

#This will run over a loop of "Include all alleles (n_to_drop=0)" and "Include only alleles present in more than two copies (n_to_drop=2)"
for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-""
	if (n_to_drop==0) n_drop_file<-"_dr_0"
	
####################################
	
for (sp in 1:length(species_names)){
	this_species<-species_names[sp]
#	setwd(paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
#	setwd(paste("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
 	setwd(paste("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
 	Spp_genind<-read.genepop(paste(substr(this_species,1,2),"_wild.gen",sep=""),ncode=3);	print(table(Spp_genind@pop))
	Spp_genpop<-genind2genpop(Spp_genind);	Spp_genind_sep<-seppop(Spp_genind)
 	max_num_trees<-max_num_trees_vect[sp]; region_makeup<-region_makeup_list[[sp]]
	#########################################3
	
	n_pops<-length(levels(Spp_genind@pop))
	n_total_indivs<- length(Spp_genind@tab[,1])
	n_ind_p_pop<-table(Spp_genind@pop)
	allele_freqs<-colSums(Spp_genpop@tab)/(n_total_indivs*2)
	num_reps<-75000
	
	list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
		
	#######################################################
	#----DETERMINE WHAT ALLELES FALL IN WHAT CATEGORIES---#
	#######################################################
	allele_cat<-get.allele.cat(Spp_genpop, region_makeup, 2, n_ind_p_pop,n_drop=n_to_drop)	
	#glob; 	glob_v_com; 	glob_com;	glob_lowfr;		glob_rare
	#reg_com_int;	loc_com_d1;		loc_com_d2;		loc_rare
	#Beware of NAs
	for (i in 1:9) alleles_existing_by_sp[sp,i]<- (sum((allele_cat[[i]])>0,na.rm=T))
	#!!!RESULT!!!- HOW MANY ALLELES IN EACH CATEGORY
	#proportion rare 
	#alleles_existing_by_sp[,5]/alleles_existing_by_sp[,1]
	
	
		cl <- makeCluster(24) # create a cluster with X cores
		registerDoParallel(cl) # register the cluster
	

		#########################################################
		#--JUST BY SAMPLING NUMBER OF TREES 				  --#
		#########################################################
			
	summ_results_tree<-array(dim=c(nrow(Spp_genind@tab)-1,11,num_scen,num_reps))
	#for (nrep in 1:num_reps) {		#if not parallel comment IN
	temp<-foreach(nrep=1:num_reps) %dopar% {	#if not parallel comment OUT
		alleles<-matrix(nrow=nrow(Spp_genind@tab)-1,ncol=length(allele_freqs))
		#For number of rows minus one		
		for (t in 2:(nrow(Spp_genind@tab)-1)){
			alleles<-colSums(Spp_genind@tab[sample(1:nrow(Spp_genind@tab), t+1),],na.rm=T)
			for (l in 1:length(allele_cat)) summ_results_tree[(t),(l+2),scen,nrep]<-sum(alleles[allele_cat[[l]]]>0, na.rm=T)
			}
			summ_results_tree[,,,nrep]	#if not parallel comment OUT
		}
	summ_results_tree[,,1,]<-abind(temp,along=3)	#if not parallel comment OUT
	save(summ_results_tree,file=paste("summ_results_tree",n_drop_file,".R",sep=""))

 stopCluster(cl)

#end of loop over species
}



		
#####################################
#	PART 2- COMPARE EX SITU AND IN SITU		#
#####################################

	#setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
	setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
	#setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
		
#This will run over a loop of "Include all alleles (n_to_drop=0)" and "Include only alleles present in more than two copies (n_to_drop=2)"
for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-""
	if (n_to_drop==0) n_drop_file<-"_dr_0"
	
	species_names<-c("Hhannerae","Hwaimeae","Masheii","Mpyramidata","Pekmanii","Psargentii","Qboyntonii","Qgeorgiana","Qoglethorpensis","Zdecumbens","Zlucayana") 
	region_makeup_list<-list(list(1:2,3:4),list(1:2,3),list(1,2),list(1:2,3:4),list(1,2),list(1:2,3),list(1:2,3,4,5:9),list(1:2,3,4,5:10),list(1:2,3,4,5:8),list(1,2),list(1:2,3))
	set_garden_p<-c(5,4,3,5,3,4,10,11,9,3,4) 
	wild_results<-matrix(nrow=length(species_names),ncol=9+1)
	alleles_existing_by_sp<-matrix(nrow=length(species_names),ncol=9)

	for (sp in 1:length(species_names)){
		this_species<-species_names[sp]
		setwd(paste("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
		#setwd(paste("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
		#setwd(paste("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/",this_species,sep=""))
		Spp_tot_genind<-read.genepop(paste(substr(this_species,1,2),"_total.gen",sep=""),ncode=3)

		#This code compares the wild to various ex situ populations or all the ex situ merged
		#Just put in garden and wild population numbers... currently all gardens are merged into one "population"
		wild_p<-unlist(region_makeup_list[[sp]]); garden_p<-set_garden_p[sp]
		n_ind_W<-table(Spp_tot_genind@pop)[wild_p];  n_ind_G<-table(Spp_tot_genind@pop)[garden_p]; 
		Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
		Spp_tot_genind_sep<-seppop(Spp_tot_genind)
		alleles_cap<-colSums(Spp_tot_genind_sep[[garden_p]]@tab,na.rm=T)

			#Allele categories based only on wild populations (can look at all wild pop'ns or only one if you want)
		allele_cat_tot<-get.allele.cat(Spp_tot_genpop[wild_p], region_makeup_list[[sp]], 2, n_ind_W, glob_only=F,n_drop=n_to_drop)
			#This goes through each allele category and divides the number captured ex situ (alleles_cap) by the number of alleles existing (allele_cat_tot)
			list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")
			
		for (i in 1:9) alleles_existing_by_sp[sp,i]<- (sum((allele_cat_tot[[i]])>0,na.rm=T))
		
		for (l in 1:length(allele_cat_tot)) wild_results[sp,l]<-round(sum(alleles_cap[allele_cat_tot[[l]]]>0)/length(allele_cat_tot[[l]]),4)
		
		wild_results[sp,10]<-n_ind_G
	}

	#setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
	setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
	#setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Oak_popgen_analyses/")
	wild_results<-cbind(species_names,wild_results)
	write.csv(wild_results,file=paste("ex_vs_in_situ",n_drop_file,".csv",sep=""))
	
} #close num to drop loop		
		


