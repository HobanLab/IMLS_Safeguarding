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

###############################################################################
# The code is arranged to first produce figures then run the stats
# The figures and tables are: 
# Table 1: Attributes of each species- composed by hand, NOT in the code
# Table 2: genetic diversity actually captured in the botanic garden collections today along with collection size
#		This is actually crated in the "sampling" R file - look for the "ex_vs_in_situ",n_drop_file,".csv" 
#		To actually make the table in the word doc did require moving numbers from the dr0 file into a table
#		from the non dr0 file in Excel, using parentheses for those dr0 numbers
# Figure 1: Genetic diversity in current collections vs. size of current collections (one point per species)
#		This is created in the code in this document- it reads in subsets of columns in Table 2
# Figure 2: genetic diversity captured on number of samples (from simulated sampling)
#		This is created in the code in this document, it uses the summ_results_tree file
# Table 3: collector gap
#		This is created in the code in this document, it uses BOTH of the aforementioned docs (ex_situ_vs and summ_results)
# The stats are:
# Stats: ANOVA of genetic diversity captured and number of samples needed ~ genus 
# Stats: regression of genetic diversity captured and number of samples needed ~ allele frequency summaries
# Stats: regression of genetic diversity captured and number of samples needed ~ FST summaries 
#
#
# Supplemental
#
# It is assumed you have already run the "sampling" code file and thus you have
# files summ_results_tree and ex_situ_vs_in_situ as output of that code 


#The following should always be run first- setwd, load libraries, set folder names and plot colors
library(adegenet);	library(hierfstat)

#Get working directory (only for sean)
which_comp<-getwd()
	if (grepl("shoban",which_comp)) prefix_wd<-"C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
	if (grepl("shoban.DE",which_comp)) prefix_wd<-"C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
	if (grepl("home",which_comp)) prefix_wd<-"/home/user/Dropbox/Projects/IN_PROGRESS/IMLS_synthesis_analysis/"
setwd(prefix_wd)
	
folder_names<-c("Hhannerae","Hwaimeae","Masheii","Mpyramidata","Pekmanii","Psargentii","Qboyntonii","Qgeorgiana","Qoglethorpensis","Zdecumbens","Zlucayana") 
species_names<-c("Hibiscus w. subsp hannerae","Hibiscus w. subsp waimeae","Magnolia ashei","Magnolia pyramidata","Pseudopheonix ekmanii","Pseudopheonix sargentii","Quercus boyntonii","Quercus georgiana","Quercus oglethorpensis", "Zamia decumbens","Zamia lucayana")

sp_colors<-c("palevioletred1","palevioletred1","limegreen","limegreen","mediumorchid3","mediumorchid3","dodgerblue3","dodgerblue3","dodgerblue3","orange","orange")	
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare")

num_reps<-75000


#####################################
#	FIGURE 1- COMPARE EX SITU AND IN SITU		
#	This figure is a simple plot of the current conservation of genetic diversity in the current collection sizes
#	One point per species, and a log relationship across them
#	The first chunk of code creates plots of all allele types and over Reduced and Full Datasets
#	The second chunk of code is only the Reduced and only the two alleles shown in figure in paper 
#####################################

for (n_to_drop in c(0,2)){
	if (n_to_drop==2) n_drop_file<-"";		if (n_to_drop==0) n_drop_file<-"_dr_0"
	
	species_names<-folder_names
	
	wild_results<-read.csv(paste("ex_vs_in_situ",n_drop_file,".csv",sep=""))[,-1]
	setwd(paste(prefix_wd,"/Fig1_ex_vs_in/",sep=""))
	pdf(file=paste("captured_ex_situ",n_drop_file,".pdf",sep=""),width=5,height=5)	
	for (i in 2:6){
	plot(wild_results[,11],as.numeric(wild_results[,i])*100,ylim=c(0,100), main="",xlab="number of plants in ex situ collections",ylab="percentage of alleles captured")
	text(as.numeric(wild_results[,11]),as.numeric(wild_results[,i])*100+5,substr(wild_results[,1],1,2),cex=1,col=sp_colors)
	#Fitting data		
		d <- as.data.frame(list(plants=wild_results[,11], gendiv=wild_results[,i]*100))
		#d<-rbind(c(0,0),d)
		mod <- lm(gendiv ~ I(log(plants)), d)		#Note: a square root relationship was tested and is ok but not as good as log
		new.df <- data.frame(plants=seq(10,250,by=1))
		out <- predict(mod, newdata = new.df)
		#plot(d,pch = 16, xlim = c(0,250), ylim = c(0,100))
		lines(unlist(c(new.df)), out, col = "grey", lwd=1.5)
		text(193,10,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black")
	}
	dev.off();		setwd("..")
	#interesting no global rare alleles for P ekmanii (considering the min freq)
} #close num to drop loop		

wild_results<-read.csv("ex_vs_in_situ.csv")[,-1]
pdf(file="captured_ex_situ_for_pub.pdf",width=9,height=5)	
par(mfrow=c(1,2))
for (i in c(2,5)){
	if (i==2) y_text<-"percentage of 'all alleles' captured"; if (i==5) y_text<-"percentage of 'low frequency alleles' captured"
	plot(wild_results[,11],as.numeric(wild_results[,i])*100,ylim=c(0,100), main="",xlab="",ylab=y_text)
	mtext(side=1,line=-2,"number of plants living in ex situ collections", outer=T)
	text(as.numeric(wild_results[,11]),as.numeric(wild_results[,i])*100+5,substr(wild_results[,1],1,2),cex=1,col=sp_colors)
		d <- as.data.frame(list(plants=wild_results[,11], gendiv=wild_results[,i]*100))
		mod <- lm(gendiv ~ I(log(plants)), d)		#Note: a square root relationship was tested and is ok but not as good as log
		new.df <- data.frame(plants=seq(10,250,by=1))
		out <- predict(mod, newdata = new.df)
		lines(unlist(c(new.df)), out, col = "grey", lwd=1.5)
		text(193,10,paste("adj R2 =",round(as.numeric(summary(mod)[9]),2)),col="black")
	}	
dev.off()
setwd("..")




##############################################################	
#		Figure 2: OVERLAY PLOTS OF SIMULATED SAMPLING, MULTIPLE SPECIES BY ALLELE
#
# This creates a single figure used in the paper, with two subpanels
# It shows the genetic diversity accumulation for all alleles
# Focusing on a 95% threshold target
# Panel A is Reduced Dataset and Panel B is Full Dataset
# Thus there is a loop over n_to_drop
############################################################

#NOTE COULD ADD GREY LINE TO THIS PLOT TO SHOW HOW LOW IT IS

min_thresh<-95; y_lower<-90	# y_lower is the lower y axis limit for the plot
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare")

setwd(paste(prefix_wd,"/Fig2_simulated_sampling_overlay_threshold/",sep=""))

pdf(file=paste("bytrees_overlay_merged_t_",min_thresh,".pdf",sep=""),width=10, height=12)
par(mfrow=c(2,1))	
for (n_to_drop in c(2,0)){
	if (n_to_drop==0) { n_drop_file<-"_dr_0" ; main_title<-"(B) Full Dataset"}
	if (n_to_drop==2) { n_drop_file<-""; main_title<-"(A) Reduced Dataset"}
	
	#calc_min will be a matrix to hold numbers for the minimum collection size thresholds for "sufficiency"
	calc_min<-matrix(nrow=length(species_names),ncol=length(list_allele_cat))
	rownames(calc_min)<-c("H. w. hannerae","H. w. waimeae","M. ashei","M. pyramidata","P. ekmanii","P. sargentii","Q. boyntonii","Q. georgiana","Q. oglethorpensis","Z. decumbens","Z. lucayana")
	colnames(calc_min)<-list_allele_cat

		i<-3	#all alleles (global) only
		sp_colors<-c("palevioletred1","palevioletred1","limegreen","limegreen","mediumorchid3","mediumorchid3","dodgerblue3","dodgerblue3","dodgerblue3","orange","orange")	#note this needs to be here because it is re-sorted below and needs to be reset

		for (sp in 1:length(folder_names)){

			this_species<-folder_names[sp]
			setwd(paste(prefix_wd,this_species,sep=""))
			load(file=paste("summ_results_tree",n_drop_file,".R",sep=""));	load(file="params_run.Rda")

			#optional if proportions- CAN DELETE THIS?
			for (n in 1:num_reps) summ_results_tree[,i,1,n]<-t(t(summ_results_tree[,i,1,n])/summ_results_tree[length(summ_results_tree[,1,1,1]),i,1,n])	
			#mean across reps
			all_mean<-apply(summ_results_tree[,i,1,1:num_reps],1,mean)*100

			if (sp==1) {
			plot(all_mean,ylim=c(y_lower,100),type="l",lwd=4,xlim=c(0,250),xlab="number of plants (using simulated sampling)",
				ylab="percentage of genetic variation", cex.lab=1.2,col=sp_colors[sp],main=main_title)
			}	else {lines(all_mean,lwd=4,col=sp_colors[sp])}
			calc_min[sp,i-2]<-min(which(all_mean>min_thresh))
		}	#end species loop
		
		#now just add text in rectangle for the figure legends with the numbers crossing 95% and line showing that threshold
		min_to_plot<-sort(calc_min[,i-2],decreasing=T)
		sp_colors<-sp_colors[order(calc_min[,i-2], decreasing=T)]
		abline(h=95, lty=2, col="darkgrey")
		rect(158,89.9,240,98.8,col=rgb(1,1,1,alpha=0.89)) ; text(200,98.5,"minimum to capture 95% of alleles")
		for (sp in 1:length(folder_names)) text(200,(89.6+.75*sp),paste(names(min_to_plot)[sp],": ",min_to_plot[sp]),col=sp_colors[sp])
		
}	#end to drop loop
dev.off();		setwd("..")


#############################################################
#		OVERLAY PLOTS OF SIMULATED SAMPLING MULTIPLE SPECIES BY ALLELE			#
#
# This will create four PDFs actually: each combination of threshold (70, 95%) and allele dropping choice (0 or 2)
# And it will go through all the allele types
# Thus it has more nested loops (n_to_drop, threshold, and allele category)
# This PDF is NOT in the publication, but is rather an expanded version of the Figure in the publication, with all results
# Some of which are in the Supplemental
############################################################

setwd(paste(prefix_wd,"/Fig2_simulated_sampling_overlay_threshold/",sep=""))
for (n_to_drop in c(0,2)){
	if (n_to_drop==0) n_drop_file<-"_dr_0" 
	if (n_to_drop==2) n_drop_file<-""
	
	#thresholds for "sufficiency"
	for (min_thresh in c(70,95)){	

	y_lower<-65; if (min_thresh==95) y_lower<-90

	pdf(file=paste("bytrees_overlay",n_drop_file,"_t_",min_thresh,".pdf",sep=""),width=10, height=6)
	calc_min<-matrix(nrow=length(species_names),ncol=length(list_allele_cat))
	rownames(calc_min)<-c("H. w. hannerae","H. w. waimeae","M. ashei","M. pyramidata","P. ekmanii","P. sargentii","Q. boyntonii","Q. georgiana","Q. oglethorpensis","Z. decumbens","Z. lucayana")
	colnames(calc_min)<-list_allele_cat

	for (i in c(3,5,6,7,8)){

		if (i==5|i==8) {x_ax_lim<-100; text_x<-80} else {x_ax_lim<-250; text_x<-200}
		sp_colors<-c("palevioletred1","palevioletred1","limegreen","limegreen","mediumorchid3","mediumorchid3","dodgerblue3","dodgerblue3","dodgerblue3","orange","orange")	

		for (sp in 1:length(folder_names)){

			this_species<-folder_names[sp]
			
			load(file=paste(prefix_wd,this_species,"/summ_results_tree",n_drop_file,".R",sep=""))

			list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_comm","loc_com_d1","loc_com_d2","loc_rare")

			#optional if proportions
			for (n in 1:num_reps) summ_results_tree[,i,1,n]<-t(t(summ_results_tree[,i,1,n])/summ_results_tree[length(summ_results_tree[,1,1,1]),i,1,n])	
			#mean across reps
			all_mean<-apply(summ_results_tree[,i,1,1:num_reps],1,mean)*100

			if (sp==1) {
			plot(all_mean,ylim=c(y_lower,100),type="l",lwd=4,xlim=c(0,x_ax_lim),xlab="number of plants",ylab="percentage of genetic variation",
				cex.lab=1.2,col=sp_colors[sp],main=paste(list_allele_cat[i-2]," alleleles"))
			}	else {lines(all_mean,lwd=4,col=sp_colors[sp])}
			calc_min[sp,i-2]<-min(which(all_mean>min_thresh))
		}	#end species loop
		
		#now just add text in rectangle and lines
		min_to_plot<-sort(calc_min[,i-2], decreasing=T)
		sp_colors<-sp_colors[order(calc_min[,i-2],decreasing=T)]
			if (min_thresh==70) {
				abline(h=70, lty=2, col="darkgrey")
				rect(158,65.9,240,97.5,col=rgb(1,1,1,alpha=0.89)) ; text(text_x,96.0,"minimum to capture 70% of alleles")
				}
			if (min_thresh==95) {
				abline(h=95, lty=2, col="darkgrey")
				rect(158,89.9,240,98.8,col=rgb(1,1,1,alpha=0.89)) ; text(text_x,98.5,"minimum to capture 95% of alleles")
				}
		for (sp in 1:length(folder_names)) {
			if (min_thresh==70) text(text_x,(65.6+2.5*sp),paste(names(min_to_plot)[sp],": ",min_to_plot[sp]),col=sp_colors[sp])
			if (min_thresh==95) text(text_x,(89.6+.75*sp),paste(names(min_to_plot)[sp],": ",min_to_plot[sp]),col=sp_colors[sp])
		}
		
		write.csv(calc_min,file=paste("min_for_",min_thresh,n_drop_file,".csv",sep=""))	#make sure this works
	}	#end allele loop
	
	dev.off()

	}	#end min thresh loop

}	#end to drop loop
setwd("..")




#####################
##	TABLE 3- COLELCTOR GAP
##  This will do two calculations: 
##	A- the number of samples to get the same % we have now (reduce sample size)
##	B- the % we could get with the samples we have now (increase genetic capture)
#####################

	allele<-3

for (n_to_drop in c(0,2)){
	if (n_to_drop==0) n_drop_file<-"_dr_0" 
	if (n_to_drop==2) n_drop_file<-""
	min_p_needed<-vector(length=11)	#min needed is the amount that current collections could be reduced to without losing diversity 

	setwd(prefix_wd)
	curr_ex_situ<-read.csv(file=paste0("ex_vs_in_situ",n_drop_file,".csv"))
	collector_gap_res<-matrix(nrow=11,ncol=7)

	for (sp in 1:length(folder_names)){

		this_species<-folder_names[sp]
		load(file=paste(prefix_wd,this_species,"/summ_results_tree",n_drop_file,".R",sep=""))

	#part A ideal number of plants
	#use the ex situ vs in situ table (read in above) and first take note of the current % conserved e.g. 71%, 95% (gd_current)
	#then we want to find the number of samples needed to reach that amount of diversity under ideal sampling (use the summ_results) -> min_p_needed
	#later we want to compare this ideal sampling minimum to the current sample size  
		
		for (n in 1:num_reps) summ_results_tree[,allele,1,n]<-t(t(summ_results_tree[,allele,1,n])/summ_results_tree[length(summ_results_tree[,1,1,1]),allele,1,n])
		all_mean<-apply(summ_results_tree[,allele,1,1:num_reps],1,mean) 		#mean across reps
		
		#take the first row that exceeds the currently protected percentage
		gd_current<-curr_ex_situ[sp,allele]
		min_p_needed[sp]<-which(all_mean>gd_current)[1]
		
	#part B ideal genetic capture 
	#take the current number of samples- what is in ex situ collections 
	#we want to find what could ideally be conserved for that number- so simply get that row number in summ_results 

		plants_current<-curr_ex_situ[sp,12]		#the number of samples ex situ
		#go to that line in the in situ results table 
		#to find the % captured, under ideal sampling, for the same number of plants currently held
		gd_ideal<-all_mean[plants_current+1]
		
		collector_gap_res[sp,]<-c(folder_names[sp],length(summ_results_tree[,1,1,1]), plants_current, 
			min_p_needed[sp]-plants_current, min_p_needed[sp]/plants_current,
			gd_ideal-gd_current, gd_ideal/gd_current)	#changed 8 March for Paper Revision, Table 3
	}
collector_gap_res<-cbind(folder_names,collector_gap_res)	#make sure this works
colnames(collector_gap_res)<-c("species", "species", "N ex situ samples", "N in situ garden samples", "reduction num garden plants", "proportional reduction garden", "diff genetic diversity from ideal", "increase in gen diversity")
write.csv(collector_gap_res,file=paste0("Table3_collector_gaps",n_drop_file,".csv"))

}


###################################
#	STATS: ANOVAS ON GENUS				
#  Here we are testing the influence of genus on three response variables
#  CURRENT EX SITU VS IN SITU: Current genetic diversity conserved ex situ in collections today
#		This is 8 ANOVAs (2 files to analyze and 4 allele categories)
#  MINIMUM TO CATCH: Minimum sampling needed to achieve a sufficiency threshold (e.g. 70%, 95%)
#		This is 16 ANOVAs (4 files to analyze and 4 allele categories)
#  COLLECTOR GAP: The collector gap in terms of reduction in collection size and improved collection capture 
#		This is 4 ANOVAs (2 files to analyze and 2 ways of viewing the collector gap)
#####################################
setwd(prefix_wd)
library(broom)
extr_r_sq<-function(aov){	#code to extract R squared value
	tidy_aov<-tidy(aov)
	sum_squares_regression <- tidy_aov$sumsq[1]; sum_squares_residuals <- tidy_aov$sumsq[2]
	R_squared <- sum_squares_regression / (sum_squares_regression + sum_squares_residuals) 
	R_squared
}

	#	CURRENT EX SITU (VS IN SITU)	#
files_for_minimum<-c("ex_vs_in_situ.csv","ex_vs_in_situ_dr_0.csv")
pval_ex_situ<-matrix(ncol=5,nrow=2); rsq_ex_situ<-matrix(ncol=5,nrow=2) 
rownames(pval_ex_situ)<-files_for_minimum; colnames(pval_ex_situ)<-c("genus","all","comm","low_freq","rare")
for(f in 1:2){
	curr_ex_situ<-read.csv(file=files_for_minimum[f])[,-1]
	curr_ex_situ[,1]<-substr(curr_ex_situ[,1],0,1)		#replace first column with genus letter (for ANOVA)
	curr_ex_situ<-curr_ex_situ[,c(1,2,4,5,6)]
	curr_ex_situ[curr_ex_situ=="NaN"]<-NA;  colnames(curr_ex_situ)<-c("genus","all","comm","low_freq","rare")
	colnames(curr_ex_situ)[1]<-"genus"
	for (i in 1:4) pval_ex_situ[f,i+1]<- as.numeric(unlist(summary(aov(curr_ex_situ[,i+1]~genus,data=curr_ex_situ))[[1]][5]))[1]
	for (i in 1:4) rsq_ex_situ[f,i+1]<- extr_r_sq(aov(curr_ex_situ[,i+1]~genus,data=curr_ex_situ))
}		

	#	MINIMUM_TO_CATCH	#
files_for_minimum<-c("70.csv","95.csv","70_dr_0.csv","95_dr_0.csv")
pval_min_to_catch<-matrix(ncol=5,nrow=4); rsq_min_to_catch<-matrix(ncol=5,nrow=4) 
rownames(pval_min_to_catch)<-files_for_minimum; colnames(pval_min_to_catch)<-c("genus","all","comm","low_freq","rare")
for(f in 1:4){
	min_to_catch<-read.csv(paste0("min_for_",files_for_minimum[f]))
	min_to_catch[,1]<-substr(min_to_catch[,1],0,1)		#replace first column with genus letter
	min_to_catch<-min_to_catch[,c(1,2,4,5,6)]
	min_to_catch[min_to_catch==Inf]<-NA; ; colnames(min_to_catch)<-c("genus","all","comm","low_freq","rare")
	colnames(min_to_catch)[1]<-"genus"
	for (i in 1:4) pval_min_to_catch[f,i+1]<- as.numeric(unlist(summary(aov(min_to_catch[,i+1]~genus,data=min_to_catch))[[1]][5]))[1]
	for (i in 1:4) rsq_min_to_catch[f,i+1]<- extr_r_sq(aov(min_to_catch[,i+1]~genus,data=min_to_catch))
}

	#	COLLECTOR GAP	#
pval_coll_gap<-matrix(ncol=2,nrow=2); rsq_coll_gap<-matrix(ncol=2,nrow=2)
rownames(pval_coll_gap)<-c("dr0","dr2"); colnames(pval_coll_gap)<-c("times_reduce","times_increase")
for (n_to_drop in c(0,2)){
	if (n_to_drop==0) {n_drop_file<-"_dr_0"; which_row<-1}
	if (n_to_drop==2) {n_drop_file<-""; which_row<-2}
	collector_gap_res<-read.csv(file=paste0("collector_gaps",n_drop_file,".csv"))
	collector_gap_res<-collector_gap_res[,-1]
	collector_gap_res[,1]<-substr(collector_gap_res[,1],0,1)		#replace first column with genus letter
	colnames(collector_gap_res)<-c("genus","n_insitu","n_exsitu","n_reduce","n_times_reduce","percent_inc","times_increase")
	pval_coll_gap[[which_row,1]]<-summary(aov(n_times_reduce~genus,data=collector_gap_res))[[1]][5][[1]][1]
	pval_coll_gap[[which_row,2]]<-summary(aov(times_increase~genus,data=collector_gap_res))[[1]][5][[1]][1]
	rsq_coll_gap[[which_row,1]]<-extr_r_sq(aov(n_times_reduce~genus,data=collector_gap_res))
	rsq_coll_gap[[which_row,2]]<-extr_r_sq(aov(times_increase~genus,data=collector_gap_res))
}

write.csv(rbind(round(pval_ex_situ,3),
	round(matrix(p.adjust(pval_ex_situ,method="BH"),
	dimnames=list(paste0("adj.",rownames(pval_ex_situ)),colnames(pval_ex_situ)),ncol=ncol(pval_ex_situ)),3)),
	file="p_values/pvals_AOV_ex_situ_current.csv")
write.csv(rbind(round(pval_min_to_catch,3),
	round(matrix(p.adjust(pval_min_to_catch,method="BH"),
	dimnames=list(paste0("adj.",rownames(pval_min_to_catch)),colnames(pval_min_to_catch)),ncol=ncol(pval_min_to_catch)),3)),
	file="p_values/pvals_AOV_min_to_reach_thresh.csv")
write.csv(rbind(round(pval_coll_gap,3),
	round(matrix(p.adjust(pval_coll_gap,method="BH"),
	dimnames=list(paste0("adj.",rownames(pval_coll_gap)),colnames(pval_coll_gap)),ncol=ncol(pval_coll_gap)),3)),
	file="p_values/pvals_AOV_collector_gap.csv")



##############################################################################
#	STATS: REGRESSIONS WITH FST 				
#	Here we test the influence of the mean, SD, max and min of each species pairwise FSTs on the genetic capture etc.
#	FST is a quantitative, continuous variable so we use linear regression 
#	As above, we have three sets of tests
#	CURRENT EX SITU VS IN SITU: Current genetic diversity conserved ex situ in collections today
#		This is 32 regressions (2 files to analyze (Full/ Red), 4 allele categories, and 4 summaries of FST)
#  MINIMUM TO CATCH: Minimum sampling needed to achieve a sufficiency threshold (e.g. 70%, 95%)
#		This is 64 regressions (2 files to analyze (Full/ Red), 4 allele categories, 2 min thresholds and 4 summaries of FST)
#  COLLECTOR GAP: The collector gap in terms of reduction in collection size and improved collection capture 
#		This is 16 regressions (2 files to analyze (Full/ Red), 2 types of collector gap, and 4 summaries of FST)
###############################################################################


setwd(prefix_wd)
#Places to put FST summary statistics
sp_summ_fst<- matrix(nrow=length(folder_names),ncol=4); rownames(sp_summ_fst)<-folder_names
colnames(sp_summ_fst)<-c("PwFst","SD","Max","Min")

##generate Fst values and input into a table
for (sp in 1:length(folder_names)){
  this_species<-folder_names[sp];	 setwd(paste("./",this_species,sep=""))
  Spp_genind<-read.genepop(paste(substr(this_species,1,2),"_wild.gen",sep=""),ncode=3);	print(table(Spp_genind@pop))
  sm_fst<-as.matrix(pairwise.fst(Spp_genind));  sm_fst[sm_fst==0]<-NA
  Qsp_pwfst<-mean(sm_fst,na.rm=T)	
  print(this_species); print(sm_fst)
  sp_summ_fst[sp,1:4]<-c(mean(Qsp_pwfst), sd(sm_fst, na.rm = T), max(sm_fst, na.rm = T), min(sm_fst, na.rm = T)) 
  setwd("..")			 #Done- double checked FSTs compare to magnolia ms (ashei) and oaks paper
 }
 

 par(mfrow=c(3,3))
results_min<-matrix(nrow=1,ncol=4)	#all have four columns because that is the number of summary stats for FST (mean,sd,max,min)
results_capture<-matrix(nrow=1,ncol=4)
results_gap<-matrix(nrow=1,ncol=4)

for (n_to_drop in c(0,2)){
	if (n_to_drop==0) n_drop_file<-"_dr_0" 
	if (n_to_drop==2) n_drop_file<-""
	
	#regressions of pairwise Fst with current ex situ genetic capture 
	gen_cap<-read.csv(paste0("ex_vs_in_situ",n_drop_file,".csv"))[,-1];	gen_cap[gen_cap=="Inf"]<-NA
	reg_pval<-matrix(nrow=5, ncol=4); reg_r2<-matrix(nrow=5, ncol=4)

	for (ss in 1:4){  #this is a loop over min, max, mean and sd
	  for (at in c(2,4:6)){   #this is the loop over each allele type
		reg_pval[at-1,ss]<-((summary(lm(gen_cap[,at] ~ sp_summ_fst[,ss]))$coefficients[,4]))[2]
		reg_r2[at-1,ss]<-summary(lm(gen_cap[,at] ~ sp_summ_fst[,ss]))$r.squared
	  }
	}
	colnames(reg_pval)<-colnames(sp_summ_fst); colnames(reg_r2)<-colnames(sp_summ_fst)
	reg_pval<-reg_pval[-2,]; reg_2<-reg_r2[-2,]
	rownames(reg_pval)<-paste(c("global", "common", "low freq", "rare"),n_drop_file)
	results_capture<-rbind(results_capture,reg_pval)

	##regression of pairwise Fst with min to reach sufficiency thresholds for allele sampling
	for (min_thresh in c(70,95)){	

		min_samp<-read.csv(paste0("min_for_",min_thresh,n_drop_file,".csv"));	min_samp[min_samp=="Inf"]<-NA
		reg_pval<-matrix(nrow=5, ncol=4); reg_r2<-matrix(nrow=5, ncol=4)

		for (ss in 1:4){  #this is a loop over min, max, mean and sd
		  for (at in c(2,4:6)){   #this is the loop over each allele type
			reg_pval[at-1,ss]<-((summary(lm(min_samp[,at] ~ sp_summ_fst[,ss]))$coefficients[,4]))[2]
			reg_r2[at-1,ss]<-summary(lm(min_samp[,at] ~ sp_summ_fst[,ss]))$r.squared
			if (reg_pval[at-1,ss]<0.05) plot(min_samp[,at], sp_summ_fst[,ss],main=paste(n_drop_file,"thresh",min_thresh,"allele type", at),
				ylab=colnames(sp_summ_fst)[ss])
		  }
		}
		colnames(reg_pval)<-colnames(sp_summ_fst); colnames(reg_r2)<-colnames(sp_summ_fst)
		reg_pval<-reg_pval[-2,]; reg_2<-reg_r2[-2,]
		rownames(reg_pval)<-paste(c("global", "common", "low freq", "rare"),paste(n_drop_file,"thresh",min_thresh))
		results_min<-rbind(results_min,reg_pval)
	}
	
	##regressions of pairwise FST for collector gap calculation
	collector_gap_res<-read.csv(file=paste0("collector_gaps",n_drop_file,".csv"))
	collector_gap_res<-collector_gap_res[,-1]
	reg_pval<-matrix(nrow=2, ncol=4); reg_r2<-matrix(nrow=2, ncol=4)
		
	for (ss in 1:4){ 
		reg_pval[1,ss]<-((summary(lm(collector_gap_res[,6] ~ sp_summ_fst[,ss]))$coefficients[,4]))[2]
		reg_r2[1,ss]<-summary(lm(collector_gap_res[,6] ~ sp_summ_fst[,ss]))$r.squared
		reg_pval[2,ss]<-((summary(lm(collector_gap_res[,8] ~ sp_summ_fst[,ss]))$coefficients[,4]))[2]
		reg_r2[2,ss]<-summary(lm(collector_gap_res[,8] ~ sp_summ_fst[,ss]))$r.squared
	}
	colnames(reg_pval)<-colnames(sp_summ_fst); colnames(reg_r2)<-colnames(sp_summ_fst)
	rownames(reg_pval)<-paste(c("times_reduce","percent_increas"),paste(n_drop_file,"thresh",min_thresh))
	results_gap<-rbind(results_gap,reg_pval)
	
}
#print raw p value results and p value correction

write.csv(rbind(round(results_capture,3),
	round(matrix(p.adjust(results_capture,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_capture)),colnames(results_capture)),ncol=ncol(results_capture)),3)),
	file="p_values/pvals_reg_fst_ex_situ_current.csv")
write.csv(rbind(round(results_min,3),
	round(matrix(p.adjust(results_min,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_min)),colnames(results_min)),ncol=ncol(results_min)),3)),
	file="p_values/pvals_reg_fst_min_to_reach_thresh.csv")
write.csv(rbind(round(results_gap,3),
	round(matrix(p.adjust(results_gap,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_gap)),colnames(results_gap)),ncol=ncol(results_gap)),3)),
	file="p_values/pvals_reg_fst_collector_gap.csv")



		
######################################
#	STATS: Regressions on properties of the Allele freq histograms	
#	Here we test the influence of properties of the allele frequency spectrum (basically proportion of rare alleles,
#	as you'd see looking at a frequency histogram) on the genetic capture etc.
#	proportion of alleles below a given frequency (0.1, 0.05, 0.01, 0.005) will be our predictor variable
#	this is a quantitative, continuous variable so we use linear regression 
#	As above, we have three sets of tests
#	CURRENT EX SITU VS IN SITU: Current genetic diversity conserved ex situ in collections today
#		This is 32 regressions (2 files to analyze (Full/ Red), 4 allele categories, and 4 summaries of allele freqs)
#  MINIMUM TO CATCH: Minimum sampling needed to achieve a sufficiency threshold (e.g. 70%, 95%)
#		This is 64 regressions (2 files to analyze (Full/ Red), 4 allele categories, 2 min thresholds and 4 summaries of allele freqs)
#  COLLECTOR GAP: The collector gap in terms of reduction in collection size and improved collection capture 
#		This is 16 regressions (2 files to analyze (Full/ Red), 2 types of collector gap, and 4 summaries of allele freqs)
#
#	NOTE: The variables for results matrices (results_gap, results_capture, and results_min) are the same as those from above code
#	It was easy to jut copy and re-use them so just be careful- these are not totally unique to this section
###############################################################################

#####################################

#This will go through each species and make allele frequency histograms as well as report proportion of alleles in bins
#below 0.05, 0.01, and 0.005 frequency, and save this table as .csv 
setwd(prefix_wd)

#This needed simply to know which populations are "wild"
region_makeup_list<-list(list(1:2,3:4),list(1:2,3),list(1,2),list(1:2,3:4),list(1,2),list(1:2,3),list(1:2,3,4,5:9),list(1:2,3,4,5:10),list(1:2,3,4,5:8),list(1,2),list(1:2,3))
#Place to store results- columns are types of allele, rows are species, elements are the proportion of alleles that fall below a freq threshold
allele_freq_counts<-matrix(nrow=length(folder_names),ncol=4); rownames(allele_freq_counts)<-folder_names
colnames(allele_freq_counts)<-c("f_0.005","f_0.01","f_0.05","f_0.1")
pdf("allele_freq_hist.pdf",height=10,width=14)

for (sp in 1:length(folder_names)){
	this_species<-folder_names[sp]
	Spp_tot_genind<-read.genepop(paste(prefix_wd,this_species,"/",substr(this_species,1,2),"_total.gen",sep=""),ncode=3)

	wild_p<-unlist(region_makeup_list[[sp]]);		n_ind_W<-sum(table(Spp_tot_genind@pop)[wild_p])
	Spp_tot_genpop<-genind2genpop(Spp_tot_genind)
	#Get the frequencies, and get the total number of alleles- this will allow calculating how many alleles are below a frequency
	wild_all_freqs<-colSums(Spp_tot_genpop[wild_p]@tab,na.rm=T)/(n_ind_W*2)
	total_num_alleles<-length(Spp_tot_genpop[wild_p]@tab[1,])
	
	#plots of allele frequency histograms
	par(fig = c(0,1,0,1))
	hist(sort(wild_all_freqs),breaks=seq(0,1,by=.01),xlim=c(0,.15),main=paste0("allele freq histogram for ",species_names[sp]), xlab="allele frequency",ylab="number alleles in that frequency"); abline(v=c(0.01,0.05),col="red")
	par(fig = c(0.4,1, 0.4, 1), new = T) 
	hist(sort(wild_all_freqs),breaks=seq(0,1,by=.02),main="",ylab="",xlab="")
	
	#calculating how many alleles are below a frequency
	allele_freq_counts[sp,1:4]<-c((sum(wild_all_freqs<0.005)/total_num_alleles),(sum(wild_all_freqs<0.01)/total_num_alleles),
		(sum(wild_all_freqs<0.05)/total_num_alleles), (sum(wild_all_freqs<0.1)/total_num_alleles))
}
dev.off()
write.csv(allele_freq_counts,"allele_freq_counts_by_spp.csv")


 par(mfrow=c(3,3))
results_min<-matrix(nrow=1,ncol=4)	#all have four columns because that is the number of summary stats of the allele freqs
results_capture<-matrix(nrow=1,ncol=4)
results_gap<-matrix(nrow=1,ncol=4)

for (n_to_drop in c(0,2)){
	if (n_to_drop==0) n_drop_file<-"_dr_0" 
	if (n_to_drop==2) n_drop_file<-""
	
	#regressions of prop alleles below thresholds with current ex situ genetic capture 
	gen_cap<-read.csv(paste0("ex_vs_in_situ",n_drop_file,".csv"))[,-1];	gen_cap[gen_cap=="Inf"]<-NA
	reg_pval<-matrix(nrow=5, ncol=4); reg_r2<-matrix(nrow=5, ncol=4)

	for (ss in 1:4){  #this is a loop over min, max, mean and sd
	  for (at in c(2,4:6)){   #this is the loop over each allele type
		reg_pval[at-1,ss]<-((summary(lm(gen_cap[,at] ~ allele_freq_counts[,ss]))$coefficients[,4]))[2]
		reg_r2[at-1,ss]<-summary(lm(gen_cap[,at] ~ allele_freq_counts[,ss]))$r.squared
	  }
	}
	colnames(reg_pval)<-colnames(allele_freq_counts); colnames(reg_r2)<-colnames(allele_freq_counts)
	reg_pval<-reg_pval[-2,]; reg_2<-reg_r2[-2,]
	rownames(reg_pval)<-paste(c("global", "common", "low freq", "rare"),n_drop_file)
	results_capture<-rbind(results_capture,reg_pval)

	##regression of prop alleles below thresholds with min to reach sufficiency thresholds for allele sampling
	for (min_thresh in c(70,95)){	

		min_samp<-read.csv(paste0("min_for_",min_thresh,n_drop_file,".csv"));	min_samp[min_samp=="Inf"]<-NA
		reg_pval<-matrix(nrow=5, ncol=4); reg_r2<-matrix(nrow=5, ncol=4)

		for (ss in 1:4){  #this is a loop over min, max, mean and sd
		  for (at in c(2,4:6)){   #this is the loop over each allele type
			reg_pval[at-1,ss]<-((summary(lm(min_samp[,at] ~ allele_freq_counts[,ss]))$coefficients[,4]))[2]
			reg_r2[at-1,ss]<-summary(lm(min_samp[,at] ~ allele_freq_counts[,ss]))$r.squared
			if (reg_pval[at-1,ss]<0.05) plot(min_samp[,at], allele_freq_counts[,ss],main=paste(n_drop_file,"thresh",min_thresh,"allele type", at),
				ylab=colnames(allele_freq_counts)[ss])
		  }
		}
		colnames(reg_pval)<-colnames(allele_freq_counts); colnames(reg_r2)<-colnames(allele_freq_counts)
		reg_pval<-reg_pval[-2,]; reg_2<-reg_r2[-2,]
		rownames(reg_pval)<-paste(c("global", "common", "low freq", "rare"),paste(n_drop_file,"thresh",min_thresh))
		results_min<-rbind(results_min,reg_pval)
	}
	
	##regressions of prop alleles below thresholds for collector gap calculation
	collector_gap_res<-read.csv(file=paste0("collector_gaps",n_drop_file,".csv"))
	collector_gap_res<-collector_gap_res[,-1]
	reg_pval<-matrix(nrow=2, ncol=4); reg_r2<-matrix(nrow=2, ncol=4)
		
	for (ss in 1:4){ 
		reg_pval[1,ss]<-((summary(lm(collector_gap_res[,6] ~ allele_freq_counts[,ss]))$coefficients[,4]))[2]
		reg_r2[1,ss]<-summary(lm(collector_gap_res[,6] ~ allele_freq_counts[,ss]))$r.squared
		reg_pval[2,ss]<-((summary(lm(collector_gap_res[,8] ~ allele_freq_counts[,ss]))$coefficients[,4]))[2]
		reg_r2[2,ss]<-summary(lm(collector_gap_res[,8] ~ allele_freq_counts[,ss]))$r.squared
	}
	colnames(reg_pval)<-colnames(allele_freq_counts); colnames(reg_r2)<-colnames(allele_freq_counts)
	rownames(reg_pval)<-paste(c("times_reduce","percent_increas"),paste(n_drop_file,"thresh",min_thresh))
	results_gap<-rbind(results_gap,reg_pval)
	
}

write.csv(rbind(round(results_capture,3),
	round(matrix(p.adjust(results_capture,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_capture)),colnames(results_capture)),ncol=ncol(results_capture)),3)),
	file="p_values/pvals_reg_af_ex_situ_current.csv")
write.csv(rbind(round(results_min,3),
	round(matrix(p.adjust(results_min,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_min)),colnames(results_min)),ncol=ncol(results_min)),3)),
	file="p_values/pvals_reg_af_min_to_reach_thresh.csv")
write.csv(rbind(round(results_gap,3),
	round(matrix(p.adjust(results_gap,method="BH"),
	dimnames=list(paste0("adj.",rownames(results_gap)),colnames(results_gap)),ncol=ncol(results_gap)),3)),
	file="p_values/pvals_reg_af_collector_gap.csv")
 round(matrix(p.adjust(results_capture,method="BH"),ncol=4,dimnames=list(rownames(results_capture),colnames(results_capture))),2)
 round(matrix(p.adjust(results_min,method="BH"),ncol=4,dimnames=list(rownames(results_min),colnames(results_min))),2)
 round(matrix(p.adjust(results_gap,method="BH"),ncol=4,dimnames=list(rownames(results_gap),colnames(results_gap))),2)
  

