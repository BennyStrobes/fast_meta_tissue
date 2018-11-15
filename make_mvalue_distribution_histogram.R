args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(grid)
library(gplots)
library(RColorBrewer)
library(MASS)
require(reshape2)
require(cowplot)


get_positional_counts <- function(num_m) {
	y_value <- rep(0,49)
	for (row_num in 1:nrow(num_m)) {
		y_value[num_m$num_m[row_num]] <- y_value[num_m$num_m[row_num]] +1 
	}
	norm_y_value = y_value/sum(y_value)
	new = c(norm_y_value,c(norm_y_value[49]))
	return(new)
}


get_positional_counts_zeros <- function(num_m) {
	y_value <- rep(0,50)
	for (row_num in 0:nrow(num_m)) {
		y_value[num_m$num_m[row_num]] <- y_value[num_m$num_m[row_num]] +1 
	}
	norm_y_value = y_value/sum(y_value)
	new = c(norm_y_value,c(norm_y_value[50]))
	return(new)
}

get_num_m_gt_thresh <- function(df,thresh) {
	counts <- rep(0,nrow(df))
	for (row_num in 1:nrow(df)) {
		ms <- df[row_num,]
		count <- 0
		for (col_num in 1:length(ms)) {
			if(!is.na(ms[col_num])) {
				if(ms[col_num] > thresh) {
					count <- count + 1
				}
			}
		}
		counts[row_num] <- count
	}
	num_m <- counts[ which(!counts == 0)]
	return(as.data.frame(num_m))
}


get_num_m_gt_thresh_include_zeros <- function(df,thresh) {
	counts <- rep(0,nrow(df))
	for (row_num in 1:nrow(df)) {
		ms <- df[row_num,]
		count <- 0
		for (col_num in 1:length(ms)) {
			if(!is.na(ms[col_num])) {
				if(ms[col_num] > thresh) {
					count <- count + 1
				}
			}
		}
		counts[row_num] <- count
	}
	num_m <- counts[ which(counts >= 0)]
	return(as.data.frame(num_m))
}

generate_cis_top_cis_trans_mvalue_histogram <-function(trans, top_cis,cis,threshold) {
	###Create a vector the length of number of trans eqtls. Where each element is number of tissues with m-value > (threshold)
	trans_num_m <- get_num_m_gt_thresh(trans,threshold)
	###Create a vector the length of number of cis eqtls. Where each element is number of tissues with m-value > (threshold)
	cis_num_m <- get_num_m_gt_thresh(cis,threshold)

	top_cis_num_m <- get_num_m_gt_thresh(top_cis, threshold)
	###Convert cis_num_m vector to a vector of length 49 (1:49). Where each element i is the number of eqtls that have i tissues with m-value > .5. (this is used to create an "open histogram")

	cis_poitional_counts <- get_positional_counts(cis_num_m)
	top_cis_positional_counts <- get_positional_counts(top_cis_num_m)
	###Plot away!
	g_plot = ggplot(trans_num_m,aes(x=num_m)) + 
	geom_histogram(aes(y=..count../sum(..count..),col="black",fill="Trans"),size=.1,binwidth=1) +     ##Plot trans histogram
	geom_step(data = as.data.frame(cis_poitional_counts),aes(x=seq_along(cis_poitional_counts)-.5,y=cis_poitional_counts,colour='red',fill="Cis"),size=.6) +  ##Plot cis open Histogram
	geom_step(data = as.data.frame(top_cis_positional_counts),aes(x=seq_along(top_cis_positional_counts)-.5,y=top_cis_positional_counts,colour='pink',fill="Cis_top"),size=.6) +  ##Plot cis open Histogram

	labs(x=paste0("Number of tissues with m-values > ", threshold),y="Fraction of total genes",fill="",colour="") +        ##Add x-axis and y-axis labels
	#labs(x = "Sample size of discovery tissue", y = "# of shared cis eQTL tissues") + 
	scale_fill_manual(values=alpha(c("Cis"="red","Cis_top"="pink","Trans" = "darkturquoise"),c(.6,1,1))) +      ##Fill in with color!
	scale_colour_manual(values=alpha(c("black" ="black","pink"="pink","red" = "red"),c(1,1,.6)))+ guides(colour=FALSE) +
	theme(panel.background=element_rect(fill = "white",colour=alpha("black",.5)),legend.text=element_text(size=8),legend.key.size=unit(.8,"cm"),legend.key.width=unit(.2,"cm"),legend.key.height=unit(.2,"cm")) +     ##Size formatting...
	coord_cartesian(xlim = c(1.35, 47.42),ylim = c(.03,.53))    #size formatting
	
	fsize=6.7;
   g_plot =  g_plot + 
                theme(axis.text.x=element_text( hjust=0.5, vjust=0.5, size=fsize),      ## text on the x ticks
                axis.text.y=element_text(size=fsize),                                          ## text on the y ticks         
                axis.title=element_text(size = 8),                                              ## x label
                axis.ticks = element_line(size=0.1),  
                strip.text = element_text(size=fsize),                                              ## margin of the plot
                panel.grid.minor = element_blank(),                                        ## remove minor grid
                panel.grid.major = element_blank(),											##Remove major grid
            	legend.justification = c(1, 1),legend.position=c(1.02,1.03))               ## move the legend inside the plot


	return(g_plot)
	#ggplot_function(g_plot,output_directory,output_file,74,53,0)
}



generate_cis_trans_mvalue_histogram <-function(trans,cis,threshold) {
	###Create a vector the length of number of trans eqtls. Where each element is number of tissues with m-value > (threshold)
	trans_num_m <- get_num_m_gt_thresh(trans,threshold)
	###Create a vector the length of number of cis eqtls. Where each element is number of tissues with m-value > (threshold)
	cis_num_m <- get_num_m_gt_thresh(cis,threshold)
	###Convert cis_num_m vector to a vector of length 49 (1:49). Where each element i is the number of eqtls that have i tissues with m-value > .5. (this is used to create an "open histogram")
	cis_poitional_counts <- get_positional_counts(cis_num_m)
	print(cis_poitional_counts)
	###Plot away!
	g_plot = ggplot(trans_num_m,aes(x=num_m)) + 
	geom_histogram(aes(y=..count../sum(..count..),col="black",fill="Trans"),size=.1,binwidth=1) +     ##Plot trans histogram
	geom_step(data = as.data.frame(cis_poitional_counts),aes(x=seq_along(cis_poitional_counts)-.5,y=cis_poitional_counts,colour='red',fill="Cis"),size=.6) +  ##Plot cis open Histogram
	labs(x=paste0("Number of tissues with m-values > ", threshold),y="Fraction of total genes",fill="",colour="") +        ##Add x-axis and y-axis labels
	#labs(x = "Sample size of discovery tissue", y = "# of shared cis eQTL tissues") + 
	scale_fill_manual(values=alpha(c("Cis"="red","Trans" = "darkturquoise"),c(.6,1))) +      ##Fill in with color!
	scale_colour_manual(values=alpha(c("black" ="black","red" = "red"),c(1,.6)))+ guides(colour=FALSE) +
	theme(panel.background=element_rect(fill = "white",colour=alpha("black",.5)),legend.text=element_text(size=8),legend.key.size=unit(.8,"cm"),legend.key.width=unit(.2,"cm"),legend.key.height=unit(.2,"cm")) +     ##Size formatting...
	coord_cartesian(xlim = c(1.35, 47.42),ylim = c(.03,.53))    #size formatting
	
	fsize=6.7;
   g_plot =  g_plot + 
                theme(axis.text.x=element_text( hjust=0.5, vjust=0.5, size=fsize),      ## text on the x ticks
                axis.text.y=element_text(size=fsize),                                          ## text on the y ticks         
                axis.title=element_text(size = 8),                                              ## x label
                axis.ticks = element_line(size=0.1),  
                strip.text = element_text(size=fsize),                                              ## margin of the plot
                panel.grid.minor = element_blank(),                                        ## remove minor grid
                panel.grid.major = element_blank(),											##Remove major grid
            	legend.justification = c(1, 1),legend.position=c(1.02,1.03))               ## move the legend inside the plot


	return(g_plot)
	#ggplot_function(g_plot,output_directory,output_file,74,53,0)
}




trans_egene_file <- args[1]
cis_egene_file <- args[2]
top_cis_egene_file <- args[3]
output_directory <- args[4]



trans_mvalues <- (read.table(trans_egene_file, header = TRUE))[,56:104]
cis_mvalues <- (read.table(cis_egene_file,header = TRUE))[,56:104]
top_cis_mvalues <- (read.table(top_cis_egene_file, header = TRUE))[,56:104]


threshold <- .5
histo <- generate_cis_top_cis_trans_mvalue_histogram(trans_mvalues, top_cis_mvalues, cis_mvalues, threshold)
ggsave(paste0(output_directory, "cis_top_cis_trans_mvalue_histogram_m_thresh_", threshold, ".png"),histo, width=10, height=5, units="cm")



threshold <- .5
# histo <- generate_cis_trans_mvalue_histogram(trans_mvalues, cis_mvalues, threshold)
# ggsave(paste0(output_directory, "cis_trans_mvalue_histogram_m_thresh_", threshold, ".png"),histo, width=10, height=5, units="cm")

threshold <- .7
# histo <- generate_cis_trans_mvalue_histogram(trans_mvalues, cis_mvalues, threshold)
# ggsave(paste0(output_directory, "cis_trans_mvalue_histogram_m_thresh_", threshold, ".png"),histo, width=10, height=5, units="cm")

threshold <- .9
# histo <- generate_cis_trans_mvalue_histogram(trans_mvalues, cis_mvalues, threshold)
# ggsave(paste0(output_directory, "cis_trans_mvalue_histogram_m_thresh_", threshold, ".png"),histo, width=10, height=5, units="cm")

