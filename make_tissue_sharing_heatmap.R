args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(grid)
library(gplots)
library(RColorBrewer)
library(MASS)
require(reshape2)
require(cowplot)

get_heatmap <- function(df,dendy,tissues) {
	cc = matrix(NA, 49, 49);
	for (i in 1:49) {
		for (j in i:49) {
    		cc[i,j] = cor(df[,i], df[,j], use='pairwise.complete.obs', method='spearman');
    		cc[j,i] =cc[i,j]
  		}
   		cc[i,i] = NA;

	}


	cc = abs(cc)
	row.names(cc) <-tissues;
	colnames(cc) <- tissues;
	my_palette <- colorRampPalette(c("#FC8D59", "#FFFFBF","#c8e7af", "#467323"))(n = 256);

	heatmap_df <- heatmap.2(t(cc), Colv=dendy,keysize=1,col=my_palette,Rowv=rev(dendy),dendrogram='none',trace="none",density.info="none",labCol=NA,key.title=NA,key.xlab=" Spearman correlation",lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(1.5,4,1), lwid=c(1.5,4));
	return(heatmap_df)
}

hierarchial_clustering <- function(df_beta_cis,tissue_file) {
	##cc is tissue by tissue correlation matrix
	cc = matrix(NA, 49, 49);
	##We fill in elements of cc using spearman correlation on only observed elements
	for (i in 1:49) {
		for (j in i:49) {
    		cc[i,j] = cor(df_beta_cis[,i], df_beta_cis[,j], use='pairwise.complete.obs', method='spearman');
    		cc[j,i] =cc[i,j]
  		}
   		cc[i,i] = 1.0;
	}
	tissues = as.matrix(read.table(tissue_file));
	row.names(cc) <- tissues;
	colnames(cc) <- tissues;
	cc = abs(cc)
	##Use 1 - correlation as distance metric
	disty = as.dist(t(1-cc))
	##Perform heirarchial aggolomerative clustering using 'average' linkage based on cis eqtls
	clust = hclust(disty,method='average')
	dendy = as.dendrogram(clust)	
	cis_heatmap <- get_heatmap(df_beta_cis,dendy,tissues)
	ordered_tissues <- rev(tissues[cis_heatmap$rowInd])
	return(rev(cis_heatmap$rowInd))
}

clustering <- function(trans,cis,out_dir,tissue_file) {
	indices <- hierarchial_clustering(cis, tissue_file)	
	tissues = (as.vector(as.matrix(read.table(tissue_file))))[indices]
	clustered_cis_file <- paste0(out_dir,"cis_clustered.txt")
	median_cis <- print_clustered_matrix_cis(indices,cis,tissue_file,clustered_cis_file)
	clustered_trans_file <- paste0(out_dir,"trans_clustered.txt")
	median_trans <- print_clustered_matrix_trans(indices,trans,tissue_file,clustered_trans_file)
	return(c(median_cis, median_trans))

}

print_clustered_matrix_cis <- function(indices,df_beta_cis,tissue_file,output_file) {
	cc = matrix(NA, 49, 49);
	corrz <- c()
	##We fill in elements of cc using spearman correlation on only observed elements
	for (i in 1:49) {
		for (j in i:49) {
    		cc[i,j] = cor(df_beta_cis[,i], df_beta_cis[,j], use='pairwise.complete.obs', method='spearman');
    		cc[j,i] =cc[i,j]
    		if (i != j) {
    			corrz <- c(corrz, cc[i,j])
    		}
  		}
   		cc[i,i] = 0.0;
	}

	max_value = max(cc)
	for (i in 1:49) {
		cc[i,i] = max_value;
	}
	tissues = as.matrix(read.table(tissue_file));
	ordered_tissues <- tissues[indices]
	ordered_tissues <- sub("-", "", ordered_tissues, fixed = TRUE)
	clustered_mat = matrix(NA,50,49);
	for (i in 1:49) {
		for (j in 1:49) {
			new_i = which(indices == i) +1;
			new_j = which(indices == j);
			if (new_j >= new_i-1) {
				clustered_mat[new_i,new_j] = cc[i,j];
			} else {
				clustered_mat[new_i,new_j] = NA;
			}
		}
		clustered_mat[1,i] = ordered_tissues[i];

	}

	row.names(clustered_mat) <- c("Row",ordered_tissues);
	write.table(clustered_mat,file=output_file,sep="\t",quote=FALSE,col.names=FALSE)
	return(median(corrz))
}

print_clustered_matrix_trans <- function(indices,df_beta_cis,tissue_file,output_file) {
	cc = matrix(NA, 49, 49);
	corrz <- c()
	##We fill in elements of cc using spearman correlation on only observed elements
	for (i in 1:49) {
		for (j in i:49) {
    		cc[i,j] = cor(df_beta_cis[,i], df_beta_cis[,j], use='pairwise.complete.obs', method='spearman');
    		cc[j,i] =cc[i,j]
    		if (i != j) {
    			corrz <- c(corrz, cc[i,j])
    		}
  		}
   		cc[i,i] = 0.0;
	}

	max_value = max(cc)
	for (i in 1:49) {
		cc[i,i] = max_value;
	}
	tissues = as.matrix(read.table(tissue_file));
	ordered_tissues <- tissues[indices]
	ordered_tissues <- sub("-", "", ordered_tissues, fixed = TRUE)
	clustered_mat = matrix(NA,50,49);
	for (i in 1:49) {
		for (j in 1:49) {
			new_i = which(indices == i) +1;
			new_j = which(indices == j);
			if (new_i -1 >= new_j) {
				clustered_mat[new_i,new_j] = cc[i,j];
			} else {
				clustered_mat[new_i,new_j] = NA;
			}
		}
		clustered_mat[1,i] = ordered_tissues[i];

	}

	row.names(clustered_mat) <- c("Row",ordered_tissues);
	write.table(clustered_mat,file=output_file,sep="\t",quote=FALSE,col.names=FALSE)
	return(median(corrz))
}

make.sharing.plot.cis = function(df, middler) {
    p = ggplot(data = df, aes(x = variable, y = Row)) +
        geom_tile(aes(fill = value), colour = 'black') + theme_bw() +
	scale_fill_gradient2(low = 'dodgerblue4',mid="white", high = 'orangered4',midpoint=middler, name = 'cis spearman\ncorrelation', na.value = "transparent") + #Midpoint taken to be median
	theme(axis.ticks = element_blank(),
			  axis.text = element_blank(),
              axis.title = element_blank(),
              #plot.margin = unit(c(-1,-1,-1,-1), "mm"),
              panel.background = element_rect(fill="transparent",colour="transparent"),
              panel.border = element_blank(),
              legend.position = "top",
              legend.title=element_text(size=8), 
              legend.text=element_text(size=8),
              legend.key.size=unit(5,"mm"),
              plot.background = element_rect(fill="transparent",colour="transparent"),
              panel.grid=element_blank())
    p$layout$clip[p$layout$name=="panel"] <- "off"
    return(p)
}

make.sharing.plot.trans = function(df, middler) {
    p = ggplot(data = df, aes(x = variable, y = Row)) +
        geom_tile(aes(fill = value), colour = 'black') + theme_bw() +
	scale_fill_gradient2(low = 'dodgerblue4',mid="white", high = 'orangered4',midpoint=middler, name = 'trans spearman\ncorrelation', na.value = "transparent") + #Midpoint taken to be median
	theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              #plot.margin = unit(c(-3,-3,-3,-3), "mm"),
              panel.background = element_rect(fill="transparent",colour="transparent"),
              panel.border = element_blank(),
              legend.position = "bottom",
              legend.title=element_text(size=8), 
              legend.text=element_text(size=8),
              legend.key.size=unit(5,"mm"),
              plot.background = element_rect(fill="transparent",colour="transparent"),
              panel.grid=element_blank())
    p$layout$clip[p$layout$name=="panel"] <- "off"
    return(p)
}

make.combined.plot = function(cis.sharing, trans.sharing,colors.vertical, colors.horizontal) {
    combined = ggdraw() +
        draw_plot(cis.sharing, 0.052,0.1860,.7,0.807) +
        draw_plot(trans.sharing, 0.052,0.02,0.7,0.807) + 
        draw_plot(colors.vertical, -.282, .155, .7, .7) + 
        draw_plot(colors.horizontal,0.02,.429, .76,.76)
        #draw_plot(legend, 0.65,0,0.35,1)
        #draw_plot(vert, .005,0.068,0.05,0.608) +
        #draw_plot(horiz, 0.001,0.065,0.602,0.05) +
        #draw_plot(cis.sharing, 0.052,0.114,0.5,0.607) +
        #draw_plot(trans.sharing, 0.052,0.02,0.5,0.607) 
    return(combined)
}

make.point.plot = function(tissuesdf, colors, vertical = TRUE){
    if (vertical) {
        p = ggplot(tissuesdf, aes(x = 1, y = -order, label = tissue_id))
    } else {
        p = ggplot(tissuesdf, aes(x = order, y = 1))
    }
    p = p + geom_point(aes(colour = tissue_id), size = .9) +
        scale_colour_manual(values = colors) + guides(colour = FALSE) + 
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())
    return(p)
}

cis_05_file_name <- args[1]
trans_file_name <- args[2]
tissue_name_input_file <- args[3]
gtex_tissue_colors_file <- args[4]
output_directory <- args[5]



# Load in data
trans <- read.table(trans_file_name, header = TRUE)
cis_05 <- read.table(cis_05_file_name,header=TRUE)
# Read in tissue colors and names
tissues = read.table(gtex_tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")


# Extract metatissue betas from matrices
beta_indices <- (105:(dim(trans)[2]))[c(TRUE,FALSE)]
trans_beta <- trans[,beta_indices]
cis_beta <- cis_05[,beta_indices]

medians <- clustering(trans_beta, cis_beta, output_directory, tissue_name_input_file)


cis.sharing.raw = read.table(paste0(output_directory,"cis_clustered.txt"), header = T, stringsAsFactors = F, sep = '\t')
trans.sharing.raw = read.table(paste0(output_directory,"trans_clustered.txt"), header = T, stringsAsFactors = F, sep = '\t')


# Melt the dataframe
cis.sharing = melt(cis.sharing.raw)
cis.sharing$variable = factor(as.character(cis.sharing$variable), levels = colnames(cis.sharing.raw)[-1])
cis.sharing$Row = factor(cis.sharing$Row, levels = rev(levels(cis.sharing$variable)))

cis.sharing$value[is.na(cis.sharing$value)] = NA; # turn NaN into NA

trans.sharing = melt(trans.sharing.raw)
trans.sharing$variable = factor(as.character(trans.sharing$variable), levels = colnames(trans.sharing.raw)[-1])
trans.sharing$Row = factor(trans.sharing$Row, levels = rev(levels(trans.sharing$variable)))

trans.sharing$value[is.na(trans.sharing$value)] = NA 


tissues$tissue_id = factor(tissues$tissue_id, levels = as.character(levels(cis.sharing$variable)))
tissues$order = as.numeric(tissues$tissue_id)
tissues = tissues[!is.na(tissues$order),]
colors = paste0("#",tissues$tissue_color_hex)
names(colors) = tissues$tissue_id


colors.vertical = make.point.plot(tissues, colors)
colors.horizontal = make.point.plot(tissues, colors, vertical = FALSE)


cis.sharing.plot = make.sharing.plot.cis(cis.sharing, medians[1])

trans.sharing.plot = make.sharing.plot.trans(trans.sharing, medians[2])

combined_plot <- make.combined.plot(cis.sharing.plot, trans.sharing.plot, colors.vertical, colors.horizontal)

ggsave(paste0(output_directory, "cis_trans_colors.png"),combined_plot, width=10, height=10, units="cm")

