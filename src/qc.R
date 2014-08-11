#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(optparse))
suppressMessages(library(compiler))
suppressMessages(library(stats))
suppressMessages(library(grDevices))
suppressMessages(library(graphics))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(gridExtra))

options(warn=-1)

qc.dataMatrix = function(input_matrix){
  
  ## TO DO: CHECK IF MATRIX IS WELL FORMATTED 
  bad_cols = 0
  data_col_idx = 5
  data_row_idx = 3
  prey_col_idx = 1
  bait_names = input_matrix[1,data_col_idx:ncol(input_matrix)]
  ip_names = colnames(input_matrix)[data_col_idx:ncol(input_matrix)]
  #prey_names = rownames(input_matrix)[data_row_idx:nrow(input_matrix)]
  data_matrix = input_matrix[data_row_idx:nrow(input_matrix),data_col_idx:ncol(input_matrix)]
  rownames(data_matrix) = input_matrix[data_row_idx:nrow(input_matrix),prey_col_idx]
  if(any(data_matrix==""))
    stop("BLANK ENTRIES FOUND IN THE DATA MATRIX.")
  data.matrix(data_matrix)  
}

qc.getIpToBaits = function(input_matrix){
  bait_row_idx = 1
  data_col_idx = 5
  df = data.frame(colnames(input_matrix)[data_col_idx:ncol(input_matrix)],t(input_matrix[bait_row_idx,data_col_idx:ncol(input_matrix)]), 0)
  colnames(df) = c("ip","bait","replicate")
  bait_reps = aggregate( ip ~ bait, FUN=length, data=df)
  
  for(b in 1:nrow(bait_reps)){
    bait_cur = bait_reps[b,]
    df[df$bait==bait_cur$bait,]$replicate = 1:bait_cur$ip 
  }
  df
}

qc.clusterHeatmap = function(data_matrix, output_file, ip_baits, font_scale){
  data_matrix_w_names = data_matrix
  colnames(data_matrix_w_names) = do.call(paste, c(ip_baits[,c('bait','ip')], sep = " "))
  cor_matrix = cor(data_matrix_w_names, use="pairwise.complete.obs", method="pearson")
  color_scale = colorRampPalette(brewer.pal(7,"Blues"))(9)
  pheatmap(cor_matrix, cluster_rows=T, cluster_cols=T, scale="none",fontsize_row=font_scale,fontsize_col=font_scale, cellwidth=font_scale, cellheight=font_scale, border_color=NA, filename=output_file, color=color_scale, breaks=1:10/10, clustering_distance_rows="correlation", clustering_distance_cols="correlation", treeheight_row=0, treeheight_col=0) 
}

qc.ipDists = function(data_matrix, ip_baits, baseName){
  tmp = melt(data_matrix, varnames=c("prey","ip"), value.name="count")
  data_long = merge(tmp, ip_baits, by="ip")
  data_long = data_long[data_long$count > 0,]
  theme_set(theme_bw(base_size = 12,base_family='Helvetica')) 
  #pdf(gsub('.txt','_proteincounts.pdf',baseName), width=10, height=(ceiling(length(unique(data_long$bait))/20)*8) )
  #pdf(gsub('.txt','_proteincounts.pdf',baseName), width=10, height=(ceiling(length(unique(data_long$bait))/20)*8))   
  p = ggplot(data_long, aes(x=ip)) + geom_bar(stat='bin') + facet_wrap(facets= ~bait, scales='free_x', ncol=10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #dev.off()
  ggsave(filename=gsub('.txt','_proteincounts.pdf',baseName), plot=p, scale=length(unique(data_long$bait))/10, limitsize=FALSE)
  
  ## solution with ggplot
  #ggplot(data_long, aes(x=log2(count), colour=factor(replicate))) + geom_density(alpha=1, size=1) + facet_wrap(facets= ~bait, scales='fixed') 
  
  ## solution to give every plot its separate legend with gridExtra library 
#   pdf(gsub('.txt','_peptidecount_dists.pdf',baseName), width=30, height=30)
#   out = by(data = data_long, INDICES = data_long$bait, FUN = function(m) {
#     m = droplevels(m)
#     #print(m)
#     bait_title = unique(m$bait)
#     m = ggplot(m, aes(x=count, colour=ip)) + geom_density(alpha=1, size=1) + ggtitle(bait_title)
#     ## + theme(legend.title=element_text(unique(m$bait)))
#   })
#   do.call(grid.arrange, out)
#   dev.off()
}

qc.NumUniquePlot <- function(ip_matrix, matrix_file){
  x <- ip_matrix[,-c(2:4)]
  baits = data.frame(id=names(x)[-1], bait=unlist(x[1,-1]), stringsAsFactors=F )
  x <- x[-c(1,2),]
  x <- melt(x, id=c('a'))
  x = merge(x, baits, by.x=c('variable'), by.y=c('id'))
  names(x) = c('ID','Prey', 'value','Bait')
  x$value = as.integer(x$value)
  x = x[!x$value==0,]
  
  bait_num = length(unique(x$Bait))
  plots_per_col = 5
  plot_width = 15
  plot_height = ((plot_width/plots_per_col) * ceiling(bait_num/plots_per_col))
  
  outfile = gsub(".txt", "_PepCount_Distributions.pdf", matrix_file)
  #plot
  
  pdf(file=outfile, width=plot_width, height=plot_height)
    p = ggplot(x, aes(x=factor(ID), y=value))
    print( p + geom_boxplot() + facet_wrap(~ Bait, scales="free", drop=T, ncol=plots_per_col) + theme(axis.text=element_text(size=9),  axis.text.x=element_text(angle=45, hjust=1))  )
  dev.off()	
}


qc.main = function(matrix_file, font_scale, cluster=T, ip_dists=T){
  ip_matrix = read.delim(matrix_file, stringsAsFactors=F)
  data_matrix = qc.dataMatrix(ip_matrix)  
  ip_baits = qc.getIpToBaits(ip_matrix)
  if(cluster){
    qc.clusterHeatmap(data_matrix, gsub('.txt','.pdf',matrix_file), ip_baits, font_scale)  
  }
  if(ip_dists){
    qc.ipDists(data_matrix, ip_baits, matrix_file)
    qc.NumUniquePlot(ip_matrix, matrix_file)
  }
}

if(!exists("PIPELINE") || PIPELINE==F){
  option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
                dest="verbose", help="Print little output"),
    make_option(c("-d", "--data_file"),
                help="data file containing matrix"),
    make_option(c("-o", "--output_file"),
                help="output file for cluster plot"),
    make_option(c("-s", "--font_scale"), default=40,
                help="scaling factor for fonts on the rows and columns")
  )
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
}
## TODO: make the following code into a unit-test 
# config = yaml.load(string=paste(readLines("tests/APMS_TEST.yml"),collapse='\n'))
# config = yaml.load(string=paste(readLines("tests/entero/APMS_ENTERO.yml"),collapse='\n'))
# qc.main(matrix_file=config$qc$matrix_file, font_scale=config$qc$cluster_font_scale, cluster=config$qc$cluster, ip_dists=config$qc$ip_distributions)
