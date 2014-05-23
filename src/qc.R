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
  prey_names = rownames(input_matrix)[data_row_idx:nrow(input_matrix)]
  data_matrix = input_matrix[data_row_idx:nrow(input_matrix),data_col_idx:ncol(input_matrix)]
  rownames(data_matrix) = input_matrix[data_row_idx:nrow(input_matrix),prey_col_idx]
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
  colnames(data_matrix_w_names) = do.call(paste, c(ip_baits, sep = " "))
  cor_matrix = cor(data_matrix_w_names, use="pairwise.complete.obs", method="pearson")
  color_scale = colorRampPalette(brewer.pal(7,"RdBu"))(9)
  pheatmap(cor_matrix, cluster_rows=T, cluster_cols=T, scale="none",fontsize_row=font_scale,fontsize_col=font_scale, cellwidth=font_scale, cellheight=font_scale, border_color=NA, filename=output_file, color=color_scale, breaks=1:10/10, clustering_distance_rows="correlation", clustering_distance_cols="correlation", treeheight_row=0, treeheight_col=0)  
}

qc.ipDists = function(data_matrix, ip_baits, baseName){
  tmp = melt(data_matrix, varnames=c("prey","ip"), value.name="count")
  data_long = merge(tmp, ip_baits, by="ip")
  data_long = data_long[data_long$count > 0,]
  theme_set(theme_bw(base_size = 18,base_family='Helvetica'))
  
  pdf(gsub('.txt','_proteincounts.pdf',baseName), width=10, height=10)
  print(ggplot(data_long, aes(x=ip)) + geom_bar(stat='bin') + facet_wrap(facets= ~bait, scales='free_x', nrow=5))
  dev.off()
  
  ## solution with ggplot
  #ggplot(data_long, aes(x=log2(count), colour=factor(replicate))) + geom_density(alpha=1, size=1) + facet_wrap(facets= ~bait, scales='fixed') 
  
  ## solution to give every plot its separate legend with gridExtra library 
  pdf(gsub('.txt','_peptidecount_dists.pdf',baseName), width=10, height=10)
  out = by(data = data_long, INDICES = data_long$bait, FUN = function(m) {
    m = droplevels(m)
    m = ggplot(m, aes(x=log2(count), colour=ip)) + geom_density(alpha=1, size=1)
  })
  do.call(grid.arrange, out)
  dev.off()
}

qc.main = function(matrix_file, font_scale, cluster=T, ip_dists=T){
  ip_matrix = read.delim(matrix_file, stringsAsFactors=F)
  data_matrix = qc.dataMatrix(ip_matrix)
  ip_baits = qc.getIpToBaits(ip_matrix)
  if(cluster){
    qc.clusterHeatmap(data_matrix, gsub('.txt','.pdf',config$qc$matrix_file), ip_baits, font_scale)  
  }
  if(ip_dists){
    qc.ipDists(data_matrix, ip_baits, config$qc$matrix_file)
  }
}

if(is.null(PIPELINE)){
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
config = yaml.load(string=paste(readLines("tests/APMS_TEST.yml"),collapse='\n'))
qc.main(matrix_file=config$qc$matrix_file, font_scale=config$qc$cluster_font_scale, cluster=config$qc$cluster, ip_dists=config$qc$ip_distributions)