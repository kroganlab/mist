#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GO.db))
suppressMessages(library(KEGG.db))
suppressMessages(library(PFAM.db))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))



# Get GO terms for CC, MF, & BP
overRepresented.getGOTerms <- function(go_ids){
  go_terms = unique(select(GO.db, keys=go_ids$GO, columns="TERM", keytype="GOID"))
  res = merge(go_ids, go_terms, by.x="GO", by.y="GOID")
  return(res)
}

# Get PFAM terms
overRepresented.getPFAMTerms <- function(go_ids){
  PF = PFAMDE
  mapped_keys = mappedkeys(PF)
  PF = as.data.frame(PF[mapped_keys])
  res = merge(go_ids, PF, by.x='PFAM', by.y='ac')  
  names(res)[grep('^de$', names(res))] = "TERM"
  return(res)
}

# get KEGG pathway names associated with PREY
overRepresented.getKEGGTerms <- function(go_ids){
  KEGG = org.Hs.egPATH
  mapped_keys = mappedkeys(KEGG)
  KEGG_Entrez = as.data.frame(KEGG[mapped_keys])
  
  KEGG = KEGGPATHID2NAME
  mapped_keys = mappedkeys(KEGG)
  KEGG_Terms = as.data.frame(KEGG[mapped_keys])
  
  KEGG = unique(merge(KEGG_Entrez, KEGG_Terms, by='path_id'))
  res = unique(merge(go_ids, KEGG, by.x='ENTREZID', by.y='gene_id'))
  names(res)[9] = "TERM"
  return(res)
}


toLongForm <- function(dat, termid, out_file){
  print("CONVERTING TO LONG FORM")
  # cycle through each term in each bait and explode it based on the uniprot_ac
  n = dim(dat)[1]
  
  prey_list = c()
  ID_list = c()
  
  for(i in 1:n){
    if(i %% 100 == 0) print(paste(i,n,sep="/"))
    r = dat[i,]
    tmp = unlist(strsplit(r$UNIPROT,","))
    prey_list = c(prey_list, tmp)
    ID_list = c(ID_list, rep(r[,termid], length(tmp)))
  }
  
  ID_to_prey = data.frame(PREY=prey_list,ID=ID_list)
  ID_to_prey = unique(ID_to_prey)
  TERMS_l = merge(dat[,c("setid",termid),], ID_to_prey, by.x=termid, by.y="ID")
  colnames(TERMS_l) =c("ID","setid","Term","PREY") 
  write.table(TERMS_l, file=out_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  TERMS_l
}
toLongForm_C = cmpfun(toLongForm)


getSampleDistros <- function (scores, unique_counts, iterations) {
  set.seed(seed=23)
  sample_params =c()
  for(cnt in unique_counts){
    sample_dist = c()
    for(i in 1:iterations){
      dummy_score = mean(sample(scores, cnt, replace=F))   
      sample_dist = c(sample_dist, dummy_score)
    }
    sample_m = mean(sample_dist)
    sample_s = sd(sample_dist)
    sample_params = rbind(sample_params,c(cnt,sample_m, sample_s))
  }
  sample_params = as.data.frame(sample_params)
  colnames(sample_params) = c("size","mean","sd")
  sample_params
}

getScoreProbabilities <- function(term_groups, sample_params) {
  p_values_adj = c()
  p_values = c()
  for(s in unique(term_groups$BAIT)){
    subset = term_groups[term_groups$BAIT==s,]
    p_tmp = c()
    
    for(term_row in 1:nrow(subset)){
      term = subset[term_row,]
      term_m = sample_params[sample_params$size == term$count,"mean"][[1]]
      term_s = sample_params[sample_params$size == term$count,"sd"][[1]]
      term_p = pnorm(term$score_avg, mean=term_m, sd=term_s, lower.tail=F)
      p_tmp = c(p_tmp,term_p)
    }
    p_values = c(p_values, p_tmp)
    p_values_adj = c(p_values_adj, p.adjust(p_tmp, method="fdr"))
  }
  data.frame(term_groups=term_groups, p_values=p_values, q_values=p_values_adj)
}

scoreMatrix <- function(term_groups_selected, HEAT_STYLE="Q") {
  
  ## recover values to build term x bait matrix
  unique_sets = unique(term_groups_selected$setid)
  unique_terms = unique(term_groups_selected$Term)
  
  if(HEAT_STYLE=="Q"){
    value.var = "q_values_lt"
    term_groups_selected = data.frame(term_groups_selected, q_values_lt = -log10(as.double(term_groups_selected$q_values)))
  }else if(HEAT_STYLE=="SCORE"){
    value.var = "score_avg"
  }
  fun.aggregate=max
  term_groups_selected_w = dcast(data=term_groups_selected, formula=Term~setid, value.var=value.var, fun.aggregate=fun.aggregate, na.rm=T)  
  rownames(term_groups_selected_w) = term_groups_selected_w[,1]
  term_groups_selected_w=term_groups_selected_w[,-1]
  term_groups_selected_w = data.matrix(term_groups_selected_w)
  term_groups_selected_w[is.na(term_groups_selected_w) | is.infinite(term_groups_selected_w)]=0
  term_groups_selected_w
}

## following code will print column labels at 45 degrees
library(grid)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


resample.main = function(termsNscores, out_file, termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T=0.7, score_name="MIST", ITERATIONS, P_T, BREAKS, LABELS){
  term_src = data.frame(termsNscores)
  
  ## select only list of proteins that are on the interactome
  if(FILTER_HICONF){
      term_src = term_src[ term_src[,score_name]>SCORE_T, ]
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HERE!
# need to include "uniprot_id" into the new mist scores or else remove them from here!!

  ## computing bait-term aggregate stats (mean/median)
  tmp = term_src[,c("BAIT","TERM","uniprot_id",score_name)]
  tmp_ag_scores = aggregate(formula(paste(score_name,"~.", sep="")), tmp[,c("BAIT","TERM",score_name)], FUN=AGGREGATE_FUN)
  tmp_ag_preys = aggregate(uniprot_id ~ ., tmp[,c("BAIT","TERM","uniprot_id")], FUN=unique)
  tmp_ag_preys$uniprot_id = gsub('c|\\(|\\)|[\"\']|\n','',tmp_ag_preys$uniprot_id)
  tmp_ag_count = aggregate(uniprot_id ~ ., tmp[,c("BAIT","TERM","uniprot_id")], FUN=length)
  term_groups = merge(tmp_ag_scores, tmp_ag_preys, by=c("BAIT","TERM"))
  term_groups = merge(term_groups, tmp_ag_count, by=c("BAIT","TERM"))
  colnames(term_groups) = c('BAIT','TERM','score_avg','preys_grouped','count')
  term_groups = term_groups[term_groups$count > 1, ] ## remove singletons
  cat(sprintf("SELECTING %s NON-SINGLETON BAIT-TERM PAIRS\n",nrow(unique(term_groups[,c("BAIT","TERM")]))))
  
  ## compute sample distributions
  SAMPLE_DISTROS = getSampleDistros(scores=termsNscores[[score_name]], unique_counts=unique(term_groups$count), iterations=ITERATIONS) 
  
  ## compute score probabilities
  term_groups_tmp = c()
  for(s in unique(term_groups$BAIT)){
    tmp = getScoreProbabilities(term_groups[term_groups$BAIT==s,], SAMPLE_DISTROS)
    term_groups_tmp = rbind(term_groups_tmp, tmp)
  }
  term_groups = term_groups_tmp
  colnames(term_groups) = c("setid","Term","score_avg","preys_grouped","count","p_values","q_values")
  write.table(term_groups, file=out_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  
  ## format into data matrix with only significant hits
  term_groups_selected_w = scoreMatrix(term_groups_selected=term_groups,HEAT_STYLE=heat_style)
  mask = apply(term_groups_selected_w,2,function(x)x>=-log10(P_T)) ## select everything which is significant (log10 base)
  colmask = unname(colSums(mask)) > 0
  rowmask = unname(rowSums(mask)) > 0
  term_groups_selected_w = term_groups_selected_w[rowmask,colmask]
  colnames(term_groups_selected_w) = gsub('_','.',colnames(term_groups_selected_w))
  
  ## format stuff for heatmap
  colors = c(colorRampPalette(brewer.pal(n = 7, name = "Reds"))(length(BREAKS)-1))
  LOWER_T = BREAKS[length(BREAKS)]
  term_groups_selected_w_display = term_groups_selected_w
  term_groups_selected_w_display[term_groups_selected_w_display>LOWER_T]=LOWER_T
  pheatmap(term_groups_selected_w_display, cellheight=10, cellwidth=10, scale="none", filename=gsub('.txt','_heatmap.pdf',out_file), fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, color = colors, breaks=BREAKS, legend_breaks=BREAKS,legend_labels=LABELS, fontfamily="Helvetica")
  term_groups_selected_w
}


# Get GO, PFAM, and KEGG terms for preys
overRepresented.main <- function(score_file, output_dir, uniprot_dir="/files/", score_name="MIST_hiv", prey_name="Prey", bait_name="Bait", score_threshold=0.7){
  cat("\tMapping GO Terms...\n")
  scores = data.table(read.delim(score_file, sep="\t", stringsAsFactors=F))
  setnames(scores, prey_name, "UNIPROT")  # change colnames to use merge.data.table for speeeeed
  setnames(scores, bait_name, "BAIT")
  ids = unique(scores[,UNIPROT])
  
  # Check if uniprot_id is in scores
  if( !any(grepl("uniprot_id", names(scores))) & any(grepl("Entry.name", names(scores))) ){
    setnames(scores, "Entry.name", "uniprot_id")
  }else if( !any(grepl("uniprot_id", names(scores))) & !any(grepl("Entry.name", names(scores))) ){
    uniprot_file = paste(uniprot_dir, "uniprot_protein_descriptions_HUMAN.txt", sep="")
    genes = tryCatch(data.table(read.delim(uniprot_file, header=TRUE, sep="\t", stringsAsFactors=F)), 
                     error = function(e){ stop(paste(uniprot_file,'not found. Please check to make sure it exists.'))})
    cat("\tNeither 'uniprot_id' nor 'Entry.name' column found. Loading id's from file.\n")
    setnames(genes, "Entry", "UNIPROT")
    scores = merge(scores, genes[,list(UNIPROT,Entry.name)], by="UNIPROT")
    setnames(scores, "Entry.name", "uniprot_id")
  }
  
  
  #~~~~~~~ GET ENRICHMENT TERMS ~~~~~~~
  # Get Entrez ID's etc for mapping to terms
  go_ids = select(org.Hs.eg.db, keys = ids, columns=c("GO", "PFAM", "SYMBOL", "UNIPROT", "ENTREZID"), keytype="UNIPROT")
  # Get GO terms for CC, MF, BP
  GO_terms = overRepresented.getGOTerms(go_ids)
  
  # GO-CC
  cat("\tFinding GO Cellular Complex terms...\n")
  out_file = paste(output_dir,"GO_CC_Terms.txt", sep="")
  goCC = merge(scores, data.table(GO_terms[GO_terms$ONTOLOGY=="CC",]) , by="UNIPROT", allow.cartesian=T, all.x=T)
  write.table(goCC, out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  # GO-BP
  cat("\tFinding GO Biological Process terms...\n")
  out_file = paste(output_dir,"GO_BP_Terms.txt", sep="")
  goBP = merge(scores, data.table(GO_terms[GO_terms$ONTOLOGY=="BP",]) , by="UNIPROT", allow.cartesian=T, all.x=T)
  write.table(goBP, out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  # GO-MF
  cat("\tFinding GO Molecular Function terms...\n")
  out_file = paste(output_dir,"GO_MF_Terms.txt", sep="")
  goMF = merge(scores, data.table(GO_terms[GO_terms$ONTOLOGY=="MF",]) , by="UNIPROT", allow.cartesian=T, all.x=T)
  write.table(goMF, out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  # PFAM
  cat("\tFinding PFAM terms...\n")
  out_file = paste(output_dir,"PFAM_Terms.txt", sep="")
  pfam = merge(scores, data.table(overRepresented.getPFAMTerms(go_ids)), by="UNIPROT", allow.cartesian=T, all.x=T)
  write.table(pfam, out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    
  # KEGG
  cat("\tFinding KEGG terms...\n")
  out_file = paste(output_dir,"KEGG_Terms.txt", sep="")
  kegg = merge(scores, overRepresented.getKEGGTerms(go_ids), by="UNIPROT", allow.cartesian=T, all.x=T)
  write.table(overRepresented.getKEGGTerms(go_ids), out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
  #~~~~~~~ CALCULATE OVER-REPRESENTATION OF PROTEINS ~~~~~~~
  cat("\tCalculating over-representation of proteins:\n")
  ## CONFIGURATION
  ITERATIONS = 1000                                                                                                     # <----
  #databases =                                                                                                          # <----
  HEAT_STYLE="Q" ##HEAT_STYLE="SCORE"
  P_T = 0.05
  BREAKS = c(0,-log10(c(P_T,P_T/10,P_T/100,P_T/1000)))
  LABELS = c(0,"  5E-2","  5E-3","< 5E-4","")
  termid="TERM"   # column name of the Terms column 
  
  # GO-CC Resampling
  cat("\t  GO-CC terms... ")
  out_file = paste(output_dir,"GO_CC_Terms_enriched.txt", sep="")
  term_groups_selected_w = resample.main(termsNscores=goCC, out_file=out_file, termid=termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T= score_threshold, score_name, ITERATIONS, P_T, BREAKS, LABELS)
  write.table(term_groups_selected_w, file=gsub('.txt','_heatmap.txt',out_file), eol='\n', sep='\t',quote=F, col.names=NA, row.names=T)

  # GO-BP Resampling
  cat("\t  GO-BP terms... ")
  out_file = paste(output_dir,"GO_BP_Terms_enriched.txt", sep="")
  term_groups_selected_w = resample.main(termsNscores=goBP, out_file=out_file, termid=termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T=score_threshold, score_name, ITERATIONS, P_T, BREAKS, LABELS)
  write.table(term_groups_selected_w, file=gsub('.txt','_heatmap.txt',out_file), eol='\n', sep='\t',quote=F, col.names=NA, row.names=T)
  
  # GO-MF Resampling
  cat("\t  GO-MF terms... ")
  out_file = paste(output_dir,"GO_MF_Terms_enriched.txt", sep="")
  term_groups_selected_w = resample.main(termsNscores=goMF, out_file=out_file, termid=termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T=score_threshold, score_name, ITERATIONS, P_T, BREAKS, LABELS)
  write.table(term_groups_selected_w, file=gsub('.txt','_heatmap.txt',out_file), eol='\n', sep='\t',quote=F, col.names=NA, row.names=T)
  
  # PFAM Resampling
  cat("\t  PFAM terms... ")
  out_file = paste(output_dir,"PFAM_Terms_enriched.txt", sep="")
  term_groups_selected_w = resample.main(termsNscores=pfam, out_file=out_file, termid=termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T=score_threshold, score_name, ITERATIONS, P_T, BREAKS, LABELS)
  write.table(term_groups_selected_w, file=gsub('.txt','_heatmap.txt',out_file), eol='\n', sep='\t',quote=F, col.names=NA, row.names=T)
  
  # KEGG Resampling
  cat("\t  KEGG terms... ")
  out_file = paste(output_dir,"KEGG_Terms_enriched.txt", sep="")
  term_groups_selected_w = resample.main(termsNscores=kegg, out_file=out_file, termid=termid, heat_style="Q", FILTER_HICONF=T, AGGREGATE_FUN=mean, SCORE_T=score_threshold, score_name, ITERATIONS, P_T, BREAKS, LABELS)
  write.table(term_groups_selected_w, file=gsub('.txt','_heatmap.txt',out_file), eol='\n', sep='\t',quote=F, col.names=NA, row.names=T)
  
  
}


# score_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/processed/kinases_data_wKEYS_MAT_ALLSCORES.txt"
# score_name = "MIST_hiv"
# output_dir = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/processed/Enrichment/"


# score_file = "~/projects/mist/tests//small/processed/preprocessed_NoC_MAT_ALLSCORES.txt"
# score_name = "MIST"
# prey_name = "Prey"
# bait_name = "Bait"
# output_dir = "~/projects/mist/tests//small/processed/Enrichment/"
# 
# 
# overRepresented.main(score_file, output_dir, scores_name, prey_name, bait_name)









