library(sqldf)
library("Biostrings")

remove_negative_baits = function(df,pattern='negative|NEGATIVE|Negative',col='BAIT'){
  df[grep(pattern,df[,col],invert=T),]
}

# filter_on_negative_scores = function(df=test, mist_hiv_neg=F, comppass_wd_neg=T, mist_new=T){
#   df_neg = df[grep('NEGATIVE|negative|negative_flag|negative_strep|control|vector|Vector',df$BAIT), ]
#   df_neg = sqldf("select PREY,max(MIST_hiv) as MIST_hiv , max(COMPPASS_WD) as COMPPASS_WD, max(newmistscore) as newmistscore from df_neg group by PREY")
#   if(mist_hiv_neg){
#     df = sqldf("select D.* from df D join df_neg DN on D.PREY=DN.PREY where D.MIST_hiv > DN.MIST_hiv group by D.BAIT, D.PREY")
#   }
#   if(comppass_wd_neg){
#     df= sqldf("select D.* from df D join df_neg DN on D.PREY=DN.PREY where D.COMPPASS_WD > DN.COMPPASS_WD group by D.BAIT, D.PREY")
#   }
#   if(mist_new){
#     df= sqldf("select D.* from df D join df_neg DN on D.PREY=DN.PREY where D.newmistscore > DN.newmistscore group by D.BAIT, D.PREY")
#   }
#   df
# }

## OLD VERSION, previous doesnt work for now
filter_on_negative_scores = function(df, mist_hiv_neg=T, comppass_wd_neg=T){
  
  if(mist_hiv_neg){
    #df = df[is.na(df$MIST_hiv_negative) | df$MIST_hiv > df$MIST_hiv_negative, ]  
    df = df[df$MIST_hiv > df$MIST_hiv_negative, ]  
  }
  if(comppass_wd_neg){
    #df = df[is.na(df$COMPPASS_WD_NEGATIVE) | df$COMPPASS_WD > df$COMPPASS_WD_NEGATIVE, ]
    df = df[df$COMPPASS_WD > df$COMPPASS_WD_NEGATIVE, ]
  }
  df
}

filter_on_negative_features = function(df=test){
  df_neg = df[grep('NEGATIVE|negative|negative_flag|negative_strep|control|vector|Vector',df$BAIT), ]
  df_neg = sqldf("select PREY, max(MIST_R) as MIST_R, max(MIST_A) as MIST_A from df_neg group by PREY")
  tmp = sqldf("select D.* from df D join df_neg DN on D.PREY=DN.PREY where D.MIST_R > DN.MIST_R or D.MIST_A > DN.MIST_A group by D.BAIT, D.PREY")
  tmp
}

addBaitRank = function(df){
  df = transform(df, BAIT_mist_rank = ave(MIST_hiv, BAIT, FUN = function(x) rank(-x, ties.method = "first")))
  df = transform(df, BAIT_comppass_rank = ave(COMPPASS_WD, BAIT, FUN = function(x) rank(-x, ties.method = "first")))
  df = transform(df, BAIT_saint_rank = ave(SAINT_AVG_P, BAIT, FUN = function(x) rank(-x, ties.method = "first")))
  df
}

simplifyColumns = function(df, cols=c('BAIT','PREY','MIST_hiv','MIST_A','MIST_R','MIST_S','COMPPASS_WD','SAINT_AVG_P','TSC_AVG','uniprot_id', 'gene_name',  'description', 'species')){
  df[, cols]
}

combineSetsOnScore = function(df1, df2, m1='m1', m2='m2', max_score='MIST_hiv'){
  df = rbind(data.frame(df1, machine=m1), data.frame(df2, machine=m2))
  ##  remove doubles based on mist score
  tmp = sqldf(sprintf("select BAIT, PREY, max(%s) as 'max_score' wfrom df group by BAIT, PREY",max_score))
  tmp2 = sqldf(sprintf("select D.* from tmp T join df D on T.BAIT = D.BAIT and T.PREY=D.PREY and T.max_score = D.%s group by BAIT, PREY",max_score))
  tmp2
}

combineSetsOnMachine = function(df1, df2, m1='m1', m2='m2', prefered_machine='m2'){
  df = rbind(data.frame(df1, machine=m1), data.frame(df2, machine=m2))
  ##  remove doubles based on mist score
  preferred = df[df$machine == prefered_machine,]
  other = df[df$machine != prefered_machine,]
  rest = sqldf("select * from other O where O.BAIT not in (select distinct(P.BAIT) from preferred P)")
  res = rbind(preferred, rest)
}

filter_on_baits = function(df, re, re_invert=F){
  filtered_re = df[grep(re, df$BAIT, invert=re_invert),]
  print(sprintf("filtered baits : %s",unique(filtered_re$BAIT)))
  filtered_re
}

threshold_on_fixed_condition = function(df, condition){
  ## select
  print(condition)
  tmp = sqldf(sprintf("select * from df D where %s",condition))
  tmp
}

df_threshold_bait_stats = function(df){
  tmp = sqldf("select BAIT, count(*) as 'interactors', type from df group by type, BAIT")
  tmp
}

add_host_interactions = function(df_T, host_ppi, restrict_to_bait=F){
  df_T = data.frame(df_T, type='APMS')
  if(restrict_to_bait){
    bait_clause = "on (D1.BAIT=D2.BAIT)"
  }else{
    bait_clause = ""
  }
  host_tmp2 = sqldf(sprintf("select H.uniprot_ac_1 as 'BAIT', H.uniprot_ac_2 as 'PREY', H.complex_name as 'description', 'host_ppi' as 'type', D1.BAIT, D2.BAIT from df_T D1 join df_T D2 %s join host_ppi H on (((H.uniprot_ac_1 = D1.PREY) and (H.uniprot_ac_2 = D2.PREY)) or ((H.uniprot_ac_1 = D2.PREY) and (H.uniprot_ac_2 = D1.PREY))) group by H.uniprot_ac_1, H.uniprot_ac_2", bait_clause))  
  ## DIRTY: make DF with the same size as host_tmp2 but with the column scheme of DF_T, set all values to NA and then copy the columns with values fron host_tmp2 we're interested in
  if(nrow(host_tmp2)>0){
    host_tmp = df_T[1:nrow(host_tmp2),]
    host_tmp[,] = NA
    host_tmp$BAIT = host_tmp2$BAIT
    host_tmp$PREY = host_tmp2$PREY
    host_tmp$description = host_tmp2$description
    host_tmp$type = host_tmp2$type
    res = rbind(df_T, host_tmp)
    res  
  }else{
    df_T
  }
}

networkStats = function(df_T_hh, df_T_hh_gold, df_2){
  apms_p_nodes = sqldf("select BAIT,PREY from df_T_hh where type='APMS' group by BAIT,PREY")
  apms_b_nodes_num = sqldf("select count(distinct(BAIT)) from df_T_hh where type='APMS'")[1,1]
  apms_p_nodes_num = nrow(apms_p_nodes)
  
  apms_bp_edges_num = sqldf("select count(*) from df_T_hh where type='APMS'")[1,1]
  apms_bp_g_edges_num = sqldf("select count(*) from df_T_hh_gold")[1,1]
  apms_pp_edges_num = sqldf("select count(*) from df_T_hh where type='host_ppi'")[1,1]
  
  median_bp_edges_num = median(sqldf("select count(*) from df_T_hh where type='APMS' group by BAIT")[,1])
  max_bp_edges_num = max(sqldf("select count(*) from df_T_hh where type='APMS' group by BAIT")[,1])
  min_bp_edges_num = min(sqldf("select count(*) from df_T_hh where type='APMS' group by BAIT")[,1])
  
  df_2_p_nodes = sqldf("select BAIT,PREY from df_2 group by BAIT,PREY")
  df_huh_293_bp_overlap = merge(apms_p_nodes, df_2_p_nodes)
  df_huh_293_p_overlap = intersect(apms_p_nodes$PREY, df_2_p_nodes$PREY)
  
  res = data.frame(apms_baits=apms_b_nodes_num,apms_preys=length(unique(apms_p_nodes$PREY)), apms_edges=apms_bp_edges_num, apms_edges_gold=apms_bp_g_edges_num, host_host_edges=apms_pp_edges_num, median_apms_edges=median_bp_edges_num, max_apms_edges=max_bp_edges_num, min_apms_edges=min_bp_edges_num, edges_overlap_huh_293=nrow(df_huh_293_bp_overlap) , prey_overlap_huh_293=length(df_huh_293_p_overlap))
  res
}


# 
# thresholdNetworkAndHostPPI = function(df, df_host, condition, nw_out, nw_stats_out){
#   df_bait_stats_total = c()
#   df_threshold = threshold_on_fixed_condition(df, condition)
#   df_threshold = data.frame(df_threshold, type='APMS')
#   df_hh = add_host_interactions(df_threshold, df_host)
#   df_bait_stats = data.frame(df_threshold_bait_stats(df_hh), condition=nw_out)
#   write.table(df_hh, file=nw_out, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
#   write.table(df_bait_stats, file=nw_stats_out, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
#   df_hh
# }

interactionsDecected = function(dataset, bait_prey_uniprot_goldset, merge_conditions="PREY_ONLY"){
  if(merge_conditions=="PREY_ONLY"){
    tmp = sqldf("select * from dataset D join bait_prey_uniprot_goldset G on D.PREY=G.PREY group by PREY")  
  }else if(merge_conditions=="BAIT_PREY"){
    tmp = merge(df_T[,c("BAIT","PREY")], bait_prey_uniprot_goldset, by=c("BAIT","PREY"), all=T)
  }
  tmp[,colnames(dataset)]
}

filterContaminants = function(contaminant_fasta, df, col='ms_uniprot_ac'){
  ## use biostrings function to read in fasta file
  fastaF = as.data.frame(readAAStringSet(filepath=contaminant_fasta))
  cont_df = data.frame(name=rownames(fastaF), seq=fastaF$x)
  ## split out the first part of the fasta header which should be the identifier -- split could be eithe space or pipe
  fastaFFirstNames = as.data.frame(apply(cont_df,1,FUN=function(x) as.character(strsplit(x[1],split="\\s|\\|")[[1]][1])))
  ## check whether the key is in uniprot format
  UniprotNames = as.data.frame(fastaFFirstNames[grep('^[A-Z]{1}[0-9]{1}[A-Z,0-9]{4}(\\-1)*',fastaFFirstNames[,1]),])
  colnames(UniprotNames)[1] = 'contaminant_key'
  df_F = df[!(df[,col] %in% UniprotNames$contaminant_key),]
  print(sprintf("FILTERED %s CONTAMINANTS", nrow(df)-nrow(df_F)))
  df_F
}
