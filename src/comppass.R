#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(reshape2))
suppressMessages(library(optparse))
suppressMessages(library(compiler))
suppressMessages(library(stats))

Comppass.cleanMatrix = function(data){
  ## get rid of unassigned columns
  data_tmp = data[,colnames(data) != ""]
  ## get rid of specificity exclusion row and meta-columns
  data_tmp = data_tmp[-1,c(-2:-4)]
  colnames(data_tmp)[1] = "Preys"
  data_tmp
}

Comppass.cnames = function(data){
  cnames = colnames(data)[colnames(data)!=""]
  cnames = cnames[c(-1:-4)]
  cnames
}

## dirty function to clean up matrix and convert it to long format
Comppass.convertMatrix = function(data, cnames){
  ## make clean version of #replicated bait names w/o R re-naming
  cnames_rep = rep(cnames, each=nrow(data))
  ## convert into long format so we can later re-convert to wide format with an aggregate function
  # print(colnames(data))
  data_l = melt(data, id=c("Preys"),  measure.vars=colnames(data)[2:ncol(data)])
  ## put the clean version of of #replicated bait names
  data_l$variable = cnames_rep
  ## make into numbers
  data_l$value = as.numeric(data_l$value)
  ## filter out 0's so that the mean gets computed correctly over # observations instead of # matrix cols   
  data_lf = data_l[data_l$value>0,]
  data_lf
}

Comppass.StatsTable = function(data_long){
  #stats_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=mean, fill=0)
  stats_tab = tryCatch( dcast(Preys ~ variable, data=data_long, fun.aggregate=mean, fill=0),
                        error = function(e){ 
                          detach("package:reshape2")
                          suppressMessages(library(reshape))
                          temp = cast(Preys ~ variable, data=data_long, fun.aggregate=mean, fill=0)
                          detach("package:reshape")
                          suppressMessages(library(reshape2))
                          retun(temp)
                        } ) #trycatch for when dataset is too large -> segfault
  rownames(stats_tab) = stats_tab[,1]
  stats_tab = as.matrix(stats_tab[,2:ncol(stats_tab)])
  stats_tab
}

Comppass.SpeciTable = function(stats_tab){
  stats_msk = stats_tab
  stats_msk[stats_msk>0]=1
  stats_msk
}

Comppass.ReproTable = function(data_long){
  repro_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=length, fill=0)
  repro_tab = tryCatch( dcast(Preys ~ variable, data=data_long, fun.aggregate=length, fill=0),
                        error = function(e){ 
                          detach("package:reshape2")
                          suppressMessages(library(reshape))
                          temp = cast(Preys ~ variable, data=data_long, fun.aggregate=length, fill=0)
                          detach("package:reshape")
                          suppressMessages(library(reshape2))
                          retun(temp)
                        } ) #trycatch for when dataset is too large -> segfault
  rownames(repro_tab) = repro_tab[,1]
  repro_tab = as.matrix(repro_tab[,2:ncol(repro_tab)])
  repro_tab
}

Comppass.Z = function(stats_tab){
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  Z = (stats_tab - m) / sd
  Z
}

Comppass.S = function(stats_tab, speci_tab){
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  S = sqrt( inv_freq * stats_tab )
  S
}

Comppass.D = function(stats_tab, speci_tab, repro_tab, normalized=F, D_T=1){ ## 95% score as normalization bar
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  tmp = stats_tab
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] * (inv_freq[i] ^ (repro_tab[i,j]))
    }
  }
  D = sqrt(tmp)
  if(normalized){
    D_threshold = quantile(D, probs=D_T)
    D / D_threshold
  }else{
    D
  }
}

Comppass.WD = function(stats_tab, speci_tab, repro_tab, normalized=F, WD_T=1){ ## 95% score as normalization bar
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  w = sd / m
  w[w<=1]=1
  
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  tmp = stats_tab
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] * (w[i] * inv_freq[i] ^ (repro_tab[i,j]))
    }
  }
  WD = sqrt(tmp)
  if(normalized){
    WD_threshold = quantile(WD, probs=WD_T)
    WD / WD_threshold
  }else{
    WD
  }
}

Comppass.WDv2 = function(stats_tab, speci_tab, repro_tab, num_reps, normalized=F, WD_T=1){ ## 95% score as normalization bar
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  w = sd / m
  w[w<=1]=1
  
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  tmp = stats_tab

  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] * (w[i] * inv_freq[i] ^ (repro_tab[i,j]/num_reps[j]))
     }
  }
  WD = sqrt(tmp)
  if(normalized){
    WD_threshold = quantile(WD, probs=WD_T)
    WD / WD_threshold
  }else{
    WD
  }
}
# 
# Comppass.Summary = function(stats_tab, Z, S, D, WD, WD_T=0){ ## WD_T takes everything over this  score
#   
#   score_table = c()
#   key_table = c()
#   
#   for(i in 1:nrow(stats_tab)){
#     for(j in 1:ncol(stats_tab)){
#       if(WD[i,j] >= WD_T){
#         bait = colnames(stats_tab)[j]
#         prey = rownames(stats_tab)[i]
#         z_score = Z[i,j]
#         s_score = S[i,j]
#         d_score = D[i,j]
#         wd_score = WD[i,j]
#         tsc = stats_tab[i,j]
#         key_table = rbind(key_table, c(bait, prey))
#         score_table = rbind(score_table, c(tsc, z_score, s_score, d_score, wd_score))
#       }
#     }
#   }
#   total_table = cbind(key_table, score_table)
#   colnames(total_table) = c("Bait","Prey","Abundance","Z","S","D","WD")
#   total_table
# }

Comppass.Summary = function(stats_tab, Z, S, D, WD, WDv2, WD_T=0){ 
  baits = colnames(stats_tab)
  
  stats_m = cbind(rownames(stats_tab),stats_tab)
  colnames(stats_m)[1]="Prey"
  stats_m = melt(stats_tab, id.vars=baits)
  
  Z_m = cbind(rownames(Z),Z)
  colnames(Z_m)[1]="Prey"
  Z_m = melt(Z, id.vars=baits)
  
  S_m = cbind(rownames(S),S)
  colnames(S_m)[1]="Prey"
  S_m = melt(S, id.vars=baits)
  
  D_m = cbind(rownames(D),D)
  colnames(D_m)[1]="Prey"
  D_m = melt(D, id.vars=baits)
  
  WD_m = cbind(rownames(WD),WD)
  colnames(WD_m)[1]="Prey"
  WD_m = melt(WD, id.vars=baits)
  
  WDv2_m = cbind(rownames(WDv2),WDv2)
  colnames(WD_m)[1]="Prey"
  WDv2_m = melt(WDv2, id.vars=baits)
  
  M = cbind(stats_m[,2],stats_m[,c(1,3)], Z_m[,3], S_m[,3], D_m[,3], WD_m[,3], WDv2_m[,3])
  colnames(M) = c("Bait","Prey","Abundance","Z","S","D","WD","WDv2")  
  M
}

Comppass.Summary = cmpfun(Comppass.Summary, options=list("optimize",3))  

Comppass.ResampledScreen = function(data_tmp){
  ## count # identified proteins in each run
  protein_cnts = apply(data.matrix(data_tmp[,2:ncol(data_tmp)]),2,function(x)sum(x>0))
  ## count # identified TSC in each run
  tsc_counts = apply(data.matrix(data_tmp[,2:ncol(data_tmp)]), 2, sum)
  ## cleanup
  Preys = data_tmp$Preys
  rownames(data_tmp) = Preys
  data_tmp = data_tmp[,-1]
  data_tmp = data.matrix(data_tmp, rownames.force=T)
  ## get TOTAL # of times prey identified
  prey_totals = apply(data_tmp, 1, sum)
  ## make random proteome 
  random_proteome = rep(rownames(data_tmp), times=prey_totals)
  random_proteome_size = length(random_proteome)
  
  ## initialize random screen with 0 observations
  random_screen = data_tmp
  random_screen[] = 0
  
  for(i in 1:ncol(random_screen)){
    
    ## make a random run
    p=0
    tsc=0
    tsc_total = tsc_counts[i]
    protein_total = protein_cnts[i]
    
    while(p < protein_total | tsc < tsc_total){
      
      ## get a random protein and add it to run
      random_protein_idx = sample(1:random_proteome_size, 1)
      random_protein = random_proteome[random_protein_idx]
     
      current_tsc_count = random_screen[random_protein,i]
      if(current_tsc_count==0){ ## first time we sample this protein
        p = p + 1 
      }
      random_screen[random_protein,i] = current_tsc_count + 1
      tsc = tsc + 1
    }  
  }
  random_screen = cbind(Preys, random_screen)
  colnames(random_screen)[1] = "Preys"
  as.data.frame(random_screen)
}
Comppass.ResampledScreen = cmpfun(Comppass.ResampledScreen, options=list("optimize",3))

Comppass.ResampledPvalues = function(data_tmp, cnames, summary){
  random_screen = Comppass.ResampledScreen(data_tmp)
  random_data_long = Comppass.convertMatrix(random_screen, cnames)
  random_stats_tab = Comppass.StatsTable(random_data_long)
  random_repro_tab = Comppass.ReproTable(random_data_long)
  random_speci_tab = Comppass.SpeciTable(random_stats_tab)
  random_Z = Comppass.Z(random_stats_tab)
  random_S = Comppass.S(random_stats_tab, random_speci_tab)
  random_D = Comppass.D(random_stats_tab, random_speci_tab, random_repro_tab)
  random_WD = Comppass.WD(random_stats_tab, random_repro_tab, random_speci_tab)
  random_Z_ecdf = ecdf(random_Z)
  random_S_ecdf = ecdf(random_S)
  random_D_ecdf = ecdf(random_D)
  random_WD_ecdf = ecdf(random_WD)
  P = c()
  P = cbind(P, sapply(summary[,"Z"], function(x) 1-random_Z_ecdf(x)))
  P = cbind(P, sapply(summary[,"S"], function(x) 1-random_S_ecdf(x)))
  P = cbind(P, sapply(summary[,"D"], function(x) 1-random_D_ecdf(x)))
  P = cbind(P, sapply(summary[,"WD"], function(x) 1-random_Z_ecdf(x)))
  colnames(P) = c("pZ","pS","pD","pWD")
  P
}

Comppass.main = function(data_file, output_file, resampling=F){
  cat("\tREADING\n")
  data = read.delim(data_file,  stringsAsFactors=F, header=T,skip=1, check.names=FALSE)
  cat("\tCONVERTING\n")
  ## convert into intermediate format
  data_tmp = Comppass.cleanMatrix(data)
  cnames = Comppass.cnames(data)
  data_long = Comppass.convertMatrix(data_tmp, cnames)
  
  cat("\tCOMPUTING STATS\n")
  ## make stats table
  stats_tab = Comppass.StatsTable(data_long)
  
  ## make reproducibility table
  repro_tab = Comppass.ReproTable(data_long)
  
  ## make specificity table
  speci_tab = Comppass.SpeciTable(stats_tab)
  
  cat("\tCOMPUTING SCORES\n")
  ## compute scores
  Z = Comppass.Z(stats_tab)
  S = Comppass.S(stats_tab, speci_tab)
  D = Comppass.D(stats_tab, speci_tab, repro_tab)
  WD = Comppass.WD(stats_tab, speci_tab, repro_tab)
  #get number of replicates
  num_reps = as.matrix(table(as.character(names(data)[-c(1:4)])))
  WDv2 = Comppass.WDv2(stats_tab, speci_tab, repro_tab, num_reps)
  
  cat("\tSUMMARIZING\n")
  ## compile all scores
  summary = Comppass.Summary(stats_tab, Z, S, D, WD, WDv2)
  
  if(resampling){
    print("RESAMPLING SCREEN")
    ## resample the screen
    P = Comppass.ResampledPvalues(data_tmp, cnames, summary)
    
    summary = cbind(summary, P)
  }
  cat("\tWRITING\n")
  ## write out
  write.table(summary, file=output_file, row.names=F, col.names=T, eol="\n", sep="\t", quote=F) 
  return(summary)
}




