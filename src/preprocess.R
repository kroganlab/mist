suppressMessages(library(optparse))
suppressMessages(library(reshape2))

# filter out contaminants
preprocess.filterContaminants = function(contaminants_file, df, prey_colname) {
	# read in contaminants file
  contaminants <- read.delim(contaminants_file, header=F, stringsAsFactors=F, sep='\t')
  contaminants <- paste(unlist(contaminants), collapse="|")
  # find and remove contaminants
  idx <- grep(contaminants, df[,prey_colname])
  df <- df[-idx,]
  cat(sprintf("\t  FILTERED %s CONTAMINANTS\n", length(idx)))
  return(df)
}

# merge data with keys
preprocess.mergeData <- function(dat, keys, id_colname){
  ids = union(unique(dat[,id_colname]), unique(keys[,id_colname]) )
  cat(paste("\t  ", length(ids), " TOTAL BAITS DETECTED\n", sep=""))
  cat(paste("\t  ", length(unique(dat[,id_colname])), "/", length(ids), " BAITS DETECTED IN DATA FILE\n", sep=""))
  if(length(unique(dat[,id_colname]))<length(ids)){
      cat(sprintf("\t  MISSING BAITS: %s \n",setdiff(ids, unique(dat[,1]))))
  }
  cat(paste("\t  ", length(unique(keys[,id_colname])), "/", length(ids), " BAITS DETECTED IN KEYS FILE\n", sep=""))
  if(length(unique(keys[,id_colname]))<length(ids)){
    cat(sprintf("\t  MISSING BAITS: %s\n",setdiff(ids, unique(keys[,id_colname]))))
  }
  
	x <- merge(dat, keys, by=id_colname)
	x <- x[,c(dim(x)[2],1:(dim(x)[2]-1))]
	return(x)
}

# remove duplicate prey entries per IP. Rare, but happens
preprocess.removeDuplicates = function(y, id_colname, prey_colname){
  #find duplicate prey entries per IP
  idx <- duplicated(y[,c(id_colname, prey_colname)])
  if(sum(idx)>0){
    cat("\t!! DUPLICATE ID/PREY PAIRS FOUND IN DATA! DUPLICATES REMOVED\n")
    y <- y[!idx,]
  }
  return(y)
}

# order the data based on it's id
preprocess.orderExperiments <- function(y, id_colname){
	idnum <- sub("(^[A-Za-z_]+|^[A-Za-z_]+-)","",y[,id_colname])	#strip out numbers
	idname <- sub("([0-9].*|-[0-9].*)","",y[,id_colname])			#strip out id characters
	tmp = strsplit(idnum,"-")
	tmp = matrix( unlist(lapply(tmp, function(x) if(length(x)>1){c(x[1], x[2])}else{c(x[1],0)} )), ncol=2, byrow=T)
	tmp = cbind(idname, tmp)
  tmp[,2] = sub("[^0-9]","",tmp[,2])  #added to remove extra characters/puctuation here after number like "rerun"
	idx = order(tmp[,1], as.numeric(tmp[,2]), as.numeric(tmp[,3]))
	return(idx)
}

# Find potential carryover and print to a file. To be used/appended to final score consolidation sheet
preprocess.findCarryover <- function(x, id_colname, prey_colname, pepcount_colname){
	# order experiments
	x <- x[preprocess.orderExperiments(x, id_colname),]
	experiments <- unique(x[,id_colname])
	baits <- unique(x$BAIT)
	tab = table(x[,prey_colname])	#for prey
	preys = as.numeric(tab)
	names(preys) = names(tab)
	
	idx <- which(x[,pepcount_colname] > 10)	# 0) find which ms_uniq_pep > 10
	to_remove <- NULL
	for( i in idx){
		# get index of next 4 experiments run after this one
		idxj <- (which(experiments == x[i,id_colname])+1):(which(experiments == x[i,id_colname])+4)
		tmp <- x[which(!is.na(match(x[,id_colname], experiments[idxj]))),]		# pull these 4 samples
		tmp <- tmp[which(tmp[,prey_colname] == x[i,prey_colname]),]
    
		# check all criteria from carryover v1
		prey = preys[names(preys) == x[i,prey_colname]]	#the number of occurrences of this prey in the data set
		be_gone <- which(tmp[,pepcount_colname] > 0 & tmp[,pepcount_colname] < (x[i, pepcount_colname]/2) & prey<length(experiments)/3)
		if(length(be_gone)){
			to_remove <- rbind(to_remove, tmp[be_gone,c(1:6)] )
		}
	}
  if(length(to_remove)>0){
	  to_remove = unique(to_remove)	#remove duplicates
	  to_remove <- to_remove[order(to_remove[,id_colname],decreasing=FALSE),]	#re-order
  }
  return(to_remove)
}

# create the data matrix that will be used as input by MiST and Saint algorithms
preprocess.createMatrix <- function(y, collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname){
  # Creating consistent names
  names(y)[grep(paste("^",id_colname,"$",sep=""), names(y))] = "id_colname"
  names(y)[grep(paste("^",prey_colname,"$",sep=""), names(y))] = "prey_colname"
  names(y)[grep(paste("^",pepcount_colname,"$",sep=""), names(y))] = "pepcount_colname"
  names(y)[grep(paste("^",mw_colname,"$",sep=""), names(y))] = "mw_colname"
  
  # collapse bait names from "collapse" file
  if( (file.info(collapse_file)$size >0) & (file.exists(collapse_file)) & (!file.info(collapse_file)$isdir) ){
    collapse <- read.delim(collapse_file, sep="\t", header=F, stringsAsFactors=FALSE)
    for( i in 1:dim(collapse)[1] ){
      y$BAIT[y$BAIT == collapse[i,1]] = collapse[i,2]
    }
  }else if( !(file.exists(collapse_file)) | (file.info(collapse_file)$isdir)){
    cat("\tCOLLAPSE FILE DOES NOT EXIST. NOT USING COLLAPSE FEATURE.\n")
  }else{
    cat("\tCOLLAPSE FILE IS EMPTY\n")
  }
  
  # remove the "remove" ip's
  if( (file.info(remove_file)$size >0) & (file.exists(remove_file)) & (!file.info(remove_file)$isdir) ){
    removals <- read.delim(remove_file, sep="\t", header=F, stringsAsFactors=FALSE)
    y.len = dim(y)[1]
    y <- y[!y$id %in% removals[,1],]
    if(dim(removals)[1]>0 & y.len == dim(y)[1])
      cat("\tWARNING: REMOVE > 0 BUT NO ENTRIES REMOVED\n")
  }else if( !(file.exists(remove_file)) | (file.info(remove_file)$isdir) ){
    cat("\tREMOVE FILE DOES NOT EXIST. NOT USING REMOVE FEATURE.\n")
  }else{
    cat("\tRemove file is empty\n")
  }
  
  # Create matrix using either "number of unique peptides" or "spectral count"
  if(!'pepcount_colname' %in% colnames(y)){
    cat(sprintf("\tPEPTIDE COUNT COLUMN %s NOT FOUND\n\t\tPLEASE CHECK DATA FILE\n",pepcount_colname))
    quit()
  }else if(!'prey_colname' %in% colnames(y)){
    cat(sprintf("\tPREY IDENTIFIER COLUMN %s NOT FOUND\n\t\tPLEASE CHECK DATA FILE\n",prey_colname))
    quit()
  }else{
    datmat <- dcast(y, prey_colname ~ id_colname + BAIT, value.var = c('pepcount_colname'), sum)
  }
  
  # get "Lengths" (molecular weights)
  preys <- unique(y[,c('prey_colname', 'mw_colname')])
  # handle multiple molecular weights for a prey protein
  if(length(preys$prey_colname) > length(unique(preys$prey_colname))){
    dup_prey = preys$prey_colname[which(duplicated(preys$prey_colname))]
    cat(sprintf("\tDIFFERENT MOLECULAR WEIGHTS DETECTED FOR THE FOLLOWING PREY: \n"))
    cat(sprintf("\t\t%s\n",dup_prey))
    cat(sprintf("\tUSING MEDIAN WEIGHT PER PROTEIN\n"))
    preys = aggregate(mw_colname~prey_colname, data=preys, median)
  }

  preys$mw_colname <- floor(as.numeric(preys$mw_colname)/110)
  datmat <- merge(preys, datmat, by.x='prey_colname', by.y='prey_colname', all.y=T)
  # add other columns for saint. (currently not used)
  datmat <- cbind(PepAtlas=1, 'PreyType/BaitCov'='N', datmat)
  datmat <- datmat[,c(3,1,4,2,5:dim(datmat)[2])]
  
  # handle exclusions
  if(file.info(exclusions_file)$size>0){
    exclusions <- unique(read.delim(exclusions_file, sep="\t", header=F, stringsAsFactors=FALSE))
    
    # if multiple instances of bait in col1, combine all of the col2 exclusions
    if( any(duplicated(exclusions[,1])) ){
      idx<-which(duplicated(exclusions[,1]))
      for(i in unique(exclusions[idx,1])){
        idx2 = which(exclusions[,1]==i)
        exclusions[idx2,2] = paste(exclusions[idx2,2], collapse="|")
      }
      exclusions = unique(exclusions)
    }
    
    ips <-unique(y[,c('id_colname','BAIT')])
    ips <- merge(ips, exclusions, by.x="BAIT", by.y="V1", all.x=TRUE)[, c(2,1,3)]
    ips <- ips[order(ips$id_colname, ips$BAIT),]
    names(ips)[3] = "PreyType/BaitCov"
    idx <- which(is.na(ips[,3]))
    ips[idx,3] = ips[idx,2]
  }else{
    cat("\tExclusions file is empty\n")
    ips <-unique(y[,c('id_colname','BAIT')])
    ips <- ips[order(ips$id_colname),]
    ips$"PreyType/BaitCov" = ips$BAIT
  }
  
  # Create headers for matrix
  corner_titles = cbind(c("a","b","c","IP"), c("a","b","c","Bait"), c("Preys","PepAtlas","Length","PreyType/BaitCov"))
  ips <- t(rbind(corner_titles, as.matrix(ips)) )
  
  return(list(ips,datmat))
  
}

preprocess.checkNames <- function(x, id_colname, prey_colname, pepcount_colname, mw_colname){
  config_names = c(id_colname, prey_colname, pepcount_colname, mw_colname) 
  if(!all(config_names %in% names(x))){
    cat(paste("\t'",config_names[!config_names %in% names(x)], "' NOT FOUND IN COLUMN NAMES. PLEASE CHECK DATA FILE.\n", sep=""))
    quit()
  }
  
}

# wrapper to filter data and merge with keys
preprocess.main <- function(data_file, keys_file, output_file, filter_data, contaminants_file, rm_co=T, collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname){
  cat("\tREADING FILES\n")
  keys=tryCatch(read.delim(keys_file, sep="\t", header=F, stringsAsFactors=FALSE), error = function(e) cat(sprintf('\tERROR reading keys from : %s\n\t  Please make sure full root path to file is used in yaml file.\n',keys_file)) )
  df=tryCatch(read.delim(data_file, sep="\t", header=T, stringsAsFactors=FALSE), error = function(e) cat(sprintf('\tERROR reading data from : %s\n\t  Please make sure full root path to file is used in yaml file.\n',data_file)))
	preprocess.checkNames(df, id_colname, prey_colname, pepcount_colname, mw_colname) 
	
  names(keys) = c(id_colname, "BAIT")
  keys$BAIT = gsub(' ', '', keys$BAIT)
  
  # quality control
  ## TO DO GIT ISSUE #1
  if(class(df[,3])=="character"){
    cat("\t!!! CHARACTERS FOUND IN unique_pep COLUMN. CONVERTING TO NUMERIC.\n")
    df[,3] = as.numeric(df[,3])
  }
	df <- df[which(df[,3] > 0 | is.na(df[,3]) | df[,3] == ""),]   # remove ms_unique_pep <= 0
  df <- preprocess.removeDuplicates(df, id_colname, prey_colname)
  
	#filter contaminants out
  cat("\tFILTERING COMMON CONTAMINANTS\n")
	if(filter_data == 1)
		df <- preprocess.filterContaminants(contaminants_file, df, prey_colname)
  
	#merge keys with data
	cat("\tMERGING KEYS WITH DATA\n")
	df <- preprocess.mergeData(df, keys, id_colname)
	df <- df[preprocess.orderExperiments(df, id_colname),]  #GENERATES WARNINGS WHEN ID# HAS CHARACTERS IN IT: FIXED
	write.table(df, output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
  
  # Remove Carryover
  if(rm_co==1){
    to_remove <- preprocess.findCarryover(df, id_colname, prey_colname, pepcount_colname)
    if(length(to_remove)>0){
      write.table(to_remove, gsub('.txt','_ToRemove.txt',output_file), eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
      # remove carryover proteins
       tmp = merge(df, data.frame(to_remove[,c(id_colname, prey_colname)], here=1), by=c(id_colname, prey_colname), all.x=TRUE) #get index of carryover proteins
      df <- tmp[which(is.na(tmp$here)),-dim(tmp)[2]]  # remove the 'here' column
      df <- df[preprocess.orderExperiments(df, id_colname),]
      output_file = gsub('.txt','_NoC.txt',output_file)
      write.table(df, output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
    }
  }
  
  # create matrix to be used by MiST
	cat("\tCONVERTING TO MATRIX\n")
  matrix_output_file = gsub('.txt','_MAT.txt',output_file)
  df_mat <- preprocess.createMatrix(df, collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname) #return a list b/c of space padding
	
  write.table(df_mat[[1]], matrix_output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="")
  write.table(df_mat[[2]], matrix_output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="", append=TRUE)
  
  return(matrix_output_file)
}

if(!exists("PIPELINE") || PIPELINE==F){
  option_list <- list(
    make_option(c("-d", "--data_file"),
                help="data file containing values"), 
    make_option(c("-k", "--keys_file"),
                help="keys file containing bait names"), 
    make_option(c("-o", "--output_file"),
                help="output file"), 
    make_option(c("-a", "--collapse_file"),
                help="collapse file"),
    make_option(c("-b", "--exclusions_file"),
                help="exclusions file"),
    make_option(c("-e", "--remove_file"),
                help="remove file"),
    make_option(c("-g", "--filter_data"),
                help="filter out common contaminants first"),
    make_option(c("-r", "--remove_carryover"),
                help="prints potential carryover to file"),
    make_option(c("-s", "--nupsc_flag"),
                help="the column name of the quantifiable value used to score data (num unique peptides or spectral count)")  )
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
}
