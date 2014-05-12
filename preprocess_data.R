#! /usr/bin/Rscript
suppressMessages(library(Biostrings))
suppressMessages(library(optparse))
suppressMessages(library(reshape2))

# filter out contaminants
filterContaminants = function(contaminant_fasta, df) {
	## use biostrings function to read in fasta file
	fastaF = as.data.frame(readAAStringSet(filepath = contaminant_fasta))
	cont_df = data.frame(name = rownames(fastaF), seq = fastaF$x)
	## split out the first part of the fasta header which should be the identifier -- split could be eithe space or pipe
	fastaFFirstNames = as.data.frame(apply(cont_df, 1, FUN = function(x) as.character(strsplit(x[1], split = "\\s|\\|")[[1]][1])))
	## check whether the key is in uniprot format
	UniprotNames = as.data.frame(fastaFFirstNames[grep("^[A-Z]{1}[0-9]{1}[A-Z,0-9]{4}(\\-1)*", fastaFFirstNames[, 1]), ])
	colnames(UniprotNames)[1] = "contaminant_key"
	df_F = df[!(df$ms_uniprot_ac %in% UniprotNames$contaminant_key),]
	print(sprintf("FILTERED %s CONTAMINANTS", nrow(df) - nrow(df_F)))
	df_F
}


# merge data with keys
mergeData <- function(dat, keys){
  ids = union(unique(dat[,1]), unique(keys[,1]) )
  print(paste("  ", length(ids), " TOTAL BAITS DETECTED", sep=""))
  print(paste("  ", length(unique(dat[,1])), "/", length(ids), " BAITS DETECTED IN DATA FILE", sep=""))
  if(length(unique(dat[,1]))<length(ids)){
      print("    MISSING BAITS: ")
      print(setdiff(ids, unique(dat[,1])))
  }
  print(paste("  ", length(unique(dat[,1])), "/", length(ids), " BAITS DETECTED IN KEYS FILE", sep=""))
  if(length(unique(dat[,1]))<length(ids)){
    print("    MISSING BAITS: ")
    print(setdiff(ids, unique(keys[,1])))
  }
  
	x <- merge(dat, keys, by="id")
	x <- x[,c(dim(x)[2],1:(dim(x)[2]-1))]
	return(x)
}


# remove duplicate prey entries per IP. Rare, but happens
removeDuplicates = function(y){
  #find duplicate prey entries per IP
  idx <- duplicated(y[,c('id','ms_uniprot_ac')])
  if(sum(idx)>0){
    print("  !! DUPLICATE ID/PREY PAIRS FOUND IN DATA! DUPLICATES REMOVED.")
    y <- y[!idx,]
  }
  return(y)
}


# order the data based on it's id
orderExperiments <- function(y){
	idnum <- sub("(^[A-Za-z]+|^[A-Za-z]+-)","",y$id)	#strip out numbers
	idname <- sub("([0-9].*|-[0-9].*)","",y$id)			#strip out id characters
	tmp = strsplit(idnum,"-")
	tmp = matrix( unlist(lapply(tmp, function(x) if(length(x)>1){c(x[1], x[2])}else{c(x[1],0)} )), ncol=2, byrow=T)
	tmp = cbind(idname, tmp)
  tmp[,2] = sub("[^0-9]","",tmp[,2])  #added to remove extra characters/puctuation here after number like "rerun"
	idx = order(tmp[,1], as.numeric(tmp[,2]), as.numeric(tmp[,3]))
	return(idx)
}

# Find potential carryover and print to a file. To be used/appended to final score consolidation sheet
findCarryover <- function(x){
	# order experiments
	x <- x[orderExperiments(x),]
	experiments <- unique(x$id)
	baits <- unique(x$BAIT)
	tab = table(x$ms_uniprot_ac)	#for prey
	preys = as.numeric(tab)
	names(preys) = names(tab)
	
	idx <- which(x$ms_num_unique_pep > 10)	# 0) find which ms_uniq_pep > 10
	to_remove <- NULL
	for( i in idx){
		# get index of next 4 experiments run after this one
		idxj <- (which(experiments == x[i,]$id)+1):(which(experiments == x[i,]$id)+4)
		
		tmp <- x[which(!is.na(match(x$id, experiments[idxj]))),]		# pull these 4 samples
		tmp <- tmp[which(tmp$ms_uniprot_ac == x$ms_uniprot_ac[i]),]

		# check all criteria from carryover v1
		prey = preys[names(preys) == x$ms_uniprot_ac[i]]	#the number of occurrences of this prey in the data set
		be_gone <- which(tmp$ms_num_unique_pep > 0 & tmp$ms_num_unique_peptide < (x$ms_num_unique_peptide[i]/2) & prey<length(experiments)/3)
		#print(be_gone)
		if(length(be_gone)){
			#print(tmp[be_gone,1:6])
			to_remove <- rbind(to_remove, tmp[be_gone,c(1:6)] )
		}
	}
	to_remove = unique(to_remove)	#remove duplicates
	to_remove <- to_remove[order(to_remove$id,decreasing=FALSE),]	#re-order
	
}


# create the data matrix that will be used as input by MiST and Saint algorithms
createMatrix <- function(y, nupsc_flag, collapse_file, exclusions_file, remove_file){

  # collapse bait names from "collapse" file
  if(file.info(collapse_file)$size >0){
    collapse <- read.delim(collapse_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    for( i in 1:dim(collapse)[1] ){
      y$BAIT[y$BAIT == collapse[i,1]] = collapse[i,2]
    }
  }else{
    print("Collapse file is empty.")
  }
  
  # remove the "remove" ip's
  if(file.info(remove_file)$size>0){
    removals <- read.delim(remove_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    y.len = dim(y)[1]
    y <- y[!y$id %in% removals[,1],]
    if(dim(removals)[1]>0 & y.len == dim(y)[1])
      print(" !! WARNING: REMOVE > 0 BUT NO ENTRIES REMOVED!")
  }else{
    print("Remove file is empty.")
  }
  
  # Create matrix using either "number of unique peptides" or "spectral count"
  if(nupsc_flag == 0 | is.null(y$ms_spectral_counts) ){
    if(length(y$ms_num_unique_peptide)==0){
      print("  !!ERROR: ms_num_unique_peptide column needs to be named 'ms_num_unique_peptide'")
    }
    datmat <- dcast(y, ms_uniprot_ac ~ id + BAIT, value.var = "ms_num_unique_peptide", sum)
    if(is.null(y$ms_spectral_counts)){
      print(" !!! SPECTRAL COUNTS NOT FOUND! Replace or label as 'ms_spectral_counts'.")
      print("     -> USING 'ms_num_unique_peptide'.")
    }
  }else{
    datmat <- dcast(y, ms_uniprot_ac ~ id + BAIT, value.var = "ms_spectral_counts", sum)
  }

  # get "Lengths" (molecular weights)
  preys <- unique(y[,c('ms_uniprot_ac','ms_protein_mw')])
  preys$ms_protein_mw <- floor(preys$ms_protein_mw/110)
  datmat <- merge(preys, datmat, by="ms_uniprot_ac", all.y=T)
  # add other columns for saint. (currently not used)
  datmat <- cbind(PepAtlas=0, 'PreyType/BaitCov'='N', datmat)
  datmat <- datmat[,c(3,1,4,2,5:dim(datmat)[2])]
  
  # handle exclusions
  if(file.info(exclusions_file)$size>0){
    exclusions <- unique(read.delim(exclusions_file, sep="\t", header=FALSE, stringsAsFactors=FALSE))
    ips <-unique(y[,c('id','BAIT')])
    ips <- merge(ips, exclusions, by.x="BAIT", by.y="V1", all.x=TRUE)[, c(2,1,3)]
    ips <- ips[order(ips$id, ips$BAIT),]
    names(ips)[3] = "PreyType/BaitCov"
    idx <- which(is.na(ips[,3]))
    ips[idx,3] = ips[idx,2]
  }else{
    print("Exclusions file is empty.")
    ips <-unique(y[,c('id','BAIT')])
    ips$"PreyType/BaitCov" = ips$BAIT
  }
  
  # Create headers for matrix
  corner_titles = cbind(c("a","b","c","IP"), c("a","b","c","Bait"), c("Preys","PepAtlas","Length","PreyType/BaitCov"))
  ips <- t(rbind(corner_titles, as.matrix(ips)) )
  #colnames(ips) = colnames(datmat) = NULL
  #datmat <- as.data.frame(rbind(ips, as.matrix(datmat)), straingsAsFactors=FALSE)
  
  return(list(ips,datmat))
  
  #############################################################################################################################################
}



# wrapper to filter data and merge with keys
main <- function(data_file, keys_file, output_file, filter_data, rm_co, nupsc_flag, collapse_file, exclusions_file, remove_file){
  print("Reading Files")
  out_file <- unlist(strsplit(output_file,"\\."))[1]		#get ouput_dir from output_file
	#out_dir <- paste(out_dir[-length(out_dir)],collapse="/")
	df <- read.delim(data_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	keys <- read.delim(keys_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	names(keys) = c("id", "BAIT")
	
  # quality control
  print("Removing decoys and prey with 0 unique peptides")
  df <- df[-grep("decoy",df[,4]),]               # remove "decoys"
	df <- df[which(df[,3] > 0 | is.na(df[,3])),]   # remove ms_unique_pep <= 0
  df <- removeDuplicates(df)
  
	#filter contaminants out
  print("FILTERING COMMON CONTAMINANTS")
	if(filter_data == 1)
		df <- filterContaminants("~/HPC/MSPipeline/files/contaminants.fasta",df)
	
	#merge keys with data
	print("MERGING KEYS WITH DATA")
	df <- mergeData(df, keys)
	df <- df[orderExperiments(df),]  #GENERATES WARNINGS WHEN ID# HAS CHARACTERS IN IT: FIXED
  outfile = paste(out_file, ".txt", sep="")
	write.table(df, outfile, eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")

  # create matrix to be used by MiST and Saint
	print("CONVERTING TO MATRIX")
  outfile = paste(out_file, "_MAT.txt", sep="")
  df <- createMatrix(df, nupsc_flag, collapse_file, exclusions_file, remove_file) #return a list b/c of space padding
	write.table(df[[1]], outfile, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="")
  write.table(df[[2]], outfile, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="", append=TRUE)
  
  
	# Remove Carryover
	#if(rm_co==1){
	#	to_remove <- findCarryover(df)
	#	write.table(to_remove, paste(out_file_temp,"_ToRemove.txt",sep=""), eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
	#	write.table(df, paste(out_file_temp,"_NoC.txt",sep=""), eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
	#}
  	
}


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
              help="use number of unique peptides or spectral count as matrix values")  
)
parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

#data_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/input/kinases_data.txt"
#keys_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/input/kinases_keys.txt"
#collapse_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/input/kinases_collapse.txt"
#exclusions_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/input/kinases_exclusions.txt"
#remove_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/input/kinases_remove.txt"
#output_file = "~/Desktop/test_matrix.txt"
#nupsc_flag = 1
#remove_carryover=0
#filter_data=1
#main(data_file, keys_file, output_file, filter_data, remove_carryover, nupsc_flag, collapse_file, exclusions_file, remove_file)



main(parsedArgs$data_file, parsedArgs$keys_file, parsedArgs$output_file, parsedArgs$filter_data, parsedArgs$remove_carryover, parsedArgs$nupsc_flag, parsedArgs$collapse_file, parsedArgs$exclusions_file, parsedArgs$remove_file)
