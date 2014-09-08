#! /usr/bin/Rscript
suppressMessages(library(optparse))
#source('src/training.R')

# simplify/convert data into workable form
mist.processMatrix <- function(x){
	info = data.frame(Ip=colnames(x)[5:dim(x)[2]], Bait=as.vector(t(x[1,5:dim(x)[2]])),Specificity_exclusion=as.vector(t(x[2,5:dim(x)[2]])), stringsAsFactors=F)
	x <- x[-c(1:2),]
	names(x)[1:3] <- c("Preys", "PepAtlas", "Length")
	row.names(x) <- x$Preys
	idx <- c(2,3,5:dim(x)[2])
	#convert columns from characters to numbers
	for(i in idx){
		x[,i] <- as.numeric(x[,i])
	}
	return(list(x, info))
}

# "normalize" the M3D variable
mist.getM3D_normalized <- function(x){
	# divide cols by their bait length
	x1 <- x[,5:dim(x)[2]]/x[,3]
	# normalize by column (of x)	--	M3D
	x1 <- scale(x1, center=FALSE, scale=colSums(x[,5:dim(x)[2]]))
	# normalize by column (of x1)	--	M3D normalized
	x1 <- scale(x1, center=FALSE, scale=colSums(x1))	
	return(x1)
}

# calculate Abundance, Reproducibility, and Specificity
mist.getMetrics <- function(x, info){
	reproducibility <- c()
	abundance <- c()
	
	# look at data per bait
	for( i in unique(info$Bait)){
		bidx <- which(info$Bait == i)
		y <- x[,bidx]
		# normalize by row by bait
		
    if(length(bidx)>1){
      y <- y/apply(y, 1, sum)
    }
		y[is.na(y)] <- 0
		
		# Reproducibility ("entropies")
		# -----------------------------
		y[y!=1 & y>0] = y[y!=1 & y>0] * log2(y[y!=1 & y>0])
		y[y==1] = (1-1e-10)*log2(1-1e-10)
		
		if(length(bidx)<2){ #if there's only one row for that bait
      y <- sum(y)
		}else{
      y <- rowSums(y)
		}
    
    if(length(bidx)!=1){
			y = y/log2(1/length(bidx))
		}else{
			y = y*0  # set reproducibility to 0 if there is only one replicate
			cat(sprintf("\t!!ONLY ONE REPLICATE USING '%s' AS A BAIT PROTEIN. ASSIGNING REPRODUCIBILITY SCORE OF 0.\n",i))
		}
		reproducibility <- cbind(reproducibility, y)
		colnames(reproducibility)[dim(reproducibility)[2]] <- i
		
		# Abundance ("averages")
		# ----------------------
    if(length(bidx)>1){
  		abundance <- cbind(abundance, rowSums(x[,bidx])/length(bidx))
    }else{ #if there's only one row for that bait
      abundance <- cbind(abundance, x[,bidx]/length(bidx))
    }
		colnames(abundance)[dim(abundance)[2]] <- i
	}
	
	# Specificity
	# ----------------------
  # Must account for specificity exclusions
	specificity <- mist.getSpecificity(abundance, info)
	#specificity <- t(apply(abundance,1,function(z) z/sum(z)))
	return(list(reproducibility, abundance, specificity))
}

# calculate specificity taking the exclusions into account
mist.getSpecificity <- function(abundance, info){
  specificity <- abundance
  
  for( i in 1:dim(abundance)[2]){
    bait <- colnames(abundance)[i]
    # find baits to exclude in exclusions list
    baits_to_exclude = unique(info[info$Bait==bait, 'Specificity_exclusion'])
    baits_to_exclude = setdiff(baits_to_exclude,bait)
    baits_to_include = setdiff(colnames(abundance),baits_to_exclude)
    
    if(length(baits_to_include)==1){
      specificity[,i] <- 1
    }else if(length(baits_to_exclude)>0 & sum(rowSums(abundance[,baits_to_include]))>0){   # there are exclusions AND all the rowsums of the data != 0
      specificity[,i] <- specificity[,i]/rowSums(abundance[,baits_to_include])
    }else if(length(baits_to_exclude)==0 & sum(rowSums(abundance))>0){    # no exclusions AND all the rowsums of the data != 0
      specificity[,i] <- specificity[,i]/rowSums(abundance)
    }else{    # rowsums of the data = 0
      specificity[,i] <- 0
    }
    ## remove NA's from dividing by 0
    specificity[is.na(specificity[,i]),i] <- 0
    ## if specificity == Inf it means all other rows were 0, so set to 1 (max.value)  
  }
  return(specificity)
}

# vectorize the metrics while keeping the names straight
mist.vectorize <- function(x){
  temp = melt(x, factorsAsStrings=TRUE)
	names(temp) <- c("Prey", "Bait", "Xscore")
	temp$Xscore <- as.numeric(temp$Xscore)
  temp$Bait <- as.character(temp$Bait)
  temp$Prey <- as.character(temp$Prey)
	return(temp[,c('Bait','Prey','Xscore')])
}

# Perform PCA analysis AS DONE IN MIST
mist.doPCA <- function(R,A,S){
	# vectorize the variables
	m <- cbind(R=R$Xscore, A=A$Xscore, S=S$Xscore )
	x <- princomp(m)	# <- shouldn't we mean scale per variable before this step?
	
	#now do some other stuff??
	scores <- -x$scores[,1]
	scores <- (scores - min(scores))/(max(scores)-min(scores))
	scores <- cbind(R=R, A=A$Xscore, S=S$Xscore, MiST=scores)
	names(scores) = c("Bait", "Prey", "Reproducibility", "Abundance", "Specificity", "MIST_self")
	return(scores)
}

mist.getSampleOccurences = function(m3d_norm, info){
  data_melted = melt(m3d_norm, value.name='abundance', varnames=c('Prey','Ip'))
  data_merged = merge(info, data_melted[data_melted$abundance>0,], by='Ip')
  ip_occurences = aggregate(data_merged[,c('Bait','Prey','Ip')], by=as.list(data_merged[,c('Bait','Prey')]), FUN=function(x)paste(unique(x),collapse=','))
  ip_occurences[,c('Bait','Prey','Ip')]
}

#scores <- cbind(scores, mist_hiv=scores$Repro*0.30853 + scores$Abundance*0.00596 + scores$Specificity*0.68551 )
##############################################################################################################

mist.main <- function(matrix_file, weights='fixed', w_R=0.30853, w_A=0.00596, w_S=0.68551, training_file, training_steps=0.1){
  dat <- read.delim(matrix_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  dat <- mist.processMatrix(dat)
  m3d_norm <- mist.getM3D_normalized(dat[[1]])
  info = dat[[2]]

  dat <- mist.getMetrics(m3d_norm, info)
  R <- mist.vectorize(dat[[1]])
  A <- mist.vectorize(dat[[2]])
  S <- mist.vectorize(dat[[3]])
  
  metrics = data.frame(Bait=A$Bait,Prey=A$Prey,Abundance=A$Xscore,Reproducibility=R$Xscore,Specificity=S$Xscore)
  ## only retain non-zero results
  metrics = metrics[metrics$Abundance>0,]
  ## for debug purposes
#   output_file <- gsub('.txt', "_MIST_METRICS.txt", matrix_file)
#   write.table(metrics, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
  
  if(weights == 'fixed'){
    mist_scores = metrics$Reproducibility*w_R + metrics$Abundance*w_A + metrics$Specificity*w_S
    results = data.frame(metrics, MIST=mist_scores)  
  }else if(weights == 'PCA'){
    results <- mist.doPCA(R,A,S)
  }else if(weights == 'training'){
    training_set = read.delim(training_file, header=F, stringsAsFactors=F)
    colnames(training_set) = c('Bait','Prey')
    output_file <- paste(dirname(matrix_file), "prediction_rates.txt", sep="/")
    training_weights = mist.train.main(metrics, training_set, output_file, training_steps)
    cat(sprintf("\tWEIGHTS BASED ON TRAINING SET:\n\t  REPRODUCIBILITY: %s\n\t  ABUNDANCE: %s\n\t  SPECIFICITY: %s\n",training_weights$R, training_weights$A, training_weights$S))
    mist_scores = metrics$Reproducibility*training_weights$R + metrics$Abundance*training_weights$A + metrics$Specificity*training_weights$S
    results = data.frame(metrics, MIST=mist_scores)  
  }else{
    print(sprintf('unrecognized MIST option: %s',weights))
  }
  
  ## per bait get all IPs the prey was found in
  ip_occurences = mist.getSampleOccurences(m3d_norm, info)
  
  if(nrow(ip_occurences) != nrow(results)){
    print(" ERROR : INCONSISTENCY BETWEEN SCORE MATRIX AND IP OCCURENCE MATRIX ")
  }
  
  results_with_samples = merge(results, ip_occurences, by=c('Bait','Prey'))
  results_with_samples
}

if(!exists("PIPELINE") || PIPELINE==F){
  option_list <- list(
    make_option(c("-d", "--data_file"),
                help="data file containing values"), 
    make_option(c("-b", "--exclusions_file"),
                help="exclusions file"),
    make_option(c("-n", "--MAT_name"),
                help="name of the MAT data file"),
    make_option(c("-s", "--self_score"),
                help="calculate MIST self score (0/1)") 
  )
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
}

## TODO: make the following code into a unit-test 
# config = yaml.load(string=paste(readLines("tests/APMS_TEST.yml"),collapse='\n'))
# config = yaml.load(string=paste(readLines("tests/entero/APMS_ENTERO.yml"),collapse='\n'))
# mist.main(matrix_file=config$mist$matrix_file, weights=config$mist$weights, w_R=config$mist$reproducibility, w_A=config$mist$abundance, w_S=config$mist$specificity, training_file=config$mist$training_file)