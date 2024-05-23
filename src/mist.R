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
mist.getMetrics <- function(x, info, standardize_specificity = NULL){
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
	# specificity <- t(apply(abundance,1,function(z) z/sum(z)))
	# Must account for specificity exclusions
	sc.list <- mist.getSpecificity(abundance, info, returnCounts = TRUE)
	specificity <- sc.list$specificity
	
	## transform specificity to FCvsBG using number of other baits
	fcVsBg <- mist.specificity2FC(specificity, sc.list$otherBaitCounts)
	
	## transform fc to a standard specificity. 18 other baits is count from original HIV study, and is default we calculate even if it won't be used
	if (is.null(standardize_specificity) || !is.integer(standardize_specificity) || standardize_specificity < 2){
	  standardize_specificity = 18
	}
	standardized_specificity <- mist.FC2Specificity(fcVsBg, otherBaitCount = standardize_specificity)	
	
  result <- list(Reproducibility = reproducibility, 
                 Abundance = abundance,
                 Specificity = specificity,
                 FcVsBg = fcVsBg)
  # special name for standardized specificity
  result[[sprintf("Specificity%d", standardize_specificity)]] <- standardized_specificity
	return(result)
}

# calculate specificity taking the exclusions into account
mist.getSpecificity <- function(abundance, info, returnCounts = FALSE){
  specificity <- abundance
  otherBaitCounts <- integer(ncol(abundance))
  names(otherBaitCounts) <- colnames(abundance)
  
  for( i in 1:dim(abundance)[2]){
    bait <- colnames(abundance)[i]
    # find baits to exclude in exclusions list
    baits_to_exclude = unique(info[info$Bait==bait, 'Specificity_exclusion'])
    baits_to_exclude = setdiff(baits_to_exclude,bait)
    baits_to_include = setdiff(colnames(abundance),baits_to_exclude)
    
    otherBaitCounts[i] <- length(baits_to_include)-1 # don't count bait_i
    
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
  }
  if (returnCounts){
    return (list(specificity = specificity, otherBaitCounts = otherBaitCounts))
  } else{
    return(specificity)
  }
}

# given matrix of specificity and number of other baits used in calculation, transform to a fold change vs average other
mist.specificity2FC <- function(specificity.mat, otherBaitCounts){
  # Calculation is (N-1) * S/(1-S), see notes
  # Here I use n instead of N-1 because its already been subtracted
  # otherBaitCounts = N-1, and n = N-1 below
  .fc <- function(s, n){
    n * s/(1-s)
  }

  # calculate per column
  sweep (specificity.mat, 2, otherBaitCounts, .fc)
}

# given matrix of fold changes vs background, and number of other baits desired, transform to a specificity
mist.FC2Specificity <- function (foldChange.mat, otherBaitCount = 18){
  # To explain...
  # Set average in other baits O_hat = 1, then abundance in B_i is same as FC. 'cause FC = B_i/O_hat
  # By definition of specificity, S = B/(B + sum_others) = B/(B + 1*otherBaitCount) =  FC/(FC + n)
  fc.mat <- foldChange.mat/(foldChange.mat + otherBaitCount)
  
  # Inf/(Inf + n) produces NA instead of 1
  fc.mat[foldChange.mat == Inf] <- 1
  
  return (fc.mat)
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

mist.main <- function(matrix_file, weights='fixed', w_R=0.30853, w_A=0.00596, w_S=0.68551, training_file, training_steps=0.1, standardize_specificity = NULL){
  dat <- read.delim(matrix_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  dat <- mist.processMatrix(dat)
  m3d_norm <- mist.getM3D_normalized(dat[[1]])
  info = dat[[2]]

  dat <- mist.getMetrics(m3d_norm, info, standardize_specificity = standardize_specificity)
  dat.vct <- lapply (dat,mist.vectorize )


  # old metrics...order may matter in some legacy code, so preserve
  R <- dat.vct$Reproducibility #mist.vectorize(dat[[1]])
  A <- dat.vct$Abundance #mist.vectorize(dat[[2]])
  S <- dat.vct$Specificity #mist.vectorize(dat[[3]])
  metrics = data.frame(Bait=A$Bait,Prey=A$Prey,Abundance=A$Xscore,
                       Reproducibility=R$Xscore,
                       Specificity=S$Xscore)
  
  # new metrics, order shouldn't matter, so use whatever order mist.getMetrics delivers
  # first delete the old ones
  for (name in c("Reproducibility", "Abundance", "Specificity"))
    dat.vct[[name]] <- NULL
  # make a table of every score available.
  newMetrics <-   data.frame(lapply(dat.vct, function(x)x$Xscore) )
  
  # combined old with new
  metrics <- cbind (metrics, newMetrics)
  
  ## only retain non-zero results
  metrics = metrics[metrics$Abundance>0,]
  ## for debug purposes
#   output_file <- gsub('.txt', "_MIST_METRICS.txt", matrix_file)
#   write.table(metrics, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
  
  
  # switch specificity based on standardize_specificity argument
  # clean up the argument
  if (is.null(standardize_specificity)){
    message ("New in V1.5: Specificity is now automatically scaled to dataset size by default. If you need consistency with previous versions, set standardize_specificity to FALSE in yml config file.")
    standardize_specificity <- 18L
  } else { 
    if (is.numeric(standardize_specificity)){
      if (standardize_specificity < 2){
        standardize_specificity <- as.logical(standardize_specificity)
      } else{
        standardize_specificity <- as.integer(standardize_specificity)
      }
    }
    if(TRUE == standardize_specificity)
      standardize_specificity <- 18L
    if (! (is.logical(standardize_specificity) | is.integer(standardize_specificity))){
      message("Unexpected value for standardize_specificity: ", standardize_specificity, " It should be TRUE/FALSE or an integer. Setting to FALSE")
      standardize_specificity <- FALSE
    }
  }
  
  if(standardize_specificity){
    .S <- metrics[[sprintf("Specificity%d", standardize_specificity)]]
    .metrics <- data.frame(Bait = metrics$Bait,
                           Prey = metrics$Prey,
                           Abundance = metrics$Abundance,
                           Reproducibility = metrics$Reproducibility,
                           Specificity = .S)    
  } else{
    .S <- metrics$Specificity
    .metrics <- metrics
  }
  
  
  if(weights == 'fixed'){
    mist_scores = metrics$Reproducibility*w_R + metrics$Abundance*w_A + .S*w_S
    results = data.frame(metrics, MIST=mist_scores)  
  }else if(weights == 'PCA'){
    results <- mist.doPCA(R,A,.S)
  }else if(weights == 'training'){
    training_set = read.delim(training_file, header=F, stringsAsFactors=F)
    colnames(training_set) = c('Bait','Prey')
    output_file <- paste(dirname(matrix_file), "prediction_rates.txt", sep="/")
    training_weights = mist.train.main(.metrics, training_set, output_file, training_steps)
    cat(sprintf("\tWEIGHTS BASED ON TRAINING SET:\n\t  REPRODUCIBILITY: %s\n\t  ABUNDANCE: %s\n\t  SPECIFICITY: %s\n",training_weights$R, training_weights$A, training_weights$S))
    mist_scores = metrics$Reproducibility*training_weights$R + metrics$Abundance*training_weights$A + .S*training_weights$S
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