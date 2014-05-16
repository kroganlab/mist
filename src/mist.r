#! /usr/bin/Rscript
suppressMessages(library(optparse))

# simplify/convert data into workable form
processMatrix <- function(x){
	baits <- x[1,5:dim(x)[2]]
	x <- x[-c(1:2),]
	names(x)[1:3] <- c("Preys", "PepAtlas", "Length")
	row.names(x) <- x$Preys
	idx <- c(2,3,5:dim(x)[2])
	#convert columns from characters to numbers
	for(i in idx){
		x[,i] <- as.numeric(x[,i])
	}
	return(list(x, baits))
}

# "normalize" the M3D variable
getM3D_normalized <- function(x){
	# divide cols by their bait length
	x1 <- x[,5:dim(x)[2]]/x[,3]
	# normalize by column (of x)	--	M3D
	x1 <- scale(x1, center=FALSE, scale=colSums(x[,5:dim(x)[2]]))
	# normalize by column (of x1)	--	M3D normalized
	x1 <- scale(x1, center=FALSE, scale=colSums(x1))	
	return(x1)
}


# calculate Abundance, Reproducibility, and Specificity
getMetrics <- function(x, baits, exclusions){
	reproducibility <- c()
	abundance <- c()
	
	# look at data per bait
	for( i in unique(baits)){
		bidx <- which(baits == i)
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
			y = y*0+1
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
	specificity <- getSpecificity(abundance, exclusions)
	#specificity <- t(apply(abundance,1,function(z) z/sum(z)))
	return(list(reproducibility, abundance, specificity))
}

# calculate specificity taking the exclusions into account
getSpecificity <- function(abundance, exclusions){
  specificity <- abundance
  for( i in 1:dim(abundance)[2]){
    bait <- colnames(abundance)[i]
    # find bait in exclusions list
    idx <- which(exclusions[,1]==bait)
    exes <- unique(unlist(strsplit(exclusions[idx,2],"\\|")))
    # find where excluded baits are in matrix
    idx <- which(!is.na(match(colnames(abundance), exes)))
    
    if(length(idx)>0 & sum(rowSums(abundance[,-idx]))>0){   # there are exclusions AND all the rowsums of the data != 0
      specificity[,i] <- specificity[,i]/rowSums(abundance[,-idx])
    }else if(length(idx)==0 & sum(rowSums(abundance))>0){    # no exclusions AND all the rowsums of the data != 0
      specificity[,i] <- specificity[,i]/rowSums(abundance)
    }else{    # rowsums of the data = 0
      specificity[,i] <- 0
    }
    # remove NA's from dividing by 0
    specificity[is.na(specificity[,i]),i] <- 0
  }
  return(specificity)
}

# vectorize the metrics while keeping the names straight
vectorize <- function(x){
	temp <- c()
	for(i in 1:dim(x)[2] ){
		tmp <- as.data.frame(cbind(x[,i], colnames(x)[i], names(x[,i])), stringsAsFactors=FALSE)
		temp <- rbind(temp,tmp)
	}
	names(temp) <- c("Xscore", "Bait", "Prey")
	temp$Xscore <- as.numeric(temp$Xscore)
	return(temp[,c(3,2,1)])
}


# Perform PCA analysis AS DONE IN MIST
doPCA <- function(R,A,S){
	# vectorize the variables
	m <- cbind(R=R$Xscore, A=A$Xscore, S=S$Xscore )
	x <- princomp(m)	# <- shouldn't we mean scale per variable before this step?
	
	#now do some other stuff??
	scores <- -x$scores[,1]
	scores <- 1 - (scores - min(scores))/(max(scores)-min(scores))
	scores <- cbind(R=R, A=A$Xscore, S=S$Xscore, MiST=scores)
	names(scores) = c("Prey", "Bait", "Reproducibility", "Abundance", "Specificity", "MIST_self")
	return(scores)
}

#scores <- cbind(scores, mist_hiv=scores$Repro*0.30853 + scores$Abundance*0.00596 + scores$Specificity*0.68551 )
##############################################################################################################

main <- function(data_file, exclusions_file, MAT_name, self_score){
  dat <- read.delim(data_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  exclusions <- read.delim(exclusions_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  
  dat <- processMatrix(dat)
  m3d_norm <- getM3D_normalized(dat[[1]])
  
  dat <- getMetrics(m3d_norm, unlist(dat[[2]]), exclusions)
  R <- vectorize(dat[[1]])
  A <- vectorize(dat[[2]])
  S <- vectorize(dat[[3]])
  
  x <- doPCA(R,A,S) #This is more or less ignored since we always use hiv weights
  x <- cbind(x, MIST_hiv=x$Repro*0.30853 + x$Abundance*0.00596 + x$Specificity*0.68551 )
  x <- x[,c('Bait','Prey','Reproducibility','Abundance','Specificity','MIST_self','MIST_hiv')]
  
  #write out MIST_hiv score/metrics
  out_file <- paste(MAT_name, "_MIST_HIV_metrics.txt", sep="")
  write.table(x[,c('Bait','Prey','Reproducibility','Abundance','Specificity')], out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
  out_file <- paste(MAT_name, "_MIST_HIV_scores.txt", sep="")
  write.table(x[,c('Bait','Prey','MIST_hiv')], out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
  
  #write out MIST_self score/metrics
  if(self_score == 1){
    out_file <- paste(MAT_name, "_MIST_SELF_metrics.txt", sep="")
    write.table(x[,c('Bait','Prey','Reproducibility','Abundance','Specificity')], out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
    out_file <- paste(MAT_name, "_MIST_SELF_scores.txt", sep="")
    write.table(x[,c('Bait','Prey','MIST_self')], out_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
  }
}



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
main(parsedArgs$data_file, parsedArgs$exclusions_file, parsedArgs$MAT_name, parsedArgs$self_score)



