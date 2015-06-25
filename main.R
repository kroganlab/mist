#! /usr/bin/Rscript --vanilla

cat(">> LOADING DEPENDENCIES...\n")
suppressMessages(library(getopt))
suppressMessages(library(yaml))

#########################
## MIST MAIN FILE #######

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'config'  , 'c', 1, "character", "configuration file in YAML format"),
  byrow=TRUE, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

PIPELINE=T

# set source directory
args <- commandArgs(trailingOnly = F) 
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))  # get the path to where main.R is located

## load all externeal files
source(paste(scriptPath,"/src/preprocess.R",sep=""))
source(paste(scriptPath,"/src/qc.R",sep=""))
source(paste(scriptPath,"/src/mist.R",sep=""))
source(paste(scriptPath,'/src/training.R',sep=""))

# check if libraries are installed
checkForLibraries <- function(required_libraries){
  pkgs = data.frame(installed.packages(), stringsAsFactors=F)
  #if the packages haven't been installed on the computer yet
  if(!all(required_libraries %in% pkgs$Package)){
    missing_packages = required_libraries[ !required_libraries %in% pkgs$Package ]
    stop( paste("The following packages need to be installed on your computer in order to run MIST:\n"), paste('\t',missing_packages,'\n', sep=""), call.=FALSE ) 
    
  }else{
    cat(">> All required R libraries accounted for.\n")
  }
}

# check for all packages reqiured
required_libraries = c('yaml','getopt','optparse','reshape2','compiler','stats','grDevices','graphics','pheatmap','RColorBrewer','ggplot2','gridExtra','MESS')
checkForLibraries(required_libraries)

# some qc to make sure config file exists and is well formatted
getConfig <- function(config_file){
  if( !file.exists( config_file )){
    stop( cat(paste("ERROR!!! The yml configuration file:\n\t",paste(getwd(),config_file,sep='/'), '\ndoes not exist. Please enter a correct path/filename.\n', sep='')), call.=F)
  }
  
  x = readLines(config_file)
  x = x[x!=""]  #remove \n\n cases (blank Lines)
  x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
  x = gsub('\t', '  ', x)
  config = paste(x, collapse='\n')
  config = tryCatch(yaml.load(config), error = function(e) { print("!!! Error loading the config file. Please make sure the file follows YAML format."); stop()} )
  
  # check formatting of config file
  config = formatConfig(config)
  return(config)
}

# chec out YAML config file and make sure it's setup properly
formatConfig = function(config){
  # replace any spaces in the colname with a "." to match how R reads in headers with spaces. 
  config$preprocess$ip_colname <- gsub(' ','.',config$preprocess$ip_colname)
  config$preprocess$prey_colname <- gsub(' ','.',config$preprocess$prey_colname)
  config$preprocess$pepcount_colname <- gsub(' ','.',config$preprocess$pepcount_colname)
  config$preprocess$mw_colname <- gsub(' ','.',config$preprocess$mw_colname)
  
  return(config)
}


main <- function(opt){
  # read and qc config file
  config = getConfig(opt$config)
  
  ##  create an outputdir if it doesn't exist 
  if(is.null(config$files$output_dir) || config$files$output_dir == '') config$files$output_dir = sprintf('%s/processed/',getwd())
  
  dir.create(config$files$output_dir, showWarnings = T)
  
  ## main switches between parts of the pipeline
  if(config$preprocess$enabled){
    cat(">> PREPROCESSING FILES\n")
    matrix_file = preprocess.main(data_file=config$files$data, keys_file=config$files$keys, output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/'), filter_data=config$preprocess$filter_contaminants, contaminants_file=config$preprocess$contaminants_file, rm_co=config$preprocess$remove_carryover, collapse_file=config$files$collapse, exclusions_file=config$files$specificity_exclusions, remove_file=config$files$remove, id_colname=config$preprocess$id_colname, prey_colname=config$preprocess$prey_colname, pepcount_colname=config$preprocess$pepcount_colname, mw_colname=config$preprocess$mw_colname)  
  }
  if(config$qc$enabled){
    cat(">> QUALITY CONTROL\n")
    if(!config$preprocess$enabled){ ## use previous data matrix instead of the one from pre-processing call 
      matrix_file = config$qc$matrix_file
    }
    qc.main(matrix_file=matrix_file, font_scale=config$qc$cluster_font_scale, cluster=config$qc$cluster, ip_dists=config$qc$ip_distributions)
  }
  if(config$mist$enabled){
    cat(">> MIST\n")
    if(!config$preprocess$enabled){ ## use previous data matrix instead of the one from pre-processing call 
      matrix_file = config$mist$matrix_file
    }
    mist.results = mist.main(matrix_file=matrix_file, weights=config$mist$weights, w_R=config$mist$reproducibility, w_A=config$mist$abundance, w_S=config$mist$specificity, training_file=config$mist$training_file, training_steps=config$mist$training_steps)
    output_file = gsub('.txt', "_MIST.txt", matrix_file)
    write.table(mist.results, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  }
  
  cat(">> SCORING FINISHED!\n")
}

if(!exists('DEBUG') || DEBUG==F) main(opt)

