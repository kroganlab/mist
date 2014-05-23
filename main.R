#! /usr/bin/Rscript --vanilla
suppressMessages(library(getopt))

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

###############################
suppressMessages(library(yaml))

PIPELINE=T

## load all externeal files
source("src/preprocess.R")
source("src/qc.R")
source("src/mist.R")

getConfig <- function(config_file){
  x = readLines(config_file)
  x = x[x!=""]  #remove \n\n cases (blank Lines)
  x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
  x = gsub('\t', '  ', x)
  config = paste(x, collapse='\n')
  config = yaml.load(config)
  return(config)
}

main <- function(opt){
  config = tryCatch(getConfig(opt$config), error = function(e) { print("!!! Error loading the config file. Please make sure the file follows YAML format."); break} )
  
  ## set the input dir if defined
  if(is.null(config$files$dir) || config$files$dir == '') config$files$dir = getwd()
  
  ##  create an outputdir if it doesn't exist 
  if(is.null(config$files$output_dir) || config$files$output_dir == '') config$files$output_dir = sprintf('%s/processed/',config$files$dir)
    
  dir.create(config$files$output_dir, showWarnings = T)
  
  ## main switches between parts of the pipeline
  if(config$preprocess$enabled){
    cat(">> PREPROCESSING FILES\n")
    matrix_file = preprocess.main(data_file=config$files$data, keys_file=config$files$keys, output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/'), filter_data=config$preprocess$filter_contaminants, contaminants_file=config$preprocess$contaminants_file, rm_co=config$preprocess$remove_carryover, collapse_file=config$files$collapse, exclusions_file=config$files$specificity_exclusions, remove_file=config$files$remove, prey_colname=config$preprocess$prey_colname, pepcount_colname=config$preprocess$pepcount_colname)  
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
    mist.main(matrix_file=matrix_file, weights=config$mist$weights, w_R=config$mist$reproducibility, w_A=config$mist$abundance, w_S=config$mist$specificity, training_file=config$mist$training_file)
  }
}

if(!exists('DEBUG') || DEBUG==F) main(opt)

