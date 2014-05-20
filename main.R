#! /usr/bin/Rscript --vanilla
suppressMessages(library(yaml))
suppressMessages(library(optparse))

PIPELINE=T

## load all externeal files
source("src/preprocess_data.R")

getConfig <- function(config_file){
  x = readLines(config_file)
  x = x[x!=""]  #remove \n\n cases (blank Lines)
  x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
  x = gsub('\t', '  ', x)
  config = paste(x, collapse='\n')
  config = yaml.load(config)
  return(config)
}

main <- function(parsedArgs){
  ## load config file
  config = tryCatch(getConfig(parsedArgs$config_file), error = function(e) { print("!!! Error loading the config file. Please make sure the file follows YAML format."); break} )
  
  ##  CREATE AN OUTOPUT DIRECTORY IF IT DOESNT EXIST YET
  #dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  dir.create(config$files$output_dir, showWarnings = T)
  
  ## main switches between parts of the pipeline
  if(config$preprocess$enabled){
    print(">> PREPROCESSING FILES")
    preprocess.main(data_file=config$files$data, keys_file=config$files$keys, output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/'), filter_data=config$general$filter_contaminants, rm_co=config$general$remove_carryover, nupsc_flag=config$general$spectral_counts, collapse_file=config$files$collapse, exclusions_file=config$files$specifity_exclusions, remove_file=config$files$remove)
  }
  
  if(config$qc$enabled){
    print(">> QUALITY CONTROL")
    qc.main(matrix_file=config$qc$matrix_file, font_scale=config$qc$cluster_font_scale, cluster=config$qc$cluster, ip_dists=config$qc$ip_distributions)
  }
  
  
}

option_list <- list(
  make_option(c("-h", "--help"),
              help="available arguments (this screen)"), 
  make_option(c("-c", "--config_file"),
              help="configuration file in YAML format")
)
parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
## TEST override
parsedArgs$config_file = paste0(getwd(),"/tests/APMS_TEST.yml")

## CALL MAIN WITH ALL ARGS
main(parsedArgs)