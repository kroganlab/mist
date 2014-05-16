#! /usr/bin/Rscript --vanilla
suppressMessages(library(yaml))

getConfig <- function(config_file){
  x = readLines(config_file)
  x = x[x!=""]  #remove \n\n cases (blank Lines)
  x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
  x = gsub('\t', '  ', x)
  config = paste(x, collapse='\n')
  config = yaml.load(config)
  return(config)
}

main <- function(config_file){
  # load config file
  config = tryCatch(getConfig(config_file), error = function(e) { print("!!! Error in loading the config file. Please make sure the file follows YAML format."); break} )
  
}

option_list <- list(
  make_option(c("-h", "--help"),
              help="available arguments (this screen)"), 
  make_option(c("-c", "--config_file"),
              help="configuration file in YAML format")
)
parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

## test override
parsedArgs$config_file = "~/projects/mist/APMS_TEMPLATE.yml"

main(parsedArgs$config_file)