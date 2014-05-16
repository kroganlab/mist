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



config_file = "~/projects/mist/APMS_TEMPLATE.yml"
#config = yaml.load_file("~/projects/mist/APMS_TEMPLATE.yml")
main(config_file)



