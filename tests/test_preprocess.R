library(testthat)

config = yaml.load(string=paste(readLines("tests/APMS_TEST.yml"),collapse='\n'))
config = yaml.load(string=paste(readLines("tests/entero/APMS_ENTERO.yml"),collapse='\n'))
