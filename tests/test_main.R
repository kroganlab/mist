library(testthat)
DEBUG=T
source('main.R')

config_small = 'tests/small/mist_small_test.yml'
parsedArgs$config_file = config_small
main(parsedArgs)