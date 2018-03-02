#!/usr/bin/env Rscript

### R script companion to species_tracker.pl

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))

start_date <- as.Date(args[1])
name <- args[2]
n <- as.integer(args[3])
rm(args)
