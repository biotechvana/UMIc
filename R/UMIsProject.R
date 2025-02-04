#!/usr/bin/Rscript

########## Loading libraries ########## 
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(stringdist)
library(pryr)

rm(list = ls())

source("R/casesWorkflows.R")
source("R/functions.R")

########## Inputs ##########

args = commandArgs(trailingOnly=TRUE)

if (length(args)<9) {
  stop("Need all arguments: input_fodler output_folder pariedData UMIlocation UMIlength sequenceLength countsCutoff UMIdistance sequenceDistance", call.=FALSE)
} else{

  #inputs folder / working directory
  inputsFolder <- args[1]
  
  #outputs folder
  outputsFolder <- args[2]
  
  #type of data - paired or single
  pairedData <- args[3] == "T"
  
  #UMI located in Read1 --> "R1"
  #UMI located in Read1 and Read2 --> "R1 & R2"
  #UMI located in Read1 --> "R2"
  UMIlocation <- args[4]
  
  #length of the UMI
  UMIlength <- as.numeric(args[5])
  
  #length of th sequence
  sequenceLength <- as.numeric(args[6])
  
  #min read counts per UMI, for initial data cleaning
  countsCutoff <- as.numeric(args[7])
  
  #max UMI distance for UMI merging
  UMIdistance <- as.numeric(args[8])
  
  #max sequence distance for UMI correction
  sequenceDistance <- as.numeric(args[9])
  
}

########## Run the appropriate scenario ##########


if (pairedData & UMIlocation == "R1"){   #case 1 -- paired data and UMI only in Read1
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    commonPart <- as.character(str_split(file1,"R1", simplify = T))
    commonPart <- commonPart[length(commonPart)]
    
    file2 <- str_replace(file1,paste0("R1",commonPart),paste0("R2",commonPart))
    
    filepath1 <- paste0(inputsFolder,"/",file1)
    filepath2 <- paste0(inputsFolder,"/",file2)
    
    pairedR1(filepath1, 
             filepath2, 
             outputsFolder, 
             UMIlength, 
             UMIdistance, 
             sequenceLength, 
             sequenceDistance, 
             countsCutoff)
    
    inputFiles <- inputFiles[str_detect(inputFiles,paste0(file1,"|",file2), negate = T)]
  }

} else if (pairedData & UMIlocation == "R2"){ #case 2 -- paired data and UMI only in Read2
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
  
    file1 <- inputFiles[1]
  
    commonPart <- as.character(str_split(file1,"R1", simplify = T))
    commonPart <- commonPart[length(commonPart)]
  
    file2 <- str_replace(file1, paste0("R1", commonPart), paste0("R2", commonPart))
  
    filepath1 <- paste0(inputsFolder,"/",file1)
    filepath2 <- paste0(inputsFolder,"/",file2)
  
    pairedR2(filepath1,
             filepath2,
             outputsFolder,
             UMIlength, 
             UMIdistance, 
             sequenceLength, 
             sequenceDistance, 
             countsCutoff)
  
  inputFiles <- inputFiles[str_detect(inputFiles,paste0(file1,"|",file2), negate = T)]
  }   
  
} else if (pairedData & UMIlocation == "R1 & R2"){   #case 2 -- paired data and UMI in Read1 and Read2
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    commonPart <- as.character(str_split(file1, "R1", simplify = T))
    commonPart <- commonPart[length(commonPart)]
    
    file2 <- str_replace(file1, 
                         paste0("R1", commonPart), 
                         paste0("R2", commonPart))
    
    filepath1 <- paste0(inputsFolder, "/", file1)
    filepath2 <- paste0(inputsFolder, "/", file2)
    
    cat(c("Files:", file1, file2, "\n"))

    pairedR1R2(filepath1, 
               filepath2, 
               outputsFolder, 
               UMIlength, 
               UMIdistance, 
               sequenceLength, 
               sequenceDistance, 
               countsCutoff)
    
    
    inputFiles <- inputFiles[str_detect(inputFiles, paste0(file1, "|", file2), negate = T)]
  }

  
} else if (!pairedData){  #case 3 -- single data
  
  inputFiles <- list.files(inputsFolder, pattern = "fastq") 
  
  while(length(inputFiles) > 0){
    
    file1 <- inputFiles[1]
    
    filepath1 <- paste0(inputsFolder,"/",file1)

    single(filepath1, 
           outputsFolder, 
           UMIlength, 
           UMIdistance, 
           sequenceLength, 
           sequenceDistance, 
           countsCutoff)
    
    inputFiles <- inputFiles[str_detect(inputFiles,file1, negate = T)]
  }
 
}  
