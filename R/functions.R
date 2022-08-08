groupingPairedR1 <- function(intermediate.table, # r1, 
                             # counts, 
                             full, 
                             quality, 
                             full2, 
                             quality2, 
                             UMIlength){
  
  intermediate.table.c1 = intermediate.table[which(intermediate.table$count == 1), ]
  intermediate.table.c2 = intermediate.table[which(intermediate.table$count != 1), ]
  
  rm(intermediate.table)
  
  result = list()
  
  if(nrow(intermediate.table.c1) > 0){
    
      f1 = full[which(full$UMI %in% intermediate.table.c1$UMI),] 
      f1 = f1[order(f1$id), ]
      
      f2 = full2[which(full2$id %in% f1$id), ]
      f2 = f2[order(f2$id), ]
      
      qr1 = quality[which(quality$id %in% f1$id), ]
      qr1 = qr1[order(qr1$id), ]
      qr1 = as.matrix(qr1[,1:(ncol(qr1)-1)])
      qr1 = qr1 + 33
      qr1 = base::apply(qr1, 1, intToUtf8)
      qr1 = as.character(qr1)
      
      qr2 = quality2[which(quality2$id %in% f1$id), ]
      qr2 = qr2[order(qr2$id), ]
      qr2 = as.matrix(qr2[,1:(ncol(qr2)-1)])
      qr2 = qr2 + 33
      qr2 = base::apply(qr2, 1, intToUtf8)
      qr2 = as.character(qr2)
      
      result.1 = data.table(UMI = f1$UMI12,
                            read1 = f1$read, quality1 = qr1, 
                            read2 = f2$read, quality2 = qr2)
      
      colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
      
      rm(f1, f2, qr1, qr2)
      
      result[[1]] = result.1
  
  }
  
  if(nrow(intermediate.table.c2) > 0){
    
    for(i in 1:nrow(intermediate.table.c2)){
      
      r1 = intermediate.table.c2$UMI[i]
      
      #reads with specific UMI
      grouping = full[which(full$UMI == r1), ]
      grouping2 <- full2[which(full2$id %in% grouping$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      #File 1
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      #File 2
      grouping2$UMI <- grouping$UMI
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      rm(grouping, grouping2, quality.1, quality.2)
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = substr(r1,1,UMIlength),
                           read1 = result1[1,1], quality1 = result1[1,2], 
                           read2 = result2[1,1], quality2 = result2[1,2])
      
      colnames(result.2) <- c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i+1]] = result.2
    }
  }
  
  result = rbindlist(result)
  
  return(result)
}

groupingPairedR2 <- function(intermediate.table,
                             full,
                             quality,
                             full2,
                             quality2,
                             UMIlength){
  
  intermediate.table.c1 = intermediate.table[which(intermediate.table$count == 1), ]
  intermediate.table.c2 = intermediate.table[which(intermediate.table$count != 1), ]
  
  result = list()
  
  if(nrow(intermediate.table.c1) > 0){
    
    f2 = full2[which(full2$UMI %in% intermediate.table.c1$UMI),] 
    f2 = f2[order(f2$id), ]
    
    f1 = full[which(full$id %in% f2$id), ]
    f1 = f1[order(f1$id), ]
    
    qr1 = quality[which(quality$id %in% f1$id), ]
    qr1 = qr1[order(qr1$id), ]
    qr1 = as.matrix(qr1[,1:(ncol(qr1)-1)])
    qr1 = qr1 + 33
    qr1 = base::apply(qr1, 1, intToUtf8)
    qr1 = as.character(qr1)
    
    qr2 = quality2[which(quality2$id %in% f1$id), ]
    qr2 = qr2[order(qr2$id), ]
    qr2 = as.matrix(qr2[,1:(ncol(qr2)-1)])
    qr2 = qr2 + 33
    qr2 = base::apply(qr2, 1, intToUtf8)
    qr2 = as.character(qr2)
    
    result.1 = data.table(UMI = f2$UMI12,
                          read1 = f1$read, quality1 = qr1, 
                          read2 = f2$read, quality2 = qr2)
    
    colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
    
    rm(f1, f2, qr1, qr2)
    
    result[[1]] = result.1
    
  }
  
  if(nrow(intermediate.table.c2) > 0){
    
    for(i in 1:nrow(intermediate.table.c2)){
      
      r1 = intermediate.table.c2$UMI[i]
      
      #reads with specific UMI
      grouping2 = full2[which(full2$UMI == r1), ]
      grouping <- full[which(full$id %in% grouping2$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      grouping2 = grouping2[order(grouping2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      #File 1
      grouping$UMI <- grouping2$UMI
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      #File 2
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      rm(grouping, grouping2, quality.1, quality.2)
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = substr(r1,1,UMIlength),
                             read1 = result1[1,1], quality1 = result1[1,2], 
                             read2 = result2[1,1], quality2 = result2[1,2])
      
      colnames(result.2) <- c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i+1]] = result.2
    }
  }
  
  result_mean = rbindlist(result)
  
  return(result_mean)
}


groupingPairedR1R2 <- function(intermediate.table, # counts, 
                               full, 
                               quality, 
                               full2, 
                               quality2){
  
  intermediate.table.c1 = intermediate.table[which(intermediate.table$count == 1), ]
  intermediate.table.c2 = intermediate.table[which(intermediate.table$count != 1), ]
  
  rm(intermediate.table)
  
  result = list()
  
  if(nrow(intermediate.table.c1) > 0){
    
    f1 = full[which(full$UMI12 %in% intermediate.table.c1$UMI12),]
    f1 = f1[order(f1$id), ]
    
    f2 = full2[which(full2$id %in% f1$id), ]
    f2 = f2[order(f2$id), ]
    
    qr1 = quality[which(quality$id %in% f1$id), ]
    qr1 = qr1[order(qr1$id), ]
    qr1 = as.matrix(qr1[,1:(ncol(qr1)-1)])
    qr1 = qr1 + 33
    qr1 = base::apply(qr1, 1, intToUtf8)
    qr1 = as.character(qr1)
    
    qr2 = quality2[which(quality2$id %in% f1$id), ]
    qr2 = qr2[order(qr2$id), ]
    qr2 = as.matrix(qr2[,1:(ncol(qr2)-1)])
    qr2 = qr2 + 33
    qr2 = base::apply(qr2, 1, intToUtf8)
    qr2 = as.character(qr2)
    
    result.1 = data.table(UMI = f1$UMI12,
                          read1 = f1$read, quality1 = qr1, 
                          read2 = f2$read, quality2 = qr2)
    
    colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
    
    rm(f1, f2, qr1, qr2)
    
    result[[1]] = result.1
  }
  
  if(nrow(intermediate.table.c2) > 0){
    for(i in 1:nrow(intermediate.table.c2)){
      
      r1 = intermediate.table.c2$UMI12[i]

      # reads with specific UMI
      grouping = full[which(full$UMI12 == r1), -"UMI12"]
      grouping2 = full2[which(full2$id %in% grouping$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      grouping2 = grouping2[order(grouping2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      # File 1
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      # File 2
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      # rm(grouping, grouping2, quality.1, quality.2)
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = r1,
                           read1 = result1[1,1], quality.1 = result1[1,2], 
                           read2 = result2[1,1], quality.2 = result2[1,2])
      
      colnames(result.2) = c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i + 2]] = result.2

    }
  }
  
  result = data.table::rbindlist(result)
  
  return(result)
  
  
  # if(counts == 1) {
  #   
  #   temp.result <- full[which(full$UMI12 == r1),]
  #   read1 <- temp.result$read
  #   read2 <- full2[which(full2$id == temp.result$id), read]
  #   
  #   quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
  #   quality.read1 <- as.numeric(quality.read1) + 33
  #   quality.read1 <- intToUtf8(quality.read1)
  #   
  #   quality.read2 <- as.character(quality2[which(quality2$id == temp.result$id), 1:(ncol(quality2)-1)])
  #   quality.read2 <- as.numeric(quality.read2) + 33
  #   quality.read2 <- base::intToUtf8(quality.read2)
  #   
  #   result <- data.table(UMI = r1,
  #                        read1 = read1, quality1 = quality.read1, 
  #                        read2 = read2, quality2 = quality.read2)
  #   
  #   colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  #   
  # } else {
  #   
  #   # reads with specific UMI
  #   grouping = full[which(full$UMI12 == r1), -"UMI12"]
  #   grouping2 = full2[which(full2$id %in% grouping$id), ]
  #   
  #   quality = quality[which(quality$id %in% grouping$id), ]
  #   quality2 = quality2[which(quality2$id %in% grouping2$id), ]
  #   
  #   grouping = grouping[order(grouping$id), ]
  #   grouping2 = grouping2[order(grouping2$id), ]
  #   
  #   quality = quality[order(quality$id), ]
  #   quality2 = quality2[order(quality2$id), ]
  #   
  #   # File 1
  #   grouping_q = cbind(grouping, quality[,-c("id")])
  #   
  #   # File 2
  #   grouping_q2 = cbind(grouping2, quality2[,-c("id")])
  #   
  #   rm(grouping, grouping2, quality, quality2)
  #   
  #   result1 <- calculationsFunction(grouping_q) 
  #   result2 <- calculationsFunction(grouping_q2)
  #   
  #   result <- data.table(UMI = r1,
  #                        read1 = result1[1,1], quality1 = result1[1,2], 
  #                        read2 = result2[1,1], quality2 = result2[1,2])
  #   
  #   colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  # }
  
}

groupingSingle <- function(intermediate.table, # r1, 
                               # counts, 
                               full, 
                               quality,  
                               UMIlength){
    
    intermediate.table.c1 = intermediate.table[which(intermediate.table$count == 1), ]
    intermediate.table.c2 = intermediate.table[which(intermediate.table$count != 1), ]
    
    rm(intermediate.table)
    
    result = list()
    
    if(nrow(intermediate.table.c1) > 0){
      
      f1 = full[which(full$UMI %in% intermediate.table.c1$UMI),] 
      f1 = f1[order(f1$id), ]
      
      qr1 = quality[which(quality$id %in% f1$id), ]
      qr1 = qr1[order(qr1$id), ]
      qr1 = as.matrix(qr1[,1:(ncol(qr1)-1)])
      qr1 = qr1 + 33
      qr1 = base::apply(qr1, 1, intToUtf8)
      qr1 = as.character(qr1)
      
      result.1 = data.table(UMI = f1$UMI12,
                            read1 = f1$read, quality1 = qr1)
      
      colnames(result.1) = c("UMI" , "read1", "quality1")
      
      rm(f1, qr1)
      
      result[[1]] = result.1
      
    }
    
    if(nrow(intermediate.table.c2) > 0){
      
      for(i in 1:nrow(intermediate.table.c2)){
        
        r1 = intermediate.table.c2$UMI[i]
        
        #reads with specific UMI
        grouping = full[which(full$UMI == r1), ]
        
        quality.1 = quality[which(quality$id %in% grouping$id), ]
        
        grouping = grouping[order(grouping$id), ]
        
        quality.1 = quality.1[order(quality.1$id), ]
        
        #File 1
        grouping_q = cbind(grouping, quality.1[,-c("id")])
        
        rm(grouping, quality.1)
        
        result1 <- calculationsFunction(grouping_q) 
        
        result.2 <- data.table(UMI = substr(r1,1,UMIlength),
                               read1 = result1[1,1], quality1 = result1[1,2])
        
        colnames(result.2) <- c("UMI" , "read1", "quality1")
        
        result[[i+1]] = result.2
      }
    }
    
    result = rbindlist(result)
    
    return(result)
  
  #if ((counts == 1)){
    
   # temp.result <- full[which(full$UMI == r1),]
   #  read1 <- temp.result$read
    
   # quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
   # quality.read1 <- as.numeric(quality.read1) + 33
   # quality.read1 <- intToUtf8(quality.read1)
    
   #  result <- data.table(UMI = substr(r1,1,UMIlength),
   #                      read1 = read1, quality1 = quality.read1)
    
  #  colnames(result) <- c("UMI" , "read1", "quality1")
    
  #}else{
  
    #print(r1)
    
    #reads with specific UMI
   # grouping = full[which(full$UMI == r1), ]
    
   # quality = quality[which(quality$id %in% grouping$id), ]
    
   # grouping = grouping[order(grouping$id), ]
    
   # quality = quality[order(quality$id), ]
    
    #File 1
   # grouping_q = cbind(grouping, quality[,-c("id")])
    
   # rm(grouping, quality)
    
   # result1 <- calculationsFunction(grouping_q) 
    
   # result <- data.table(UMI = substr(r1,1,UMIlength),
    #                     read1 = result1[1,1], quality1 = result1[1,2])
    
   # colnames(result) <- c("UMI" , "read1", "quality1")
 # }
  
  #return(result)
}

groupingFinalPairedR1 <- function(newUMIs, # r1, 
                                  full, 
                                  quality, 
                                  full2, 
                                  quality2, 
                                  first_consensus, 
                                  UMIlength){
  
  newUMIs.1 = newUMIs[which(str_length(newUMIs$UMI) == UMIlength), ]
  newUMIs.2 = newUMIs[which(str_length(newUMIs$UMI) != UMIlength), ]
  
  rm(newUMIs)
  
  result = list()
  
  if(nrow(newUMIs.1) > 0){
    
    result.1 = first_consensus[which(first_consensus$UMI %in% newUMIs.1$UMI),]
    colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
    
    result[[1]] = result.1
    
  }
  
  if(nrow(newUMIs.2) > 0){
    
    for(i in 1:nrow(newUMIs.2)){
      
      r1 = newUMIs.2$UMI[i]
      
      #reads with specific UMI
      grouping = full[str_detect(full$UMI,as.character(r1)), ]
      grouping2 <- full2[which(full2$id %in% grouping$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      grouping2 = grouping2[order(grouping2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      #File 1
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      #File 2
      grouping2$UMI <- grouping$UMI
      
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      rm(grouping, grouping2, quality.1, quality.2)
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = substr(r1,1,12),
                           read1 = result1[1,1], quality1 = result1[1,2], 
                           read2 = result2[1,1], quality2 = result2[1,2])
      
      colnames(result.2) <- c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i+2]] = result.2
    }
    
  }
  
  result = rbindlist(result)
  
  return(result)
}

groupingFinalPairedR2 <- function(newUMIs,
                                  full,
                                  quality,
                                  full2,
                                  quality2,
                                  result_mean,
                                  UMIlength){
  newUMIs.1 = newUMIs[which(str_length(newUMIs$UMI) == UMIlength), ]
  newUMIs.2 = newUMIs[which(str_length(newUMIs$UMI) != UMIlength), ]
  
  rm(newUMIs)
  
  result = list()
  
  if(nrow(newUMIs.1) > 0){
    
    result.1 = result_mean[which(result_mean$UMI %in% newUMIs.1$UMI),]
    colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
    
    result[[1]] = result.1
    
  }
  
  if(nrow(newUMIs.2) > 0){
    
    for(i in 1:nrow(newUMIs.2)){
      
      r1 = newUMIs.2$UMI[i]
      
      #reads with specific UMI
      grouping = full[str_detect(full$UMI,as.character(r1)), ]
      grouping2 <- full2[which(full2$id %in% grouping$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      grouping2 = grouping2[order(grouping2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      #File 1
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      #File 2
      grouping2$UMI <- grouping$UMI
      
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      rm(grouping, grouping2, quality.1, quality.2)
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = substr(r1,1,12),
                             read1 = result1[1,1], quality1 = result1[1,2], 
                             read2 = result2[1,1], quality2 = result2[1,2])
      
      colnames(result.2) <- c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i+2]] = result.2
    }
    
  }
  
  consensus_mean = rbindlist(result)
  
  return(consensus_mean)
}

groupingFinalPairedR1R2 <- function(newUMIs, # r1, 
                                    full, 
                                    quality, 
                                    full2, 
                                    quality2, 
                                    first_consensus, 
                                    UMIlength){
  
  newUMIs.1 = newUMIs[which(str_length(newUMIs$UMI) == 2*UMIlength), ]
  newUMIs.2 = newUMIs[which(str_length(newUMIs$UMI) != 2*UMIlength), ]
  
  rm(newUMIs)
  
  result = list()
  
  if(nrow(newUMIs.1) > 0){
    result.1 = first_consensus[which(first_consensus$UMI %in% newUMIs.1$UMI),]
    colnames(result.1) = c("UMI" , "read1", "quality1", "read2", "quality2")
    
    result[[1]] = result.1
  }
  
  if(nrow(newUMIs.2) > 0){
    for(i in 1:nrow(newUMIs.2)){
      
      r1 = newUMIs.2$UMI[i]
      
      # reads with specific UMI
      grouping = full[str_detect(full$UMI12, r1), -"UMI12"]
      grouping2 = full2[which(full2$id %in% grouping$id), ]
      
      quality.1 = quality[which(quality$id %in% grouping$id), ]
      quality.2 = quality2[which(quality2$id %in% grouping2$id), ]
      
      grouping = grouping[order(grouping$id), ]
      grouping2 = grouping2[order(grouping2$id), ]
      
      quality.1 = quality.1[order(quality.1$id), ]
      quality.2 = quality.2[order(quality.2$id), ]
      
      # File 1
      grouping_q = cbind(grouping, quality.1[,-c("id")])
      
      # File 2
      grouping_q2 = cbind(grouping2, quality.2[,-c("id")])
      
      result1 <- calculationsFunction(grouping_q) 
      result2 <- calculationsFunction(grouping_q2)
      
      result.2 <- data.table(UMI = substr(r1,1,2*UMIlength),
                           read1 = result1[1,1], quality.1 = result1[1,2], 
                           read2 = result2[1,1], quality.2 = result2[1,2])
      
      colnames(result.2) <- c("UMI" , "read1", "quality1", "read2", "quality2")
      
      result[[i + 1]] = result.2
    }
  }
  
  result = rbindlist(result)
  
  return(result)
  
  # if (nchar(r1) == 2*UMIlength){
  #   
  #   result <- first_consensus[which(first_consensus$UMI == r1),]
  #   colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  #   
  # } else {
  #   
  #   # reads with specific UMI
  #   grouping = full[str_detect(full$UMI12, r1), -"UMI12"]
  #   grouping2 <- full2[which(full2$id %in% grouping$id), ]
  #   
  #   quality = quality[which(quality$id %in% grouping$id), ]
  #   quality2 = quality2[which(quality2$id %in% grouping2$id), ]
  #   
  #   grouping = grouping[order(grouping$id), ]
  #   grouping2 = grouping2[order(grouping2$id), ]
  #   quality = quality[order(quality$id), ]
  #   quality2 = quality2[order(quality2$id), ]
  #   
  #   # File 1
  #   grouping_q = cbind(grouping, quality[,-c("id")])
  #   
  #   # File 2
  #   grouping_q2 = cbind(grouping2, quality2[,-c("id")])
  #   
  #   rm(grouping, grouping2, quality, quality2)
  #   
  #   result1 <- calculationsFunction(grouping_q) 
  #   result2 <- calculationsFunction(grouping_q2)
  #   
  #   result <- data.table(UMI = substr(r1,1,2*UMIlength),
  #                        read1 = result1[1,1], quality1 = result1[1,2], 
  #                        read2 = result2[1,1], quality2 = result2[1,2])
  #   
  #   colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  # }
}

groupingFinalSingle <-function(newUMIs, # r1, 
                                    full, 
                                    quality, 
                                    first_consensus, 
                                    UMIlength){
    
    newUMIs.1 = newUMIs[which(str_length(newUMIs$UMI) == UMIlength), ]
    newUMIs.2 = newUMIs[which(str_length(newUMIs$UMI) != UMIlength), ]
    
    rm(newUMIs)
    
    result = list()
    
    if(nrow(newUMIs.1) > 0){
      
      result.1 = first_consensus[which(first_consensus$UMI %in% newUMIs.1$UMI),]
      colnames(result.1) = c("UMI" , "read1", "quality1")
      
      result[[1]] = result.1
      
    }
    
    if(nrow(newUMIs.2) > 0){
      
      for(i in 1:nrow(newUMIs.2)){
        
        r1 = newUMIs.2$UMI[i]
        
        #reads with specific UMI
        grouping = full[str_detect(full$UMI,as.character(r1)), ]
        
        quality.1 = quality[which(quality$id %in% grouping$id), ]
        
        grouping = grouping[order(grouping$id), ]
        
        quality.1 = quality.1[order(quality.1$id), ]
        
        #File 1
        grouping_q = cbind(grouping, quality.1[,-c("id")])
        
        rm(grouping, quality.1)
        
        result1 <- calculationsFunction(grouping_q) 
        
        result.2 <- data.table(UMI = substr(r1,1,12),
                               read1 = result1[1,1], quality1 = result1[1,2])
        
        colnames(result.2) <- c("UMI" , "read1", "quality1")
        
        result[[i+2]] = result.2
      }
      
    }
    
    result = rbindlist(result)
    
    return(result)
  
}

UMIcorrectionPairedR1 <- function(intermediate.table,
                                  first_consensus, 
                                  sequenceDistance, 
                                  UMIdistance,
                                  outputsFolder){
  dist.calc <- 0
  
  uniqueUMIs <- c()
  IDs_1 <- c()
  IDs_2 <- c()
  counts <- c()
  
  #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  while(nrow(intermediate.table) > 1){
    
    best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI[1]),]
    list.best <- best$UMI[1]
    list.counts <- intermediate.table$count[1]
    
    
    temp.intermediate = intermediate.table[2:nrow(intermediate.table), ]
    
    base_dist <- stringdist::stringdistmatrix(a = best$UMI[1], 
                                              b = temp.intermediate$UMI, 
                                              method = "hamming")[1,]
    
    dist.calc <- dist.calc + length(base_dist)
    
    who = which(base_dist <= UMIdistance)
    
    if(length(who) > 0){
      
      temp.intermediate = temp.intermediate[who,]
      
      temp_read1 = first_consensus[which(first_consensus$UMI %in% temp.intermediate$UMI),]$read1
      temp_read2 = first_consensus[which(first_consensus$UMI %in% temp.intermediate$UMI),]$read2
      
      dist1 = stringdist::stringdistmatrix(a = best$read1,
                                           b = temp_read1,
                                           method = "hamming")[1,]
      
      dist2 = stringdist::stringdistmatrix(a = best$read2,
                                           b = temp_read2,
                                           method = "hamming")[1,]
      
      dist.calc <- dist.calc + (2*length(dist1))
      
      who = which((dist1 <= sequenceDistance) & (dist2 <= sequenceDistance))
      
      if(length(who) > 0){
        
        temp.intermediate = temp.intermediate[who,]
        list.best = c(list.best, temp.intermediate$UMI)
        list.counts = c(list.counts, temp.intermediate$count)
        
        list.best <- paste(list.best, collapse = "|")
        list.counts <- paste(list.counts, collapse = "|")
        
      }
      
    }
    
    
    # iterations <- c(2:nrow(intermediate.table))
    # 
    # for (i in iterations){
    #   
    #   base_dist <- stringDist(c(best$UMI[1], intermediate.table$UMI[i]), method = "hamming")
    #   
    #   if (base_dist <= UMIdistance){
    #     
    #     temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read1]
    #     dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
    #     
    #     temp_read2 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read2]
    #     dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
    #     
    # 
    #     if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
    #       list.best <- paste0(list.best,"|",intermediate.table$UMI[i])
    #       list.counts <- paste0(list.counts,"|",intermediate.table$count[i])
    #     }
    #   }
    # }
    
    IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2, intermediate.table$ID2[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs, list.best)
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI, 
                                                        as.character(list.best), 
                                                        negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI[1])
    
  }
  
  newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, ID2 = IDs_2, Counts = counts))
  
  line = paste0("Number of Hamming distances calculated: ", dist.calc)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  return(newUMIs)
}

UMIcorrectionPairedR2 <- function(intermediate.table,
                                  result_mean,
                                  sequenceDistance,
                                  UMIdistance,
                                  outputsFolder){
  dist.calc <- 0
  
  uniqueUMIs <- c()
  IDs_1 <- c()
  IDs_2 <- c()
  counts <- c()
  
  #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  while(nrow(intermediate.table) > 1){
    
    best <- result_mean[which(result_mean$UMI == intermediate.table$UMI[1]),]
    list.best <- best$UMI[1]
    list.counts <- intermediate.table$count[1]
    
    
    temp.intermediate = intermediate.table[2:nrow(intermediate.table), ]
    
    base_dist <- stringdist::stringdistmatrix(a = best$UMI[1], 
                                              b = temp.intermediate$UMI, 
                                              method = "hamming")[1,]
    
    dist.calc <- dist.calc + length(base_dist)
    
    who = which(base_dist <= UMIdistance)
    
    if(length(who) > 0){
      
      temp.intermediate = temp.intermediate[who,]
      
      temp_read1 = result_mean[which(result_mean$UMI %in% temp.intermediate$UMI),]$read1
      temp_read2 = result_mean[which(result_mean$UMI %in% temp.intermediate$UMI),]$read2
      
      dist1 = stringdist::stringdistmatrix(a = best$read1,
                                           b = temp_read1,
                                           method = "hamming")[1,]
      
      dist2 = stringdist::stringdistmatrix(a = best$read2,
                                           b = temp_read2,
                                           method = "hamming")[1,]
      
      dist.calc <- dist.calc + (2*length(dist1))
      
      who = which((dist1 <= sequenceDistance) & (dist2 <= sequenceDistance))
      
      if(length(who) > 0){
        
        temp.intermediate = temp.intermediate[who,]
        list.best = c(list.best, temp.intermediate$UMI)
        list.counts = c(list.counts, temp.intermediate$count)
        
        list.best <- paste(list.best, collapse = "|")
        list.counts <- paste(list.counts, collapse = "|")
        
      }
      
    }
    
    
    # iterations <- c(2:nrow(intermediate.table))
    # 
    # for (i in iterations){
    #   
    #   base_dist <- stringDist(c(best$UMI[1], intermediate.table$UMI[i]), method = "hamming")
    #   
    #   if (base_dist <= UMIdistance){
    #     
    #     temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read1]
    #     dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
    #     
    #     temp_read2 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read2]
    #     dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
    #     
    # 
    #     if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
    #       list.best <- paste0(list.best,"|",intermediate.table$UMI[i])
    #       list.counts <- paste0(list.counts,"|",intermediate.table$count[i])
    #     }
    #   }
    # }
    
    IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2, intermediate.table$ID2[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs, list.best)
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI, 
                                                        as.character(list.best), 
                                                        negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI[1])
    
  }
  
  newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, ID2 = IDs_2, Counts = counts))
  
  line = paste0("Number of Hamming distances calculated: ", dist.calc)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  return(newUMIs)
}

UMIcorrectionPairedR1R2 <- function(intermediate.table, 
                                    first_consensus, 
                                    sequenceDistance, 
                                    UMIdistance, 
                                    UMIlength,
                                    outputsFolder){
  dist.calc <- 0
  uniqueUMIs <- c()
  IDs_1 <- c()
  IDs_2 <- c()
  counts <- c()
  
  # while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  
  while(nrow(intermediate.table) > 1){
    
    best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI12[1]),]
    list.best.UMI <- best$UMI[1]
    list.counts <- intermediate.table$count[1]

    # iterations <- c(2:nrow(intermediate.table))
    
    best.UMI1 <- substring(intermediate.table$UMI12[1], 1, UMIlength)
    best.UMI2 <- substring(intermediate.table$UMI12[1], UMIlength + 1, 2 * UMIlength)
    
    temp.intermediate = intermediate.table[2:nrow(intermediate.table),]
    
    temp.UMI1 <- substring(temp.intermediate$UMI12, 1, UMIlength)
    temp.UMI2 <- substring(temp.intermediate$UMI12, UMIlength + 1, 2 * UMIlength)
    
    base_dist1 = stringdist::stringdistmatrix(a = best.UMI1,
                                              b = temp.UMI1,
                                              method = "hamming")[1,]
    
    base_dist2 = stringdist::stringdistmatrix(a = best.UMI2,
                                              b = temp.UMI2,
                                              method = "hamming")[1,]
    dist.calc <- dist.calc +length(base_dist1) +length(base_dist2)
    
    
    who = which((base_dist1 <= UMIdistance) & (base_dist2 <= UMIdistance))
    
    if(length(who) > 0){
      temp.intermediate = temp.intermediate[who,]
      
      temp_read1 = first_consensus[which(first_consensus$UMI %in% temp.intermediate$UMI12),]$read1
      temp_read2 = first_consensus[which(first_consensus$UMI %in% temp.intermediate$UMI12),]$read2
      
      dist1 = stringdist::stringdistmatrix(a = best$read1[1], 
                                           b = temp_read1, 
                                           method = "hamming")[1,]
      
      dist2 = stringdist::stringdistmatrix(a = best$read2[1], 
                                           b = temp_read2, 
                                           method = "hamming")[1,]
      dist.calc <- dist.calc +length(dist1) +length(dist2)
      
      who = which((dist1 <= sequenceDistance) & (dist2 <= sequenceDistance))
      
      if(length(who) > 0){
        temp.intermediate = temp.intermediate[who,]
        list.best.UMI = c(list.best.UMI, temp.intermediate$UMI12)
        list.counts = c(list.counts, temp.intermediate$count)
        
        list.best.UMI = paste(list.best.UMI, collapse = "|")
        list.counts = paste(list.counts, collapse = "|")
      }
    }
    
    # for (i in 2:nrow(intermediate.table)){
    #   
    #   best.UMI1 <- substring(intermediate.table$UMI12[1], 1, UMIlength)
    #   temp.UMI1 <- substring(intermediate.table$UMI12[i], 1, UMIlength)
    # 
    #   best.UMI2 <- substring(intermediate.table$UMI12[1], UMIlength + 1, 2 * UMIlength)
    #   temp.UMI2 <- substring(intermediate.table$UMI12[i], UMIlength + 1, 2 * UMIlength)
    # 
    #   base_dist1 <- Biostrings::stringDist(c(best.UMI1, temp.UMI1), method = "hamming")
    #   base_dist2 <- Biostrings::stringDist(c(best.UMI2, temp.UMI2), method = "hamming")
    # 
    #   if ((base_dist1 <= UMIdistance) & (base_dist2 <= UMIdistance)){
    #    
    #     temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI12[i]), read1]
    #     dist1 <- Biostrings::stringDist(c(best$read1[1], temp_read1), method = "hamming") 
    #     
    #     temp_read2 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI12[i]), read2]
    #     dist2 <- Biostrings::stringDist(c(best$read2[1], temp_read2), method = "hamming") 
    #     
    #     if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
    #     
    #       list.best.UMI <- paste0(list.best.UMI, "|", intermediate.table$UMI12[i])
    #       list.counts <- paste0(list.counts, "|", intermediate.table$count[i])
    #       
    #     }
    #   }
    #   
    # }
    
    IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2, intermediate.table$ID2[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs, list.best.UMI)
  
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI12, as.character(list.best.UMI), 
                                                        negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 = append(IDs_1, intermediate.table$ID1[1])
    IDs_2 = append(IDs_2, intermediate.table$ID2[1])
    counts = append(counts, intermediate.table$count[1])
    uniqueUMIs = append(uniqueUMIs, intermediate.table$UMI12[1])
    
  }
  
  newUMIs = as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, ID2 = IDs_2, Counts = counts))
  
  line = paste0("Number of Hamming distances calculated: ", dist.calc)
  write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
  
  return(newUMIs)
}

UMIcorrectionSingle<- function(intermediate.table,
                                    first_consensus, 
                                    sequenceDistance, 
                                    UMIdistance,
                                    outputsFolder){
    dist.calc <- 0
  
    uniqueUMIs <- c()
    IDs_1 <- c()
    counts <- c()
    
    #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
    while(nrow(intermediate.table) > 1){
      
      best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI[1]),]
      list.best <- best$UMI[1]
      list.counts <- intermediate.table$count[1]
      
      
      temp.intermediate = intermediate.table[2:nrow(intermediate.table), ]
      
      base_dist <- stringdist::stringdistmatrix(a = best$UMI[1], 
                                                b = temp.intermediate$UMI, 
                                                method = "hamming")[1,]
      
      dist.calc <- dist.calc + length(base_dist)
      
      who = which(base_dist <= UMIdistance)
      
      if(length(who) > 0){
        
        temp.intermediate = temp.intermediate[who,]
        
        temp_read1 = first_consensus[which(first_consensus$UMI %in% temp.intermediate$UMI),]$read1
        
        dist1 = stringdist::stringdistmatrix(a = best$read1,
                                             b = temp_read1,
                                             method = "hamming")[1,]
        
        dist.calc <- dist.calc + (length(dist1))
        who = which(dist1 <= sequenceDistance) 
        
        if(length(who) > 0){
          
          temp.intermediate = temp.intermediate[who,]
          list.best = c(list.best, temp.intermediate$UMI)
          list.counts = c(list.counts, temp.intermediate$count)
          
          list.best <- paste(list.best, collapse = "|")
          list.counts <- paste(list.counts, collapse = "|")
          
        }
        
      }
      
      
      IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
      counts <- append(counts, list.counts)
      uniqueUMIs <- append(uniqueUMIs, list.best)
      intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI, 
                                                          as.character(list.best), 
                                                          negate = T), ]
      
    }
    
    if (!is.null(intermediate.table)){
      
      IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
      counts <- append(counts, intermediate.table$count[1])
      uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI[1])
      
    }
    
    newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1,  Counts = counts))
    
    line = paste0("Number of Hamming distances calculated: ", dist.calc)
    write(line, paste0(outputsFolder,"/extra_info.txt"), append = T)
    return(newUMIs)  
  
}



one.run.calculationsFunction <- function(one.base){
  
  # setDT(exampleone.base)
  one.base = one.base[,.(mean = mean(V2), 
                         count = .N, 
                         perc_qual = mean(V2) / 93*100,
                         perc_count = .N/nrow(one.base)*100), by = list(V1)]
  
  one.base$criterion = rowMeans(one.base[,c("perc_qual", "perc_count")])
  
  #####!!! what happens if we have the same maximum criterion for different V1 !!!##########
  
  one.base <- one.base[which(one.base$criterion == max(one.base$criterion)), ]
  
  return(data.table(V1 = as.character(one.base[which.max(one.base$perc_qual),]$V1), 
                    mean = round(one.base[which.max(one.base$perc_qual),]$mean)))
  
  # return(data.table(V1 = as.character(one.base[which.max(one.base$criterion), ]$V1)[1], 
  #                   mean = round(one.base[which.max(one.base$criterion), ]$mean)[1]))
  
}

calculationsFunction <- function(grouping_q){
  
  # list with all nts
  cons = str_split(grouping_q$read, pattern = "", simplify = TRUE)
  cons = as.list(as.data.table(cons))
  
  # list with all qualities
  grouping_q = as.list(grouping_q[,4:ncol(grouping_q)])
  
  # merge lists 
  grouping_q = base::Map(data.table, cons, grouping_q)
  
  rm(cons)
  
  ##
  cons_corr = lapply(grouping_q, one.run.calculationsFunction)
  cons_corr = rbindlist(cons_corr)
  
  #join again in one final sequence
  consensus = str_c(cons_corr$V1, collapse = "")
  meanQuality = as.numeric(cons_corr$mean) + 33 
  meanQuality = intToUtf8(meanQuality)
  result = data.table(seq = consensus, qual = meanQuality)
  
  return(result)
  
}
