


##### csbtxtread function #####
# Generic function to read in population structure from a text file produced with the csb2txt program from the EBT tools.

csbtxtread <- function(csbtxtfile, envnames, statenames, stateindex=NULL, indexonly=FALSE){
  
  ### function input
  # csbtxtfile: path to a txt file created with the csb2txt function from the ebt tools with .txt as extention.
  # envnames: optional array with the names of the environmental variables in the file
  # statenames: optional array with the names of the individual state variables in the file
  # stateindex: index of the states to read in.
  # indexonly: if TRUE, only a summary of the states in the file is given.
  
  ## note that runtime of this function scales with the size of the document and the number of states to read
  
  # check packages
  require(tools)
  require(readr)
  
  # quick check csbtxtfile
  if(!file.exists(csbtxtfile)){stop("The csbtxt-file does not exist.")}
  if(!identical(file_ext(csbtxtfile), "txt")){ stop("The csbtxt-file should be a .txt file.")}
  
  message("Reading in data")
  
  # read data file
  textdata <- suppressWarnings(read_table(file=csbtxtfile, col_names=as.character(rep(NA, max(length(envnames), length(statenames)+1))), guess_max=300, progress = show_progress(), show_col_types=FALSE, skip_empty_rows = FALSE))
  #textdata <- read.table(csbtxtfile, header=FALSE, fill=TRUE, blank.lines.skip = FALSE)
  
  # check which lines are empty
  nalines <-  which(rowSums(is.na(textdata))==ncol(textdata))
  
  # calculate the start and end of each state based on the empty lines
  datalayout <- data.frame(startindex = c(1, nalines[which(nalines[2:length(nalines)] - nalines[1:(length(nalines)-1)] == 1)+1]+1),
                           endindex = c(nalines[which(nalines[2:length(nalines)] - nalines[1:(length(nalines)-1)] == 1)+1]-2, nrow(textdata)))
  
  datalayout$index <- 1:nrow(datalayout)
  
  # correct for possible empty lines at the end
  datalayout <- datalayout[which(datalayout$startindex < nrow(textdata)),]
  
  # extract dimentions from first structure
  nenv <- sum(!is.na(textdata[1,]))
  nstate <- sum(!is.na(textdata[3,]))
  nstruc <- nrow(datalayout)
  
  message("Number population structures: ", nstruc)
  message("Number of environmental variables: ", nenv)
  message("Number of individual state variables: ", nstate-1)
  
  # print summary if only indexes are requested
  if(indexonly){
    for(i in 1:nstruc){
      message("State", i, "has", datalayout$endindex[i]-datalayout$startindex[i]-2, "cohorts. First environmental variable value is", textdata[datalayout$startindex[i],1] )
    }
    return()
  }
  
  # check values of provided indexes
  if(length(stateindex)>0){
    if(!is.integer(stateindex)) stop("Provided state index should be an integers")
    stateindex <- stateindex[which(stateindex > 0 & stateindex <= nstruc)]
    if(length(stateindex)==0)stop(paste("provided state indexes should be between 0 and", nstruc))
    datalayout=datalayout[stateindex,]
  }
  
  # check provided names
  if(!is.null(envnames) & length(envnames) != nenv){
    warning(paste("Envnames ignored because wrong number of names was provided."))
    envnames <- NULL
  }
  if(!is.null(statenames) & length(statenames) != nstate - 1 ){
    warning(paste("Statenames ignorded because wrong number of names was provided"))
    statenames <- NULL
  }

  # separate the structure blocks
  struct <- apply(X=datalayout, MARGIN=1, FUN=function(x){
    
    # first line is the environmental variables
    envdata <- textdata[as.integer(x[1]),]
    envdata <- envdata[which(!is.na(envdata))]
    row.names(envdata) <- NULL
    if(!is.null(envnames) ){names(envdata) <- envnames}
    
    # extract data on structure
    allstrucdata <- textdata[(as.integer(x[1])+2):as.integer(x[2]),]
    
    # extract start and ends of every population
    struclayout <- data.frame(startindex = c(1, which(rowSums(is.na(allstrucdata))==ncol(allstrucdata))+1), endindex = c(which(rowSums(is.na(allstrucdata))==ncol(allstrucdata))-1, nrow(allstrucdata) ) )
    
    # if more than one population, split data into multiple populations
    if(nrow(struclayout)>1){
      strucdata <- apply(X=struclayout, MARGIN=1, FUN=function(y){
        
        # extract data from one population
        popdata <- as.data.frame(allstrucdata[as.integer(y[1]):as.integer(y[2]),])
        popdata <- popdata[,which(colSums(is.na(popdata))<nrow(popdata))]
        
        # set names of states
        row.names(popdata) <- NULL
        if(!is.null(statenames) ){colnames(popdata) <- c("density", statenames)}
        else{colnames(popdata) <- c("density", names(popdata)[1:(ncol(popdata)-1)])}
        
        return(popdata)
      })
    }
    else{
      strucdata <- allstrucdata[,which(colSums(is.na(allstrucdata))<nrow(allstrucdata))]
      row.names(strucdata) <- NULL
      if(!is.null(statenames) ){colnames(strucdata) <- c("density", statenames)}
      else{colnames(strucdata) <- c("density", names(strucdata)[1:(ncol(strucdata)-1)])}
      strucdata <- list(strucdata)
    }
    totaldata <- c(list(envdata), strucdata)
    names(totaldata) <- c("env", paste0("pop",1:nrow(struclayout)))
    
    return(totaldata)
  } )
  names(struct) <- paste0("struc", datalayout$index)
  return(struct)
}

