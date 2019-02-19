ValidateSubmission <- function(submission.files = submission.files, nrow=nrow, ncol=ncol){
  valid.header <- c("Nucleotide_position", "Variant", "P-value_MO", "Standard_deviation", "P-value_MO+WT", "Standard_deviation", "Functional_effect", "Confidence:", "Comments")
  
  files.path <- submission.files
  for (i in 1:length(files.path)) {
    infile <- files.path[i]
    cat('checking file', infile, '\n')
    #check if file exists
    if (!file.exists(infile)) stop("File ",infile," does not exist!")
    # and if file can be read
    if (file.access(infile,mode=4) < 0) stop("No permission to read file ",infile,"!")
    
    #read the lines in the file
    lines <- scan(infile,what=character(),sep="\n")
    if (length(lines) != nrow) stop("File does not have the required 39 lines")
    
    #break lines into data fields
    fields <- strsplit(lines,"\t")
    #check if each line has six fields
    if (any(sapply(fields,length) != ncol)) stop("Must have 8 tab-delimited columns including Nucleotide position and Variant!")
    #Comments column is not considered since none used it
    #Comments are optional

    
    #check if the header is valid
    header <- fields[[1]]
    if (any(header != valid.header)) stop("Invalid table headers!")
    
    #extract data columns for the remaining lines
    body <- fields[-1]
    body.cols <- lapply(1:8,function(i) sapply(body,`[[`,i))
    
    
    #check if the first column contains valid allele descriptors, where ^ marks start and $ marks end
    if (any(regexpr("^\\w\\d+\\w(,\\w\\d+\\w)*$",body.cols[[1]]) != 1)) {
      stop("First column must contain valid alleles!")
    }
    
    #check if the second (variant) column contains valid allele descriptors
    if (any(regexpr("^\\w\\d+\\w(,\\w\\d+\\w)*$",body.cols[[2]]) != 1)) {
      stop("Second column must contain valid alleles!")
    }
    
    #check if the third column "P-value_MO" is numerical values or stars
    if (!all(sapply(body.cols[[3]],function(x) x == "*" || !is.na(as.numeric(x))))) {
      stop("Third column must be numeric or \"*\"!")
    }
    
    #check if the fourth column "Standard_deviation" is numerical values or stars
    if (!all(sapply(body.cols[[4]],function(x) x == "*" || (!is.na(as.numeric(x))) ))) {
      stop("Fourth column must be numeric or \"*\"!")
    }
    
    #check if the fifth column "P-value_MO+WT" is numerical values or stars
    if (!all(sapply(body.cols[[5]],function(x) x == "*" || (!is.na(as.numeric(x))) ))) {
      stop("Fifth column must be numeric or \"*\"!")
    }
    
    #check if the sixth column "Standard_deviation" is numerical values or stars
    if (!all(sapply(body.cols[[6]],function(x) x == "*" || (!is.na(as.numeric(x))) ))) {
      stop("Sixth column must be numeric or \"*\"!")
    }
    
    #check if the seventh column "Functional_effect" is numerical values or stars
    if (!all(sapply(body.cols[[7]],function(x) x == "*" || (!is.na(as.numeric(x))) ))) {
      system(paste('cat', infile, '| sed "s/benign/0/g" | sed "s/hypomorphic/1/g" | sed "s/pathogenic/2/g" > corrected_file'))
      system(paste('cat', infile, '| sed "s/Benign/0/g" | sed "s/Hypomorphic/1/g" | sed "s/Pathogenic/2/g" > corrected_file')) # for submission 6
      stop("Seventh column must be numeric or \"*\"!")
    }
    
    #check if the eigth column "Confidence:" is numerical values or stars
    if (!all(sapply(body.cols[[8]],function(x) x == "*" || (!is.na(as.numeric(x))) ))) {
      stop("eigth column must be numeric or \"*\"!")
    }
    
    cat('File format check passed!', '\n')
  }

}

validateClasses <- function(submission.files = submission.files)
{
  threshold <- 0.5
  for (i in 1:length(submission.files)) {
    print(paste("Check submission: ", submission.files[i]))
    current.submission <- read.table(submission.files[i], sep = '\t', header = TRUE) #file with experimental data
    p.MO <- current.submission$P.value_MO
    p.MO.WT <- current.submission$P.value_MO.WT
    class <- current.submission$Functional_effect
    LOF.condition <- which(p.MO > threshold & p.MO.WT <= threshold)
    pLOF.condition <- which(p.MO <= threshold & p.MO.WT <= threshold)
    benign.condition <- which(p.MO <= threshold & p.MO.WT > threshold)
    other.condition <- which(p.MO > threshold & p.MO.WT > threshold)
    if(length(which(class[LOF.condition] != "2")) != 0)
      print(paste("Check class 2:", paste(which(class[LOF.condition] != "2"), collapse = " ")))
    if(length(which(class[pLOF.condition] != "1")) != 0)
      print(paste("Check class 1:", paste(which(class[pLOF.condition] != "1"), collapse = " ")))
    if(length(which(class[benign.condition] != "0")) != 0)
      print(paste("Check class 0:", paste(which(class[benign.condition] != "0"), collapse = " ")))
    if(length(other.condition) != 0)
      print(paste("Check other condition:", paste(other.condition, collapse = " ")))
  }
}


