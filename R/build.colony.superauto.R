#' This function builds an automatic .DAT file with parameters that I (Jenna) am
#' using to test some simulated data. It should NOT be used to build .DAT files
#' for other projects, as it will NOT prompt for all of the options.
#' This is purely for my convenience and to speed up testing!
#'
#'
#' @param wd The directory where the generated file will be placed. The default
#' is the current working directory.
#' @param name The desired filename for the Colony input file (the default is
#' Colony2.DAT).
#' @param delim What is the delimiter for the input files (default is that they
#' are whitespace delimited).
#' @return A text file is produced. This file can be used by Colony2 as an
#' input file.
#' @author Owen R. Jones
#' @seealso \code{\link{run.colony}}
#' @references Wang, J. (2004) Sibship reconstruction from genetic data with
#' typing errors.  Genetics 166: 1963-1979.
#' @keywords manip
#' @export

build.colony.superauto <- function(wd=getwd(), 
                                   name,
                                   datasetname,
                                   delim = "",
                                   sample_size,
                                   num_loci,
                                   error_rates_path,
                                   genotypes_path,
                                   exclusion_path
                                   ){
  
  colonyfile = NULL
  
  #Functions used here
  is.whole <- function(a) {floor(a) == a}
  
  #######################################################
  #  ! C, Dataset name, Length<51
  #######################################################
  colonyfile$datasetname = datasetname
  write(paste(colonyfile$datasetname, "! C, Dataset name, Length<51"), name, append = FALSE)
  
  #######################################################
  #  ! C, Main output file name, Length<21
  #######################################################
  colonyfile$outfile = colonyfile$datasetname
  write(paste(colonyfile$outfile, "! C, Main output file name, Length<21"), name, append = TRUE)
  
  
  #######################################################
  #  ! C, Note to the project
  #######################################################
  
  colonyfile$note = ""
  
  #######################################################
  #  ! I, Number of offspring in the sample
  #######################################################
  colonyfile$n.offspring = sample_size
  write(paste(colonyfile$n.offspring, "! I, Number of offspring in the sample"), name, append = TRUE)
    
    if(length(colonyfile$n.offspring) != 0){
      #Whole number warning
      if(is.whole(colonyfile$n.offspring) == FALSE){
        utils::flush.console()
        colonyfile = colonyfile[which(names(colonyfile) != "n.offspring")]
        warning("The number of offspring must be a whole number!\n", immediate. = TRUE)
      }
    }
  
  
  
  #######################################################
  #  ! I, Number of loci
  #######################################################
  colonyfile$n.loci = num_loci
  write(paste(colonyfile$n.loci, "! I, Number of loci"), name, append = TRUE)
    
    if(length(colonyfile$n.loci) != 0){
      #Whole number warning
      if(is.whole(colonyfile$n.loci) == FALSE){
        utils::flush.console()
        colonyfile = colonyfile[which(names(colonyfile) != "n.loci")]
        warning("The number of loci must be a whole number!\n", immediate. = TRUE)
      }
    }
  
  
  #######################################################
  #  ! I, Seed for random number generator
  #######################################################
  colonyfile$rseed = as.numeric(123)
  write(paste(colonyfile$rseed, "! I, Seed for random number generator"), name, append = TRUE)
  
  #######################################################
  #  ! B, 0/1=Not updating/updating allele frequency
  #######################################################
  colonyfile$updateallelefreq = 0
  write(paste(colonyfile$updateallelefreq, "! B, 0/1=Not updating/updating allele frequency"), name, append = TRUE)
  
  
  #######################################################
  #  ! 2/1=Dioecious/Monoecious
  #######################################################
  colonyfile$diomonoecy <- 2
  write(paste(colonyfile$diomonoecy, "! 2/1=Dioecious/Monoecious"), name, append = TRUE)
  
  
  #######################################################
  #  ! 0/1=No inbreeding/inbreeding
  #######################################################
  colonyfile$inbreeding <- 0
  write(paste(colonyfile$inbreeding, "! 0/1=No inbreeding/inbreeding"), name, append = TRUE)
  
  
  #######################################################
  #  ! B, 0/1=Diploid species/HaploDiploid species
  #######################################################
  colonyfile$ploidy <- 1
  write(paste(colonyfile$ploidy, "! B, 0/1=Diploid species/HaploDiploid species"), name, append = TRUE)
  
  #######################################################
  #  ! B, 0/1=Polygamy/Monogamy for males & females
  #######################################################
  colonyfile$malepolygamy <- 1
  colonyfile$femalepolygamy <- 1
  write(paste(colonyfile$malepolygamy, colonyfile$femalepolygamy, "! B, 0/1=Polygamy/Monogamy for males & females"), name, append = TRUE)
  
  #######################################################
  #  ! 0/1=No clone inference/clone inference
  #######################################################
  colonyfile$cloneinference <- 0
  write(paste(colonyfile$cloneinference, "! 0/1=No clone inference/clone inference"), name, append = TRUE)
  
  #######################################################
  #  ! 0/1=No scaling/sibship size scaling
  #######################################################
  colonyfile$sibscaling <- 0
  write(paste(colonyfile$sibscaling, "! 0/1=No scaling/sibship size scaling"), name, append = TRUE)
  
  #######################################################
  #  ! B, R, R : Use sibship prior, Y/N=1/0. If Yes, give mean paternal, maternal sibship size
  #######################################################
  colonyfile$sibship.prior <- 1
  
  if(colonyfile$sibship.prior==0){
    colonyfile$sibship.prior.paternal = 0
    colonyfile$sibship.prior.maternal = 0
    write(paste(colonyfile$sibship.prior, colonyfile$sibship.prior.paternal, colonyfile$sibship.prior.maternal, "! B, R, R : Use sibship prior, Y/N=1/0. If Yes, give mean paternal, maternal sibship size"), name, append = TRUE)
  }else{
    colonyfile$sibship.prior.paternal = as.numeric(1)
    colonyfile$sibship.prior.maternal = as.numeric(1)
    
    write(paste(colonyfile$sibship.prior, colonyfile$sibship.prior.paternal, colonyfile$sibship.prior.maternal, "! B, R, R : Use sibship prior, Y/N=1/0. If Yes, give mean paternal, maternal sibship size"), name, append = TRUE)
  }
  
  #######################################################
  #  ! B, 0/1=Unknown/Known population allele frequency
  #######################################################
  colonyfile$knownAFreq <- 0
  
  write(paste(colonyfile$knownAFreq, "! B, 0/1=Unknown/Known population allele frequency\n"), name, append = TRUE)
  
  #######################################################
  #  ! I, Number of runs
  #######################################################
  colonyfile$n.runs = as.numeric(1)
  write(paste(colonyfile$n.runs, "! I, Number of runs"), name, append = TRUE)
  
  if(length(colonyfile$n.runs) != 0){
    #Whole number warning
    if(is.whole(colonyfile$n.runs) == FALSE){
      utils::flush.console()
      colonyfile = colonyfile[which(names(colonyfile) != "n.runs")]
      warning("The number of runs must be a whole number!\n", immediate. = TRUE)
    }
  }
  
  
  #######################################################
  #  ! I, Length of Run (1, 2, 3) = (Short, Medium, Long)
  #######################################################
  colonyfile$runlength <- 3
  write(paste(colonyfile$runlength, "! I, Length of Run (1, 2, 3) = (Short, Medium, Long)"), name, append = TRUE)
  
  
  #######################################################
  #  ! B, 0/1=Monitor method by Iterate#/Time in second
  #######################################################
  colonyfile$monitortype <- 1
  write(paste(colonyfile$monitortype, "! B, 0/1=Monitor method by Iterate#/Time in second"), name, append = TRUE)
  
  #######################################################
  #  ! I, Monitor interval in Iterate#/Seconds
  #######################################################
  colonyfile$interval = as.numeric(60)
  write(paste(format(colonyfile$interval, scientific = FALSE), "! I, Monitor interval in Iterate#/Seconds"), name, append = TRUE)
  
  
  #######################################################
  #  ! B, 0/1=Other platform/Windows execution
  #######################################################
  colonyfile$sys <- 0
  write(paste(colonyfile$sys, "! B, 0/1=Other platform/Windows execution"), name, append = TRUE)
  
  #######################################################
  #  ! 1/0=Full-likelihood/pair-likelihood score method
  #######################################################
  colonyfile$likelihood.method <- 1
  write(paste(colonyfile$likelihood.method, "! 1/0=Full-likelihood/pairwise-likelihood score method"), name, append = TRUE)
  
  #######################################################
  #  ! 1/2/3=low/medium/high precision
  #######################################################
  colonyfile$precision <- 3
  write(paste(colonyfile$precision, "! 1/2/3=low/medium/high precision"), name, append = TRUE)
  write("\n", name, append = TRUE)
  
  #######################################################
  #Marker file import
  #######################################################
  
  #Give the path to the marker types and error rate file. This should be a file with a number of columns equal to the number of markers used.
  #There should be 4 rows, 1) marker ID, 2) marker type, 3) marker specific allelic dropout rate, 4) marker specific other typing error rate.
  
  colonyfile$MarkerPATH = error_rates_path
  colonyfile$Markers = utils::read.table(colonyfile$MarkerPATH, header = FALSE, colClasses = c("character"), sep = delim)
    
    if(colonyfile$n.loci != dim(colonyfile$Markers)[2]){
      colonyfile = colonyfile[which(names(colonyfile) != "MarkerPATH")]
      warning(paste("The number of defined loci ", "(",  colonyfile$n.loci, ") does not equal the number of markers provided in the file selected (", dim(colonyfile$Markers)[2], ").\n\n", sep = ""), immediate. = TRUE)
    }
  
  
  colonyfile$Markers[, 1 + dim(colonyfile$Markers)[2]] = c("!Marker IDs", "!Marker types, 0/1=Codominant/Dominant", "!Marker-specific allelic dropout rate", "!Other marker-specific typing-error rate")
  utils::write.table(colonyfile$Markers, name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write("\n", name, append = TRUE)
  
  #######################################################
  #Offspring genotype import
  #######################################################
  
  #Give the path to the offpring ID and genotype file
  #This should have a first column giving the ID, then 2 columns for each locus (1 for each allele at that locus), at least for diploid species.
  #Therefore, with 4 loci, there should be 9 columns.
  
  colonyfile$OSGenotypePATH = genotypes_path
  colonyfile$Offspring = utils::read.table(colonyfile$OSGenotypePATH, header = FALSE, colClasses = c("character"), sep = delim)
    
    if(colonyfile$n.offspring != dim(colonyfile$Offspring)[1]){
      colonyfile = colonyfile[which(names(colonyfile) != "OSGenotypePATH")]
      utils::flush.console()
      warning(paste("The number of defined offspring ", "(", colonyfile$n.offspring, ") does not equal the number of offspring provided in the file selected (", dim(colonyfile$Offspring)[1], ").\n\n", sep = ""), immediate. = TRUE)
    }
    
    fileloci = (dim(colonyfile$Offspring)[2] - 1) / 2
    
    if(colonyfile$ploidy == 0){
      if((colonyfile$n.loci) != fileloci){
        colonyfile = colonyfile[which(names(colonyfile) != "OSGenotypePATH")]
        utils::flush.console()
        warning(paste("The number of defined loci ", "(", colonyfile$n.loci, ") does not appear to equal the number of loci provided in the file selected (", fileloci, ").\n\n", sep = ""), immediate. = TRUE)
      }
    }
  
  
  colonyfile$Offspring[, 1 + dim(colonyfile$Offspring)[2]] = c("!Offspring ID and genotypes", rep("", dim(colonyfile$Offspring)[1] - 1))
  utils::write.table(colonyfile$Offspring, name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write("", name, append = TRUE)
  
  ######################################################
  #Sampling of candidate parents
  #######################################################
  #FATHERS - probability of inclusion in candidate set
  #######################################################
  colonyfile$fatherprob = as.numeric(0)
  
  if(colonyfile$fatherprob > 1){
    utils::flush.console()
    cat("Probabilities must be less than or equal to 1.\n")
    colonyfile = colonyfile[which(names(colonyfile) != "fatherprob")]
  }
  
  
  #######################################################
  #FATHERS - number of candidate fathers
  #######################################################
  colonyfile$n.father = as.numeric(0)
  
  #######################################################
  #FATHERS - Import candidate FATHERS file
  #######################################################
  
  colonyfile$fathersPATH = NA
  colonyfile$fathers = matrix(nrow = 1, ncol = 1)
  
  
  #######################################################
  #MOTHERS - probability of inclusion in candidate set
  #######################################################
  colonyfile$motherprob = as.numeric(0)
  
  #######################################################
  #MOTHERS - Number of candidate mothers
  #######################################################
  colonyfile$n.mother = as.numeric(0)
  
  #######################################################
  #MOTHERS - Import candidate MOTHERS
  #######################################################
  colonyfile$mothersPATH = NA
  colonyfile$mothers = matrix(nrow=1, ncol = 1)
  
  
  write(paste(colonyfile$fatherprob, colonyfile$motherprob, "!Probabilities that the father and mother of an offspring included in candidates"), name, append = TRUE)
  
  write(paste(colonyfile$n.father, colonyfile$n.mother, "!Numbers of candidate males and females"), name, append = TRUE)
  
  write("", name, append = TRUE)
  
  if(colonyfile$n.father != 0){
    colonyfile$fathers[, 1 + dim(colonyfile$fathers)[2]] = c("!Candidate M ID and genotypes", rep("", dim(colonyfile$fathers)[1] - 1))
    utils::write.table(colonyfile$fathers, name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  write("", name, append = TRUE)
  
  if(colonyfile$n.mother != 0){
    colonyfile$mothers[, 1 + dim(colonyfile$mothers)[2]] = c("!Candidate F ID and genotypes", rep("", dim(colonyfile$mothers)[1] - 1))
    utils::write.table(colonyfile$mothers, name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  write("", name, append = TRUE)
  
  #######################################################
  #Define known PATERNAL dyads and sibships
  #######################################################
  colonyfile$n.known.paternities.and.sibships = as.numeric(0)
  colonyfile$paternal.dyads = NA
  colonyfile$paternal.sibships = NA
  
  
  
  #######################################################
  #Define known MATERNAL dyads and sibships
  #######################################################
  
  colonyfile$n.known.maternities.and.sibships = as.numeric(0)
  colonyfile$maternal.dyads = NA
  colonyfile$maternal.sibships = NA
  
  #Paternal Dyads
  if(is.na(colonyfile$paternal.dyads)){
    utils::write.table("0 0 !Number of known paternities", name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    utils::write.table(paste(dim(colonyfile$paternal.dyads)[1], "!Number of known paternities"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    colonyfile$paternal.dyads = colonyfile$paternal.dyads[2:1]
    colonyfile$paternal.dyads[, 1+dim(colonyfile$paternal.dyads)[2]] = c("!IDs of known offspring-father dyad", rep("", dim(colonyfile$paternal.dyads)[1] - 1))
    utils::write.table(colonyfile$paternal.dyads, name, append = TRUE, quote = FALSE, na = " ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  #Maternal Dyads
  if(is.na(colonyfile$maternal.dyads)){
    utils::write.table("0 0 !Number of known maternities", name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    utils::write.table(paste(dim(colonyfile$maternal.dyads)[1], "!Number of known maternities"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    colonyfile$maternal.dyads = colonyfile$maternal.dyads[2:1]
    colonyfile$maternal.dyads[, 1 + dim(colonyfile$maternal.dyads)[2]] = c("!IDs of known offspring-mother dyad", rep("", dim(colonyfile$maternal.dyads)[1] - 1))
    utils::write.table(colonyfile$maternal.dyads, name, append = TRUE, quote = FALSE, na = " ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  #Paternal sibships
  if(is.na(colonyfile$paternal.sibships)){
    utils::write.table("0 !Number of known paternal sibships with unknown fathers", name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    utils::write.table(paste(dim(colonyfile$paternal.sibships)[1], "!Number of known paternal sibships with unknown fathers "), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.table(colonyfile$paternal.sibships, name, append = TRUE, quote = FALSE, na = " ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  #Maternal sibships
  if(is.na(colonyfile$maternal.sibships)){
    utils::write.table("0  !Number of known maternal sibships with unknown mothers", name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    utils::write.table(paste(dim(colonyfile$maternal.sibships)[1], "!Number of known maternal sibships with unknown mothers "), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.table(colonyfile$maternal.sibships, name, append = TRUE, quote = FALSE, na=" ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  
  
  #######################################################
  #Define excluded PATERNITIES
  #######################################################
  
  colonyfile$n.excluded.paternities = as.numeric(0)
    
    if(length(colonyfile$n.excluded.paternities) != 0){
      #Whole number warning
      if(is.whole(colonyfile$n.excluded.paternities) == FALSE){
        utils::flush.console()
        colonyfile = colonyfile[which(names(colonyfile) != "n.excluded.paternities")]
        warning("The number of excluded paternities must be a whole number!\n", immediate. = TRUE)
      }
    }
  
  
  if(colonyfile$n.excluded.paternities > 0){
    
    
    #Get the path, and delimiter, to the file...
    while(length(colonyfile$excluded.paternities.PATH) == 0){
      cat("Provide the path to the excluded PATERNITY file.\n\n\n")
      utils::flush.console()
      colonyfile$excluded.paternities.PATH = file.choose()
      
      
      #Read in the data...
      colonyfile$excluded.paternities = utils::read.table(colonyfile$excluded.paternities.PATH, header = FALSE,
                                                          sep=delim, colClasses=c("character"), fill = TRUE, flush = TRUE, na.strings="")
      
      
      if(colonyfile$n.excluded.paternities > 0){#Write ExcludedPaternity.txt file
        temp1 = as.data.frame(colonyfile$excluded.paternities)
        names(temp1) = c("OffspringID", paste("ExcludedFatherID", 1:(dim(temp1)[2] - 1), sep=""))
        utils::write.table(temp1, "ExcludedPaternity.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
      }
      
      #Check the data
      if(colonyfile$n.excluded.paternities != dim(colonyfile$excluded.paternities)[1]){
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.paternities.PATH")]
        utils::flush.console()
        warning(paste("The number of defined excluded paternities ", "(", colonyfile$n.excluded.paternities, ") does not equal the number provided in the file selected (", dim(colonyfile$excluded.paternities)[1], ").\n\n", sep=""), immediate. = TRUE)
      }
      
      
      #Further checks
      #if this is true, then all offspring in the dyad file are present in the offspring genotype file
      if(sum(colonyfile$excluded.paternities$V1 %in% colonyfile$Offspring[, 1]) == length(colonyfile$excluded.paternities$V1)){
        #Do nothing.
      }else{
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.paternities.PATH")]
        utils::flush.console()
        warning(paste("Offspring in excluded paternities file are not present in the offspring genotype data:", paste(colonyfile$excluded.paternities$V1[which(colonyfile$excluded.paternities$V1%in%colonyfile$fathers[, 1] == FALSE)], collapse=", ")), immediate. = TRUE)
      }
      
      os = stats::na.omit(as.vector(as.matrix(colonyfile$excluded.paternities[, 2:dim(colonyfile$excluded.paternities)[2]])))
      
      if(sum(os%in%colonyfile$fathers[, 1]) == length(os)){
        #Do nothing
      }else{
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.paternities.PATH")]
        utils::flush.console()
        warning(paste("Fathers in excluded paternities file are not present in the fathers genotype data:", paste(os[which(os%in%colonyfile$fathers[, 1] == FALSE)], collapse=", ")), immediate. = TRUE)
      }
    }
    
    csum = NULL
    for (i in 1:dim(colonyfile$excluded.paternities)[1]){
      csum[i] = length(colonyfile$excluded.paternities[i, ][!is.na(colonyfile$excluded.paternities[i, ])])
    }
    
    csum = csum - 1
    
    utils::write.table(paste(colonyfile$n.excluded.paternities, "!Number of offspring with known excluded paternity"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    colonyfile$excluded.paternities = cbind(as.character(colonyfile$excluded.paternities[, 1]), csum, colonyfile$excluded.paternities[, 2:dim(colonyfile$excluded.paternities)[2]])
    
    colonyfile$excluded.paternities[, 1 + dim(colonyfile$excluded.paternities)[2]] = c("!Offspring ID, number of excluded males, the IDs of excluded males", rep("", dim(colonyfile$excluded.paternities)[1] - 1))
    
    utils::write.table(colonyfile$excluded.paternities, name, append = TRUE, quote = FALSE, na=" ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    #If there are no excluded paternities
    utils::write.table(paste(colonyfile$n.excluded.paternities, " !Number of offspring with known excluded paternity"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  #######################################################
  #Define excluded MATERNITIES
  #######################################################
  
  colonyfile$n.excluded.maternities = as.numeric(0)
    
    if(length(colonyfile$n.excluded.maternities) != 0){
      #Whole number warning
      if(is.whole(colonyfile$n.excluded.maternities) == FALSE){
        utils::flush.console()
        colonyfile = colonyfile[which(names(colonyfile) != "n.excluded.maternities")]
        warning("The number of excluded maternities must be a whole number!\n", immediate. = TRUE)
      }
    }
  
  
  if(colonyfile$n.excluded.maternities > 0){
    #Get the path, and delimiter, to the file...
    while(length(colonyfile$excluded.maternities.PATH) == 0){
      cat("Provide the path to the excluded MATERNITY file.\n\n\n")
      utils::flush.console()
      colonyfile$excluded.maternities.PATH = file.choose()
      
      #Read in the data...
      colonyfile$excluded.maternities = utils::read.table(colonyfile$excluded.maternities.PATH, header = FALSE,
                                                          sep=delim, colClasses=c("character"), fill = TRUE, flush = TRUE, na.strings="")
      
      #Check the data
      if(colonyfile$n.excluded.maternities != dim(colonyfile$excluded.maternities)[1]){
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.maternities.PATH")]
        utils::flush.console()
        warning(paste("The number of defined excluded maternities ", "(", colonyfile$n.excluded.maternities, ") does not equal the number provided in the file selected (", dim(colonyfile$excluded.maternities)[1], ").\n\n", sep=""), immediate. = TRUE)
      }
      
      if(colonyfile$n.excluded.maternities > 0){#Write ExcludedMaternity.txt file
        temp1 = as.data.frame(colonyfile$excluded.maternities)
        names(temp1) = c("OffspringID", paste("ExcludedMotherID", 1:(dim(temp1)[2] - 1), sep = ""))
        
        utils::write.table(temp1, "ExcludedMaternity.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
      }
      
      #Futher checks
      #if this is true, then all offspring in the dyad file are present in the offspring genotype file
      if(sum(colonyfile$excluded.maternities$V1%in%colonyfile$Offspring[, 1]) == length(colonyfile$excluded.maternities$V1)){
        #Do nothing.
      }else{
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.maternities.PATH")]
        utils::flush.console()
        warning(paste("Offspring in excluded maternities file are not present in the offspring genotype data:", paste(colonyfile$excluded.maternities$V1[which(colonyfile$excluded.maternities$V1%in%colonyfile$mothers[, 1] == FALSE)], collapse=", ")), immediate. = TRUE)
      }
      
      os = stats::na.omit(as.vector(as.matrix(colonyfile$excluded.maternities[, 2:dim(colonyfile$excluded.maternities)[2]])))
      
      if(sum(os%in%colonyfile$mothers[, 1]) == length(os)){
        #Do nothing
      }else{
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.maternities.PATH")]
        utils::flush.console()
        warning(paste("Mothers in excluded maternities file are not present in the mothers genotype data:", paste(os[which(os%in%colonyfile$mothers[, 1] == FALSE)], collapse=", ")), immediate. = TRUE)
      }
    }
    
    csum = NULL
    for (i in 1:dim(colonyfile$excluded.maternities)[1]){
      csum[i] = length(colonyfile$excluded.maternities[i, ][!is.na(colonyfile$excluded.maternities[i, ])])
    }
    
    csum = csum - 1
    
    utils::write.table(paste(colonyfile$n.excluded.maternities, "!Number of offspring with known excluded maternity"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    colonyfile$excluded.maternities = cbind(as.character(colonyfile$excluded.maternities[, 1]), csum, colonyfile$excluded.maternities[, 2:dim(colonyfile$excluded.maternities)[2]])
    
    colonyfile$excluded.maternities[, 1 + dim(colonyfile$excluded.maternities)[2]] = c("!Offspring ID, number of excluded females, the IDs of excluded females", rep("", dim(colonyfile$excluded.maternities)[1] - 1))
    
    utils::write.table(colonyfile$excluded.maternities, name, append = TRUE, quote = FALSE, na=" ", row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }else{
    #If there are no excluded maternities
    utils::write.table(paste(colonyfile$n.excluded.maternities, " !Number of offspring with known excluded maternity"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  
  #######################################################
  #Define EXCLUDED PATERNAL sibships
  #######################################################
  
  colonyfile$excluded.paternal.sibships.PATH = exclusion_path
  
  if(!is.null(colonyfile$excluded.paternal.sibships.PATH)){
    
      
      #Read in the data...
      colonyfile$excluded.paternal.sibships = utils::read.table(colonyfile$excluded.paternal.sibships.PATH, header = FALSE,
                                                                sep=delim, colClasses=c("character"), fill = TRUE, flush = TRUE, na.strings="")
      
      colonyfile$n.excluded.paternal.sibships = dim(colonyfile$excluded.paternal.sibships)[1]
      
      #Further checks - do excluded sibs appear in offspring file
      os = stats::na.omit(as.vector(as.matrix(colonyfile$excluded.paternal.sibships[, 2:dim(colonyfile$excluded.paternal.sibships)[2]])))
      
      if(!sum(os%in%colonyfile$Offspring[, 1]) == length(os)){
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.paternal.sibships.PATH")]
        utils::flush.console()
        warning(paste("Offspring in excluded paternal sibships file are not present in the offspring genotype data:", paste(os[which(os%in%colonyfile$Offspring[, 1] == FALSE)], collapse=", ")), immediate. = TRUE)
      }
    
    
    colonyfile$excluded.paternal.sibships[, 1 + dim(colonyfile$excluded.paternal.sibships)[2]] = c(rep("", dim(colonyfile$excluded.paternal.sibships)[1]))
    csum = NULL
    for (i in 1:dim(colonyfile$excluded.paternal.sibships)[1]){
      csum[i] = length(colonyfile$excluded.paternal.sibships[i, ][!is.na(colonyfile$excluded.paternal.sibships[i, ])])
    }
    csum = csum - 2
    
    colonyfile$excluded.paternal.sibships = cbind(colonyfile$excluded.paternal.sibships[, 1], csum, colonyfile$excluded.paternal.sibships[, 2:ncol(colonyfile$excluded.paternal.sibships)])
    
    utils::write.table(paste(colonyfile$n.excluded.paternal.sibships, "!Number of offspring with known excluded paternal sibships"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    #utils::write.table(colonyfile$excluded.paternal.sibships, name, append = TRUE, quote = FALSE, na=" ", row.names = FALSE, col.names = FALSE)
    #write("", name, append = TRUE)
    
    lines = apply(colonyfile$excluded.paternal.sibships, 1, function(row) {
      # remove NAs
      vals = row[!is.na(row)]
      paste(vals, collapse = " ")
    })
    con <- file(name, open = "a")  # 'a' = append mode
    writeLines(unlist(lines), con)
    close(con)
    
    
    
  }else{
    #If there are no excluded sibships
    utils::write.table(paste(colonyfile$n.excluded.paternal.sibships, " !Number of offspring with known excluded paternal sibships"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  
  #######################################################
  #Define EXCLUDED MATERNAL sibships
  ######################################################
  
  colonyfile$excluded.maternal.sibships.PATH = exclusion_path
  
  if(!is.null(colonyfile$excluded.maternal.sibships.PATH)){
    
    #Get the path, and delimiter, to the file...
    
      #Read in the data...
      colonyfile$excluded.maternal.sibships = utils::read.table(colonyfile$excluded.maternal.sibships.PATH, header = FALSE,
                                                                sep=delim, colClasses = c("character"), fill = TRUE, flush = TRUE, na.strings="")
      colonyfile$n.excluded.maternal.sibships = dim(colonyfile$excluded.maternal.sibships)[1]
      
      
      #Further checks - do excluded sibs appear in offspring file
      os = stats::na.omit(as.vector(as.matrix(colonyfile$excluded.maternal.sibships[, 2:dim(colonyfile$excluded.maternal.sibships)[2]])))
      
      if(!sum(os%in%colonyfile$Offspring[, 1]) == length(os)){
        colonyfile = colonyfile[which(names(colonyfile) != "excluded.maternal.sibships.PATH")]
        utils::flush.console()
        warning(paste("Offspring in excluded maternal sibships file are not present in the offspring genotype data:", paste(os[which(os%in%colonyfile$Offspring[, 1] ==  FALSE)], collapse=", ")), immediate. = TRUE)
      }
    
    
    colonyfile$excluded.maternal.sibships[, 1 + dim(colonyfile$excluded.maternal.sibships)[2]] = c(rep("", dim(colonyfile$excluded.maternal.sibships)[1]))
    csum = NULL
    for (i in 1:dim(colonyfile$excluded.maternal.sibships)[1]){
      csum[i] = length(colonyfile$excluded.maternal.sibships[i, ][!is.na(colonyfile$excluded.maternal.sibships[i, ])])
    }
    csum = csum - 2
    
    colonyfile$excluded.maternal.sibships = cbind(colonyfile$excluded.maternal.sibships[, 1], csum, colonyfile$excluded.maternal.sibships[, 2:ncol(colonyfile$excluded.maternal.sibships)])
    
    utils::write.table(paste(colonyfile$n.excluded.maternal.sibships, "!Number of offspring with known excluded maternal sibships"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    lines = apply(colonyfile$excluded.maternal.sibships, 1, function(row) {
      # remove NAs
      vals = row[!is.na(row)]
      paste(vals, collapse = " ")
    })
    con <- file(name, open = "a")  # 'a' = append mode
    writeLines(unlist(lines), con)
    close(con)
  }else{
    #If there are no excluded sibships
    utils::write.table(paste(colonyfile$n.excluded.maternal.sibships, " !Number of offspring with known excluded maternal sibships"), name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write("", name, append = TRUE)
  }
  
  
  #######################################################
  #Other outputs
  #######################################################
  
  #MarkerTypeErrorRate.txt
  temp1 = as.data.frame(colonyfile$Markers)
  temp1 = temp1[, 1:dim(temp1)[2] - 1]
  names(temp1) = paste("Locus-", 1:dim(temp1)[2], sep="")
  utils::write.table(temp1, "MarkerTypeErrorRate.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  #AlleleFrequency.txt
  if(colonyfile$knownAFreq == 1){
    temp1 = as.data.frame(colonyfile$allele.frequency)
    temp1 = temp1[, 1:dim(temp1)[2] - 1]
    names(temp1) = paste(paste("Locus-", rep(1:(dim(temp1)[2] / 2), each = 2), sep = ""), rep(1:2, (dim(temp1)[2]/2)), sep=".")
    utils::write.table(temp1, "MarkerTypeErrorRate.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
  }
  
  #OffspringGenotype.txt
  temp1 = as.data.frame(colonyfile$Offspring)
  temp1 = temp1[, 1:dim(temp1)[2] - 1]
  n = (dim(temp1)[2])-1
  names(temp1) = c("Offspring", paste(paste("Marker", rep(1:(n / 2), each = 2), sep=""), rep(1:2, (n/2)), sep = "-"))
  utils::write.table(temp1, "OffspringGenotype.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  #MaleGenotype.txt
  if(colonyfile$n.father != 0){
    temp1 = as.data.frame(colonyfile$fathers)
    temp1 = temp1[, 1:dim(temp1)[2]-1]
    n = (dim(temp1)[2])-1
    names(temp1) = c("Male", paste(paste("Marker", rep(1:(n / 2), each = 2), sep=""), rep(1:2, (n / 2)), sep = "-"))
    utils::write.table(temp1, "MaleGenotype.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
  }
  
  if(colonyfile$n.mother != 0){
    #FemaleGenotype.txt
    temp1 = as.data.frame(colonyfile$mothers)
    temp1 = temp1[, 1:dim(temp1)[2] - 1]
    n = (dim(temp1)[2]) - 1
    names(temp1) = c("Female", paste(paste("Marker", rep(1:(n / 2), each=2), sep = ""), rep(1:2, (n / 2)), sep = "-"))
    utils::write.table(temp1, "FemaleGenotype.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)
  }
  
  #KnownPaternalDyads.txt
  if(!is.na(colonyfile$paternal.dyads)){
    temp1 = colonyfile$paternal.dyads[c(2, 1)]
    names(temp1) = c("OffspringID", "Father")
    utils::write.table(temp1, "KnownPaternalDyads.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)}
  
  #KnownMaternalDyads.txt
  if(!is.na(colonyfile$maternal.dyads)){
    temp1 = colonyfile$maternal.dyads[c(2, 1)]
    names(temp1) = c("OffspringID", "Mother")
    utils::write.table(temp1, "KnownMaternalDyads.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)}
  
  
  #KnownPaternity.txt - see earlier
  #KnownMaternity.txt - see earlier
  
  #ExcludedPaternity.txt - see earlier
  #ExcludedMaternity.txt - see earlier
  
  #Produce summary information text file.
  utils::write.table(paste("Output file path & name : ", wd, name, "\n",
                           "Number of loci : ", colonyfile$n.loci, "\n",
                           "Number of offspring in the sample : ", colonyfile$n.offspring, "\n",
                           "Number of male candidates : ", colonyfile$n.father, "\n",
                           "Number of female candidates : ", colonyfile$n.mother, "\n",
                           "Number of known paternal sibships : ", colonyfile$n.paternal.sibs.or.paternities, "\n",
                           "Number of known maternal sibships : ", colonyfile$n.maternal.sibs.or.maternities, "\n",
                           "Number of offspring with excluded fathers : ", colonyfile$n.excluded.paternities, "\n",
                           "Number of offspring with excluded mothers : ", colonyfile$n.excluded.maternities, "\n",
                           "Male mating system : ", if(colonyfile$malepolygamy == 0){"Polygamous"}else{"Monogamous"}, "\n",
                           "Female mating system : ", if(colonyfile$femalepolygamy == 0){"Polygamous"}else{"Monogamous"}, "\n",
                           "Dioecious/Monoecious : ", if(colonyfile$diomonoecy == 1){"Monoecious"}else{"Dioecious"}, "\n",
                           "Number of threads : ", "nA", "\n",
                           "Number of Excluded Paternal Sibships : ", colonyfile$n.excluded.paternal.sibships, "\n",
                           "Number of Excluded Maternal Sibships : ", colonyfile$n.excluded.paternal.sibships, "\n",
                           "Seed for random number generator : ", colonyfile$rseed, "\n",
                           "Allele frequency : ", if(colonyfile$updateallelefreq == 1){"Updating by accounting for the inferred relationship"}else{"No updating by accounting for the inferred relationship"}, "\n",
                           "Species : ", if(colonyfile$ploidy == 1){"HaploDiploid"}else{"Diploid"}, "\n",
                           "Known population allele frequency : ", if(colonyfile$knownAFreq == 1){"Yes"}else{"No"}, "\n",
                           "Number of run : ", colonyfile$n.runs, "\n",
                           "Length of run : ", if(colonyfile$runlength == 1){"Small"}else{if(colonyfile$runlength == 2){"Medium"}else{"Long"}}, "\n",
                           "Monitor intermiediate results by : ", if(colonyfile$monitortype == 1){paste("Every", colonyfile$interval, "seconds")}else{paste("Every", colonyfile$interval, "iterations")}, "\n",
                           "Project data input produced : ", date(), "\n",
                           "NOTE to the Project: ", colonyfile$note, "\n", sep=""), "ProjectInformation.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  cat("Finished!")
  cat(paste("Your file is called", name, "and is placed in", wd, "...\n\n\n"))
  
  #This could be useful at some point.
  #return(colonyfile)
}
