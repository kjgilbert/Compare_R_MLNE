

# want to get the moment estimate of migration from the mlne runs I did where I know mig and the starting p and q
# see if the error variances as calculated from the known quantities (because I gave it known ps's and q's) differes significantly from the variance in the estimate of migration



# estimate migration from the mlne moment method:
orig.mlne.moment <- function(input){
		B_T <- NULL
		B_0 <- NULL
		A_0 <- NULL

	dat <- readLines(input)
	focal_t0 <- dat[11]
	focal_current <- dat[13]
	source_t0 <- dat[16]
	num.loci <- 1
	for(i in 1:num.loci){
		temp.focal_t0 <- as.numeric(strsplit(focal_t0[i], split=" ")[[1]]) / sum(as.numeric(strsplit(focal_t0[i], split=" ")[[1]]))
		temp.focal_t1 <- as.numeric(strsplit(focal_current[i], split=" ")[[1]]) / sum(as.numeric(strsplit(focal_current[i], split=" ")[[1]]))
		temp.source_t0 <- as.numeric(strsplit(source_t0[i], split=" ")[[1]]) / sum(as.numeric(strsplit(source_t0[i], split=" ")[[1]]))
		B_T <- c(B_T, temp.focal_t1)
		B_0 <- c(B_0, temp.focal_t0)
		A_0 <- c(A_0, temp.source_t0)
	}
	x_ab <- A_0 - B_0
	x_bt <- B_T - B_0

	w_l <- (x_ab^2)  #*(p_b0*q_b0 * (((1*S_bt))) )
	F_l <- x_bt/x_ab
	W <- sum(w_l)

	return(mig.rate <- (1/W)*sum(w_l*F_l, na.rm=TRUE))	# Fortran eqn: MTM2=1.-ABS(1.-WF/WW)**(1./nt)
}






read.nemo.write.mlne <- function (currentdata, ancientdata, ninds, nsub, maxall=256, nsamps=2, mlne.gens="1", maxNe=1000, IBD=FALSE, na.s = c("0", "00", "000", "0000", "00000", "000000", "NA")) 
{
    x <- scan(currentdata, n = 4) #read first line of file, gives number of loci, alleles, etc
    #first line tells me:  13 patches, number of loci + 4, number of alleles, digits per genotype
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames <- scan(currentdata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames <- c("Pop", lnames)  #add first column heading name to be "pop"
    
    dat <- scan(currentdata, skip = nloc + 5, what = character(), na.strings = na.s)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    dat <- data.frame(matrix(dat, ncol = nloc + 5, byrow = TRUE))  #genotype data in matrix format now
    dat <- dat[,-((length(dat[1,])-3):length(dat[1,]))]  #remove last 4 columns of nemo data
    names(dat) <- lnames  #add column names

#split genotypes into alleles of 3 digits
    for(i in 2:(nloc+1)){  #go from first column of genotype data to the last column
    	if(i %% 100==0) print(i)	
    	ans <- matrix(NA, ncol=2, nrow=length(dat[,i]))
		for(j in 1:length(dat[,i])){
			test <- strsplit(as.character(dat[,i]), "")[[j]]
			ans[j,1] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[1]
			ans[j,2] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[2]
			}  #end for loop going through individual genotype column, splitting genotypes
		dat[,i] <- ans	
    	}  #end for loop going through genotype columns
    
### REPEAT FOR SECOND DATASET
x <- scan(ancientdata, n = 4) #read first line of file, gives number of loci, alleles, etc
    nloc <- x[2]-4  #get number of loci (-4 because nemo adds 4 extra rows where loci are)
    lnames2 <- scan(ancientdata, what = character(), skip = 1, nlines = nloc)  #locus names
    lnames2 <- c("Pop", lnames2)  #add first column heading name to be "pop"
    dat2 <- scan(ancientdata, skip = nloc + 5, what = character(), na.strings = na.s)  #read genotype data
    	#adding what=character reads data as characters instead of numbers and keeps the leading zeroes
    dat2 <- data.frame(matrix(dat2, ncol = nloc + 5, byrow = TRUE))  #genotype data in matrix format now
    dat2 <- dat2[,-((length(dat2[1,])-3):length(dat2[1,]))]  #remove last 4 columns of nemo data
    names(dat2) <- lnames2  #add column names
    
#split genotypes into alleles of 3 digits
    for(i in 2:(nloc+1)){  #go from first column of genotype data to the last column
    	if(i %% 100==0) print(i)	
    	ans <- matrix(NA, ncol=2, nrow=length(dat2[,i]))
		for(j in 1:length(dat2[,i])){
			test <- strsplit(as.character(dat2[,i]), "")[[j]]
			ans[j,1] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[1]
			ans[j,2] <- paste0(test[c(1,4)], test[c(2,5)], test[c(3,6)])[2]
			}  #end for loop going through individual genotype column, splitting genotypes
		dat2[,i] <- ans
    	}  #end for loop going through genotype columns
    	
finalcurrentdat_presubset <- dat
finalancientdat_presubset <- dat2
dat <- NULL
dat2 <- NULL

finalcurrentdat_1 <- finalcurrentdat_presubset[finalcurrentdat_presubset$Pop=="1",]
finalancientdat_1 <- finalancientdat_presubset[finalancientdat_presubset$Pop=="1",]

finalcurrentdat_2 <- finalcurrentdat_presubset[finalcurrentdat_presubset$Pop=="2",]
finalancientdat_2 <- finalancientdat_presubset[finalancientdat_presubset$Pop=="2",]

number.of.loci <- nloc

## ##  Nf = 1       source = 2		

## PUT IT IN THE PROPER ORDER SO THAT THE INPUT FILE COMES OUT WITH focal t0 followed by focal t1 followed by source t0 datasets
#	I did this at the end when the datasets are combined and printed to the output

mlne.currdat <- finalcurrentdat_1
mlne.ancdat <- finalancientdat_1
mlne.srcdat <- finalancientdat_2

#get allele counts for current dataset
	mfulloutcurr <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqscurr <- as.data.frame(table(mlne.currdat[,i]))
		moutputcurr <- cbind(i, names(mlne.currdat[i]), mfreqscurr)
		mfulloutcurr <- rbind(mfulloutcurr, moutputcurr)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutcurr) <- c("Number", "Locus", "Allele", "Count")
## REPEAT FOR SECOND TIME POINT
	mfulloutanc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqsanc <- as.data.frame(table(mlne.ancdat[,i]))
		moutputanc <- cbind(i, names(mlne.ancdat[i]), mfreqsanc)
		mfulloutanc <- rbind(mfulloutanc, moutputanc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutanc) <- c("Number", "Locus", "Allele", "Count")	
## REPEAT FOR SOURCE POP
	mfulloutsrc <- NULL
	for(i in 2:(number.of.loci+1)){
		mfreqssrc <- as.data.frame(table(mlne.srcdat[,i]))
		moutputsrc <- cbind(i, names(mlne.srcdat[i]), mfreqssrc)
		mfulloutsrc <- rbind(mfulloutsrc, moutputsrc)
	}
#now have a table of alleles and frequencies per locus
	colnames(mfulloutsrc) <- c("Number", "Locus", "Allele", "Count")
#now need to put zeros where some alleles are not there
mlne.fulllist <- NULL
allele.list <- NULL
allelescurr <- list()
allelesanc <- list()
allelessrc <- list()
#for(m in 2:(number.of.loci+1)){
m <- 2
	mtempallelescurr <- subset(mfulloutcurr, mfulloutcurr$Number==m)
	mtempallelesanc <- subset(mfulloutanc, mfulloutanc$Number==m)
	mtempallelessrc <- subset(mfulloutsrc, mfulloutsrc$Number==m)
	mpotalleles <- as.matrix(cbind(seq(0, maxall), rep(0, maxall+1)))
	mfpotallelescurr <- as.matrix(cbind(as.factor(mpotalleles[,1]), mpotalleles[,2]))
	mfpotallelesanc <- mfpotallelescurr
	mfpotallelessrc <- mfpotallelescurr
	#time point 1, put in allele counts
	for(k in 1:length(mfpotallelescurr[,1])){
		for(i in 1:length(mtempallelescurr$Allele)){
			if(mfpotallelescurr[k,1]==as.numeric(as.character(mtempallelescurr$Allele[i]))){mfpotallelescurr[k,2] <- mtempallelescurr$Count[i]}
	}}
	#time point 2, put in allele counts
	for(k in 1:length(mfpotallelesanc[,1])){
		for(i in 1:length(mtempallelesanc$Allele)){
			if(mfpotallelesanc[k,1]==as.numeric(as.character(mtempallelesanc$Allele[i]))){mfpotallelesanc[k,2] <- mtempallelesanc$Count[i]}
	}}
	#time point 2 source pop, put in allele counts
	for(k in 1:length(mfpotallelessrc[,1])){
		for(i in 1:length(mtempallelessrc$Allele)){
			if(mfpotallelessrc[k,1]==as.numeric(as.character(mtempallelessrc$Allele[i]))){mfpotallelessrc[k,2] <- mtempallelessrc$Count[i]}
	}}
	#delete rows where both time points have zeros
	bb <- NULL
	for(zz in 1:length(mfpotallelescurr[,1])){
		if(mfpotallelescurr[zz,2]==0 && mfpotallelesanc[zz,2]==0 && mfpotallelessrc[zz,2]==0){
			aa <- zz
			bb <- c(bb,aa)
	}}
	if(is.null(bb)==FALSE){
		mfinalalleles1 <- mfpotallelescurr[-bb,]
		mfinalalleles2 <- mfpotallelesanc[-bb,]
		mfinalalleles3 <- mfpotallelessrc[-bb,]
	}else{mfinalalleles1 <- mfpotallelescurr ; mfinalalleles2 <- mfpotallelesanc ; mfinalalleles3 <- mfpotallelessrc}
tempalls1 <- paste(if(length(mfinalalleles1)>2){paste(mfinalalleles1[,2])}else{paste(mfinalalleles1[2])}, collapse=" ")
tempalls2 <- paste(if(length(mfinalalleles2)>2){paste(mfinalalleles2[,2])}else{paste(mfinalalleles2[2])}, collapse=" ")
tempalls3 <- paste(if(length(mfinalalleles3)>2){paste(mfinalalleles3[,2])}else{paste(mfinalalleles3[2])}, collapse=" ")
allelescurr <- c(allelescurr, list(tempalls1))
allelesanc <- c(allelesanc, list(tempalls2))
allelessrc <- c(allelessrc, list(tempalls3))
allcounts <- paste(if(length(mfinalalleles1)>2){paste(length(mfinalalleles1[,2]))}else{paste(length(mfinalalleles1[2]))}, collapse=" ")
allele.list <- c(allele.list, allcounts)
#}  #end for loop

#now there are 4 lists, allele.list with each locus's allele count, allelescurr with each element of the list being an allele count for a locus at the current time sample, allelesanc for the ancient sample allele counts per locus, and allelessrc for the ancient source pop sample allele counts per locus
#so loop through and paste out each element of the lists together
for(f in 1:length(allelescurr)){
		mtemplist <- paste(c(allelescurr[f], ""))
		mlne.fulllist <- c(mlne.fulllist, mtemplist)
	}
mlne.fulllist <- paste(c(unlist(allelesanc),"", unlist(allelescurr),"","", unlist(allelessrc)))
between <- paste(c("0", mlne.gens), collapse=",")

mlne_input <- paste(c(paste("1"), paste("0"), paste(maxNe), paste("2"), paste("0"), paste(number.of.loci), paste(allele.list, collapse=","), paste(nsamps), between, "", mlne.fulllist, "", paste("1"), sep="\n"))  #add 0, number of time points, and number of loci at start of file

	writeLines(mlne_input, paste("~/Desktop/NewMLNeSims/Compare_R_MLNE/NemoSims_TestVariance/MLNeInputs/MLNeIN_", tail(unlist(strsplit(currentdata, "/")),n=1), sep=""))
} 



