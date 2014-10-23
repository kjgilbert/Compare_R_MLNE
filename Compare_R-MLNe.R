#### WITH THE ACTUAL MLNE INPUT FILES

input <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/slingshot_Oct8_MLNeIN_1.4_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
(ans <- mlne.moment(input)) 

# equation 8: m = 1-t_root(1 - (1/W) * sum over all l of w_l*F_l)
#	equals for t=1: (1/W)*sum(w_l*F_l)


mlne.moment <- function(input){
		B_T <- NULL
		B_0 <- NULL
		A_0 <- NULL

	dat <- readLines(input)
	focal_t0 <- dat[11:50]
	focal_current <- dat[52:91]
	source_t0 <- dat[94:133]
	#max.loci <- dat[7]	#don't want these anymore
	#max.loci.list <- strsplit(max.loci, split=",")
	num.loci <- 40

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

	w_l <- x_ab^2
	F_l <- x_bt/x_ab
	W <- sum(w_l)

	return(mig.rate <- (1/W)*sum(w_l*F_l, na.rm=TRUE))	# Fortran eqn: MTM2=1.-ABS(1.-WF/WW)**(1./nt)
}


### I had time points mixed up in the input, now R matches MLNe output, but the mlne outputs are now originally wrong, since I gave the input in the wrong order, so the TRUE ANSWER SHOULD BE THE FIRST VALUES I HAD FROM R which are larger migration rates, so we haven't solved the problem of why there are non-zero migration rates for these cases
      
