#### WITH THE ACTUAL MLNE INPUT FILES

input <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_1.4_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
(ans <- mlne.moment(input)) 


# equation 8: m = 1-t_root(1 - (1/W) * sum over all l of w_l*F_l)
#	equals for t=1: (1/W)*sum(w_l*F_l)

mlne.moment <- function(input){
		B_T <- NULL
		B_0 <- NULL
		A_0 <- NULL

	dat <- readLines(input)
	focal_current <- dat[11:50]
	focal_t0 <- dat[52:91]
	source_t0 <- dat[94:133]
	#max.loci <- dat[7]
	#max.loci.list <- strsplit(max.loci, split=",")
	
	for(i in 1:40){
		temp.focal_t1 <- as.numeric(strsplit(focal_current[i], split=" ")[[1]]) / sum(as.numeric(strsplit(focal_current[i], split=" ")[[1]]))
		temp.focal_t0 <- as.numeric(strsplit(focal_t0[i], split=" ")[[1]]) / sum(as.numeric(strsplit(focal_t0[i], split=" ")[[1]]))
		temp.source_t0 <- as.numeric(strsplit(source_t0[i], split=" ")[[1]]) / sum(as.numeric(strsplit(source_t0[i], split=" ")[[1]]))
		B_T <- c(B_T, temp.focal_t1[1])
		B_0 <- c(B_0, temp.focal_t0[1])
		A_0 <- c(A_0, temp.source_t0[1])
	}
	x_ab <- A_0 - B_0
	x_bt <- B_T - B_0

	w_l <- x_ab^2
	F_l <- x_bt/x_ab
	W <- sum(w_l)

	return(mig.rate <- (1/W)*sum(w_l*F_l, na.rm=TRUE))	# Fortran eqn: MTM2=1.-ABS(1.-WF/WW)**(1./nt)
}




#input <- "~/Desktop/NewMLNeSims/CalculateAlleleFreqs/MLNeIN_5.6_NoMig_mig500_meta0.1_gen16_NoMig_t2_sampled_1_1.dat"
#(ans1 <- mlne.moment(input))



in1 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_1.4_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in2 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_1.4-14_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in3 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_2.4-14_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in4 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_2.7_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in5 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_3.4-14_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in6 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_3.11_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in7 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_5.4.6-14_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"
in8 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_5.6_NoMig_mig500_meta0.1_gen17_NoMig_2alleles_sampled_1_1.dat"

in9 <- "~/Desktop/NewMLNeSims/Compare_R_MLNE/Oct8_MLNeIN_1.4-14_NoMig_mig500_meta0.01_gen167_NoMig_2alleles_sampled_1_1.dat"

(ans1 <- mlne.moment(in1))	# 1.4		0.1		R says m = 0.1362797, 	MLNe says m = 0.0740
(ans2 <- mlne.moment(in2))	# 1.4-14	0.1		R says m = 0.1635086, 	MLNe says m = 0.0768
(ans3 <- mlne.moment(in3))	# 2.4-14	0.1		R says m = -0.09536082, 	MLNe says m = 0.3
(ans4 <- mlne.moment(in4))	# 2.7		0.1		R says m = 0.03055741, 	MLNe says m = 
(ans5 <- mlne.moment(in5))	# 3.4-14	0.1		R says m = 0.3163315, 	MLNe says m = 
(ans6 <- mlne.moment(in6))	# 3.11		0.1		R says m = 0.2825987, 	MLNe says m = 
(ans7 <- mlne.moment(in7))	# 5.4,6-14	0.1		R says m = 0.2090635, 	MLNe says m = 
(ans8 <- mlne.moment(in8))	# 5.6		0.1		R says m = 0.1514543, 	MLNe says m = 


(ans9 <- mlne.moment(in9))	# 2.4-14	0.01 rep 1		R says m = 0.04111848, 	MLNe says m = 0.0103

















#	MTM=0.
#      NN=0
#      DO KK=1,NSAMB-1
#        KK1=KK+1
#        IF(M_ESTIMATE==1)THEN           !doing stuff here for moment migration estimate
#          DBT=0.
#          WNE=0.
#          DAB0=0.
#          WW=0.
#          WF=0.
#        ELSE
#          FC=0.
#          NLOCI_FC=0
#        ENDIF
#        DO I=1,NTLOCI
#          B0=NSAMP(I,KK)/REAL(NS2(I,KK))			# dividing by total?
#          BT=NSAMP(I,KK1)/REAL(NS2(I,KK1))
#          IF(M_ESTIMATE==1)THEN
#            A0=NSAMP(I,0)/REAL(NS2(I,0))
#            XT=abs(A0-B0)      
#            DAB0=DAB0-XT*XT
#            IF(XT>TINY)THEN
#              WI=(A0-B0)**2
#              WW=WW+WI
#              FI=(BT-B0)/(A0-B0)
#              WF=WF+WI*FI
#            ENDIF
#            PQ=(B0*(1.-B0)+BT*(1.-BT))*0.5
#            WNE=WNE+PQ
#            DBT=DBT+(BT-B0)**2-B0*(1-B0)/NS2(I,KK)-BT*(1.-BT)/NS2(I,KK1)
#          ELSE  !only if not estimating m
#            XT=BT+B0
#            IF(XT>=TINY.AND.XT<2.)THEN
#              FC=FC+2.*(BT-B0)**2*(1./XT+1./(2.-XT))-1./NS2(I,KK)-1./NS2(I,KK1)
#              NLOCI_FC=NLOCI_FC+1
#            ENDIF
#          ENDIF
#        END DO
#        NT=NGEN(KK1)-NGEN(KK)
#        IF(M_ESTIMATE==1)THEN
#          IF(WW<=0.)CYCLE
#          NN=NN+1
#          MTM2=1.-ABS(1.-WF/WW)**(1./nt)
#          MTM=MTM+MTM2
#          IF(ABS(MTM2)>TINY)THEN
#            MTNE2=(DBT+DAB0*(1.-(1.-MTM2)**NT)**2)/(WNE*  &
#                 (1.-(1.-MTM2)**(2*NT))/(MTM2*(2.-MTM2)))
#          ELSE
#            MTNE2=DBT/(WNE*NT)
#          ENDIF
#        ELSE
#          MTNE2=FC/NLOCI_FC/REAL(NT)
#          NN=NN+1
#        ENDIF
#        IF(MTNE2>0.)MTNE=MTNE+MTNE2
#      ENDDO
#      IF(NN==0)THEN
#        MTNE=-1.
#        MTM=-10.
#      ELSE
#        IF(MTNE>TINY)THEN
#          MTNE=NN/MTNE*.5
#        ELSE
#          MTNE=-1.
#        ENDIF
#        IF(M_ESTIMATE==1) MTM=MTM/NN
#      ENDIF
#      END SUBROUTINE MTEST
      
