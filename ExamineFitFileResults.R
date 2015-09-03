library(data.table)
#library(animation)	# only if you want to use animation
library(fields)

# functions used later to be able to sort by cross sections of the landscape
make.coords <- function(cell.num){
	if(cell.num %% dim.y == 0){x <- 40}else{x <- cell.num %% dim.y}
	if(x==40){ y <- cell.num %/% dim.y}else{y <- cell.num %/% dim.y + 1}
	return(matrix(c(x,y), nrow=1, dimnames=list(c("coords"), c("lengthwise_coord", "cross_section"))))
}
add.coords <- function(data){
	temp <- do.call(rbind, lapply(data$pop, FUN=make.coords))
	dat <- cbind(data, temp)
	return(dat)
}

dim.y <- 40
dim.x <- 2000

# these are the landscape values over space (there is no cross-sectional variation)
landscape.b0 <- seq(0,0, length.out=2000)
landscape.b0375 <- seq(-37.5,37.5, length.out=2000)
landscape.b375 <- seq(-375,375, length.out=2000)

# make them into matrices to match when we later compare to phenotypic values in each cell
mat.landscape.b0 <- do.call("rbind", replicate(dim.y, landscape.b0, simplify=FALSE))
mat.landscape.b0375 <- do.call("rbind", replicate(dim.y, landscape.b0375, simplify=FALSE))
mat.landscape.b375 <- do.call("rbind", replicate(dim.y, landscape.b375, simplify=FALSE))





### YOU NEED TO SET YOUR OWN WORKING DIRECTORY, AND CAN CHOOSE TO LOOK AT DIFFERENT FILES IF YOU LIKE

setwd("~/Documents/My_Documents/UBC/Research/ExpansionLoad_GroupProject/Westgrid_Runs/Aug28_NewRuns/Results/U0p1_K6_fec3p5_cell50_sig100_b375")


# the files to be analyzed in order of increasing generation time, THESE MUST BE IN YOUR WD
filenames <- c(
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_14000_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_14500_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_14990_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_14995_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15000_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15010_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15020_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15030_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15040_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15050_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15060_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15070_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15080_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15090_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15100_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15110_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15120_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15130_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15140_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15150_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15160_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15170_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15180_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15190_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15200_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15210_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15220_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15230_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15240_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15250_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15260_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15270_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15280_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15290_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15300_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15350_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15400_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15450_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15500_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15550_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15600_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15650_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15700_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15750_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15800_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15850_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15900_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_15950_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16000_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16050_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16100_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16150_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16200_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16250_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16300_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16350_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16400_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16450_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16500_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16550_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16600_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16650_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16700_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16750_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16800_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16850_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16900_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_16950_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17000_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17050_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17100_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17150_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17200_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17250_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17300_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17350_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17400_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17450_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17500_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17550_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17600_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17650_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17700_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17750_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17800_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17850_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17900_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_17950_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18000_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18050_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18100_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18150_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18200_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18250_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18300_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18350_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18400_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18450_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18500_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18550_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18600_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18650_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18700_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18750_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18800_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18850_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18900_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_18950_1.fit",
"Sep1_U0p1_K6_fec3p5_cell50_sig100_b375_19000_1.fit"
)




# make a 1-D plot (saved in wd) of average fitness over space at each generation 
png("1D_ExpansionFitness.png", height=(length(filenames)*1.75), width=8, units="in", res=100)
par(mfrow=c(length(filenames),1), mar=c(2.5,2.5,1,1))
for(i in 1:length(filenames)){
#	dev.hold()	# if you want to animate instead
	gen <- tail(strsplit(filenames[i], split="_")[[1]])[5]	# the generation we are looking at with this file
	#read in data
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	# add coords (see top of file)
	dat <- add.coords(dat)
	# take the mean per cross section
	agg.dat.mean <- aggregate(.~ cross_section, data=dat, mean)
	# plot it over time
	plot(1:length(agg.dat.mean$fitness_trait_2), agg.dat.mean$fitness_trait_2, xlim=c(0,2000), ylim=c(0.5,1), type="l", xlab="Space", ylab="Fitness of Landscape Cross Section")
	text(c(1700, 1900),0.6, paste(c("generation", gen)))
	print(i)
#	ani.pause(0.1)	# if you want to animate instead; larger number makes it slower
}	
dev.off()
#_________________________________________________________________________________________________________________________#





# make a 2-D plot (saved in wd) of fitness over space at each generation 
png("2D_ExpansionFitness.png", height=8, width=(length(filenames)*1.75), units="in", res=100)
par(mfrow=c(1,length(filenames)), mar=c(1,1,2.5,2.5))

# YOU NEED TO SET THIS because in each case the scale of image.plot will adjust to fit the dataset at that generation, but to make the plots comparable across time points, set this value to the lowest fitness that is experienced in any of the plots
min <- 0.45

for(i in 1:length(filenames)){
	gen <- tail(strsplit(filenames[i], split="_")[[1]])[5]	# the generation we are looking at with this file
	#read in data
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	# get mean fitness per cell (cells are identified by "pop")
	agg.dat.mean <- aggregate(.~pop, data=dat, mean)
	nas <- rep(NA,80000)
	nas[agg.dat.mean$pop] <- agg.dat.mean$fitness_trait_2		# *** NEEDS TO BE QUANTI TRAIT
	mat.results <- matrix(nas, nrow=dim.y, ncol=dim.x)
	mat.results[1,1] <- min; mat.results[1,2] <- 1				# ***** ARTIFICIALLY SETTING so multiple plots have same scale
	image(z=matrix(0,nrow=dim.y,ncol=dim.x),x=1:dim.y,y=1:dim.x,col=topo.colors(1),ylab="",xlab="",cex.axis=0.5,main=gen)
	image.plot(z=mat.results, x=1:dim.y, y=1:dim.x, legend.width=4, col=heat.colors(20), main="", add=TRUE)
	print(i)
}	
dev.off()
#_________________________________________________________________________________________________________________________#





# look at genetic variance and mean genotypic value
png("GeneticVarianceAndMean.png", height=(length(filenames)*2.5), width=12, units="in", res=100)
par(mfrow=c(length(filenames),2), mar=c(2,4,1,1))

for(i in 1:length(filenames)){
	gen <- tail(strsplit(filenames[i], split="_")[[1]])[5] 
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	# add coords (see top of file)
	dat <- add.coords(dat)
	agg.dat.mean <- aggregate(.~cross_section, data=dat, mean)	# mean genotypic value per cross section
	agg.dat.var <- aggregate(.~cross_section, data=dat, var)	# variance in genotypic value per cross section
	plot(1:length(agg.dat.mean$geno), agg.dat.mean$geno, xlim=c(0,2000), ylim=c(-0.5,0.5), type="l", xlab="Space", ylab="Mean Genetic Value of Landscape Cross Section")
	text(1900,0.2, gen)
	points(1:2000, landscape.b0, col="grey", type="l")
	plot(1:length(agg.dat.var$geno), agg.dat.var$geno, xlim=c(0,2000), ylim=c(0,0.5), col="red", type="l", xlab="Space", ylab="Genetic Variance of Landscape Cross Section")
	text(1900,0.5, gen)
	print(i)
}
dev.off()
#_________________________________________________________________________________________________________________________#





# make a 2-D plot (saved in wd) of number of individuals per cell at each generation 
png("2D_PopSize.png", height=8, width=(length(filenames)*1.75), units="in", res=100)
par(mfrow=c(1,length(filenames)), mar=c(1,1,2.5,2.5))

for(i in 1:length(filenames)){
	gen <- tail(strsplit(filenames[i], split="_")[[1]])[5]	# the generation we are looking at with this file
	#read in data
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	# get mean fitness per cell (cells are identified by "pop")
	agg.dat.num.inds <- aggregate(.~pop, data=dat, length)
	nas <- rep(NA,80000)
	nas[agg.dat.num.inds$pop] <- agg.dat.num.inds$geno
	mat.results <- matrix(nas, nrow=dim.y, ncol=dim.x)
	mat.results[1,1] <- 0; mat.results[1,2] <- 1				# ***** ARTIFICIALLY SETTING so multiple plots have same scale
	image(z=matrix(0,nrow=dim.y,ncol=dim.x),x=1:dim.y,y=1:dim.x,col=terrain.colors(1),ylab="",xlab="",cex.axis=0.5,main=gen)
	image.plot(z=mat.results, x=1:dim.y, y=1:dim.x, legend.width=4, nlevels=7, main="", add=TRUE)
}	
dev.off()
#_________________________________________________________________________________________________________________________#






# make a 1-D plot (saved in wd) of the difference between mean cross section phenotype and landscape optimum
png("1D_Pheno-Landscape.png", height=(length(filenames)*2.5), width=12, units="in", res=100)
par(mfrow=c(length(filenames),2), mar=c(2,4,1,1))

# UNCOMMENT THE CORRECT ONE
landscape <- landscape.b0; pylim=c(-1,1)
#landscape <- landscape.b0375; pylim=c(-40,40)
#landscape <- landscape.b375; pylim=c(-400,400)

for(i in 1:length(filenames)){
	gen <- tail(strsplit(filenames[i], split="_")[[1]])[5] 
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	# add coords (see top of file)
	dat <- add.coords(dat)
	agg.dat.mean <- aggregate(.~cross_section, data=dat, mean)	# mean genotypic value per cross section
	plot(1:2000, landscape, type="l", col="darkgray", ylab="Mean Phenotype", ylim=pylim)
	text(1900,0.5, gen, cex=0.9)
	points(1:length(agg.dat.mean$pheno), agg.dat.mean$pheno, col="darkred", type="l")
	plot(agg.dat.mean$pheno - landscape[1:length(agg.dat.mean$pheno)], xlim=c(0,2000), ylim=c(-0.5,0.5), type="l", ylab="Mean Phenotype - Landscape Optimum")
	abline(h=0, col="darkgray")
	print(i)
}
dev.off()
#_________________________________________________________________________________________________________________________#




# make a 2-D plot (saved in wd) of the difference between mean cell phenotype and landscape optimum
png("2D_Pheno-Landscape.png", height=6, width=(length(filenames)*2+2), units="in", res=100)
par(mar=c(1,1,1,2.5), mfrow=c(1,length(filenames)+1))

# UNCOMMENT THE CORRECT ONE
landscape <- mat.landscape.b0
#landscape <- mat.landscape.b0375
#landscape <- mat.landscape.b375

# SET to keep scale constant across plots; increase or decrease accordingly if scales don't match or are too extreme
min <- -3.5
max <- 2.5

#plot the landscape first then the phenotype difference results
image.plot(z= landscape, x=1:dim.y, y=1:dim.x, legend.width=4, col=heat.colors(2000), main="LANDSCAPE")

for(i in 1:length(filenames)){
	#read in data
	dat <- fread(filenames[i], header=TRUE, sep=" ")
	#get mean pheno per cell
	agg.dat.mean <- aggregate(.~pop, data=dat, mean)
	# put into matrix format
	nas <- rep(NA,80000)
	nas[agg.dat.mean$pop] <- agg.dat.mean$pheno
	mat.results <- matrix(nas, nrow=dim.y, ncol=dim.x)
	diffs <- mat.results - landscape
	mat.results[1,1] <- min; mat.results[1,2] <- max			# ***** ARTIFICIALLY SETTING so multiple plots have same scale
	image(z=matrix(0,nrow=dim.y,ncol=dim.x),x=1:dim.y,y=1:dim.x,col=terrain.colors(1),ylab="",xlab="",cex.axis=0.5,main=gen)
	image.plot(z=mat.results, x=1:dim.y, y=1:dim.x, legend.width=4, col=heat.colors(100), main="", add=TRUE)
	print(i)
}
dev.off()

