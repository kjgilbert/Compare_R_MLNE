#############################################################################################################################
#
#  code to write the dispersal matrix for a given number of patches and dispersal kernel to go into Nemo
#  19 Jan 2015
#
#############################################################################################################################



# FOR A LANDSCAPE OF A GIVEN SIZE

#patches.x <- 1000	# number of patches going horizontally in the cylinder, the range will expand in this direction
#patches.y <- 200		# number of patches lined up vertically in the cylinder
	# mostly irrelevant for figuring out dispersal kernel, except don't want it to loop entirely around the cylinder
	# cut in half when creating the landscape because I want symmetry so the edges match up, the landscape will be made then mirrored

# how big in 1-D is a patch (meters here)
cell <- 50
# landscape size
land.x <- 1000	# 1000 km = 1,000,000 m
land.y <- 500		# 50 km = 50,000 m

cells.x <- land.x/cell
cells.y <- land.y/cell

(total.cells <- cells.x * cells.y)



# CREATE THE DISPERSAL MATRIX

# dispersal kernel from 1 patch to surrounding patches
#	want it to be Gaussian and vary to having more LDD (long distance dispersal)
dist.mean <- 0		# mean should stay at zero as dispersal is centered around natal patch
dist.sd <- 150		# this is sigma, one standard deviation, in meters

# scale the landscape, and want to cut off at 4 sigma, so how many cells (in 1-D) are there within 4 sigma:
units <- dist.sd/cell		# number of cells per one sigma
four.sigma.units <- 4*units	# number of cells per four sigma IN ONE DIRECTION

center.cell <- four.sigma.units + 1
cells.in.1D.kernel <- four.sigma.units + four.sigma.units + 1
		# total cells going across the one-dimensional kernel are those one either side plus the central/natal one

# make the 1-D kernel
kernel.1d <- rep(NA, cells.in.1D.kernel)
k <- 1
for(i in -four.sigma.units:four.sigma.units){
	left.boundary <- (i-0.5)*cell
	right.boundary <- (i+0.5)*cell
	kernel.1d[k] <- pnorm(q=right.boundary, sd=dist.sd) - pnorm(q=left.boundary, sd=dist.sd)
	k <- k+1
}
	# correct to have the halves because the center cell in a sense can be thought of as zero space since not migrating leaves you there
	# still need to renormalize because only have the cumulative within 4 sigma
normalized.kernel.1d <- kernel.1d/sum(kernel.1d)	# normalize the probabilities to sum to 1
plot(1:cells.in.1D.kernel, normalized.kernel.1d)


########################################################################################################################
########################################################################################################################
#
#	for summing 2 distributions to create kurtosis:
#
# do the above again for a second distribution, sum the 1-D kernels and normalize
#
#	*** NOTE ***
#	WHEN SUMMING THE KERNELS, MAKE SURE TO ADD AFTER CENTERING AROUND 0
#
# dispersal kernel from 1 patch to surrounding patches
#	want it to be Gaussian and vary to having more LDD (long distance dispersal)
second.dist.mean <- 0		# mean should stay at zero as dispersal is centered around natal patch
second.dist.sd <- 1000		# this is sigma, one standard deviation, in meters

# scale the landscape, and want to cut off at 4 sigma, so how many cells (in 1-D) are there within 4 sigma:
second.units <- second.dist.sd/cell		# number of cells per one sigma
second.four.sigma.units <- 4*second.units	# number of cells per four sigma IN ONE DIRECTION

second.center.cell <- second.four.sigma.units + 1
second.cells.in.1D.kernel <- second.four.sigma.units + second.four.sigma.units + 1
		# total cells going across the one-dimensional kernel are those one either side plus the central/natal one

# make the 1-D kernel
second.kernel.1d <- rep(NA, second.cells.in.1D.kernel)
k <- 1
for(i in -second.four.sigma.units:second.four.sigma.units){
	left.boundary <- (i-0.5)*cell
	right.boundary <- (i+0.5)*cell
	second.kernel.1d[k] <- pnorm(q=right.boundary, sd=second.dist.sd) - pnorm(q=left.boundary, sd=second.dist.sd)
	k <- k+1
}
	# correct to have the halves because the center cell in a sense can be thought of as zero space since not migrating leaves you there
	# still need to renormalize because only have the cumulative within 4 sigma
second.normalized.kernel.1d <- second.kernel.1d/sum(second.kernel.1d)	# normalize the probabilities to sum to 1

diff.in.length <- abs(length(normalized.kernel.1d) - length(second.normalized.kernel.1d))
each.end.diff <- diff.in.length/2
# which is the shorter one:
if((length(normalized.kernel.1d) - length(second.normalized.kernel.1d)) < 0){	# then first kernel is shorter, add to it's length
	second.shorter <- FALSE
	long.normalized.kernel.1d <- c(rep(0, each.end.diff), normalized.kernel.1d, rep(0, each.end.diff))
	plot(-second.four.sigma.units:second.four.sigma.units, second.normalized.kernel.1d, ylim=c(0,0.5), col="red")
	points(-second.four.sigma.units:second.four.sigma.units, long.normalized.kernel.1d, col="blue")
}
if((length(normalized.kernel.1d) - length(second.normalized.kernel.1d)) > 0){	# then second kernel is shorter, add to it's length
	second.shorter <- TRUE
	second.normalized.kernel.1d <- c(rep(0, each.end.diff), second.normalized.kernel.1d, rep(0, each.end.diff))
	plot(-four.sigma.units:four.sigma.units, normalized.kernel.1d, ylim=c(0,0.5), col="red")
	points(-four.sigma.units:four.sigma.units, second.normalized.kernel.1d, col="blue")
}

prop.first.kernel <- 0.5
prop.second.kernel <- 1 - prop.first.kernel
summed.1d.kernels <- (long.normalized.kernel.1d * prop.first.kernel) + (second.normalized.kernel.1d * prop.second.kernel)

# just for the plot, to know which x values to use
if(second.shorter==FALSE){
	points(-second.four.sigma.units:second.four.sigma.units, summed.1d.kernels, pch=20, cex=0.75)
}else{
	points(-four.sigma.units:four.sigma.units, summed.1d.kernels, pch=20, cex=0.75)
}

########################################################################################################################
########################################################################################################################


# IF USING COMBINATION OF 2 KERNELS, UNCODE FOLLOWING LINE
#
# normalized.kernel.1d <- summed.1d.kernels

#normalized.kernel.1d <- kernel.1d/sum(kernel.1d)	# normalize the probabilities to sum to 1
plot(1:cells.in.1D.kernel, normalized.kernel.1d)

# multiply it by itself to create the 2-D
horizontal.to.multiply <- matrix(NA, nrow=length(normalized.kernel.1d), ncol=length(normalized.kernel.1d))
vertical.to.multiply <- matrix(NA, nrow=length(normalized.kernel.1d), ncol=length(normalized.kernel.1d))
for(j in 1:length(normalized.kernel.1d)){
	horizontal.to.multiply[j,] <- normalized.kernel.1d
	vertical.to.multiply[,j] <- normalized.kernel.1d
}

multiplied.kernels <- horizontal.to.multiply * vertical.to.multiply


# plot the multiplied matrix
library(graphics)
contour(multiplied.kernels, asp=1, nlevels= four.sigma.units)

# find the cutoff value for distance travelled, i.e. which contours don't make the full circle around because they travel farther than the central column/row?
lower.cutoff <- multiplied.kernels[1,center.cell]
multiplied.kernels[multiplied.kernels < lower.cutoff] <- 0	# replace values lower than the cutoff with zero
# plot the multiplied matrix that's been cut off to circular distances
contour(multiplied.kernels, asp=1, nlevels= four.sigma.units)

# restandardize so all sums to 1
restandardized.multiplied.kernels <- multiplied.kernels/sum(multiplied.kernels)
contour(restandardized.multiplied.kernels, asp=1, nlevels= four.sigma.units)
disp.kernel <- restandardized.multiplied.kernels

middle.of.kernel <- center.cell
	# ceiling(length(disp.kernel[,1])/2)	# works because there are an odd number for the length/width
max.dist.include.natal <- middle.of.kernel	# the maximum number of cells an individual might migrate including 1 as staying in natal patch
max.dist.exclude.natal <- middle.of.kernel - 1


#############################################################################################################################
#
#	NOW WE HAVE THE DISPERSAL KERNEL
#
#	GO AHEAD AND MAKE THE DISPERSAL MATRIX
#
#############################################################################################################################



# dispersal matrix needs to be patch number x patch number with columns summing to 1
#patches.x
#patches.y
total.patches <- patches.x * patches.y

disp.mat <- matrix(0, nrow= total.patches, ncol= total.patches) 	
	# dispersal matrix is patch number by patch number
	# it is all zeros so I can simply fill in the non-zero values
	
# make an identifying number for all patches, go down each column before moving on to the next column
patch.number.ids <- matrix(0, nrow= patches.y, ncol= patches.x)
for(i in 1: (patches.x * patches.y)){
	patch.number.ids[i] <- i
}



# a function to mirror image a matrix
rotate90 <- function(x) t(apply(x, 2, rev))
mirror <- function(x) rotate90(rotate90(x))


# FOR A GIVEN PATCH, WHAT PATCHES ARE THERE TO CHOOSE FROM FOR MIGRATION

column.length <- patches.y
row.length <- patches.x

for(i in 1: patches.y){		# row
	for(j in 1: patches.x){	# column

		focal.patch.id <- patch.number.ids[i,j]
		potential.rows.vertical <- (i - max.dist.exclude.natal):(i + max.dist.exclude.natal)
		potential.columns.horizontal <- (j - max.dist.exclude.natal):(j + max.dist.exclude.natal)

		# check within bounds for rows
		if(TRUE %in% (potential.rows.vertical > column.length)){	# at the bottom of the cylinder, loop to top
			how.many.over <- table((potential.rows.vertical > column.length))[[2]]	
				# number of trues, i.e. those that are greater than the max column height
			potential.rows.vertical <- c((i - max.dist.exclude.natal):(column.length), 1:how.many.over)
		}
		if(TRUE %in% (potential.rows.vertical < 1)){	# at the top of the cylinder, loop to bottom
			how.many.under <- table((potential.rows.vertical < 1))[[2]]	
				# number of trues, i.e. those that are greater than the max column height
			potential.rows.vertical <- c((column.length-(how.many.under-1)):column.length, 1:(i + max.dist.exclude.natal))
		}
		# check within bounds for columns; DO NOT want these to loop to the other end of the x spectrum
		if(TRUE %in% (potential.columns.horizontal > row.length)){	# at the right end of the cylinder, don't loop
			how.many.over <- table((potential.columns.horizontal > row.length))[[2]]	
				# number of trues, i.e. those that are greater than the right end of the cylinder
			potential.columns.horizontal <- (j - max.dist.exclude.natal):(row.length)
		}
		if(TRUE %in% (potential.columns.horizontal < 1)){	# at the left end of the cylinder, don't loop 
			how.many.under <- table((potential.columns.horizontal < 1))[[2]]	
				# number of trues, i.e. those that are less than the left end of the cylinder
			potential.columns.horizontal <- 1:(j + max.dist.exclude.natal)
		}

		# take those columns and rows and create a list of potential patches
		potential.patches.matching.kernel <- patch.number.ids[potential.rows.vertical, potential.columns.horizontal]




		if(length(potential.patches.matching.kernel) < length(disp.kernel)){
			print("Patch dispersal kernel hits the left or right border!")
			# when at the left or right end of the cylinder, need to add in the extra migration - make boundaries reflecting
			
			# how many columns should be in the kernel total:
			kernel.length <- max.dist.include.natal + max.dist.exclude.natal
			
			
			# are we at the left or right edge?
			
			# left edge
			if(TRUE %in% (1: patches.y %in% potential.patches.matching.kernel)){
				# how many columns do we have / are we missing in the kernel
				short.kernel <- length(potential.patches.matching.kernel[1,])
				cols.missing <- kernel.length - short.kernel
				
				# missing that many columns, and in this case on the left
				# so mirror those 2 columns and add them on to the dispersal kernel existing leftmost columns
				
				modified.disp.kernel <- disp.kernel[, -c(1:cols.missing)]
				add.on.missing <- disp.kernel[, c(1:cols.missing)]
				if(cols.missing > 1){
					# mirror the dispersal kernel bit being added back on
					add.on.missing <- mirror(add.on.missing)
					for(k in 1:cols.missing){
						modified.disp.kernel[,k] <- modified.disp.kernel[,k] + add.on.missing[,k]
					}
				}else{
					modified.disp.kernel[,1] <- modified.disp.kernel[,1] + add.on.missing
				}
			}
			
			
			# right edge
			if(TRUE %in% (((total.patches-patches.y): total.patches) %in% potential.patches.matching.kernel)){
				# how many columns do we have / are we missing in the kernel
				short.kernel <- length(potential.patches.matching.kernel[1,])
				cols.missing <- kernel.length - short.kernel
				
				# missing that many columns, and in this case on the right
				# so mirror those columns if there are more than 1 and add them on to the dispersal kernel existing rightmost columns
				
				modified.disp.kernel <- disp.kernel[, -c(kernel.length:(kernel.length - (cols.missing-1)))]
				add.on.missing <- disp.kernel[, c(kernel.length:(kernel.length - (cols.missing-1)))]	# this ALREADY mirrored any columns greater than 1
				
				# harder on the right side to line up matrices for adding, so just add zeros to the existing reflected probabilities that are being added on
				empty.cols.to.add <- short.kernel - cols.missing
				
				# harder on the right side to line up matrices for adding, so just add zeros to the existing reflected probabilities that are being added on
				for(kk in 1: empty.cols.to.add){
					add.on.missing  <- cbind(rep(0, kernel.length), add.on.missing)
				}
				modified.disp.kernel <- modified.disp.kernel + add.on.missing				
			}
			# put the dispersal probability values into the dispersal matrix
			disp.mat[focal.patch.id, potential.patches.matching.kernel] <- modified.disp.kernel
			
		}else{
			#this is the simple case where values can just be inserted because we are not at the left or right boundary
			
			disp.mat[focal.patch.id, potential.patches.matching.kernel] <- disp.kernel	
				# for the focal patch we are on, that is the row number in the final dispersal matrix for nemo
				# insert into that row the values of where that patch sends migrants
		}
	}
}



setwd("~/Documents/My_Documents/UBC/Research/Range_SimulationProject/R_Code")

write(disp.mat, file="DispMatrix.txt", ncolumns=(total.patches), sep=", ")
	# this adds in the commas and maintains the columns, need the brackets another way


# can't get the brackets in properly without it writing them as a row instead of a column.  Quicker fix to just paste file into excel and concatenate the brackets in that way




#############################################################################################################################
#############################################################################################################################


# PATCH CAPACITY MATRIX

# also write the patch capacity matrix

# carrying capacity
K <- 10

# where to initiate the pop

# total number of pops
total.patches
# those patches are arrayed by:
# vertical number of patches
patches.y
# horizontal number of patches
patches.x

# width to initiate the pops over
initial.x <- 5

# that width creates this number of pops which are ordered correctly at the beginning of the list because of how the dispersal matriz has been constructed
initial.num.pops <- initial.x * patches.y


# FOR BURN IN
burnin <- paste(c("{", paste(rep(K, initial.num.pops), collapse=", "), ",", paste(rep(0, (total.patches - initial.num.pops)), collapse=", "), "}"), collapse=" ")

# FOR EXPANSION
expansion <- paste(c("{", paste(rep(K, total.patches), collapse=", "), "}"), collapse=" ")








#############################################################################################################################










#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#temp <- as.matrix(cbind.data.frame(as.matrix(rep("{", total.patches)), disp.mat, as.matrix(rep("}", total.patches))))
#temp2 <- disp.mat
#temp3 <- cbind(NA, disp.mat, NA)
#temp2$end <- rep("}", total.patches)
#new.mat <- read.table("DispMatrix.txt")
#test <- cbind(rep("{", total.patches), new.mat, rep("}", total.patches))
#write(as.matrix(test), file="test_DispMatrix.txt", ncolumns=(total.patches+2))


# forward migration, each row needs to sum to 1
#for(i in 1: total.patches){		# go through one row, i.e. one patch
#	for(j in 1: total.patches){	# proceed through that row's columns, i.e. one patch
		# patch 1 to patch 1
		# 2 to 2
		# 3 to 3 .... etc
#		temp.sum <- 0
#		if(i == j) {
#			disp.mat[i,j] <- disp.kernel[middle.of.kernel, middle.of.kernel] # we are in the natal patch, give it that dispersal probability
#			temp.sum <- temp.sum + disp.kernel[middle.of.kernel, middle.of.kernel]	# want to make sure we stop putting values in when we've reached the sum of 1
#		}
#		if(i != j){		# we are not in the natal patch, figure out where we are in the kernel
#			diff <- j-i	# how many patches are we from the natal patch in that row?
#			if(diff <= max.dist.include.natal){			
#			}
#		}
#	}
#}




# not sure if I need this? patch ids jsut within the dispersal kernel
#kernel.ids <- matrix(0, nrow=length(disp.kernel[,1]), ncol=length(disp.kernel[1,]))
#for(i in 1:length(disp.kernel)){
#	kernel.ids[i] <- i
#}

# don't need the array because things seem to stay lined up without it, so can just paste them in, saving the code below just in case though

#potential.patches.matching.kernel.with.disp.array <- c(as.vector(potential.patches.matching.kernel), as.vector(disp.kernel))
#dim(potential.patches.matching.kernel.with.disp.array) <- c((max.dist.include.natal+ max.dist.exclude.natal), (max.dist.include.natal+ max.dist.exclude.natal), 2)

# NOW HAVE AN ARRAY TO CALL FROM for patch ID and it's probability of being dispersed to by the focal patch



