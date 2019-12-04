
# This example R code takes two files containing the genotype data and geographic locations of sampled individuals, and sorts
# the rows of the file containing the geographical locations such that the ith rows of the two files corresponds to the ith sampled individual
# (i.e., the ith row corresponds to the same sampled individual). The code then creates a new text file of the geographic locations with the 
# sorted rows. 
# Before running this code, the working directory should be set to the one containing the two files.
# The example uses the two files "Pigeon.gen" and "Pigeon.coordinates", corresponding to data for rock pigeons in Singapore (Tang et al., 2018).

# Read in the file of genotype data and extract the names of the sampled individuals.
# The first column in the file of genotype data should contain the names of the sampled individuals.
library(adegenet)
gen_data <- read.genepop("Pigeon.gen", ncode=3)
tt <- dimnames(gen_data@tab)
tq <- tt[1]
qq <- as.data.frame(tq)

# Read in the file of geographical locations and sort the rows such that the ith row corresponds to the same individual as the ith
# row of the file of geneotype data. 
# The first column in the file of geographical locations should contain the names of the sampled individuals, and should have the column name "X".
sample.points <- read.table(file = "Pigeon.coordinates", sep = '\t', header = T, fill = T)
names(qq) <- c("X")
rs <- merge(qq, sample.points, by = "X")

# Write the file of geographical locations with sorted rows to a new text file.
write.table(rs, "Pigeon2.coordinates", sep = "\t", row.names = FALSE)
