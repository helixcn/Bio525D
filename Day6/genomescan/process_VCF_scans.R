
#source("http://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")


library(VariantAnnotation)
library(GenomicFeatures)
library(Biostrings)

setwd("~/Documents/My_Documents/UBC/Classes/BIOL525_Bioinformatics/Bio525D/Day6/genomescan")

#read in the data, don't worry about the warnings (comes from lazy maintenance of the VCF header somewhere along the way)
# too big a file   #seqcap <- readVcf(file = "snps.seqcap.vcf.reduce", genome = "scaffold_1000.fasta.contig")
seqcap <- readVcf(file = "snps.gbs.vcf", genome = "scaffold_1000.fasta.contig")

head (seqcap)

info(seqcap)[1:5,]


#show the contents of each genotype entry:
 geno(header(seqcap))



#extract the genotype calls:
GTseqcap <- geno(seqcap)$GT

#find all rows that have a "2" somewhere in them, as these are triallelic sites:
ind1 <- apply (GTseqcap,1,function(x){length(grep("2",x))})

#reduce the genotype table to remove triallelics
GTseqcap <- GTseqcap[ind1 == 0,]

#extract the genotype qualities and remove triallelics:
GQseqcap <- geno(seqcap)$GQ 
GQseqcap <- GQseqcap[ind1 == 0,]

#extract the depths of coverage and remove triallelics:
DPseqcap <- geno(seqcap)$DP
DPseqcap <- DPseqcap[ind1 == 0,]



#convert genotypes to 0,1,2 format, it's easier to work with numeric values:
GTseqcap_f2 <- as.matrix (GTseqcap)
GTseqcap_f2[GTseqcap_f2 == "0/0"] <- 0
GTseqcap_f2[GTseqcap_f2 == "0/1"] <- 1
GTseqcap_f2[GTseqcap_f2 == "1/1"] <- 2
GTseqcap_f2[GTseqcap_f2 == "."] <- NA


#this is a dirty hack to make it read the genotypes in as numerical values, if anyone knows a cleaner way, please let me know. I've just always done it this way!
write.table (GTseqcap_f2, "trash1.txt")
GTseqcap_f2 <- read.table( "trash1.txt")


min_DP <- 5
min_GQ <- 20


#filter based on minimum depth of coverage and genotype quality:

filter1 <- function (input, filter, cutoff){
	
	filter_mask <- filter
	filter_mask[filter_mask < cutoff] <- NA
	filter_mask[filter_mask >= cutoff] <- 1
	
	output <- input * filter_mask	
	
}


GTseqcap_filt1 <- filter1 (GTseqcap_f2, GQseqcap, min_GQ)
GTseqcap_filt2 <- filter1 (GTseqcap_filt1, DPseqcap, min_DP)

#Do some checks now to make sure that the filters worked as you expect them to! Are there other filters you should be applying? (based on yesterday's lecture?)



#check the number of individuals with genotypes called:
ind_seqcap_good <- rowSums (is.na (GTseqcap_filt2)==F)

#now plot a histogram the number of individuals with good genotypes called (cutting out SNPs where fewer than 10 had been called)
hist (ind_seqcap_good[ind_seqcap_good > 10], breaks=50, xlab="Number of individuals with that many good genotype calls", ylab="# SNPs")


#filter, retaining the SNPs that have more than 400 individuals with good genotype calls
GTseqcap_good <- GTseqcap_filt2 [(ind_seqcap_good > 400),]
	# got 704 snps left?


#Calculate the allele frequency at each SNP over the entire population.
GTseqcap_mask <- GTseqcap_good 
GTseqcap_mask[is.na (GTseqcap_mask) == F] <- 2

seqcap_freq <- rowSums (GTseqcap_good,na.rm = T) / rowSums (GTseqcap_mask, na.rm = T)



#calculate nucleotide diversity (pi), from Hohenlohe (2011):

allele1 <- rowSums (GTseqcap_good, na.rm = T)
all <- rowSums (GTseqcap_mask, na.rm = T)
allele2 <- all - allele1

pi <- 1 - ((choose (allele1,2) + choose (allele2,2)) / choose (all,2))

#standard way you think of expected heterozygosity:
exp_het <- 1 - (seqcap_freq^2 + (1 - seqcap_freq)^2)

#these should be identical:
cor.test (exp_het, pi)




#extract all of the positions and contig names:

# a "clean" way to extract info using the accessor functions
pos <- start (ranges(seqcap))
#except now you'd have to cut it down to the same size as the Genotypes table...

# a "dirty" way to extract info by brute force R tricks:
#sub <- unlist(strsplit (rownames (GTseqcap_good), split = ":"))
sub <- strsplit (rownames (GTseqcap_good), split = ":")
contig <- sapply (sub,"[[",1)
sub2 <- strsplit (sapply (sub,"[[",2),split = "_")
pos <- sapply (sub2,"[[",1)


info1 <- cbind (contig,pos,exp_het)


#get one contig out:
sub1 <- info1[info1[,1] == "tscaffold953_243504_255332",]
xlims <- range (as.numeric (sub1[,2]))
par (mar = c (5,5,5,5))
plot (sub1[,2],sub1[,3], type = "l", xlab = "position", ylab = "expected heterozygosity", xlim = xlims)

#extract the same contig from DPseqcap:
indgrep <- grep ("tscaffold953_243504_255332", rownames (DPseqcap))
sub_depth <- DPseqcap[indgrep,]
the_depths <- rowMeans (sub_depth,na.rm = T)

#extract the position of each SNP from the name
sub_split <- strsplit (rownames (sub_depth),split = ":")
sub_split2 <- strsplit (sapply (sub_split,"[[",2),split = "_")
sub_pos <- sapply (sub_split2,"[[",1)

#add depth to the plot:
par (new = TRUE)
plot (sub_pos,the_depths, type = "l", col = "red",xlim = xlims, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis (4)

#observed heterozygosity per individual (heterozygotes have a genotype of "1":
obs_het <- apply (GTseqcap_good,2,function(x){length(which (x == 1)) / length (is.na (x) == F)})

sub_names <- strsplit (colnames (GTseqcap_good), split = "A")
pop_names <- sapply (sub_names,"[[",1)

res_het <- cbind (pop_names, obs_het)

geog <- read.table ("geographic_data.txt", T)

res_combined <- as.matrix (merge (res_het, geog, by.x = "pop_names",by.y = "Population"))

########


#read in the XTX results
xtx <- read.table ("xtx_values.txt",T, comment.char = "&")

hist (log10(xtx$XTX_rank))
summary (xtx$XTX_rank)

hist (table (as.character (xtx$gcontig)))

#calculate the value of the 99th quantile
q1 <- quantile (xtx$XTX_rank, 0.99)

#subset the outliers
outliers <- xtx[xtx$XTX_rank > q1,]
hist (table (as.character (outliers$gcontig)))


#a simple non-parametric test of the number of outliers per contig with outliers, if the same number of SNPs were drawn at random from the genome to be outliers
results <- array (NA, 1000)
res_tab <- NULL

for (i in 1:1000){
	
	samp1 <- as.character (sample (xtx[,1],size = nrow (outliers), replace = F))		
	results[i] <- mean (table (samp1))	
	
	tab1 <- as.data.frame (table (table (samp1)))
	
	res_tab <- rbind (res_tab, tab1)
	
}

#the mean of the null distribution
mean (results)


#calculate the average number of outliers per contig with outliers
mean (table (as.character (outliers$gcontig)))


sub1 <- xtx[xtx[,1] == "C31224742_1_1547",]
plot (sub1[,2],sub1[,3], type = "l", xlab = "position", ylab = "XTX")
arrows (-1000,mean(xtx[,3]),1000000,mean(xtx[,3]),col = "red")
arrows (-1000,q1,1000000,q1,col = "purple")




