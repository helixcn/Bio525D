library(rbamtools)  #library for align tools?

setwd("/home/darwin/Kim/Exercises_Day2/Alignments_Data")

reader <- bamReader("cp_350PE_01Err_bowtie2_sorted.bam", idx=TRUE) # read in the sorted bam file from bowtie2
getRefData(reader)  # gives the name of the sequence in the reference?

align <- getNextAlign(reader)
failedQC(align) # did the thing it aligned to fail quality check? false = no
pcrORopt_duplicate(align) #
name(align)
flag(align)
refID(align)
position(align)
mapQuality(align)
cigarData(align)
nCigar(align)
mateRefID(align)

count <- bamCountAll(reader, verbose=TRUE)
count

coords <- c(0,0,15000)
range <- bamRange(reader, coords)
countNucs(range)
ncs <- nucStats(reader)
ncs


xlim <- c(1,100000)
coords <- c(0,xlim[1], xlim[2])
range <- bamRange(reader, coords)
ad <- alignDepth(range)
mean(getDepth(ad))
plotAlignDepth(ad)


reader <- bamReader("cp_350PE_01Err_bwa_sorted.bam", idx=TRUE) # read in the sorted bam file from bowtie2
getRefData(reader)  # gives the name of the sequence in the reference?

align <- getNextAlign(reader)
failedQC(align) # did the thing it aligned to fail quality check? false = no
pcrORopt_duplicate(align) #
name(align)
flag(align)
refID(align)
position(align)
mapQuality(align)
cigarData(align)
nCigar(align)
mateRefID(align)

count <- bamCountAll(reader, verbose=TRUE)
count

coords <- c(0,0,15000)
range <- bamRange(reader, coords)
countNucs(range)
ncs <- nucStats(reader)
ncs


xlim <- c(1,100000)
coords <- c(0,xlim[1], xlim[2])
range <- bamRange(reader, coords)
ad <- alignDepth(range)
mean(getDepth(ad))
plotAlignDepth(ad)

