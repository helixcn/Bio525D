library(rbamtools)  #library for align tools?

#setwd("/home/darwin/Kim/Exercises_Day2/Alignments_Data")
setwd("~/Documents/My_Documents/UBC/Classes/BIOL525_Bioinformatics/Bio525D/My_Exercises/Jan21_AlignmentExercises")

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


reader2 <- bamReader("cp_350PE_01Err_bwa_sorted.bam", idx=TRUE) # read in the sorted bam file from bowtie2
getRefData(reader2)  # gives the name of the sequence in the reference?

align2 <- getNextAlign(reader2)
failedQC(align2) # did the thing it aligned to fail quality check? false = no
pcrORopt_duplicate(align2) #
name(align2)
flag(align2)
refID(align2)
position(align2)
mapQuality(align2)
cigarData(align2)
nCigar(align2)
mateRefID(align2)

count2 <- bamCountAll(reader2, verbose=TRUE)
count2

coords2 <- c(0,0,15000)
range2 <- bamRange(reader2, coords2)
countNucs(range2)
ncs <- nucStats(reader2)
ncs


xlim2 <- c(1,100000)
coords2 <- c(0,xlim2[1], xlim2[2])
range2 <- bamRange(reader2, coords2)
ad2 <- alignDepth(range2)
mean(getDepth(ad2))
plotAlignDepth(ad2)

