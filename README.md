## RAD-SSR_extraction_scripts
## 9/8/2018 CURRENTLY BEING UPDATED
## Come back soon for a full reproducable example

## libraries
library(gdata)
library(plyr)
library(zoo)

## pipeline plan
### 1. Define the repeat motifs you would like to search for, simple bialleleic repeats are a good place to start.
### 2. Search for RAD tags that contain these motifs
### 3. Filter the selected tags to retain the best candidates. For example, those present in many individuals and without repeats at the beginning and end of the read etc.
### 4. Generate 'primers' for those RAD tags, in the form of the first 15-25 bases, to allow us to search for repeat motifs in the whole read library.
### 5. Search for selected reads across all individuals.
### 6. Run a script to convert the repeat motifs in raw reads into factors, and to call heterozygotes or homozygotes depending on the number of matching reads.
### 7. Convert this into desired file format e.g. .stru
### 8. Apply any biases, if you wish, for example selecting highly variable loci in the same manner that PCR-SSRs are selected.


#### files you will need
#### 1. sample files in .fq format
#### 2. catalog_tags.tsv file derived from Stacks software, however many other file formats can be accomodated if necessary


######## IN UNIX ########
# prepare a folder containing all samples
cp /data/all_of_your_samples* /data/rawreads.fq


# find tags of each different motif type in the catalog_tags.tsv file
# can repeat this as necessary for different motif types
grep 'ATATATATAT\|TATATATATA' catalog_tags.tsv > catalog_AT_tags1.txt 
grep 'CTCTCTCTCT\|TCTCTCTCTC' catalog_tags.tsv > catalog_CT_tags1.txt 
grep 'ACACACACAC\|CACACACACA' catalog_tags.tsv > catalog_AC_tags1.txt 
grep 'GCGCGCGCGC\|CGCGCGCGCG' catalog_tags.tsv > catalog_GC_tags1.txt 
grep 'GAGAGAGAGA\|AGAGAGAGAG' catalog_tags.tsv > catalog_GA_tags1.txt 
grep 'GTGTGTGTGT\|GTGTGTGTGT' catalog_tags.tsv > catalog_GT_tags1.txt 

# alternative approach for tri-nucleotide repeats,
grep 'ATTATTATTATTATT\|AATAATAATAATAAT\|ATAATAATAATAATA\|TAATAATAATAATAA\|TTATTATTATTATTA' catalog_tags.tsv > catalog_ATT_tags.txt
grep 'ACCACCACCACCACC\|AACAACAACAACAAC\|ACAACAACAACAACA\|CAACAACAACAACAA\|CCACCACCACCACCA' catalog_tags.tsv > catalog_ACC_tags.txt
grep 'CTTCTTCTTCTTCTT\|CCTCCTCCTCCTCCT\|CTCCTCCTCCTCCTC\|TCCTCCTCCTCCTCC\|TTCTTCTTCTTCTTC' catalog_tags.tsv > catalog_CTT_tags.txt
grep 'GGCGGCGGCGGCGGC\|GCCGCCGCCGCCGCC\|GCGGCGGCGGCGGCG\|CGCCGCCGCCGCCGC\|CCGCCGCCGCCGCCG' catalog_tags.tsv > catalog_GGC_tags.txt
grep 'ATCATCATCATCATC\|ACTACTACTACTACT\|TCATCATCATCATCA\|TACTACTACTACTAC\|CATCATCATCATCAT\|CTACTACTACTACTA' catalog_tags.tsv > catalog_ATC_tags.txt
grep 'ATGATGATGATGATG\|GATGATGATGATGAT\|TAGTAGTAGTAGTAG\|AGTAGTAGTAGTAGT\|TGATGATGATGATGA\|GTAGTAGTAGTAGTA' catalog_tags.tsv > catalog_ATG_tags.txt
grep 'AAGAAGAAGAAGAAG\|AGAAGAAGAAGAAGA\|GAAGAAGAAGAAGAA\|AGGAGGAGGAGGAGG\|GGAGGAGGAGGAGGA\|GAGGAGGAGGAGGAG' catalog_tags.tsv > catalog_AAG_tags.txt
grep 'GCAGCAGCAGCAGCA\|CAGCAGCAGCAGCAG\|AGCAGCAGCAGCAGC\|CGACGACGACGACGA\|ACGACGACGACGACG\|GACGACGACGACGAC' catalog_tags.tsv > catalog_GCA_tags.txt


######## IN R ########

# combine all candidate tags into one file to permit easy filtering.
# here, we have retained various columns of infromation, but most importantly is
# the full tag sequence and number of samples in which it is present
combined_motifs <- rbind(count_AC_test[,c(4,5,9:10)],
		count_AT_test[,c(4,5,9:10)],
		count_CT_test[,c(4,5,9:10)],
		count_GA_test[,c(4,5,9:10)],
		count_GC_test[,c(4,5,9:10)],
		count_GT_test[,c(4,5,9:10)])

# bind numbers
count <- nchar(as.character(combined_motifs$V9))
combined_motifs <- cbind(combined_motifs,count)

# subset those with > n inividuals
combined_motifs_retained <- subset(combined_motifs, count>500)

# trimming:
# here we trip the last few bases which often have poorer quality
# you could also consider trimming any sequence that has SSR repeat type motifs at the beginning or end of the sequence
trimmed_seqs <- substring(combined_motifs_retained$V10, 1,79)
combined_motifs_retained$V10 <- trimmed_seqs

# generating the 'primers' file
# here we subset the first 15-25 charachters, allowing us to seach for these motifs in the raw reads.
# the more charachters we retain for the primer, the more specific it is, but it may not identify 
# sequences where the variable SSR is early in the read. Experimenting with this parameter is important.
primers <- substr(cols_needed_count_sub2[,4],1,25)
write.table(primers, "primers_for_agrep.txt")


######## IN UNIX ########
# use agrep to search all samples for specific tag
# for example... using these 'primers'...
./agrep/agrep -1 GGAGAGAACGCTGGAGTTCATTGGC sample* > /data/results/loc1.txt
./agrep/agrep -1 AAAAGAGATGGTGAGAGAATCATAC sample* > /data/results/loc2.txt
./agrep/agrep -1 ATCTCTCTCTCCAAAATGTGAGAGG sample* > /data/results/loc3.txt
./agrep/agrep -1 ATTTTCAAATCCAGAGAGAGAGCTT sample* > /data/results/loc4.txt
./agrep/agrep -1 CAACAGCAGCAGCAGCAGATTCAAC sample* > /data/results/loc5.txt
./agrep/agrep -1 GTATGTTCTTTTCTTTATGTTTTTG sample* > /data/results/loc6.txt


# the resulting files e.g. loc1.txt will contain every read for every individual with that primer sequence, 
# this then allows us to assess variation and data quality. We start with a function:

modifyFile <- function(in_table) {
  repeat_matches <- in_table
  repeat_matches <- repeat_matches["V2"]

  results2 <- list()
  rep <- grep("*", repeat_matches[1:length(repeat_matches$V2),]) #NA
  results2[rep] <- c("NA")

  #make it a factor
  results2[rep] <- as.numeric(repeat_matches$V2)
  
  #names <- substr(repeat_matches[1:length(repeat_matches$V2),], 1, 12)
  in_table[] <- lapply(in_table, as.character)
  names <- in_table$V1

  list <- unlist(results2)
  
  data <- cbind(names,list)
  data <- as.data.frame(data)
  names(data)
  
  # NOW WE NEED TO GET FACTOR LIST OF NAMES, AND ADD LENGTHS TO THAT LIST
  count <- count(data)
  
  # remove rows with three or less counts
  count_filtered<-count[!(count$freq<3),]
  
  # this takes the most common two reads
  count_ordered <- count_filtered[with(count_filtered, order(names, -freq)), ]
  obj <- split(count_ordered, with(count_ordered, names), drop = F)

  # Then we remove the second read if the difference in the number of reads
  # is greater than N.
  N <- 10
  filter2 <- lapply(obj,function(x)if(nrow(x)>1){
    if((x[1,3]/x[2,3])>N){
      x[1,]
    }else{
      x[1:2,] 
    }
  }else{
    x
  })
  
  # REPEAT with filter2
  lines <- lapply(filter2,function(x)x[1:2,])
  
  # this replaces na values in homozygotes
  filter3 <- lapply(lines,function(x)if(row.names(x[2,])=="NA"){
    x <- na.locf(x)
  }else{
    x 
  })
  
  # make into one table again
  combined <- do.call("rbind", filter3)
  
  # drop the frequency column
  combined_final  <- combined[,1:2]
  
  out_table <- combined_final
  # for info: # NA=Null allele
  # <NA> = missing data
  
  
  return(out_table)
}


##########################################
# to test a single file
##########################################

loc5test <- read.table("~/data/loc5.txt", quote="\"")
loc5test<-loc5test[,2]
modifyFile(loc5test)

##########################################
# to run on all loci files
##########################################

name_list #use ls() or list.files to make a list of file names

setwd("~/data/locs")
for (i in 1:nrow(name_list)) {
  
  file_name <- name_list[i,1]
  
  file_name <- as.character(file_name)
  
  file <- read.table(file_name)#, stringsAsFactors=F)
  
  out_file <- modifyFile(file)
  
  colnames(out_file) <- c("IND", file_name)
  
  out_name <- paste(file_name, "csv", sep=".")
  
 write.table(out_file, out_name, col.names=T,row.names=FALSE,sep=",",quote=FALSE)
}

# then combine together
files <- list.files("~/locs_out")
combined_files <- lapply(files, read.csv)

# minor step to make all the names unique
for(i in 1:length(data)){
seqs <- rep(1:2,(length(data[[i]]$IND)/2))
names <- data[[i]]$IND
names.1.2 <- paste(names,seqs,sep=".")
data[[i]]$IND <- c(names.1.2)
}


# create a list of names
conc <- unique(c(as.character(data[[1]]$IND), #need to fix this to get all names
                 as.character(data[[2]]$IND),
                 as.character(data[[3]]$IND),
                 as.character(data[[4]]$IND),
                 as.character(data[[5]]$IND),
                 as.character(data[[6]]$IND),
                 as.character(data[[7]]$IND),
                 as.character(data[[8]]$IND),
                 as.character(data[[9]]$IND),
                 as.character(data[[10]]$IND),
                 as.character(data[[11]]$IND),
                 as.character(data[[12]]$IND),
                 as.character(data[[13]]$IND),
                 as.character(data[[14]]$IND),
                 as.character(data[[15]]$IND),
                 as.character(data[[16]]$IND)
                 ))
				 
#conc <- unique(c(as.character(data[[1]]$IND), as.character(data[[2]]$IND), as.character(data[[3]]$IND),as.character(data[[4]]$IND),as.character(data[[5]]$IND)))
conc <- as.data.frame(conc)
colnames(conc) <- c("IND")

out <- list()
for(i in 1:length(data)){
merged <- merge(conc,data[i], by='IND',all=T)
subset1 <- merged[!merged$IND == "NA.1", ]
subset2 <- subset1[!subset1$IND == "NA.2", ]
out[[i]] <- subset2
}


df <- do.call("cbind", out)
name_col <- df[,1]
rep_data <- df[, seq(2, ncol(df), by = 2)]
final_data <- cbind(name_col,rep_data)

write.table(final_data, "~/Data_out/final_table_of_factors.csv", col.names=T,row.names=FALSE,sep=",",quote=FALSE)



#####################################################
# then various ascertainment biases can be applied
#####################################################


# for example 

# count NAs in column
na_count <- sapply(final_data, function(x) sum(length(which(is.na(x)))))
na_count <- data.frame(na_count)
vec <- na_count < 180
length(which(vec==T))
cols_keep <- which(vec==T)

final_data2 <- subset(final_data,,cols_keep)
ncol(final_data2)


# OR

# Thinking about number of different allele values, 
# if this is very high it suggests a sequencing error
unique <- rapply(final_data2, function(x) sum(length(unique(x))))
unique <- data.frame(unique)
vec <- unique >30
length(which(vec==T))
summary(unique)
