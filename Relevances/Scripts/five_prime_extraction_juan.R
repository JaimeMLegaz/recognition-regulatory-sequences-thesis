library(tidyverse)
library(rtracklayer)

# GTF file can be dowloaded from Ensembl: https://www.ensembl.org/info/data/ftp/index.html (This page links to the latest version of GTF ensembl). We use rtracklayer import function to read the gtf file into a proper human readable format
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.94.gtf")
# keep only standard chromosomes (1:22,X,Y,MT)
gtf = keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode="coarse")
gtf.df=as.data.frame(gtf,stringsAsFactor=F)

# Select only protein-coding genes and transcripts
gtf.df.pc = gtf.df %>% dplyr::filter(gene_biotype %in% "protein_coding", transcript_biotype %in% "protein_coding")

# Select TSL level 1
# TSL coulumn has scoring for the reliability of the transcripts based on multiple evidence. I select only the highly confident transcripts
gtf.df.pc$transcript_support_level = gsub("\\s*\\([^\\)]+\\)","",as.numeric(gtf.df.pc$transcript_support_level))
gtf.df.pc.tsl1 = gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)

# extract 5' UTRs
five_prime = gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "five_prime_utr")

# Collapsing the 5'UTRs among the transcripts for each gene
# The data in the file is at transcript level rather than gene. So I merge multiple transcripts of a gene into one.

# this makes a list of Granges objects, with each element is Granges for a gene and contains data about it's transcripts
five_prime_grList = makeGRangesListFromDataFrame(five_prime, split.field = "gene_id", names.field = "transcript_id")
# we merge the transcript coordinates using the reduce function.
five_prime_collapse = IRanges::reduce(five_prime_grList, with.revmap=TRUE) %>% as.data.frame() %>%
  mutate(elements_collapsed = lengths(revmap), five_prime_utr_id = paste(group_name,seqnames, start, end, strand, elements_collapsed, sep=":"))




###############################################
###### Extracting DNA sequence ################
###############################################

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
options(stringsAsFactors = FALSE)

# converting to Granges object
five_prime_gr = makeGRangesFromDataFrame(five_prime_collapse, keep.extra.columns = TRUE)

# adding 'chr' in front of seqnames
newStyle <- mapSeqlevels(seqlevels(five_prime_gr), "UCSC")
five_prime_gr <- renameSeqlevels(five_prime_gr, newStyle)

# Extract sequences using the package BSgenome
data = data.frame(id = five_prime_gr$five_prime_utr_id, 
                  X5utr = BSgenome::getSeq(Hsapiens, five_prime_gr, as.character = TRUE)
                  )

write.csv(data,file="Data/UTR/PositiveUTRSequences_untreated", row.names = FALSE, quote = FALSE)
write.csv(NegSequences,file="Data/UTR/NegativeUTRSequences_untreated", row.names = FALSE, quote = FALSE)

## FALTA SABER CUÁLES SON POSITIVOS Y CUÁLES NEGATIVOS









