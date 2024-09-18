library(stringr)
library(Biostrings)
library(UniProt.ws)
# Benjamin Parker, PhD
# 8.1.24
# Retreive and tile proteins for phage libraries based on GOTerms.
# Tiles have a specifiable length and overlap between each.
# Uses a modified FASTA of the disorder score of the human proteome to exclude tiles which are too ordered.
# Requires for input (see variables below).
#   1. Codon table as a csv; correlate amino acid with codon. 
#   2. TSV of aiupred-scored human proteome; each amino acid is replaced with a number from 0 (ordered) to 9 (disordered).
#

# Codon table for reverse translation.
codon_table<-read.csv('/home/weisslab/Documents/scriptcore/codon_tables/e_coli.csv')

# AIpred scored human proteome
aipred_raw<-read.csv('/home/weisslab/Documents/scriptcore/human_proteome/hprot_aiupred.tsv',sep='\t',header=FALSE,colClasses="character")

# Set working directory
setwd('/home/weisslab/Documents/Ben/Notebook/240821-tile_phage')

# Tile length in aa
tile_size=24

#Tile overlap in aa
tile_overlap=12

# AIupredcutoff
aiupred_cutoff=0.4

# Restriction site to remove
rsite_to_remove='AAGCTT' #HindIII for now
# Sequence to replace +0 through +2 frame with
rsite_replace0='AAACTG'
rsite_replace1='ATCTTT'
rsite_replace2='AAGCCT'
# CUSTOM GOTERM FILTER (based on Eric discussion, 10.5.23)
ndr_goterms=c(
  # "GO:0032465 AND organism_id:9606"#, # regulate cytokinesis
  # "(GO:0032465 OR GO:0035329) AND organism_id:9606"# # regulation cytokinesis and hippo pathway
 # "(GO:0000910 OR GO:0035329) AND organism_id:9606"#, # cytokinesis and hippo pathway
  # "GO:0000902 AND organism_id:9606"#, # cell morphogenesis
  # GO:0006893  AND organism_id:9606" # Section+golgi to PM *****************
   # "GO:0046903 AND organism_id:9606" # Section 8575 phage
  # "GO:0005856 AND organism_id:9606" # Cytoskeleton
  # "GO:0003924  AND organism_id:9606"  # GTPase activity #ERIC WANTS THIS!!!!!!!!!
  # " GO:0051020  AND organism_id:9606"  # GTPase binding 
  
  # cytokinesis and hippo pathway #nd GTPase activity 
   "(GO:0000910 OR GO:0035329 OR GO:0003924 ) AND organism_id:9606" #5394 phage
  
 # " GO:0007010  AND organism_id:9606"  # Cytoskeletal organization #ERIC WANTS THIS!!!!!!!!!
 # " GO:0008092  AND organism_id:9606"  # Cytoskeleton protein binding 
 
  # "(GO:0030011 OR GO:0030010) AND organism_id:9606" #maintenance and establishment of cell polarity
  
)
###########################################################################################################3
#AAGCCT

aipred <- as.data.frame(aipred_raw[!aipred_raw$V1=="",]$V1)
aipred$seqscore <-aipred_raw[!aipred_raw$V1==""|!aipred_raw$V2=="",]$V2
names(aipred) <-c("accession","seqscore")

# Reference human proteome
hprot_raw<-read.csv('/home/weisslab/Documents/scriptcore/human_proteome/hprot.tsv',sep='\t',header=FALSE,colClasses="character")
names(hprot_raw) <-c("accession","sequence")


#ndr_goterms=c("GO:0032465 AND organism_id:9606")

# Secretion GO:0046903
# Hippo 
# Cytoskeleton GO:0005856
# Combine STK38 interactors and Mob2 interactors

# Determine sequences/uniprot ids in 4x chunks
# Note limit of 100k entries per query
# TEST <-mapUniProt(
#   from="UniProtKB_AC-ID",
#   to="UniProtKB-Swiss-Prot",
#   query=MergeTide_ToMap[1:length(MergeTide_ToMap)],
#   columns=c("accession","sequence","id","cc_subcellular_location","go_id")
# )

Sequences_to_tile_raw <-queryUniProt(
  query=ndr_goterms,
  fields=c("accession"),
  collapse=" OR "

)
#Sequences_to_tile_raw <-queryUniProt(
#  query=c("GO:0000902 AND organism_id:9606"),
#  fields=c("accession")
  
#)

names(Sequences_to_tile_raw) <-c("accession")

# Merge with "canonical" human proteome.
Sequences_to_tile <- merge(Sequences_to_tile_raw,hprot_raw)
# Merge sequences withaipred of human proteome.
Sequences_to_tile <- merge(Sequences_to_tile,aipred)






# Empty array of tile sequences.
tile_sequence <- c()
# Empty array of tileaipred scores
tile_disoscore <- c()
# Empty array of tile names.
tile_name <- c()
# "Pointer" to keep track of where in the protein sequence we are.
na_aiupred_tiles=0 #number of tiles aiupred doesn't score.
total_tiles=0
tiles_without_disoscores=0
# Iterate across proteins, tiling each one.
for(protein in 1:nrow(Sequences_to_tile)){
  # print(protein)

# Master sequence to break into tiles.
SEQUENCE<-Sequences_to_tile[protein,]$sequence
DISORDER<-Sequences_to_tile[protein,]$seqscore
# Even tiles are tiles with the proper length. Remainder is if not divisible
# by tile length.

# "Pointer" used to iteratively calculate where in the protein sequence we are and what the tile limits should be.
pointer=0
tile_number=1
while(pointer<str_length(SEQUENCE)){
  # First tile. If pointer is at beginning, don't tile.
  if(pointer-tile_overlap<0){
    pointer=pointer+tile_size+1
  }else if(pointer-tile_overlap+tile_size<=str_length(SEQUENCE)){
    # If the pointer hasn't gone by the end of the protein sequence, it is divisible.
    # If it has, create a separate sequence which just represents the end of the peptide without tiling.
   
      pointer=pointer-tile_overlap+tile_size
  # Otherwise, you have reached the end.
  }else{
    pointer=str_length(SEQUENCE)+1
    
  }
  # print(pointer)


  
  # tile disorder score along tile:
  tile_disorder_along=substr(DISORDER,pointer-tile_size,pointer-1)
  
  tile_disorder=mean(as.numeric(str_split_1(tile_disorder_along,"")))
  if(DISORDER==""){
    na_aiupred_tiles=na_aiupred_tiles+1
  }
  total_tiles=total_tiles+1
  # Output name , disorderand sequence if sufficiently disordered
  # Remove phage which are below aiupred_cutoff.
  if(tile_disorder>=aiupred_cutoff*10 && !is.na(tile_disorder)){
    tile_name<-c(tile_name,paste(Sequences_to_tile[protein,]$accession,tile_number,sep="_"))
    tile_sequence<-c(tile_sequence,substr(SEQUENCE,pointer-tile_size,pointer-1))
    tile_disoscore<-c(tile_disoscore,tile_disorder)
  

  }else if(is.na(tile_disorder)){
    tiles_without_disoscores=tiles_without_disoscores+1
  }
  # Increment tile number
  tile_number=tile_number+1

  # Increment tile number for naming phage sequences
    # # Output name and sequence
    # tile_name<-c(tile_name,paste(Sequences_to_tile[protein,]$accession,tile_number,sep="_"))
    # tile_sequence<-c(tile_sequence,substr(SEQUENCE,pointer-tile_size,pointer-1))
    # tile_disoscore<-c(tile_disoscore,tile_disorder)

  
}
}


tiled_phage_protein_withdups<-as.data.frame(tile_name)
tiled_phage_protein_withdups$sequence <-tile_sequence
tiled_phage_protein_withdups$disoscore<-tile_disoscore
# Remove duplicate phage.
tiled_phage_protein <- tiled_phage_protein_withdups[!duplicated(tiled_phage_protein_withdups$sequence),]
#duplicated(tiled_phage_protein$sequence)




# tiled_phage$seqscore <-tile_disoscore
filename=gsub(" ","_",ndr_goterms)
filename=gsub(":","",filename)
filename=gsub("\\(","",filename)
filename=gsub("\\)","",filename)
write.table(tiled_phage_protein,paste(filename,"protein.tsv",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

# For DNA sequences.
tile_sequence_dna<-c()
# Convert to DNA.
for(protein in 1:nrow(tiled_phage_protein)){
  
  # Replace AA with codons from codon table.
  dna_sequence <- c()
  for(amino_acid_index in 1:str_length(tiled_phage_protein$sequence[protein])){
    
    # Which AA to replace?
    amino_acid_to_replace=substr(tiled_phage_protein$sequence[protein],amino_acid_index,amino_acid_index)
    # Where in the codon table is it?
    codon_index=which(amino_acid_to_replace==codon_table$amino_acid)
    # print(codon_index)
    dna_sequence <- c(dna_sequence,codon_table[codon_index,1])
  
  }

  # Collapse into single DNA sequence.
  dna_sequence <-paste(dna_sequence,collapse="")
  #print(dna_sequence)
  tile_sequence_dna<-c(tile_sequence_dna,dna_sequence)
   
}
print(protein)
tiled_phage_dna<-as.data.frame(tiled_phage_protein$tile_name)
tiled_phage_dna$sequence <-tile_sequence_dna
names(tiled_phage_dna) <-c('tile_name','dna_sequence')

# Find restriction sites to remove.
# contains_rsite<- tiled_phage_dna[grepl(rsite_to_remove,tiled_phage_dna$sequence),]
# Will contain edited sequences of contains_rsite.
rsite_replaced_phage<-tiled_phage_dna
# Test to make sure restriciton site swapped sequences still translate
TEST<-c()
for(index_with_rsite in 1:nrow(rsite_replaced_phage)){
  # Sequence of phage DNA with restriction site.
  seq_with_rsite=rsite_replaced_phage$dna_sequence[index_with_rsite]
  # print(seq_with_rsite)
  # Where does the restriction site start?
  start_of_rsite=as.numeric(str_locate(seq_with_rsite,rsite_to_remove)[,1])-1
  # print(start_of_rsite)
  # Case if frame = 0
  if(!is.na(start_of_rsite)){
    if((start_of_rsite)%%3==0)
    {
      # Frame is +0
      rsite_replaced_phage$dna_sequence[index_with_rsite]=gsub(rsite_to_remove,rsite_replace0,seq_with_rsite) 
      # print(sub("AAGCTT","AAACTT",seq_with_rsite) )
      print(index_with_rsite)
      
    }else if((start_of_rsite+1)%%3==0)
    {
      # Frame is +1
      rsite_replaced_phage$dna_sequence[index_with_rsite]=gsub(rsite_to_remove,rsite_replace1,seq_with_rsite) 
      #print(sub("AAGCTT","AAACTT",seq_with_rsite) )
     TEST<-c(TEST,index_with_rsite)
    }else if((start_of_rsite+2)%%3==0)
    {
      # Frame is +2
      rsite_replaced_phage$dna_sequence[index_with_rsite]=gsub(rsite_to_remove,rsite_replace2,seq_with_rsite) 

    }
  }

}
names(rsite_replaced_phage) <-c('tile_name','sequence')

write.table(rsite_replaced_phage,paste(filename,"dna.tsv",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

# Convert to DNA string
rsite_replaced_phage_asdna <- as.data.frame(rsite_replaced_phage$tile_name)
rsite_replaced_phage_asdna$sequence <-DNAString(rsite_replaced_phage$sequence)


TEST <-rsite_replaced_phage[!rsite_replaced_phage$sequence %in% tiled_phage_dna$sequence,]

# +0 frame
# AAG CTT		
# Convert to AAA CTT
# +1 frame
# nAA GCT T	Convert to AAACTG
# GCT (ala) is in frame; convert GCC (ala)
# nAA GCC T
#   aAA	K	-> a AG
#   cAA	Q -> c AG
#   gAA E -> g AG
#   tAA * not used
# +2 frame Convert to ATCTTT
# A AGC TTn
# AGC (ser) is in frame; convert to TCT (ser)
# A TCT TTn	
#   TTc: F
#   TTt: F
#   TTa: L
#   TTg: L
# 
#   

# 
# phage_ord_vec <- rep(NA,times=length(tiled_phage$sequence))
# for(i in 1:length(tiled_phage$sequence)){
#   phage_vector <- str_locate(Sequences_to_tile$sequence[i],tiled_phage$sequence[i])
#   phage_ord_seq_raw<-strsplit(substr(MergeTideS$seqdisoscore[i],phage_vector[1],phage_vector[2]),split="")
#   phage_ord_seq <- as.numeric(unlist(phage_ord_seq_raw))
#   phage_ord_vec[i] <-mean(phage_ord_seq)
# }
