library(dplyr)
library(plyr)
library(ggplot2)
# library(idpr)
library(stringr)
library(UniProt.ws)
library(tidyr)
library("data.table")
# library(textreadr)
# compare-peptides.R

# Benjamin Parker, PhD
# Graphs of phage

# phage_exp_id = "B022525"
# phage_exp_id = "B011525"
# phage_exp_id = "B102523"
# phage_exp_id = "B101123"
# phage_exp_id = "B072423"
# phage_exp_id = "B091218"
phage_exp_id = "B080218"
# setwd(paste('/media/weisslab/Deep_Seq/Phage Display/B022525-PE150/15nt_fwd'))
# setwd(paste('/media/weisslab/Deep_Seq/Phage Display/B011525-PE150/15nt_fwd'))
# setwd(paste('/media/weisslab/Deep_Seq/Phage Display/B102523-PE150/30nt'))

# setwd(paste('/media/weisslab/Deep_Seq/Phage Display/B101123-PE150/30nt'))
# setwd(paste('/media/weisslab/Deep_Seq/Phage Display/B072423-T7-NM2-LM1-GST')) 
# setwd("/media/weisslab/Deep_Seq/Phage Display/230306-new_work/30nt/B091218_dataset")
setwd("/media/weisslab/Deep_Seq/Phage Display/230306-new_work/30nt/B080218_dataset")

# # # Phage librarys PL1 and later
# t7_reference <-read.csv('/home/weisslab/Documents/scriptcore/pl1_all_peptides.tsv',sep='\t',header=FALSE)
# names(t7_reference)<-c('name','sequence')
# t7_reference$accession <-gsub('_[0-9]*','',t7_reference$name)
# 

# # # Elledge library
t7_reference <-read.csv('/home/weisslab/Documents/scriptcore/t7_pep_all_peptides.csv',sep=',',header=FALSE)
names(t7_reference)<-c('name','sequence','accession','protsequence')


# Load in table from Zabin/Moses with the motifs/ereferences.
reference_property_table<-read.csv2("~/Documents/Ben/Notebook/250512-Phage_high_dimension_test_again/elife-46883-table.csv",sep=',')


# Protein code
# N=Ndr-Mob2
# # L=Lats-Mob1
# protein_code="L"
 protein_code="N"
enrich_stage=3
##################### GUTS!! ###########################
###################################!!###################################!!#################################





# load reference library
reflib<- read.csv("Input_processing/Library.csv")


# Load panned
N1<- read.csv(paste(protein_code,enrich_stage,"_processing/",protein_code,enrich_stage,"-1.csv",sep=""))
names(N1) <-c("name","n1")
N2<- read.csv(paste(protein_code,enrich_stage,"_processing/",protein_code,enrich_stage,"-2.csv",sep=""))
names(N2) <-c("name","n2")
N3<- read.csv(paste(protein_code,enrich_stage,"_processing/",protein_code,enrich_stage,"-3.csv",sep=""))
names(N3) <-c("name","n3")

# Load control librarys
G1<- read.csv(paste("G",enrich_stage,"_processing/G",enrich_stage,"-1.csv",sep=""))
names(G1) <-c("name","g1")
G2<- read.csv(paste("G",enrich_stage,"_processing/G",enrich_stage,"-2.csv",sep=""))
names(G2) <-c("name","g2")
G3<- read.csv(paste("G",enrich_stage,"_processing/G",enrich_stage,"-3.csv",sep=""))
names(G3) <-c("name","g3")





# merge libraries together
N <- rbind.fill(
  N1,
  N2,
  N3,
  reflib)
G <- rbind.fill(
  G1,
  G2,
  G3,
  reflib)

# Replace NA with 0 for enrichment
N[is.na(N)] <-  0
G[is.na(G)] <-  0

# Sum together identical peptides into a separate data frame.
Ntide_raw <- N %>% group_by(name) %>% summarize_at(c("n1","n2","n3","n"),sum)
names(Ntide_raw) <-c("name","n1","n2","n3","input")

Gtide_raw <- G %>% group_by(name) %>% summarize_at(c("g1","g2","g3","n"),sum)
names(Gtide_raw) <-c("name","g1","g2","g3","input")




# 7.4.24 Moved the NA/0 to 1/31 to below Ntide/Gtide

Ntide_raw[Ntide_raw==0] <- 0
Gtide_raw[Gtide_raw==0] <-  0



Ntide<-Ntide_raw
Gtide<-Gtide_raw



names(Ntide) <-c("name","n1","n2","n3","input")
names(Gtide) <-c("name","g1","g2","g3","input")

# Map refseq ids in MergeTideE to uniprot ids
t7_uniprots <-as.data.frame(t7_reference$name)
t7_uniprots$accession <-t7_reference$accession
t7_uniprots$sequence <-t7_reference$sequence
names(t7_uniprots) <- c('name','accession','sequence')

#################### HIGH DIMENSIONALITY
#################### Iput dataset
# Merge input dataset (ItideE) with output, then find the phage sequences.
# Ntide <-Gtide
MOD_SOURCE1<-as.data.frame(Ntide$name)
MOD_SOURCE1$output<-apply(Ntide[,2:4],1,median) #4.28.25 mess with
MOD_SOURCE1$variance<-apply(Ntide[,2:4],1,var) #4.28.25 mess with
MOD_SOURCE1$output_difference<-apply(Ntide[,2:4],1,max)-apply(Ntide[,2:4],1,median)#4.28.25 mess with



names(MOD_SOURCE1) <-c('name','median','variance','output_difference')
MOD_SOURCE1<-merge(MOD_SOURCE1,reflib,by='name',all=TRUE)
# MOD_SOURCE1<-merge(MOD_SOURCE1,t7_uniprots,by='name',all=TRUE)
MOD_SOURCE1<-merge(MOD_SOURCE1,t7_uniprots,by='name')
names(MOD_SOURCE1) <-c('name','median','variance','output_difference','input','accession','sequence')



# ############ Pull sequence from uniprot; break into 100K entries since that is max allowed ##################
# MOD_SOURCE1_ToMap_raw <- merge(MOD_SOURCE1,t7_uniprots)$accession
# ### Retrive ALL phage reference data.
# MOD_SOURCE1_ToMap_raw <-t7_reference$accession
# ###
# MOD_SOURCE1_ToMap <- MOD_SOURCE1_ToMap_raw[!MOD_SOURCE1_ToMap_raw == ""]
# 
# MaxUnis=100000 #Block length of maximum uniprots to pull using mapUniProt at one time.
# MOD_SOURCE1MapLength=length(MOD_SOURCE1_ToMap)
# MaxBlocks=MOD_SOURCE1MapLength%/%MaxUnis
# 
# # End of block sequence for Uniprot pulling
# EndOfBlockSequence=0
# StartBlock=1 #Start here 
# MOD_SOURCE1_AccID<-data.frame(matrix(ncol=6,nrow=0)) #empty data frame with Uniprot entries
# names(MOD_SOURCE1_AccID) <-c('accession','accessionTo','protsequence','protname',"location","goid")
# 
# 
# while(EndOfBlockSequence==0){
#   # StartBlock=UniBlock*MaxUnis+1 # Starting range for library
#   if(StartBlock+MaxUnis>MOD_SOURCE1MapLength){ # Have we reached the end of the array?
#     EndBlock=MOD_SOURCE1MapLength
#     # If you've reached the end of the array, don't do any more loops.
#     # print(EndBlock)
#     EndOfBlockSequence=1
#     # Otherwise, iterate over the next block.
#   }else{
#     EndBlock=StartBlock+MaxUnis-1
#   }
#   
#   # Pull block from uniprot
#   MOD_SOURCE1_AccID_block <-mapUniProt(
#     from="UniProtKB_AC-ID",
#     to="UniProtKB-Swiss-Prot",
#     query=MOD_SOURCE1_ToMap[StartBlock:EndBlock],
#     columns=c("accession","sequence","id","cc_subcellular_location","go_id")
#   )
#   
#   names(MOD_SOURCE1_AccID_block) <-c('accession','accessionTo','protsequence','protname',"location","goid")
#   print(length(MOD_SOURCE1_AccID_block$accession))
#   
#   # Iteratively merge pulled uniprot entries with existant data frame.
#   MOD_SOURCE1_AccID<-rbind(MOD_SOURCE1_AccID,MOD_SOURCE1_AccID_block)
#   # Move startblock again.
#   StartBlock=StartBlock+MaxUnis
#   
#   
# }
# 
# names(MOD_SOURCE1_AccID) <-c('accession','accessionTo','protsequence','protname',"location","goid")
# # # 
# 
# 
# 
# MOD_SOURCE1<-merge(MOD_SOURCE1,MOD_SOURCE1_AccID,by="accession")
# # RAWLIBRARY<-merge(t7_reference,MergeTide_AccID,by="accession")
# 
# names(MOD_SOURCE1) <-c('accession','name','median','variance','output_difference','input','sequence','accession','protseq','protname','go1')
# 
# 
# 
# MOD_SOURCE1_ID<-as.data.frame(MOD_SOURCE1_ID_Acc$Entry)
# MOD_SOURCE1_ID$prosteq<-MOD_SOURCE1_ID_Acc$Sequence
# names(MOD_SOURCE1_ID) <-c('accession','protseq')
# MOD_SOURCE1<-merge(MOD_SOURCE1,MOD_SOURCE1_ID,by='accession')
# 
############ Rest of the initial processing ##################
MOD_SOURCE1[is.na(MOD_SOURCE1)] <-0 #1/263123573
MOD_SOURCE1$output_ppm<-MOD_SOURCE1$median/sum(MOD_SOURCE1$median)
MOD_SOURCE1$input_ppm<-MOD_SOURCE1$input/sum(MOD_SOURCE1$input)
#Calculate PPM.

#### For curation of peptide  properties
HMOD_PROP <-data.frame(matrix(ncol=3,nrow=0))
names(HMOD_PROP)<-c('property','input_value','output_value')


##################### HIGH DIMENSIONALITY MODULES FROM Zarin et al.
#################### Module 1: S content #################### 
### Ser content from 


HDMOD1_AA_S <-c(
  sum(str_count(MOD_SOURCE1$sequence,"S")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"S")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[1,] <-c('AA_S',HDMOD1_AA_S[1],HDMOD1_AA_S[2])



#################### Module 2: P content #################### 
### Pro content from Zarin.


HDMOD2_AA_P <-c(
  sum(str_count(MOD_SOURCE1$sequence,"P")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"P")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[2,] <-c('AA_P',HDMOD2_AA_P[1],HDMOD2_AA_P[2])

#################### Module 3: T content #################### 
### Thrr content from Zarin.


HDMOD3_AA_T <-c(
  sum(str_count(MOD_SOURCE1$sequence,"T")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"T")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[3,] <-c('AA_T',HDMOD3_AA_T[1],HDMOD3_AA_T[2])


#################### Module 4: A content #################### 
### Ala content from Zarin.


HDMOD4_AA_A <-c(
  sum(str_count(MOD_SOURCE1$sequence,"A")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"A")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[4,] <-c('AA_A',HDMOD4_AA_A[1],HDMOD4_AA_A[2])



#################### Module 5: H content #################### 
### His content from Zarin.


HDMOD5_AA_H <-c(
  sum(str_count(MOD_SOURCE1$sequence,"H")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"H")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[5,] <-c('AA_H',HDMOD5_AA_H[1],HDMOD5_AA_H[2])



#################### Module 6: Q content #################### 
### Gln content from Zarin.


HDMOD6_AA_Q <-c(
  sum(str_count(MOD_SOURCE1$sequence,"Q")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"Q")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[6,] <-c('AA_Q',HDMOD6_AA_Q[1],HDMOD6_AA_Q[2])



#################### Module 7: N content #################### 
### Asn content from Zarin.


HDMOD7_AA_N <-c(
  sum(str_count(MOD_SOURCE1$sequence,"N")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"N")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[7,] <-c('AA_N',HDMOD7_AA_N[1],HDMOD7_AA_N[2])


#################### Module 8: G content #################### 
### Gly content from Zarin.


HDMOD8_AA_G <-c(
  sum(str_count(MOD_SOURCE1$sequence,"G")*MOD_SOURCE1$input_ppm),
  
  sum(str_count(MOD_SOURCE1$sequence,"G")*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[8,] <-c('AA_G',HDMOD8_AA_G[1],HDMOD8_AA_G[2])



#################### Module 9: Kappa ####################
### Separation of positive vs negative charge residues.

# Main equations:
#
#
#
#
# Blob size
g=6 # Das and pappu calculates with g=5 and 6, but g=6 is evenly compatable with phage squence of 24aa.

# function to calculate charge fraction sigma
calculate_sigma <- function(sequence){
  pos_charges<-str_count(sequence,'[KR]')/str_length(sequence)
  neg_charges<-str_count(sequence,'[DE]')/str_length(sequence)

  sigma_output<-(pos_charges-neg_charges)^2/(pos_charges+neg_charges)
  # print(sigma_output)
  if(!is.na(sigma_output)){
    return(sigma_output)
  }else{
    return(0)
  }

}
# Return the "ideal" sequence needed for delta max. Charged residues
# are maximally separated in this "ideal" peptide.
return_max_seq <- function(sequence){
  pos_charges<-str_count(sequence,'[KR]')
  neg_charges<-str_count(sequence,'[DE]')
  remain_aa_num<-str_length(sequence)+-pos_charges+-neg_charges
  pos_seq<-paste(rep("R",each=pos_charges),collapse='')
  neg_seq<-paste(rep("E",each=neg_charges),collapse='')
  remain_seq<-paste(rep("X",each=remain_aa_num),collapse='')
  max_seq <-paste(pos_seq,remain_seq,neg_seq,sep="")

  return(max_seq)

}

# return_max_seq("KKPEKEPKEKEK")
# Delta list of phages.


delta <- c()
delta_max <- c()

# Interate through phage in tiles of blob size g. Calculate sigma.
for(i in 1:length(MOD_SOURCE1$sequence)){
# for(i in 629){                     #DEBUG
  # Phage squence of interest.
  phage_to_calculate <-MOD_SOURCE1$sequence[i]
   # Maximum charge separation for phage. Needed to calculate delta max.
  max_phage_to_calculate<-return_max_seq(phage_to_calculate)
  # phage_to_calculate <-TEST                     #DEBUG
  # Total sigma for peptide.
  sigma_total<-calculate_sigma(phage_to_calculate)

  # Create phage to calculate f
  # print(phage_to_calculate)#DEBUG



  # List of sigma values to sum accorss blobs.Top of the delta equation
  sigma_blob_list<-c()
  # Calculate total number of blobs.
  # total_blobs=0 #not needed; can just use sigma_blob_list length.
  # for (j in 1:str_length(phage_to_calculate)){
  ########### Calculate delta.   ###########
  j<-1
  while(j<=str_length(phage_to_calculate)){

    phage_blob_j<-substr(phage_to_calculate,j,j+g-1)
    sigma_blob_j<-calculate_sigma(phage_blob_j)
    sigma_blob_list<-c(sigma_blob_list,(sigma_blob_j-sigma_total)^2)

    # Increment J accorss phage.

    j<-j+g
    # print(phage_blob_j)#DEBUG
    # print(sigma_blob_list)#DEBUG
    # print(j+g)#DEBUG

  }
  # print(sigma_blob_list)#DEBUG


  delta_blob_j<-mean(sigma_blob_list)
  delta<-c(delta,delta_blob_j)
  # print(i)#DEBUG


  ########### Calculate delta_max   ###########
  sigma_blob_list_max<-c()
  k<-1
  while(k<=str_length(max_phage_to_calculate)){

    phage_blob_k<-substr(max_phage_to_calculate,k,k+g-1)
    sigma_blob_k<-calculate_sigma(phage_blob_k)
    sigma_blob_list_max<-c(sigma_blob_list_max,(sigma_blob_k-sigma_total)^2)

    # Increment J accorss phage.

    k<-k+g
    # print(phage_blob_k)#DEBUG
    # print(sigma_blob_list_max)#DEBUG
    # print(k+g)#DEBUG
  }
  # print(sigma_blob_list)#DEBUG


  delta_blob_k_max<-mean(sigma_blob_list_max)
  delta_max<-c(delta_max,delta_blob_k_max)
  # print(delta_blob_k_max)#DEBUG
  # print(i)#DEBUG
  # print(delta_blob_j)#DEBUG

}


# Return this as a data frame (local) to ensure the calculations are done right.
MOD_SOURCE_kappa <- MOD_SOURCE1
MOD_SOURCE_kappa$delta <- delta
MOD_SOURCE_kappa$delta_max <- delta_max
MOD_SOURCE_kappa$kappa <-MOD_SOURCE_kappa$delta/MOD_SOURCE_kappa$delta_max
MOD_SOURCE_kappa$kappa[is.na(MOD_SOURCE_kappa$kappa)]<-0


# Kappa weighted for phage abundance.
HDMOD9_kappa <-c(
  sum(MOD_SOURCE_kappa$kappa*MOD_SOURCE1$input_ppm),
  sum(MOD_SOURCE_kappa$kappa*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[9,] <-c('AA_kappa',HDMOD9_kappa[1],HDMOD9_kappa[2])


# return_max_seq("PVPASMFAPEPSSPGAARAAAAAA") #DEBUG
# calculate_sigma("PVPASMFAPEPSSPGAARAAAAAA") #DEBUG
# calculate_sigma("RXXXXXXXXXXXXXXXXXXXXXXE") #DEBUG




# #################### Module 10: Omega####################
# ### Separation of KDEP vs non-KDEP.e charge residues. See module 9.
# # See Martin 2016.

# Main equations:
#
#
#KRDEP
#
# Blob size
# Das and pappu calculates with g=5 and 6, but g=6 is evenly compatable with phage squence of 24aa.
g=6 # Das and pappu calculates with g=5 and 6, but g=6 is evenly compatable with phage squence of 24aa.

# function to calculate charge fraction sigma
calculate_sigma <- function(sequence){
  pos_charges<-str_count(sequence,'[KRDEP]')/str_length(sequence)
  neg_charges<-str_count(sequence,'[^KRDEP]')/str_length(sequence)

  sigma_output<-(pos_charges-neg_charges)^2/(pos_charges+neg_charges)
  # print(sigma_output)
  if(!is.na(sigma_output)){
    return(sigma_output)
  }else{
    return(0)
  }

}
# Return the "ideal" sequence needed for delta max. Charged residues
# are maximally separated in this "ideal" peptide.
return_max_seq <- function(sequence){
  pos_charges<-str_count(sequence,'[KR]')
  neg_charges<-str_count(sequence,'[DE]')
  remain_aa_num<-str_length(sequence)+-pos_charges+-neg_charges
  pos_seq<-paste(rep("R",each=pos_charges),collapse='')
  neg_seq<-paste(rep("E",each=neg_charges),collapse='')
  remain_seq<-paste(rep("X",each=remain_aa_num),collapse='')
  max_seq <-paste(pos_seq,remain_seq,neg_seq,sep="")

  return(max_seq)

}

# return_max_seq("KKPEKEPKEKEK")
# Delta list of phages.


delta <- c()
delta_max <- c()

# Interate through phage in tiles of blob size g. Calculate sigma.
for(i in 1:length(MOD_SOURCE1$sequence)){
  # for(i in 629){                     #DEBUG
  # Phage squence of interest.
  phage_to_calculate <-MOD_SOURCE1$sequence[i]
  # Maximum charge separation for phage. Needed to calculate delta max.
  max_phage_to_calculate<-return_max_seq(phage_to_calculate)
  # phage_to_calculate <-TEST                     #DEBUG
  # Total sigma for peptide.
  sigma_total<-calculate_sigma(phage_to_calculate)

  # Create phage to calculate f
  # print(phage_to_calculate)#DEBUG



  # List of sigma values to sum accorss blobs.Top of the delta equation
  sigma_blob_list<-c()
  # Calculate total number of blobs.
  # total_blobs=0 #not needed; can just use sigma_blob_list length.
  # for (j in 1:str_length(phage_to_calculate)){
  ########### Calculate delta.   ###########
  j<-1
  while(j<=str_length(phage_to_calculate)){

    phage_blob_j<-substr(phage_to_calculate,j,j+g-1)
    sigma_blob_j<-calculate_sigma(phage_blob_j)
    sigma_blob_list<-c(sigma_blob_list,(sigma_blob_j-sigma_total)^2)

    # Increment J accorss phage.

    j<-j+g
    # print(phage_blob_j)#DEBUG
    # print(sigma_blob_list)#DEBUG
    # print(j+g)#DEBUG

  }
  # print(sigma_blob_list)#DEBUG


  delta_blob_j<-mean(sigma_blob_list)
  delta<-c(delta,delta_blob_j)
  # print(i)#DEBUG


  ########### Calculate delta_max   ###########
  sigma_blob_list_max<-c()
  k<-1
  while(k<=str_length(max_phage_to_calculate)){

    phage_blob_k<-substr(max_phage_to_calculate,k,k+g-1)
    sigma_blob_k<-calculate_sigma(phage_blob_k)
    sigma_blob_list_max<-c(sigma_blob_list_max,(sigma_blob_k-sigma_total)^2)

    # Increment J accorss phage.

    k<-k+g
    # print(phage_blob_k)#DEBUG
    # print(sigma_blob_list_max)#DEBUG
    # print(k+g)#DEBUG
  }
  # print(sigma_blob_list)#DEBUG


  delta_blob_k_max<-mean(sigma_blob_list_max)
  delta_max<-c(delta_max,delta_blob_k_max)
  # print(delta_blob_k_max)#DEBUG
  # print(i)#DEBUG
  # print(delta_blob_j)#DEBUG

}


# Return this as a data frame (local) to ensure the calculations are done right.
MOD_SOURCE_omega <- MOD_SOURCE1
MOD_SOURCE_omega$delta <- delta
MOD_SOURCE_omega$delta_max <- delta_max
MOD_SOURCE_omega$omega <-MOD_SOURCE_omega$delta/MOD_SOURCE_omega$delta_max
MOD_SOURCE_omega$omega[is.na(MOD_SOURCE_omega$omega)]<-0


# omega weighted for phage abundance.


HDMOD10_omega <-c(
  sum(MOD_SOURCE_omega$omega*MOD_SOURCE1$input_ppm),
  sum(MOD_SOURCE_omega$omega*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[10,] <-c('AA_omega',HDMOD10_omega[1],HDMOD10_omega[2])





return_max_seq("PVPASMFAPEPSSPGAARAAAAAA")
calculate_sigma("PVPASMFAPEPSSPGAARAAAAAA")
calculate_sigma("RXXXXXXXXXXXXXXXXXXXXXXE")



#################### Module 11: Fraction of Charged Residues. #################### 
### Fraction of residues which are charged.
# See Mao 2010
MOD_SOURCE1_FCR <-MOD_SOURCE1

MOD_SOURCE1_FCR$FCR<-(str_count(MOD_SOURCE1$sequence,'[DE]')+str_count(MOD_SOURCE1$sequence,'[KR]'))/str_length(MOD_SOURCE1$sequence)


HMOD11_FCR <-c(
  sum(MOD_SOURCE1_FCR$FCR*MOD_SOURCE1_FCR$input_ppm),
  
  sum(MOD_SOURCE1_FCR$FCR*MOD_SOURCE1_FCR$output_ppm)
)
HMOD_PROP[11,] <-c('FCR',HMOD11_FCR[1],HMOD11_FCR[2])




#################### Module 12: Net charge per residue  #################### 
### Fraction of residues which are charged.
# See Mao 2010
MOD_SOURCE1_NCPR <-MOD_SOURCE1

MOD_SOURCE1_NCPR$NCPR<-(-str_count(MOD_SOURCE1_NCPR$sequence,'[DE]')+str_count(MOD_SOURCE1_NCPR$sequence,'[KR]'))/str_length(MOD_SOURCE1_NCPR$sequence)


HMOD12_NCPR <-c(
  sum(MOD_SOURCE1_NCPR$NCPR*MOD_SOURCE1_NCPR$input_ppm),
  
  sum(MOD_SOURCE1_NCPR$NCPR*MOD_SOURCE1_NCPR$output_ppm)
)
HMOD_PROP[12,] <-c('NCPR',HMOD12_NCPR[1],HMOD12_NCPR[2])

#################### Module 13: Net charge  #################### 
### Fraction of residues which are charged.
MOD_SOURCE1_net_charge <-MOD_SOURCE1

MOD_SOURCE1_net_charge$net_charge<-(-str_count(MOD_SOURCE1_NCPR$sequence,'[DE]')+str_count(MOD_SOURCE1_NCPR$sequence,'[KR]'))


HMOD13_net_charge <-c(
  sum(MOD_SOURCE1_net_charge$net_charge*MOD_SOURCE1$input_ppm),
  
  sum(MOD_SOURCE1_net_charge$net_charge*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[13,] <-c('net_charge',HMOD13_net_charge[1],HMOD13_net_charge[2])


#################### Module 14: Net charge with MAPK phosphomotifs.  #################### 
### Fraction of residues which are charged when assuming MAPK sites (S/T)P are phosphorylated. See Zarin 2017.
MOD_SOURCE1_net_charge_P <-MOD_SOURCE1

MOD_SOURCE1_net_charge_P$net_charge_P<-(
  -str_count(MOD_SOURCE1_net_charge_P$sequence,'[DE]')+
    str_count(MOD_SOURCE1_net_charge_P$sequence,'[KR]')+
    -str_count(MOD_SOURCE1_net_charge_P$sequence,'[ST]P')*2
  )


HMOD14_net_charge_P <-c(
  sum(MOD_SOURCE1_net_charge_P$net_charge_P*MOD_SOURCE1$input_ppm),
  
  sum(MOD_SOURCE1_net_charge_P$net_charge_P*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[14,] <-c('net_charge_P',HMOD14_net_charge_P[1],HMOD14_net_charge_P[2])

#################### Module 15: Sequence charge decoration.  #################### 
# Don't feel like doing this one right now.

#################### Module 16: R/K Ratio.  #################### 
# Vernon 2018. Arginine to lysine ratio.


MOD_SOURCE1_RK_ratio <-MOD_SOURCE1

MOD_SOURCE1_RK_ratio$RK_ratio<-(
  (str_count(MOD_SOURCE1_RK_ratio$sequence,'R')+1)/
    (str_count(MOD_SOURCE1_RK_ratio$sequence,'K')+1) 

)
HMOD16_RK_ratio <-c(
  sum(MOD_SOURCE1_RK_ratio$RK_ratio*MOD_SOURCE1$input_ppm),
  
  sum(MOD_SOURCE1_RK_ratio$RK_ratio*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[16,] <-c('RK_ratio',HMOD16_RK_ratio[1],HMOD16_RK_ratio[2])

#################### Module 17: E/D Ratio.  #################### 
# Glutamate to aspartic acid ratio.

MOD_SOURCE1_ED_ratio <-MOD_SOURCE1

MOD_SOURCE1_ED_ratio$ED_ratio<-(
  (str_count(MOD_SOURCE1_ED_ratio$sequence,'E')+1)/
    (str_count(MOD_SOURCE1_ED_ratio$sequence,'D')+1) 
  
)
HMOD17_ED_ratio <-c(
  sum(MOD_SOURCE1_ED_ratio$ED_ratio*MOD_SOURCE1$input_ppm),
  
  sum(MOD_SOURCE1_ED_ratio$ED_ratio*MOD_SOURCE1$output_ppm)
)
HMOD_PROP[17,] <-c('ED_ratio',HMOD17_ED_ratio[1],HMOD17_ED_ratio[2])


#################### Module 18: Linear motifs (Zarin 2019 IDs 18-54).  #################### 
# S[IVLMH]E[IVPFMLYAQR]GR.
# See Dinkel 2018.
motifs_of_interest <-as.data.frame(reference_property_table$Regular.expression..regex.[18:54])
motifs_of_interest$ID<-reference_property_table$ID[18:54]
names(motifs_of_interest)<-c('regex','ID')

# Data frame for motif of interest
MOD_SOURCE1_linear_motifs <-MOD_SOURCE1

#Iterate through motifs to screen, creating new row for each.
for (i in 1:nrow(motifs_of_interest)){
# for (i in 20){#DEBUG
  # num_of_motifs_raw <-sum(str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i]))
  MOD_SOURCE1_linear_motifs$num_of_motifs <-str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i])
# 
#     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input_ppm))#DEBUG
#     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$output_ppm))#DEBUG
#     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median))#DEBUG
#     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input))#DEBUG
#     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median)/sum(MOD_SOURCE1_linear_motifs$median))#DEBUG
#     
#     
    
    
  HMOD18to54_linear_motifs<-c(
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$input_ppm),
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$output_ppm)
    
  )
  # Write to final array. ID is offset to properly capture the name of the motif.
  HMOD_PROP[17+i,] <-c(motifs_of_interest$ID[i],HMOD18to54_linear_motifs[1],HMOD18to54_linear_motifs[2])
}

#








#################### Module 19: Physiochemical properties (Zarin 2019 IDs 55-65 which regex is sufficient).  #################### 
# Physiochemical properties involving regex. Others such as Kyte-Doolittle scale will be done separately.
# See Dinkel 2018.
motifs_of_interest <-as.data.frame(reference_property_table$Regular.expression..regex.[55:65])
motifs_of_interest$ID<-reference_property_table$ID[55:65]
names(motifs_of_interest)<-c('regex','ID')

# Data frame for motif of interest
MOD_SOURCE1_linear_motifs <-MOD_SOURCE1

#Iterate through motifs to screen, creating new row for each.
for (i in 1:nrow(motifs_of_interest)){
  # for (i in 20){#DEBUG
  # num_of_motifs_raw <-sum(str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i]))
  MOD_SOURCE1_linear_motifs$num_of_motifs <-str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i])
  # 
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input_ppm))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$output_ppm))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median)/sum(MOD_SOURCE1_linear_motifs$median))#DEBUG
  #     
  #     
  
  
  HMOD55to65_linear_motifs<-c(
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$input_ppm),
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$output_ppm)
    
  )
  # Write to final array. ID is offset to properly capture the name of the motif.
  HMOD_PROP[54+i,] <-c(motifs_of_interest$ID[i],HMOD55to65_linear_motifs[1],HMOD55to65_linear_motifs[2])
}

#################### Module 20: Repeats and complexity (Zarin 2019 IDs 66-82 which regex is sufficient).  #################### 
# Repeats and complexity properties involving regex.
# See Dinkel 2018.
motifs_of_interest <-as.data.frame(reference_property_table$Regular.expression..regex.[66:82])
motifs_of_interest$ID<-reference_property_table$ID[66:82]
names(motifs_of_interest)<-c('regex','ID')

# Data frame for motif of interest
MOD_SOURCE1_linear_motifs <-MOD_SOURCE1

#Iterate through motifs to screen, creating new row for each.
for (i in 1:nrow(motifs_of_interest)){
  # for (i in 20){#DEBUG
  # num_of_motifs_raw <-sum(str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i]))
  MOD_SOURCE1_linear_motifs$num_of_motifs <-str_count(MOD_SOURCE1_linear_motifs$sequence,motifs_of_interest$regex[i])
  # 
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input_ppm))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$output_ppm))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$input))#DEBUG
  #     print( sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1_linear_motifs$median)/sum(MOD_SOURCE1_linear_motifs$median))#DEBUG
  #     
  #     
  
  
  HMOD66to82_linear_motifs<-c(
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$input_ppm),
    sum(MOD_SOURCE1_linear_motifs$num_of_motifs*MOD_SOURCE1$output_ppm)
    
  )
  # Write to final array. ID is offset to properly capture the name of the motif.
  HMOD_PROP[65+i,] <-c(motifs_of_interest$ID[i],HMOD66to82_linear_motifs[1],HMOD66to82_linear_motifs[2])
}


#################### Graphs of properties#################### 

HMOD_PROP_forgraph<-HMOD_PROP[which(!is.na(HMOD_PROP$input_value)),]
HMOD_PROP_forgraph$enrichment <-as.numeric(HMOD_PROP_forgraph$output_value)/as.numeric(HMOD_PROP_forgraph$input_value)
HMOD_PROP_forgraph$rownum <-c(1:nrow(HMOD_PROP_forgraph))


ggplot(HMOD_PROP_forgraph,aes(reorder(property,rownum),enrichment))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 4))


ggsave(paste(phage_exp_id,protein_code,"HDIMENSION_ENRICH.png",sep="_"))
# ggsave()
#################### Graphs2#################### 


TEST<-MOD_SOURCE1[!MOD_SOURCE1$input==0,]
# TEST[TEST==0] <-  1/31

TEST$enrichment <-TEST$median/TEST$input
TEST <- TEST[order(TEST$enrichment,decreasing=TRUE),]
TEST<-dplyr::mutate(TEST, ID = row_number())


# Is the difference between two values?
# TEST$value_test <- ifelse(between(TEST$output_difference,0,TEST$median*3),1,-1)

# Is the sdeviatn between two values?
TEST$value_test <- ifelse(between(sqrt(TEST$variance),0,TEST$median*3),1,-1)


# TEST$value_test <- ifelse(str_detect(TEST$sequence,'[S]{2,}'),1,-1)


increment<-c(0)
for (i in 1:nrow(TEST)){
  increment<-c(increment,increment[i-1]+TEST$value_test[i])
}

TEST$value_accumulation<-increment



# ggplot(TEST,aes(x=ID,y=sdstuff))+geom_point(stat="identity")
# ggplot(TEST,aes(x=ID,y=log2(enrichment)))+geom_point(stat="identity")
ggplot(TEST,aes(x=(ID),y=(value_accumulation)))+geom_line(stat="identity")
# ggplot(TEST,aes(x=median,y=(output_difference)))+geom_point(stat="identity")








############################# REFERENCES ############################# 

#Zarin, Taraneh, Bob Strome, Alex N Nguyen Ba, Simon Alberti, Julie D Forman-Kay, and Alan M Moses. “Proteome-Wide Signatures of Function in Highly Diverged Intrinsically Disordered Regions.” eLife 8 (n.d.): e46883. https://doi.org/10.7554/eLife.46883.
#
# Das RK, Pappu RV. Conformations of intrinsically disordered 
# proteins are influenced by linear sequence distributions of oppositely 
# charged residues. Proc Natl Acad Sci U S A. 2013 Aug 13;110(33):13392-7. 
# doi: 10.1073/pnas.1304749110. Epub 2013 Jul 30. PMID: 23901099; PMCID: PMC3746876.

# Martin EW, Holehouse AS, Grace CR, Hughes A, Pappu RV, Mittag T. Sequence Determinants of the Conformational Properties of an Intrinsically Disordered Protein Prior to and upon Multisite Phosphorylation. J Am Chem Soc. 2016 Nov 30;138(47):15323-15335. doi: 10.1021/jacs.6b10272. Epub 2016 Nov 17. PMID: 27807972; PMCID: PMC5675102.

# For other references, see references section from:
#Zarin, Taraneh, Bob Strome, Alex N Nguyen Ba, Simon Alberti, Julie D Forman-Kay, and Alan M Moses. “Proteome-Wide Signatures of Function in Highly Diverged Intrinsically Disordered Regions.” eLife 8 (n.d.): e46883. https://doi.org/10.7554/eLife.46883.



# Create random background
# Library data. Include.
# Make libraries in paralllellllllllll and poooooooooooool them
# Visualization of where are the "enriched" stuff. 
# Set of interest is 10 genes.
# Total is 5000 genes.
# Rank "score" (most enriched to much depleted.)
# Add score of 1/number in list in set of interest.
# Or, subtract score of 1/total number of inputs.
# Biostat squid
# What is the tendancy of things that are calling "enriched" to be present






###################### JUNKYARD

# # locations of positive charges in each phage sequence.
# pos_charge_locations<-str_locate_all(MOD_SOURCE1$sequence,'[KR]')
# neg_charge_locations<-str_locate_all(MOD_SOURCE1$sequence,'[DE]')
# 
# 
# # Testing. See
# # pos_charge_locations <-str_locate_all("KEKEKEKEKEKEKEKEKEKEKE",'[KR]')
# # pos_charge_locations <-str_locate_all("KKKKKKKKKKKEEEEEEEEEEE",'[KR]')
# 
# # Iterate through positive charge locations (one per phage)
# average_location_gap_pos_charge<-c()
# for(i in 1:length(pos_charge_locations)){
#   # Iterate through sublist of positive charge locations per phage.
#   # Divide by two, since two columns (start and end). You only care about
#   # the start in this case.
#   # print(i)
#   
#   # Calculate charge asymmetry from paper Das and Pappu. mean(str_length(MOD_SOURCE1$sequence)) is the total number of residues in peptide.
#   sigma<-(length(pos_charge_locations[[i]])/mean(str_length(MOD_SOURCE1$sequence))-length(neg_charge_locations[[i]])/mean(str_length(MOD_SOURCE1$sequence)))^2/
#     (length(pos_charge_locations[[i]])/mean(str_length(MOD_SOURCE1$sequence))+length(neg_charge_locations[[i]])/mean(str_length(MOD_SOURCE1$sequence)))
#   # print(sigma)
#   
#   
#   
#   set_of_locations<-c()
#   for(j in 1:nrow(pos_charge_locations[[i]]) ){
#     # Ensure you don't go beyond the end of residue locations.
#     if(j+1<=nrow(pos_charge_locations[[i]])){
#       difference_marker<-abs(pos_charge_locations[[i]][j,1]-pos_charge_locations[[i]][j+1,1])
#       set_of_locations<-c(set_of_locations,difference_marker)
#     }
#   }
#   set_of_locations[is.null(set_of_locations)]<-0
#   
#   
#   
#   average_location_gap_pos_charge<-c(average_location_gap_pos_charge,mean(set_of_locations))
#   # print(TEST)
# }
# # Iterate through negative charge locations (one per phage)
# average_location_gap_neg_charge<-c()
# for(i in 1:length(neg_charge_locations)){
#   # Iterate through sublist of positive charge locations per phage.
#   # Divide by two, since two columns (start and end). You only care about
#   # the start in this case.
#   # print(i)
#   set_of_locations<-c()
#   for(j in 1:nrow(neg_charge_locations[[i]]) ){
#     # Ensure you don't go beyond the end of residue locations.
#     if(j+1<=nrow(neg_charge_locations[[i]])){
#       difference_marker<-abs(neg_charge_locations[[i]][j,1]-neg_charge_locations[[i]][j+1,1])
#       set_of_locations<-c(set_of_locations,difference_marker)
#     }
#   }
#   set_of_locations[is.null(set_of_locations)]<-0
#   average_location_gap_neg_charge<-c(average_location_gap_neg_charge,mean(set_of_locations))
#   # print(TEST)
# }
# 
# 
# # Attach to data frame.
# MOD_SOURCE_kappa <- MOD_SOURCE1
# 
# HDMOD9_kappa <-c(HDMMOD_SOURCE1_IDOD9_kappa <-c(
#   sum(str_count(MOD_SOURCE1$sequence,"G")*MOD_SOURCE1$input)/sum(MOD_SOURCE1$input),
#   
#   sum(str_count(MOD_SOURCE1$sequence,"G")*MOD_SOURCE1$output)/sum(MOD_SOURCE1$output)
# )
# HMOD_PROP[9,] <-c('AA_G',HDMOD8_AA_G[1],HDMOD8_AA_G[2])
# 
