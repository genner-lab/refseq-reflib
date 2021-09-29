#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-libs.R"))


# get args
option_list <- list( 
    make_option(c("-s","--seed"), type="numeric"),
    make_option(c("-p","--primer"), type="character")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# testing
#opt <- NULL
#opt$seed <- 42
#opt$primer <- "tele02"

writeLines("\n...\nObtaining taxonomy from GBIF\n")
Sys.sleep(3)

# read in the version
version <- readLines(here("temp/RELEASE_NUMBER"))

# read in the mito refseq fasta
mito.cat <- suppressMessages(suppressWarnings(read_csv(file=here("temp/refseq.mitochondrion.cat.haplotypes.csv")))) %>% 
    mutate(taxid=as.character(taxid))

# create a db from GBIF
writeLines("\n...\nObtaining taxonomic names database from GBIF (may take up to an hour).\n")
writeLines(paste("\n...\nLocal taxalight database location:\n",tl_dir()))

# create a taxalight database
if(!file.exists(paste0(tl_dir(),"/gbif/",format(Sys.Date(),"%Y"),"/data.mdb"))){
    tl_create("gbif")
    #tl_create("ncbi")
}

# annotate mito.cat with gbif
mito.cat.annotated <- tl(unique(pull(mito.cat,scientificName)),"gbif") %>%
    as_tibble() %>% 
    mutate(scientificName=paste(genus,specificEpithet)) %>%
    distinct(kingdom,phylum,class,order,family,genus,scientificName) %>%
    dplyr::right_join(mito.cat,by=c("scientificName","genus")) %>%
    filter(!is.na(kingdom) & !is.na(phylum) & !is.na(class) & !is.na(order) & !is.na(family) & !is.na(genus)) %>%
    #filter(is.na(kingdom) | is.na(phylum) | is.na(class) | is.na(order) | is.na(family) | is.na(genus)) %>%
    mutate(refseqVersion=version) %>%
    arrange(kingdom,phylum,class,order,family,genus,scientificName) %>% 
    mutate(label=paste0(accession,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",scientificName)) %>% 
    mutate(label=str_replace_all(label," ","_"))


# ncbi provider
#tl(paste("NCBI",pull(mito.cat,taxid),sep=":"),"ncbi") %>%
#    as_tibble() %>% 
#    distinct(taxonID,kingdom,phylum,class,order,family,genus,scientificName) %>%
#    rename(taxid=taxonID) %>%
#    mutate(taxid=str_replace_all(taxid,"NCBI:","")) %>%
#    select(-scientificName,-genus) %>%
#    dplyr::right_join(mito.cat,by="taxid") %>%

# write out the table
mito.cat.annotated %>% select(refseqVersion,accession,kingdom,phylum,class,order,family,genus,scientificName,taxid,length,nHaps,nucleotides) %>%
    write_csv(here("references",paste0("refseq",version,"-annotated-",opt$primer,".csv")))

# convert to fasta
refseq.fas.all <- mito.cat.annotated %>% tab2fas(seqcol="nucleotides",namecol="label")

# make a genus only table
set.seed(opt$seed)
refseq.fas.genus <- mito.cat.annotated %>% 
    group_by(genus) %>% 
    dplyr::slice_sample(n=1) %>% 
    ungroup() %>% 
    arrange(kingdom,phylum,class,order,family,genus,scientificName) %>% 
    tab2fas(seqcol="nucleotides",namecol="label")

# write out
write.FASTA(refseq.fas.all,file=here("references",paste0("refseq",version,"-annotated-",opt$primer,".fasta")))

# write out
write.FASTA(refseq.fas.genus,file=here("references",paste0("refseq",version,"-annotated-genera-",opt$primer,".fasta")))

# report
writeLines("\n...\nCleaned RefSeq reference libraries written to 'references'.\n")
