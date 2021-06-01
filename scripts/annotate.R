#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character"),
    make_option(c("-l","--lib"), type="character"),
    make_option(c("-f","--lenfwd"), type="numeric"),
    make_option(c("-r","--lenrev"), type="numeric")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# for testing
#opt <- NULL
#opt$primer <- "tele02"
#opt$lib <- "lib3"
#opt$lenfwd <- 18
#opt$lenrev <- 20

# read in the mito refseq fasta
#mt <- read.FASTA(here("temp/refseq.mitochondrion.genomic.fna"))

# create a db from GBIF
td_create(provider="gbif",dbdir=here("temp"),overwrite=TRUE)

# connect to db
td_connect(dbdir=here("temp"))

# pull out all the db then subset and remove dups
gbif.taxonomy <- taxa_tbl(provider="gbif") %>% 
    collect() %>% 
    filter(scientificName %in% pull(mito.cat,genus)) %>% 
    select(kingdom,phylum,class,order,family,genus) %>% 
    distinct() %>% 
    add_count(genus) %>% 
    filter(n==1) %>% 
    arrange(kingdom,phylum,class,order,family,genus) %>% 
    mutate(classified=TRUE)

# annotate the mito.cat
mito.cat.annotated <- mito.cat %>% 
    left_join(gbif.taxonomy,by="genus") %>% 
    filter(classified==TRUE) %>% 
    filter(kingdom=="Animalia") %>% 
    arrange(kingdom,phylum,class,order,family,genus,scientificName) %>% 
    mutate(label=paste0(accession,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",scientificName)) %>% 
    mutate(label=str_replace_all(label," ","_"))

# drop the unwanted seqs
mt.sub <- mt[names(mt) %in% pull(mito.cat.annotated,accession)]

# rename with taxonomy
names(mt.sub) <- pull(mito.cat.annotated,label)[match(names(mt.sub),pull(mito.cat.annotated,accession))]

# write out
write.FASTA(mt.sub,file=here("temp/refseq-annotated.fasta"))

# disconnect from db
td_disconnect()
