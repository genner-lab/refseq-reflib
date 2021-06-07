#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-s","--seed"), type="character")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

writeLines("\n...\nObtaining taxonomy from GBIF\n")
Sys.sleep(3)

# read in the version
version <- readLines(here("temp/RELEASE_NUMBER"))

# read in the mito refseq fasta
mito.cat <- suppressMessages(suppressWarnings(read_csv(file=here("temp/refseq.mitochondrion.cat.haplotypes.csv"))))

# create a db from GBIF
td_create(provider="gbif",dbdir=here("temp"),overwrite=TRUE,schema="dwc")

# connect to db
td_connect(dbdir=here("temp"))

# pull out all the db then subset genera
gbif.taxonomy <- taxa_tbl(provider="gbif",schema="dwc") %>% 
    collect() %>% 
    filter(taxonRank=="genus" & taxonomicStatus=="accepted")

# filter on our list of genera
# where duplicates, select the animal
gbif.taxonomy.filtered <- gbif.taxonomy %>%
    filter(scientificName %in% unique(pull(mito.cat,genus))) %>% 
    select(kingdom,phylum,class,order,family,genus) %>% 
    distinct() %>% 
    add_count(genus) %>% 
    filter(n==1 | n>1 & kingdom=="Animalia") %>%
    select(-n) %>%
    arrange(kingdom,phylum,class,order,family,genus) %>% 
    mutate(classified=TRUE)

# annotate the mito.cat
mito.cat.annotated <- mito.cat %>% 
    left_join(gbif.taxonomy.filtered,by="genus") %>% 
    filter(classified==TRUE) %>%
    mutate(refseqVersion=version) %>%
    arrange(kingdom,phylum,class,order,family,genus,scientificName) %>% 
    mutate(label=paste0(accession,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",scientificName)) %>% 
    mutate(label=str_replace_all(label," ","_"))


# write out the table
mito.cat.annotated %>% select(refseqVersion,accession,kingdom,phylum,class,order,family,genus,scientificName,taxid,length,nHaps,nMatches,matchTax,nucleotides) %>%
    write_csv(here("references/refseq-annotated.csv"))

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
write.FASTA(refseq.fas.all,file=here("references/refseq-annotated.fasta"))

# write out
write.FASTA(refseq.fas.genus,file=here("references/refseq-annotated-genera.fasta"))

# disconnect from db
td_disconnect()

# report
writeLines("\n...\nCleaned RefSeq reference libraries written to 'references'\n")
