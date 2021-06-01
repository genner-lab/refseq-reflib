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


# subset the marker of interest from the mtDNA dump
mt.sub <- run_hmmer3(dir=here("temp"), infile="refseq.mitochondrion.genomic.fna", prefix="12s.taberlet.noprimers", evalue="10", coords="env")


#write.FASTA(mt.sub,file=here("temp/refseq.mitochondrion.genomic.SUB.fna"))


# read in the mito refseq catalogue
# first download with `taxonomic-assignment.sh`
mito.cat <- read_tsv(here("temp/refseq.mitochondrion.cat.tsv"),guess_max=999999,col_names=c("taxid","scientificName","accession","dir","status","len"))



# remove the long and short sequences leaving animal mitogenomes 
#mt <- mt[lapply(mt,length) < 21000 & lapply(mt,length) > 12000]

# get accesssions from names
accs <- names(mt.sub)

# filter the catalogue down
mito.cat %<>% filter(accession %in% accs)

# add the genus from sciName
mito.cat %<>% mutate(genus=str_split_fixed(scientificName," ",2)[,1])


# taxonomic dereplication of sequences

# add sequences

mt.sub.df <- tibble(accession=names(mt.sub),nucleotides=mapply(paste,collapse="",as.character(mt.sub),USE.NAMES=FALSE)) %>% mutate(length=str_length(nucleotides))

# merge

mito.cat.nucs <- mito.cat %>% left_join(mt.sub.df)# %>% rename(dbid=accession,sciNameValid=scientificName)


mito.cat.haps <- haps2fas(mito.cat.nucs)

mito.cat.haps %>% filter(nHaps>1)
mito.cat.haps %>% filter(scientificName=="Rattus norvegicus")


mito.cat.haps %>% count(scientificName) %>% filter(n>1)
mito.cat.nucs %>% count(scientificName) %>% filter(n>1)

mito.cat.nucs %>% filter(scientificName=="Lophius piscatorius")
mito.cat.nucs %>% filter(length<100)


derp <- read.FASTA(here("temp/refseq.mitochondrion.genomic.DEREP.fna"))

missing <- setdiff(pull(mito.cat.nucs,accession),names(derp))

mito.cat.nucs %>% filter(accession %in% missing) %>% arrange(scientificName) %>% print(n=Inf)

mito.cat.nucs %>% write_csv(file="temp/mito.cat.nucs.csv")

###
options(width=130)


vsearch --derep_fulllength refseq.mitochondrion.genomic.fna --minuniquesize 1 --fasta_width 0 --output refseq.mitochondrion.genomic.DEREP.fna

refseq.mitochondrion.genomic.fna