#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# report
writeLines("\n...\nExtracting nucleotides with hidden Markov models\n")
Sys.sleep(3)

# for testing
#opt <- NULL
#opt$primer <- "tele02"

# make prefix
if(opt$primer=="tele02") {
    prefix <- "12s.taberlet.noprimers"
} else if(opt$primer=="elas02" | opt$primer=="mifish-u" | opt$primer=="mifish-u-mod") {
    prefix <- "12s.miya.noprimers"
} else {
    stop("Primers must be 'tele02', 'elas02', 'mifish-u', or 'mifish-u-mod'")
}

# subset the marker of interest from the mtDNA dump
mt.sub <- run_hmmer3(dir=here("temp"), infile="refseq.mitochondrion.genomic.fna", prefix=prefix, evalue="10", coords="env")

# report
writeLines("\n...\nFiltering and dereplicating\n")
Sys.sleep(3)

# read in the mito refseq catalogue
mito.cat <- suppressMessages(suppressWarnings(read_tsv(here("temp/refseq.mitochondrion.cat.tsv"),guess_max=999999,col_names=c("taxid","scientificName","accession","dir","status","len"))))

# get accesssions from names
accs <- names(mt.sub)

# filter the catalogue down
mito.cat %<>% filter(accession %in% accs)

# get only species and subspecies - filter out "sp.", "x" etc  
mito.cat %<>% filter(grepl("(^[A-Z][a-z]+[[:space:]][a-z]+$)|(^[A-Z][a-z]+[[:space:]][a-z]+[[:space:]][a-z]+$)",scientificName))

# convert subspecies to species
mito.cat %<>% mutate(scientificName=paste(str_split_fixed(scientificName," ",3)[,1],str_split_fixed(scientificName," ",3)[,2]))

# add sequence data to a df
mt.sub.df <- tibble(accession=names(mt.sub),nucleotides=mapply(paste,collapse="",as.character(mt.sub),USE.NAMES=FALSE)) %>% mutate(length=str_length(nucleotides))

# merge with catalog
mito.cat.nucs <- mito.cat %>% left_join(mt.sub.df,by="accession")

# collapse by haplotypes
mito.cat.haps <- haps2fas(mito.cat.nucs)

# add the genus from sciName
mito.cat.haps %<>% mutate(genus=str_split_fixed(scientificName," ",2)[,1])

# write out
mito.cat.haps %>% write_csv(file=here("temp/refseq.mitochondrion.cat.haplotypes.csv"))

# report
writeLines("\n...\nDereplicated sequences written to 'temp/refseq.mitochondrion.cat.haplotypes.csv'.\n")
