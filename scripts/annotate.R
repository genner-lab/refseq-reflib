#!/usr/bin/env Rscript

# load libs/funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-s","--seed"), type="numeric"),
    make_option(c("-p","--primer"), type="character")
    )
# set args
opt <- optparse::parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# testing
#opt <- NULL
#opt$seed <- 42
#opt$primer <- "tele02"

writeLines("\n...\nObtaining taxonomy from GBIF\n")
Sys.sleep(3)

# read in the version
version <- readLines(here("temp/RELEASE_NUMBER"))

# read in the mito refseq fasta
mito.cat <- suppressMessages(suppressWarnings(readr::read_csv(file=here::here("temp/refseq.mitochondrion.cat.haplotypes.csv")))) |> 
    dplyr::mutate(taxid=as.character(taxid))

# create a db from GBIF
writeLines("\n...\nObtaining taxonomic names database from GBIF (may take up to an hour).\n")
writeLines(paste("\n...\nLocal taxadb database location:\n",taxadb::taxadb_dir()))

taxadb::td_create(provider="gbif",overwrite=TRUE)

# get ids 
gbif.ids <- taxadb::get_ids(names=unique(dplyr::pull(mito.cat,scientificName)),provider="gbif")

# get taxonomy
gbif.tax <- gbif.ids |> taxadb::filter_id(provider="gbif") 

mito.cat.annotated <- gbif.tax |>
    dplyr::distinct(kingdom,phylum,class,order,family,genus,scientificName) |> 
    dplyr::right_join(mito.cat,by=c("scientificName","genus")) |>
    dplyr::filter(!is.na(kingdom) & !is.na(phylum) & !is.na(class) & !is.na(order) & !is.na(family) & !is.na(genus)) |>
    #filter(is.na(kingdom) | is.na(phylum) | is.na(class) | is.na(order) | is.na(family) | is.na(genus)) |>
    dplyr::mutate(refseqVersion=version) |>
    dplyr::arrange(kingdom,phylum,class,order,family,genus,scientificName) |>
    dplyr::mutate(label=paste0(accession,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",scientificName)) |> 
    dplyr::mutate(label=stringr::str_replace_all(label," ","_"))

# write out the table
mito.cat.annotated |> dplyr::select(refseqVersion,accession,kingdom,phylum,class,order,family,genus,scientificName,taxid,length,nHaps,nucleotides) |>
    readr::write_csv(here("references",paste0("refseq",version,"-annotated-",opt$primer,".csv")))

# convert to fasta
refseq.fas.all <- mito.cat.annotated |> tab2fas(seqcol="nucleotides",namecol="label")

# make a genus only table
set.seed(opt$seed)
refseq.fas.genus <- mito.cat.annotated |> 
    dplyr::group_by(genus) |> 
    dplyr::slice_sample(n=1) |> 
    dplyr::ungroup() |> 
    dplyr::arrange(kingdom,phylum,class,order,family,genus,scientificName) |> 
    tab2fas(seqcol="nucleotides",namecol="label")

# write out
ape::write.FASTA(refseq.fas.all,file=here::here("references",paste0("refseq",version,"-annotated-",opt$primer,".fasta")))

# write out
ape::write.FASTA(refseq.fas.genus,file=here::here("references",paste0("refseq",version,"-annotated-genera-",opt$primer,".fasta")))

# report
writeLines("\n...\nCleaned RefSeq reference libraries written to 'references'.\n")
