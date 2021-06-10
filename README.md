[![DOI](https://zenodo.org/badge/xxx.svg)](https://zenodo.org/badge/latestdoi/xxx)

# refseq-reflib
A pipeline written in bash and R to access and process mitochondrial reference library data from [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/). Currently supported 12S metabarcode markers are: 'tele02', 'mifish-u', 'elas02', and 'mifish-u-mod'.

### Setup

- Install [hmmer](http://hmmer.org/) v3.1b2:
```
sudo apt install hmmer
```

- Make sure system utilities including wget, curl, md5sum, gzip are available.

- Clone the repository:
```
git clone https://github.com/genner-lab/refseq-reflib.git
```

- Change directory:
```
cd refseq-reflib
```

- Make temp and output directories:
```
mkdir temp references
```

- Obtain correct R packages:
```
Rscript -e "renv::restore()"
```

### Check RefSeq version

```
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
```

### Download latest RefSeq

- May take more than 10 minutes to download and unpack the full catalogue, depending on internet connection speed.

```
scripts/download.sh
```

### Extract sequences

- Use hidden Markov models to extract homologous DNA sequences of the fragment of interest.

- The '-p' flag is the primer set. Currently supported sets for 12S are 'tele02', 'mifish-u', 'elas02', and 'mifish-u-mod'.

- The script removes all hybrid (e.g. '_Cyprinus carpio_ x _Carassius auratus_') and undescribed taxa (e.g. '_Corydoras_ sp.').

- The script converts all subspecies to species (e.g. '_Thunnus thynnus thynnus_' becomes '_Thunnus thynnus_').

- The script then taxonomically dereplicates the sequences, meaning that within species all duplicated haplotypes will be removed, but different species can share haplotypes.

```
scripts/extract.R -p tele02
```

### Annotate taxonomy

- Annotate the mtDNA data with taxonomic information from [GBIF](https://www.gbif.org/). May take 10 minutes or more to retrieve the taxonomic database, depending on internet connection speed.

- The script writes out three files into 'references': (a) the annotated RefSeq sequences in tabular format ('refseqVERSION-annotated-PRIMER.csv'); (b) annotated RefSeq fasta sequences formatted for the sintax algorithm employed in vsearch ('refseqVERSION-annotated-PRIMER.fasta'); and (c) the same sintax fasta sequences, but for one random species selected per genus ('refseqVERSION-annotated-genera-PRIMER.fasta'). 

- The '-s' flag is the random number seed for the per genus subset.

- The '-p' flag is the primer set. Currently supported sets for 12S are 'tele02', 'mifish-u', 'elas02', and 'mifish-u-mod'.

```
scripts/annotate.R -s 42 -p tele02
```

### Clean up (optional)

- The intermediate files for the taxonomic database are quite large, and it may be required to remove these to save disk space.

```
rm temp/duckdb
```
