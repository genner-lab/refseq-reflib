# refseq-reflib
Access and process mitochondrial data from NCBI RefSeq

### make temp directory
mkdir temp references

### check refseq version
curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER

### download
scripts/download.sh

# filter out and clean
scripts/assemble.R -p tele02

# add taxonomy
scripts/annotate.R -s 42
