#!/usr/bin/env/ sh

#################################################### 
#################################################### 


# DOWNLOAD REFSEQ
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER -P temp
head temp/RELEASE_NUMBER
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/README -P temp
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz -P temp
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz -P temp
# check the release number
file="temp/RELEASE_NUMBER"
gbv=$(cat "$file")
echo "Downloaded RefSeq version ""$gbv"
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release"$gbv".catalog.gz -P temp
gzip -cd temp/RefSeq-release"$gbv".catalog.gz | grep "mitochondrion" > temp/refseq.mitochondrion.cat.tsv
rm temp/RefSeq-release"$gbv".catalog.gz

# unzip
gzip -d temp/mitochondrion.1.1.genomic.fna.gz 
gzip -d temp/mitochondrion.2.1.genomic.fna.gz 

# join
cat temp/mitochondrion.1.1.genomic.fna temp/mitochondrion.2.1.genomic.fna > temp/refseq.mitochondrion.genomic.fna
rm temp/mitochondrion.1.1.genomic.fna temp/mitochondrion.2.1.genomic.fna

# clean
sed -i -e 's/ .*//g' temp/refseq.mitochondrion.genomic.fna

# check 
head temp/refseq.mitochondrion.genomic.fna
