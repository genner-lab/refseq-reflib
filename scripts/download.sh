#!/usr/bin/env sh

#################################################### 
#################################################### 

# download refseq version
wget -q https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER -P temp
file="temp/RELEASE_NUMBER"
gbv=$(cat "$file")
echo "Downloading RefSeq version ""$gbv"" ..."
sleep 3

# download readme and DNA seqs
wget -q https://ftp.ncbi.nlm.nih.gov/refseq/release/README -P temp
wget -q https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz -P temp
wget -q https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz -P temp

# download 
wget -q https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release"$gbv".catalog.gz -P temp
printf "...\nGetting catalogue md5sum\n"
md5sum temp/RefSeq-release"$gbv".catalog.gz
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

printf "...\nRefSeq DNA and catalogue written to 'temp'.\n"
