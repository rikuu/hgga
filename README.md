# HGGA

Hierarchical Guided Genome Assembler

## Getting started

```sh
# Compile minimap2, kermit and HGGA (requires gcc and zlib)
git clone https://github.com/rikuu/hgga && cd hgga
git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)
git clone --recursive https://github.com/rikuu/kermit && (cd kermit && make)

# Compile Racon
git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../..

# Download E. coli K-12 sequence
wget -O NC_000913.fa 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&db=nuccore&dopt=fasta&val=556503834'

# Download PacBio E. coli sample reads
wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
ln -s selfSampleData/pacbio_filtered.fastq reads.fq

# Color
minimap2/minimap2 -x map-pb NC_000913.fa reads.fq > align.paf
kermit/kermit-color align.paf kermit/misc/ecoli.txt > color.cf

# Assembly (requires python3 and pysam)
./hgga.py reads.fq color.cf
```
