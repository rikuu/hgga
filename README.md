# HGGA

Hierarchical Guided Genome Assembler (HGGA)

## Getting started

```sh
# Compile minimap2, miniasm and HGGA (requires gcc and zlib)
git clone https://github.com/rikuu/hgga && cd hgga
git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)

# Compile Racon
git clone https://github.com/lbcb-sci/racon && cd racon && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make && cd ../..

# Download E. coli K-12 sequence
wget -O NC_000913.fa 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&db=nuccore&dopt=fasta&val=556503834'

# Download PacBio E. coli sample reads
wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
ln -s selfSampleData/pacbio_filtered.fastq reads.fq

# Download example linkage map
wget -O ecoli.txt https://raw.githubusercontent.com/rikuu/kermit/master/misc/map.txt

# Color
minimap2/minimap2 -x map-pb NC_000913.fa reads.fq > align.paf
./closest_color.py align.paf ecoli.txt > color.cf

# Assembly (requires python3 and pysam)
./hgga.py reads.fq color.cf
```
