# MetaBMF: An Ultra-fast Algorithm for Large-scale Reference-free Metagenomic Studies

MetaBMF is an ultra-fast algorithm for large-scale reference-free metagenomic studies. The pipeline outputs all binned species in multiple metagenomic samples and their estimated relative abundances.  A test dataset with 6 microbial species is used to illustrate the MetaBMF pipeline. The data is available at https://figshare.com/articles/test-data_zip/3491951

## Some Dependencies

* [Ray Assembler](http://denovoassembler.sourceforge.net/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://www.htslib.org/)

## Download MetaBMF and the Testing Dataset.
```
mkdir example
cd example
git clone https://github.com/didi10384/MetaBMF
wget https://ndownloader.figshare.com/files/5523092
unzip 5523092

```
## Create file to store sample names.
```
cd ./test_data
for f in *.fasta; do
g=$(echo $f |gawk '{gsub(/.*[/]|.fasta/, "", $0)} 1')
echo -e "$PWD/$g" >> read_list.txt
done

```
## Extract Features

Extract useful features used for binning and estimating relative abundances, including assembled contigs and mapped reads counts matrix.

```
cd ..
mkdir metabmf_work
./MetaBMF/MetaBMF.sh -a -c ./test_data/ray/Contigs.fasta -o ./metabmf_work -s ./test_data/read_list.txt

```
Remarks:
-a denotes that the file type is fasta file (-q for fastq file); -c is the option for the path of assembled contigs; -o is the option for the path of working directory and -s is the option for the file including the list of sample names for the single-end reads (-p for paired-end reads).

More detailed options for feature extraction are listed below

```
-o      Output directory
-c      The path to the contigs fasta file
-s      the list of single-end sample names
-p      the list of paired-ends sample names
-q      For fastq files
-a      For fasta files

```
## Binning & Estimating Relative Abundances

Options for binning & estimating relative abundances
```
-h help documentation.
-r Specify the minimum contig length, contigs shorter than this value will not be included. Default is 1000.
-m Specify the MetaBMF's installation directory.
-w Specify the working directory of current dataset.
-l Specify the minimum number of clusters.
-u Specify the maximum number of clusters.
-e Specify the increment of the number of clusters from bin_min to bin_max.

```
Install R packages
```
Rscript ./MetaBMF/R/install_R_packages.R
```
Get the binning results and relative-abundance estimation

```
Rscript ./MetaBMF/R/MetaBMF.R -m "./MetaBMF" -w "./metabmf_work" -l 2 -u 15 -e 1
```

