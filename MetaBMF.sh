#!/bin/bash
HELPDOC=$( cat <<EOF
Usage:
    bash `basename $0` [options] -c <contigs> -o <outdir> -p <reads-list>
Options:
	-o      Output directory
	-c      The path to the contigs fasta file
    -s      list of sample names for single-end reads
    -p      list of sample names for paired-ends reads
    -q      For fastq files
    -a      For fasta files
######## Options for binning
	-n Number of threads: Specify the number of CPU cores used for parallel computing. When there is a large number of contigs, it is recommended to set multiple CPU cores to accelerate the computation. The default number is 1.
	-l Specify the minimum number of clusters. The default is 2.
	-u Specify the maximum number of clusters. The default is 0, which will let the algorithm sets the maximum number of cluster automatically.
	-e Specify the increment of the number of clusters from bic_min to bic_max. The default is 1.
	-t Specify the threshold for setting the initial value. It is recommended to set this number smaller(0.01-0.1) when the number of samples is less than $10$ and larger (0.1-0.2) when the number of samples is larger than $10$. The default value is set to 0.1.
	-i Specify how many percent contigs are used to set the initial value of the algorithm. The default value is 1, which means that all the contigs is used to find initial value of the algorithm. The number can be set to a smaller one, when there are a very large number of contigs.
	-r Specify the minimum contig length, contigs shorter than this value will not be included. Default is 500.
	-c If the value is "T", output the plot of BIC scores. The default is "F".
	-m The value is "1" for the simple metagenomic community. The value is "2" for the complex metagenomic community.
    -h  This help documentation.
EOF
)


###Set default value for parameters
num_thread=1
bic_min=2
bic_max=0
bic_step=1
thred=0.1
ini_prop=1
min_ctg_len=500
bic_plot=F
bic_method=2


#Parsing options
while getopts "o:c:s:p:qahn:l:u:e:t:i:r:c:m:" opt; do
    case $opt in
    	o)
			metagen_work=$OPTARG
			;;
		c)
			ctg_file=$OPTARG
			;;
        s)
            PIRED=false
            reads_list_file=$OPTARG
            ;;
        p)
            PIRED=true
            reads_list_file=$OPTARG
            ;;
        q)
            BOWTIEOPT="-q"
            EXT="fastq"
            ;;
        a)
            BOWTIEOPT="-f"
            EXT="fasta"
            ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;

###parsing option for binnig
		n)
			num_thread=$OPTARG
			;;
		l) 
			bic_min=$OPTARG
			;;
		u)
			bic_max=$OPTARG
			;;
		e)
			bic_step=$OPTARG
			;;
		t)
			thred=$OPTARG
			;;
		i)
			ini_prop=$OPTARG
			;;
		r)
			min_ctg_len=$OPTARG
			;;
		c)
			bic_plot=$OPTARG
			;;
		m)
			bic_method=$OPTARG
			;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "$HELPDOC"
            exit 1
            ;;
        :)
    		echo "Option -$OPTARG requires an argument." >&2
    		exit 1
    		;;
    esac
done


####Check the softwere versions
if ! [ -x "$(command -v bowtie2)" ]; then
  echo 'Error: Bowtie2 is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v samtools)" ]; then
  echo 'Error: Samtools is not installed.' >&2
  exit 1
fi


####Get the path for the source code
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  TARGET="$(readlink "$SOURCE")"
  if [[ $TARGET == /* ]]; then
    echo "SOURCE '$SOURCE' is an absolute symlink to '$TARGET'"
    SOURCE="$TARGET"
  else
    DIR="$( dirname "$SOURCE" )"
    echo "SOURCE '$SOURCE' is a relative symlink to '$TARGET' (relative to '$DIR')"
    SOURCE="$DIR/$TARGET" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  fi
done
echo "SOURCE is '$SOURCE'"
metagen_path="$( dirname "$SOURCE" )"
echo "MetaGen is located at "$metagen_path


###build bowtie2 index

if [ -f $ctg_file ]
then
	echo "$ctg_file found."
else
	echo "$ctg_file not found."
fi


if [ -d "$metagen_work/contigs" ]; then 
    echo $metagen_work/contigs exists
else
    mkdir $metagen_work/contigs
fi
cp $ctg_file $metagen_work/contigs/
bowtie2-build --threads ${num_thread} $metagen_work/contigs/Contigs.fasta $metagen_work/contigs/contigs-ref

###Alignment using Bowtie2
REF=$metagen_work/contigs/contigs-ref
OUTDIR=$metagen_work/map

readarray -t read_list < $reads_list_file


if [ -d "$OUTDIR" ]; then 
    echo $OUTDIR exists
else
    mkdir $OUTDIR
fi


bowtie2_align(){
    sam_name=$1
    OUTDIR=$2
    PIRED=$3
    EXT=$4
    BOWTIEOPT=$5
    REF=$6
    echo $@
    if [ "$PIRED" == "true" ];
    then
        Q1=${sam_name}_1.$EXT
        Q2=${sam_name}_2.$EXT
        folder=$(basename $sam_name)
        if [ -d "$OUTDIR/$folder" ]; then 
            echo $OUTDIR/$folder exists
        else
            mkdir $OUTDIR/$folder
        fi
        bowtie2 $BOWTIEOPT -x $REF -1 $Q1 -2 $Q2 -S $OUTDIR/$folder/1.sam
        samtools view -Sb  $OUTDIR/$folder/1.sam > $OUTDIR/$folder/1.bam
        samtools sort -o $OUTDIR/$folder/1-sorted.bam $OUTDIR/$folder/1.bam 
        samtools index $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/1-out
        samtools idxstats $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/count.dat
    fi


    if [ "$PIRED" == "false" ];
    then
        Q0=$sam_name.$EXT
        folder=$(basename $sam_name)
        if [ -d "$OUTDIR/$folder" ]; then
            echo $OUTDIR/$folder exists
        else
            mkdir $OUTDIR/$folder
        fi
        bowtie2 $BOWTIEOPT -x $REF -U $Q0 -S $OUTDIR/$folder/1.sam
        samtools view -Sb  $OUTDIR/$folder/1.sam > $OUTDIR/$folder/1.bam
        samtools sort -o $OUTDIR/$folder/1-sorted.bam $OUTDIR/$folder/1.bam
        samtools index $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/1-out
        samtools idxstats $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/count.dat
    fi
}


export -f bowtie2_align

for i in "${read_list[@]}"
do
    printf "%s\0%s\0%s\0%s\0%s\0%s\0" "$i" "$OUTDIR" "$PIRED" "$EXT" "$BOWTIEOPT" "$REF"
done  | xargs -0 -n 6 -P $num_thread bash -c 'bowtie2_align "$@"' -- 


# if [ "$PIRED" == "true" ];
# then
#     let i=0
#     while((${#read_list[@]} > i));
#     do
#         sam_name=${read_list[i++]}
#         Q1=$sam_name.$EXT
#         Q2=$sam_name.$EXT
#         folder=$(basename $sam_name)
#         if [ -d "$OUTDIR/$folder" ]; then 
#             echo $OUTDIR/$folder exists
#         else
#             mkdir $OUTDIR/$folder
#         fi
#         bowtie2 $BOWTIEOPT -x $REF -1 $Q1 -2 $Q2 -S $OUTDIR/$folder/1.sam
#         samtools view -Sb  $OUTDIR/$folder/1.sam > $OUTDIR/$folder/1.bam
#         samtools sort -o $OUTDIR/$folder/1-sorted.bam $OUTDIR/$folder/1.bam 
#         samtools index $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/1-out
#         samtools idxstats $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/count.dat
#     done
# fi

# if [ "$PIRED" == "false" ];
# then
#     let i=0
#     while((${#read_list[@]} > i));
#     do
#         sam_name=${read_list[i++]}
#         Q0=$sam_name.$EXT
#         folder=$(basename $sam_name)
#         if [ -d "$OUTDIR/$folder" ]; then
#             echo $OUTDIR/$folder exists
#         else
#             mkdir $OUTDIR/$folder
#         fi
#         bowtie2 $BOWTIEOPT -x $REF -U $Q0 -S $OUTDIR/$folder/1.sam
#         samtools view -Sb  $OUTDIR/$folder/1.sam > $OUTDIR/$folder/1.bam
#         samtools sort -o $OUTDIR/$folder/1-sorted.bam $OUTDIR/$folder/1.bam
#         samtools index $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/1-out
#         samtools idxstats $OUTDIR/$folder/1-sorted.bam > $OUTDIR/$folder/count.dat
#     done
# fi

# echo $folder":Alignment Finished."


###parsing the reads name from files
echo "Begin to extract the read counts matrix"
# bash $metagen_path/script/combine-counts.sh -s $metagen_work
ALIGN=$metagen_work/map
OUTDIR=$metagen_work/output

if [ -d "$OUTDIR" ]; then 
    echo $OUTDIR exists
else
    mkdir $OUTDIR
fi


for i in "${read_list[@]}"
do  
    folder=$(basename $i)
    echo $folder
    if [ -d "$metagen_work/map/$folder" ] 
    then
        echo $i
        cat $metagen_work/map/$folder/count.dat|awk '{print $3}' > $OUTDIR/$folder.temp
    else
        echo File not exist
    fi
done


if [ -f $OUTDIR/count-map.tsv ]; then 
    rm $OUTDIR/count-map.tsv
fi


touch $OUTDIR/count-map.tsv
for i in "${read_list[@]}"
do  
    folder=$(basename $i)
    echo -n -e $folder'\t' >> $OUTDIR/count-map.tsv
done

sed -i -e '$a\' $OUTDIR/count-map.tsv


EXPANDED=()
for E in "${read_list[@]}"; do
    folder=$(basename $E)
    EXPANDED+=($OUTDIR/"${folder}.temp")
done
# echo "${EXPANDED[@]}"
paste ${EXPANDED[@]} | column -s $'\t' -t >> $OUTDIR/count-map.tsv


rm $OUTDIR/*.temp
sed -i '$ d' $OUTDIR/count-map.tsv 



###extract the number of reads for each sample
echo "Begin to extract the read counts for each sample"

echo "filenames_1,filenames_2,readcount_1" > $metagen_work/reads_info.txt

if [ "$PIRED" == "false" ];
then 
    let i=0
    while((${#read_list[@]} > i));
    do
        sam_name=${read_list[i++]}
        case "$EXT" in 
        fastq) 
            count=$(($(wc -l < ${sam_name}.fastq)/4))
            echo $sam_name,$sam_name,$count >> $metagen_work/reads_info.txt
            ;;
        fasta)
            count=`grep -c '>' ${sam_name}.fasta`
            echo -e $sam_name,$sam_name,$count >> $metagen_work/reads_info.txt
            ;;
        esac
    done
fi

if [ "$PIRED" == "true" ];
then 
    let i=0
    while((${#read_list[@]} > i));
    do
        sam_name=${read_list[i++]}
        # echo ${f/_1.fastq/_2.fastq}
        case "$EXT" in 
        fastq) 
            if [ -f "${sam_name}_1.fastq" ];
            then 
                count=$(($(wc -l < ${sam_name}_1.fastq)/4))
                echo ${sam_name}_1.fastq,${sam_name}_2.fastq,$((count*2)) >> $metagen_work/reads_info.txt
            else
                echo Reads files are not paired.
            fi
            ;;
        fasta)
            if [ -f "${sam_name}_1.fasta" ];
            then
                count=`grep -c '>' ${sam_name}_1.fasta`
                echo $f,${f/_1.fasta/_2.fasta},$((count*2)) >> $metagen_work/reads_info.txt
            else
                echo Reads files are not paired.
            fi 
            ;;
        esac
    done
fi


