#!/bin/bash
set -eu

# Usage function to display help message
usage() {
    echo "Usage: $0 -g <genome_fasta> -s <strandness> -l <sample_list> -p <threads> -gff <annotation_gtf> -w <work_dir>"
    echo "  -g: Path to the reference genome .fasta or .fa file"
    echo "  -s: Strandness (SS or normal)"
    echo "  -l: Path to the sample list file (.txt or .tsv)"
    echo "  -p: Number of threads to use"
    echo "  -gff: Path to the annotation .gtf file"
    echo "  -w: Path to the working directory for storing intermediate files"
}

# Initialize variables
REFERENCE=""
STRANDNESS=""
SAMPLE_LIST=""
THREADS=1  # 默认线程数为 1
ANNOTATION_GTF=""
WORK_DIR=""

# Define the data directory path (non-input data and scripts are here)
DATA_DIR="./data"

# Parse options using getopts
while getopts "g:s:l:p:gff:w:" opt; do
    case $opt in
        g)
            REFERENCE="$OPTARG"
            ;;
        s)
            STRANDNESS="$OPTARG"
            ;;
        l)
            SAMPLE_LIST="$OPTARG"
            ;;
        p)
            THREADS="$OPTARG"
            ;;
        gff)
            ANNOTATION_GTF="$OPTARG"
            ;;
        w)
            WORK_DIR="$OPTARG"
            ;;
        \?)
            usage
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$REFERENCE" ]]; then
    echo >&2 "error: reference genome file not specified"
    usage
    exit 1
fi

if [[ -z "$STRANDNESS" ]]; then
    echo >&2 "error: strandness not specified"
    usage
    exit 1
fi

if [[ -z "$SAMPLE_LIST" ]]; then
    echo >&2 "error: sample list file not specified"
    usage
    exit 1
fi

if [[ -z "$ANNOTATION_GTF" ]]; then
    echo >&2 "error: annotation GTF file not specified"
    usage
    exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
    echo >&2 "error: working directory not specified"
    usage
    exit 1
fi

# Create working directory if it doesn't exist
mkdir -p "$WORK_DIR"
# Check if the annotation file is GFF3 or GTF format and convert to GTF if necessary
if [[ "$ANNOTATION_GTF" == *.gff || "$ANNOTATION_GTF" == *.gff3 ]]; then
    echo "The provided annotation file is GFF3. Converting to GTF format..."
    GTF_CONVERTED="$WORK_DIR/annotation_converted.gtf"
    gffread "$ANNOTATION_GTF" -T -o "$GTF_CONVERTED"
    ANNOTATION_GTF="$GTF_CONVERTED"
    echo "Conversion to GTF completed. Using $ANNOTATION_GTF for downstream analysis."
elif [[ "$ANNOTATION_GTF" == *.gtf ]]; then
    echo "The provided annotation file is already in GTF format."
else
    echo >&2 "error: unsupported annotation file format. Please provide a GFF3 or GTF file."
    usage
    exit 1
fi

# Step 1: Build HiSat2 index
INDEX_FILE="$WORK_DIR/genome.index.1.ht2"
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "Building HiSat2 index..."
    hisat2-build -p "$THREADS" "$REFERENCE" "$WORK_DIR/genome.index"
    echo "HiSat2 index has been successfully built in $WORK_DIR"
else
    echo "HiSat2 index already exists, skipping this step."
fi

# Step 2: Run HiSat2 for each sample
for line in $(cat "$SAMPLE_LIST"); do
    # 使用read命令将每一行的两个列（_1和_2的reads路径）分别存储到变量中
    read -r READ1 READ2 <<< "$line"

    # 使用basename命令提取样本名称，比如 sample1_1.fastq.gz 提取为 sample1
    sample=$(basename "$READ1" | sed 's/_1.*//')

    # 定义 SAM 文件的输出路径
    SAM_FILE="$WORK_DIR/${sample}.sam"

    if [[ ! -f "$SAM_FILE" ]]; then
        if [[ "$STRANDNESS" == "SS" ]]; then
            echo "Running HiSat2 with strandness SS (RF) for $sample..."
            hisat2 --new-summary --rna-strandness RF -p "$THREADS" -x "$WORK_DIR/genome.index" \
                -1 "$READ1" -2 "$READ2" -S "$SAM_FILE"
        else
            echo "Running HiSat2 with no strand specificity for $sample..."
            hisat2 --new-summary -p "$THREADS" -x "$WORK_DIR/genome.index" \
                -1 "$READ1" -2 "$READ2" -S "$SAM_FILE"
        fi
    else
        echo "SAM file for $sample already exists, skipping this step."
    fi
done


# Step 3: Sort SAM files into BAM
while read -r READ1 READ2; do
    # 使用basename命令提取样本名称，比如 sample1_1.fastq.gz 提取为 sample1
    sample=$(basename "$READ1" | sed 's/_1.*//')

    # 定义 BAM 文件的输出路径
    BAM_FILE="$WORK_DIR/${sample}.bam"

    if [[ ! -f "$BAM_FILE" ]]; then
        echo "Sorting SAM file into BAM format for $sample..."
        samtools sort -@ "$THREADS" -o "$BAM_FILE" "$WORK_DIR/${sample}.sam"
    else
        echo "BAM file for $sample already exists, skipping this step."
    fi
done < "$SAMPLE_LIST"

# Step 4: Run StringTie for transcript assembly
while read -r READ1 READ2; do
    # 使用basename命令提取样本名称，比如 sample1_1.fastq.gz 提取为 sample1
    sample=$(basename "$READ1" | sed 's/_1.*//')

    # 定义 GTF 文件的输出路径
    GTF_FILE="$WORK_DIR/${sample}.gtf"

    if [[ ! -f "$GTF_FILE" ]]; then
        echo "Running StringTie for transcript assembly for $sample..."
        if [[ "$STRANDNESS" == "SS" ]]; then
            stringtie -p "$THREADS" --rf -G "$ANNOTATION_GTF" -o "$GTF_FILE" "$WORK_DIR/${sample}.bam"
        else
            stringtie -p "$THREADS" -G "$ANNOTATION_GTF" -o "$GTF_FILE" "$WORK_DIR/${sample}.bam"
        fi
    else
        echo "GTF file for $sample already exists, skipping this step."
    fi
done < "$SAMPLE_LIST"


# Step 5: Merge transcripts using StringTie
MERGED_GTF="$WORK_DIR/candidate_transcript.gtf"
if [[ ! -f "$MERGED_GTF" ]]; then
    echo "Merging transcripts with StringTie..."
    stringtie --merge -p "$THREADS" -o "$MERGED_GTF" -G "$ANNOTATION_GTF" "$WORK_DIR"/*.gtf
else
    echo "Merged transcript file already exists, skipping this step."
fi

# Step 6: Extract transcript sequences with gffread
TRANSCRIPT_FASTA="$WORK_DIR/candidate_transcript.fasta"
if [[ ! -f "$TRANSCRIPT_FASTA" ]]; then
    echo "Extracting transcript sequences with gffread..."
    gffread -w "$TRANSCRIPT_FASTA" -g "$REFERENCE" "$MERGED_GTF"
else
    echo "Transcript FASTA file already exists, skipping this step."
fi

# Step 7: Run FEELnc_filter to identify candidate lncRNAs
CANDIDATE_LNCRNA_GTF="$WORK_DIR/candidate_lncRNA.gtf"
if [[ ! -f "$CANDIDATE_LNCRNA_GTF" ]]; then
    echo "Running FEELnc_filter to filter lncRNAs..."
    FEELnc_filter.pl -i "$MERGED_GTF" -a "$ANNOTATION_GTF" --monoex=-1 -s 200 -p "$THREADS" > "$CANDIDATE_LNCRNA_GTF"
else
    echo "Candidate lncRNA GTF file already exists, skipping this step."
fi

# Step 8: Run LncFinder R script
LNCFINDER_RESULTS="$WORK_DIR/plant-lncFinder.txt"
if [[ ! -f "$LNCFINDER_RESULTS" ]]; then
    echo "Running LncFinder analysis..."
    Rscript --vanilla "$DATA_DIR/lncfinder.R" -m "$DATA_DIR/training_mRNA.fasta" -l "$DATA_DIR/training_lncRNA.fasta" -o "$DATA_DIR/Plant_model.rda" -c "$TRANSCRIPT_FASTA" -p "$THREADS" -r "$LNCFINDER_RESULTS"
else
    echo "LncFinder results already exist, skipping this step."
fi

# Step 9: Run CPAT analysis
CPAT_OUTPUT="$WORK_DIR/CPAT_plant.output"
if [[ ! -f "$CPAT_OUTPUT" ]]; then
    echo "Running CPAT analysis..."
    source activate py27
    python2 "$DATA_DIR/cpat.py" -x "$DATA_DIR/Plant_Hexamer.tsv" -d "$DATA_DIR/Plant.logit.RData" -g "$TRANSCRIPT_FASTA" -o "$CPAT_OUTPUT"
else
    echo "CPAT output already exists, skipping this step."
fi

# Step 10: Run Diamond blastx search against UniProt
UNIPROT_OUTPUT="$WORK_DIR/uniprotoutput.txt"
if [[ ! -f "$UNIPROT_OUTPUT" ]]; then
    echo "Running Diamond blastx against UniProt..."
    diamond blastx -d "$DATA_DIR/uniprot_out" -q "$TRANSCRIPT_FASTA" -o "$UNIPROT_OUTPUT" -p "$THREADS" -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
else
    echo "UniProt output already exists, skipping this step."
fi

# Step 11: Perform intersection of results
FINAL_LNCRNA_RESULTS="$WORK_DIR/final_lncRNA_results.txt"
if [[ ! -f "$FINAL_LNCRNA_RESULTS" ]]; then
    echo "Running R script for intersection analysis..."
    Rscript --vanilla "$DATA_DIR/insersection.sh" "$WORK_DIR/candidate_lncRNA.txt" "$WORK_DIR/candidate_lncRNA.txt" "$CPAT_OUTPUT" "$LNCFINDER_RESULTS" "$UNIPROT_OUTPUT" > "$FINAL_LNCRNA_RESULTS"
else
    echo "Final lncRNA intersection results already exist, skipping this step."
fi

# Step 12: Filter the final lncRNAs from GTF file
LNC_RNA_GTF="$WORK_DIR/lncRNA.gtf"
if [[ ! -f "$LNC_RNA_GTF" ]]; then
    echo "Filtering final lncRNAs from GTF file..."
    grep -Fwf "$FINAL_LNCRNA_RESULTS" "$MERGED_GTF" > "$LNC_RNA_GTF"
else
    echo "Filtered lncRNA GTF already exists, skipping this step."
fi

# Step 13: Classify lncRNAs using FEELnc_classifier
LNC_RNA_CLASSES="$WORK_DIR/lncRNA_classes.txt"
if [[ ! -f "$LNC_RNA_CLASSES" ]]; then
    echo "Classifying lncRNAs using FEELnc_classifier..."
    FEELnc_classifier.pl -i "$LNC_RNA_GTF" -a "$ANNOTATION_GTF" > "$LNC_RNA_CLASSES"
else
    echo "LncRNA classification already exists, skipping this step."
fi

# Step 14: Extract specific lncRNA types
if [[ ! -f "$WORK_DIR/LncRNA_antisense_exonic.txt" ]]; then
    echo "Classifying and extracting specific lncRNA types..."
    awk -F '\t' '{if($1==1 && $6 == "antisense" && $10 == "exonic") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_antisense_exonic.txt"
    awk -F '\t' '{if($1==1 && $6 == "sense" && $10 == "intronic") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_intronic.txt"
    awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_upstream.txt"
    awk -F '\t' '{if($1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "downstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_downstream.txt"
    awk -F '\t' '{if($1==1 && $7 == "intergenic" && $8 > 2000) {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_intergenic.txt"
    awk -F '\t' '{if($1==1 && $6 == "antisense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_Bidirectional.txt"
else
    echo "Specific lncRNA types already extracted, skipping this step."
fi

echo "All processes completed! All files saved in $WORK_DIR."
