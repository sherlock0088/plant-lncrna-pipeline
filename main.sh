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
THREADS=1  # Default to 1 thread
ANNOTATION_GTF=""
WORK_DIR=""

# Data directory path (non-input data and scripts)
DATA_DIR="./data"

# Parse options using getopts
while getopts "g:s:l:p:gff:w:" opt; do
    case $opt in
        g) REFERENCE="$OPTARG" ;;
        s) STRANDNESS="$OPTARG" ;;
        l) SAMPLE_LIST="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        gff) ANNOTATION_GTF="$OPTARG" ;;
        w) WORK_DIR="$OPTARG" ;;
        \?) usage; exit 1 ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$REFERENCE" || -z "$STRANDNESS" || -z "$SAMPLE_LIST" || -z "$ANNOTATION_GTF" || -z "$WORK_DIR" ]]; then
    echo "Error: Missing required argument(s)"
    usage
    exit 1
fi

# Create working directory if it doesn't exist
mkdir -p "$WORK_DIR"

# Step 1: Convert GFF3 to GTF if necessary
if [[ "$ANNOTATION_GTF" == *.gff || "$ANNOTATION_GTF" == *.gff3 ]]; then
    echo "Converting GFF3 to GTF format..."
    GTF_CONVERTED="$WORK_DIR/annotation_converted.gtf"
    gffread "$ANNOTATION_GTF" -T -o "$GTF_CONVERTED"
    ANNOTATION_GTF="$GTF_CONVERTED"
fi

# Step 2: Build HiSat2 index if not present
INDEX_FILE="$WORK_DIR/genome.index.1.ht2"
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "Building HiSat2 index..."
    hisat2-build -p "$THREADS" "$REFERENCE" "$WORK_DIR/genome.index"
fi

# Step 3: Align reads with HiSat2 for each sample
while read -r READ1 READ2; do
    sample=$(basename "$READ1" | sed 's/_1.*//')
    SAM_FILE="$WORK_DIR/${sample}.sam"
    if [[ ! -f "$SAM_FILE" ]]; then
        echo "Running HiSat2 for $sample..."
        if [[ "$STRANDNESS" == "SS" ]]; then
            hisat2 --rna-strandness RF -p "$THREADS" -x "$WORK_DIR/genome.index" -1 "$READ1" -2 "$READ2" -S "$SAM_FILE"
        else
            hisat2 -p "$THREADS" -x "$WORK_DIR/genome.index" -1 "$READ1" -2 "$READ2" -S "$SAM_FILE"
        fi
    fi
done < "$SAMPLE_LIST"

# Step 4: Sort SAM files into BAM
while read -r READ1 READ2; do
    sample=$(basename "$READ1" | sed 's/_1.*//')
    BAM_FILE="$WORK_DIR/${sample}.bam"
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "Sorting SAM file into BAM for $sample..."
        samtools sort -@ "$THREADS" -o "$BAM_FILE" "$WORK_DIR/${sample}.sam"
    fi
done < "$SAMPLE_LIST"

# Step 5: Assemble transcripts with StringTie
while read -r READ1 READ2; do
    sample=$(basename "$READ1" | sed 's/_1.*//')
    GTF_FILE="$WORK_DIR/${sample}.gtf"
    if [[ ! -f "$GTF_FILE" ]]; then
        echo "Running StringTie for $sample..."
        stringtie -p "$THREADS" -G "$ANNOTATION_GTF" -o "$GTF_FILE" "$WORK_DIR/${sample}.bam"
    fi
done < "$SAMPLE_LIST"

# Step 6: Merge transcripts
MERGED_GTF="$WORK_DIR/candidate_transcript.gtf"
if [[ ! -f "$MERGED_GTF" ]]; then
    echo "Merging transcripts with StringTie..."
    stringtie --merge -p "$THREADS" -G "$ANNOTATION_GTF" -o "$MERGED_GTF" "$WORK_DIR"/*.gtf
fi

# Step 7: Extract transcript sequences
TRANSCRIPT_FASTA="$WORK_DIR/candidate_transcript.fasta"
if [[ ! -f "$TRANSCRIPT_FASTA" ]]; then
    echo "Extracting transcript sequences with gffread..."
    gffread -w "$TRANSCRIPT_FASTA" -g "$REFERENCE" "$MERGED_GTF"
fi

# Step 8: Filter candidate lncRNAs using FEELnc
CANDIDATE_LNCRNA_GTF="$WORK_DIR/candidate_lncRNA.gtf"
if [[ ! -f "$CANDIDATE_LNCRNA_GTF" ]]; then
    echo "Filtering lncRNAs with FEELnc..."
    FEELnc_filter.pl -i "$MERGED_GTF" -a "$ANNOTATION_GTF" --monoex=-1 -s 200 -p "$THREADS" > "$CANDIDATE_LNCRNA_GTF"
fi

# Step 9: Run LncFinder
LNCFINDER_RESULTS="$WORK_DIR/plant-lncFinder.txt"
if [[ ! -f "$LNCFINDER_RESULTS" ]]; then
    echo "Running LncFinder..."
    Rscript --vanilla "$DATA_DIR/lncfinder.R" -m "$DATA_DIR/training_mRNA.fasta" -l "$DATA_DIR/training_lncRNA.fasta" -o "$DATA_DIR/Plant_model.rda" -c "$TRANSCRIPT_FASTA" -p "$THREADS" -r "$LNCFINDER_RESULTS"
fi

# Step 10: Run CPAT
CPAT_OUTPUT="$WORK_DIR/CPAT_plant.output"
if [[ ! -f "$CPAT_OUTPUT" ]]; then
    echo "Running CPAT..."
    python2.7 -m cpat -x "$DATA_DIR/Plant_Hexamer.tsv" -d "$DATA_DIR/Plant.logit.RData" -g "$TRANSCRIPT_FASTA" -o "$CPAT_OUTPUT"
fi

# Step 11: Diamond blastx search against UniProt
UNIPROT_OUTPUT="$WORK_DIR/uniprotoutput.txt"
if [[ ! -f "$UNIPROT_OUTPUT" ]]; then
    echo "Running Diamond blastx..."
    diamond blastx -d "$DATA_DIR/uniprot_out" -q "$TRANSCRIPT_FASTA" -o "$UNIPROT_OUTPUT" -p "$THREADS"
fi

# Step 12: Intersection of results
FINAL_LNCRNA_RESULTS="$WORK_DIR/final_lncRNA_results.txt"
if [[ ! -f "$FINAL_LNCRNA_RESULTS" ]]; then
    echo "Performing intersection analysis..."
    Rscript --vanilla "$DATA_DIR/insersection.sh" "$WORK_DIR/candidate_lncRNA.txt" "$CPAT_OUTPUT" "$LNCFINDER_RESULTS" "$UNIPROT_OUTPUT" > "$FINAL_LNCRNA_RESULTS"
fi

# Step 13: Filter final lncRNAs
LNC_RNA_GTF="$WORK_DIR/lncRNA.gtf"
if [[ ! -f "$LNC_RNA_GTF" ]]; then
    echo "Filtering final lncRNAs from GTF file..."
    grep -Fwf "$FINAL_LNCRNA_RESULTS" "$MERGED_GTF" > "$LNC_RNA_GTF"
fi

# Step 14: Classify lncRNAs using FEELnc
LNC_RNA_CLASSES="$WORK_DIR/lncRNA_classes.txt"
if [[ ! -f "$LNC_RNA_CLASSES" ]]; then
    echo "Classifying lncRNAs using FEELnc..."
    FEELnc_classifier.pl -i "$LNC_RNA_GTF" -a "$ANNOTATION_GTF" > "$LNC_RNA_CLASSES"
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
