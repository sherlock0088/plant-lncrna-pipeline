#!/bin/bash
set -eu  # 添加严格模式，防止未处理的错误

# Usage function to display help message
usage() {
    echo "Usage: $0 -g <genome_fasta> -s <strandness> -l <sample_list> -p <threads> -f <annotation_gtf> -w <work_dir>"
    echo "  -g: Path to the reference genome .fasta or .fa file"
    echo "  -s: Strandness (SS or normal)"
    echo "  -l: Path to the sample list file (.txt or .tsv)"
    echo "  -p: Number of threads to use"
    echo "  -f: Path to the annotation .gtf file"
    echo "  -w: Path to the working directory for storing intermediate files"
}

# Initialize variables
REFERENCE=""
STRANDNESS=""
SAMPLE_LIST=""
THREADS=1  # Default to 1 thread
ANNOTATION_GTF=""
WORK_DIR=""

# Parse options using getopts
while getopts "g:s:l:p:f:w:" opt; do
    case $opt in
        g) REFERENCE="$OPTARG" ;;
        s) STRANDNESS="$OPTARG" ;;
        l) SAMPLE_LIST="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        f) ANNOTATION_GTF="$OPTARG" ;;
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

# 创建日志文件
LOG_FILE="$WORK_DIR/pipeline.log"
exec > >(tee -i "$LOG_FILE") 2>&1  # 将输出和错误重定向到日志文件

# Create working directory if it doesn't exist
mkdir -p "$WORK_DIR"

# Step 1: Convert GFF3 to GTF if necessary
GTF_CONVERTED="$WORK_DIR/annotation_converted.gtf"
if [[ ! -f "$GTF_CONVERTED" ]]; then
    if [[ "$ANNOTATION_GTF" == *.gff || "$ANNOTATION_GTF" == *.gff3 ]]; then
        echo "Converting GFF3 to GTF format..."
        gffread "$ANNOTATION_GTF" -T -o "$GTF_CONVERTED"
        if [[ $? -ne 0 ]]; then
            echo "Error: GFF3 to GTF conversion failed."
            exit 1
        fi
        ANNOTATION_GTF="$GTF_CONVERTED"
    fi
else
    ANNOTATION_GTF="$GTF_CONVERTED"
    echo "Using existed GTF file."
fi


# Step 2: Build HiSat2 index if not present
INDEX_FILE="$WORK_DIR/genome.index.1.ht2"
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "Building HiSat2 index..."
    hisat2-build -p "$THREADS" "$REFERENCE" "$WORK_DIR/genome.index"
    if [[ $? -ne 0 ]]; then
        echo "Error: HiSat2 index building failed."
        exit 1
    fi
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
        if [[ $? -ne 0 ]]; then
            echo "Error: HiSat2 alignment failed for sample $sample."
            exit 1
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
        if [[ $? -ne 0 ]]; then
            echo "Error: SAM to BAM sorting failed for sample $sample."
            exit 1
        fi
    fi
done < "$SAMPLE_LIST"

# Step 5: Assemble transcripts with StringTie
while read -r READ1 READ2; do
    sample=$(basename "$READ1" | sed 's/_1.*//')
    GTF_FILE="$WORK_DIR/${sample}.gtf"
    if [[ ! -f "$GTF_FILE" ]]; then
        echo "Running StringTie for $sample..."
        stringtie -p "$THREADS" -G "$ANNOTATION_GTF" -o "$GTF_FILE" "$WORK_DIR/${sample}.bam"
        if [[ $? -ne 0 ]]; then
            echo "Error: StringTie failed for sample $sample."
            exit 1
        fi
    fi
done < "$SAMPLE_LIST"

# Step 6: Merge transcripts
MERGED_GTF="$WORK_DIR/candidate_transcript.gtf"
if [[ ! -f "$MERGED_GTF" ]]; then
    echo "Merging transcripts with StringTie..."
    stringtie --merge -p "$THREADS" -G "$ANNOTATION_GTF" -o "$MERGED_GTF" "$WORK_DIR"/*.gtf
    if [[ $? -ne 0 ]]; then
        echo "Error: Transcript merging failed."
        exit 1
    fi
fi

# Step 7: Extract transcript sequences
TRANSCRIPT_FASTA="$WORK_DIR/candidate_transcript.fasta"
TRANSCRIPT_FASTA_TXT="$WORK_DIR/candidate_transcript.txt"
if [[ ! -f "$TRANSCRIPT_FASTA_TXT" ]]; then
    echo "Extracting transcript sequences with gffread..."
    gffread -w "$TRANSCRIPT_FASTA" -g "$REFERENCE" "$MERGED_GTF"
    if [[ $? -ne 0 ]]; then
        echo "Error: Transcript sequence extraction failed."
        exit 1
    fi
    grep '>' "$TRANSCRIPT_FASTA" | awk '{print $1}' | sed 's/>//g' | sort -u > $TRANSCRIPT_FASTA_TXT
    if [[ $? -ne 0 ]]; then
        echo "Error: Transcript ID extraction from FASTA failed."
        exit 1
    fi
fi

# Step 8: Filter candidate lncRNAs using FEELnc
CANDIDATE_LNCRNA_GTF="$WORK_DIR/candidate_lncRNA.gtf"
CANDIDATE_LNCRNA_TXT="$WORK_DIR/candidate_lncRNA_FEELnc.txt"
if [[ ! -f "$CANDIDATE_LNCRNA_TXT" ]]; then
    echo "Filtering lncRNAs with FEELnc..."
    FEELnc_filter.pl -i "$MERGED_GTF" -a "$ANNOTATION_GTF" --monoex=-1 -s 200 -p "$THREADS" > "$CANDIDATE_LNCRNA_GTF"
    if [[ $? -ne 0 ]]; then
        echo "Error: FEELnc filtering failed."
        exit 1
    fi
    cut -d ";" -f 2 "$CANDIDATE_LNCRNA_GTF" | sed 's/ transcript_id //g' | sed 's/"//g' | sort -u > $CANDIDATE_LNCRNA_TXT
    if [[ $? -ne 0 ]]; then
        echo "Error: Candidate lncRNA extraction failed."
        exit 1
    fi
fi

# Step 9: Run LncFinder
LNCFINDER_RESULTS="$WORK_DIR/plant-lncFinder.txt"
if [[ ! -f "$LNCFINDER_RESULTS" ]]; then
    echo "Running LncFinder..."
    Rscript --vanilla "$DATA_DIR/lncfinder.R" -m "$DATA_DIR/training_mRNA.fasta" -l "$DATA_DIR/training_lncRNA.fasta" -o "$DATA_DIR/Plant_model.rda" -c "$TRANSCRIPT_FASTA" -p "$THREADS" -r "$LNCFINDER_RESULTS"
    if [[ $? -ne 0 ]]; then
        echo "Error: LncFinder failed."
        exit 1
    fi
fi

# Step 10: Run CPAT
CPAT_OUTPUT="$WORK_DIR/CPAT_plant.output"
if [[ ! -f "$CPAT_OUTPUT" ]]; then
    echo "Running CPAT..."
    cpat.py -x "$DATA_DIR/Plant_Hexamer.tsv" -d "$DATA_DIR/Plant.logit.RData" -g "$TRANSCRIPT_FASTA" -o "$CPAT_OUTPUT"
    if [[ $? -ne 0 ]]; then
        echo "Error: CPAT failed."
        exit 1
    fi
fi

# Step 11: Diamond blastx search against UniProt
UNIPROT_OUTPUT="$WORK_DIR/uniprotoutput.txt"
if [[ ! -f "$UNIPROT_OUTPUT" ]]; then
    echo "Running Diamond blastx..."
    diamond blastx -d "$DATA_DIR/uniprot_out" -q "$TRANSCRIPT_FASTA" -o "$UNIPROT_OUTPUT" -p "$THREADS"
    if [[ $? -ne 0 ]]; then
        echo "Error: Diamond blastx search failed."
        exit 1
    fi
fi

# Step 12: Intersection of results
FINAL_LNCRNA_RESULTS="final_lncRNA_results.txt"
if [[ ! -f "$FINAL_LNCRNA_RESULTS" ]]; then
    echo "Performing intersection analysis..."
    Rscript --vanilla "$DATA_DIR/insersection.sh" "$TRANSCRIPT_FASTA_TXT" "$CANDIDATE_LNCRNA_TXT" "$CPAT_OUTPUT" "$LNCFINDER_RESULTS" "$UNIPROT_OUTPUT"
    if [[ $? -ne 0 ]]; then
        echo "Error: Intersection analysis failed."
        exit 1
    fi
fi

# Step 13: Filter final lncRNAs
LNC_RNA_GTF="$WORK_DIR/lncRNA.gtf"
if [[ ! -f "$LNC_RNA_GTF" ]]; then
    echo "Filtering final lncRNAs from GTF file..."
    grep -Fwf "$FINAL_LNCRNA_RESULTS" "$MERGED_GTF" > "$LNC_RNA_GTF"
    if [[ $? -ne 0 ]]; then
        echo "Error: Filtering final lncRNAs failed."
        exit 1
    fi
fi

# Step 14: Classify lncRNAs using FEELnc
LNC_RNA_CLASSES="$WORK_DIR/lncRNA_classes.txt"
if [[ ! -f "$LNC_RNA_CLASSES" ]]; then
    echo "Classifying lncRNAs using FEELnc..."
    FEELnc_classifier.pl -i "$LNC_RNA_GTF" -a "$ANNOTATION_GTF" > "$LNC_RNA_CLASSES"
    if [[ $? -ne 0 ]]; then
        echo "Error: FEELnc classification failed."
        exit 1
    fi
fi

# Step 15: Extract specific lncRNA types
if [[ ! -f "$WORK_DIR/LncRNA_antisense_exonic.txt" ]]; then
    echo "Classifying and extracting specific lncRNA types..."
    awk -F '\t' '{if(NF >= 10 && $1==1 && $6 == "antisense" && $10 == "exonic") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_antisense_exonic.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting antisense exonic lncRNAs failed."
        exit 1
    fi

    awk -F '\t' '{if(NF >= 10 && $1==1 && $6 == "sense" && $10 == "intronic") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_intronic.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting intronic lncRNAs failed."
        exit 1
    fi

    awk -F '\t' '{if(NF >= 10 && $1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_upstream.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting upstream lncRNAs failed."
        exit 1
    fi

    awk -F '\t' '{if(NF >= 10 && $1==1 && $6 == "sense" && $7 == "intergenic" && $8 <= 2000 && $10 == "downstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_downstream.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting downstream lncRNAs failed."
        exit 1
    fi

    awk -F '\t' '{if(NF >= 10 && $1==1 && $7 == "intergenic" && $8 > 2000) {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_intergenic.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting intergenic lncRNAs failed."
        exit 1
    fi

    awk -F '\t' '{if(NF >= 10 && $1==1 && $6 == "antisense" && $7 == "intergenic" && $8 <= 2000 && $10 == "upstream") {print $0}}' "$LNC_RNA_CLASSES" > "$WORK_DIR/LncRNA_Bidirectional.txt"
    if [[ $? -ne 0 ]]; then
        echo "Error: Extracting bidirectional lncRNAs failed."
        exit 1
    fi
else
    echo "Specific lncRNA types already extracted, skipping this step."
fi

echo "All processes completed! All files saved in $WORK_DIR."

