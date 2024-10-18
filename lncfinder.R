# Load necessary libraries
library(LncFinder)
library(seqinr)
library(optparse)

# Define command line options
option_list <- list(
    make_option(c("-m", "--mRNA"), type = "character", default = NULL, 
                help = "Path to mRNA FASTA file", metavar = "character"),
    make_option(c("-l", "--lncRNA"), type = "character", default = NULL, 
                help = "Path to lncRNA FASTA file", metavar = "character"),
    make_option(c("-o", "--model"), type = "character", default = NULL, 
                help = "Path to pre-trained SVM model RDS file", metavar = "character"),
    make_option(c("-c", "--candidates"), type = "character", default = "candidate_transcript.fasta",
                help = "Path to candidate transcript FASTA file (default: candidate_transcript.fasta)", metavar = "character"),
    make_option(c("-p", "--cores"), type = "integer", default = 2, 
                help = "Number of cores for parallel processing (default: 2)", metavar = "integer"),
    make_option(c("-r", "--result"), type = "character", default = "plant-lncFinder.txt", 
                help = "Path to output results file (default: plant-lncFinder.txt)", metavar = "character")
)

# Parse the command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required options are provided
if (is.null(opt$mRNA) || is.null(opt$lncRNA) || is.null(opt$model)) {
    print_help(opt_parser)
    stop("Error: Missing required arguments. mRNA, lncRNA, and model files must be specified.", call. = FALSE)
}

# Load mRNA and lncRNA sequences
mRNA <- seqinr::read.fasta(file = opt$mRNA)
lncRNA <- seqinr::read.fasta(file = opt$lncRNA)

# Create frequency features from the provided mRNA and lncRNA sequences
frequencies <- make_frequencies(cds.seq = mRNA, lncRNA.seq = lncRNA, SS.features = FALSE, 
                                cds.format = "DNA", lnc.format = "DNA", check.cds = TRUE, ignore.illegal = TRUE)

# Load the pre-trained SVM model
plant <- readRDS(opt$model)

# Load candidate transcript sequences
Seqs <- seqinr::read.fasta(file = opt$candidates)

# Run LncFinder on the candidate transcripts
Plant_results <- LncFinder::lnc_finder(Seqs, SS.features = FALSE, format = "DNA", 
                                       frequencies.file = frequencies, svm.model = plant, parallel.cores = opt$cores)

# Write the results to a file
write.table(Plant_results, file = opt$result, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Print message when finished
cat("Results saved to", opt$result, "\n")
