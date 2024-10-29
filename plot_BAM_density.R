################################################################################
# Plot BAM density for the entire genome
#
# Author: Vladimir Bajić
# Date: 2024-10-29
#
# Description:
#   - Plot pdf showing BAM density for the entire genome
#   - Allows filtering contigs smaller than defined threshold
#
# Usage:
#
# To see help message
#   Rscript --vanilla plot_BAM_density.R --help
#
# To plot with defaults
#   Rscript --vanilla plot_BAM_density.R -f input.fasta -b input.bam -o out_path
#
################################################################################


# Library ----------------------------------------------------------------------
suppressMessages(library(karyoploteR))
suppressMessages(library(Biostrings))
library(optparse)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-f", "--fasta_file"), type = "character", help = "Path to the fasta file", metavar = "character"),
    make_option(c("-b", "--bam_file"), type = "character", help = "Path to the BAM file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix (including output path)", metavar = "character"),
    make_option(c("-s", "--window_size"), type = "numeric", help = "The size of the windows for wich the density is computed using karyoploteR::kpPlotBAMDensity [Default is set to 1000 bp]", metavar = "numeric", default = 1000),
    make_option(c("-m", "--min_contig_size"), type = "numeric", help = "Plot only contigs bigger than specified value of min_contig_size [Default is set to 10k]", metavar = "numeric", default = 10000),
    make_option(c("-x", "--plot_width"), type = "numeric", help = "Plot width [Default 40]", metavar = "numeric", default = 40),
    make_option(c("-y", "--plot_height"), type = "numeric", help = "Plot height [Default 80]", metavar = "numeric", default = 80)
)

# Parsing options --------------------------------------------------------------
opt_parser <- OptionParser(option_list = option_list, add_help_option = TRUE)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$fasta_file)) {
    print_help(opt_parser)
    stop("Provide path to the fasta file.")
}

if (is.null(opt$bam_file)) {
    print_help(opt_parser)
    stop("Provide path to the BAM file.")
}

if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output path must be provided.")
}

# ------------------------------------------------------------------------------

# Define function to plot BAM density ------------------------------------------
plot_BAM_density <- function(fasta_file, bam_file, out_plot, window_size, min_contig_size, plot_width, plot_height) {
    # Load your FASTA file and extract sequence lengths
    fasta <- readDNAStringSet(fasta_file)

    # Create a GRanges object with chromosome names and lengths
    seq_lengths <- seqlengths(fasta)
    custom_genome <- GRanges(
        seqnames = names(seq_lengths),
        ranges = IRanges(start = 1, end = seq_lengths)
    )

    # Filter custom_genome to include only ranges with length greater than 10,000 bp
    filtered_genome <- custom_genome[width(custom_genome) > min_contig_size]

    # Initialize the karyoplot with the custom genome
    pdf(
        file = paste0(out_plot, "_BAM_density.pdf"),
        width = plot_width,
        height = plot_height
    )

    kp <- plotKaryotype(genome = filtered_genome, chromosomes = "all")
    kp <- kpAddBaseNumbers(kp, tick.dist = 10000, add.units = TRUE)
    kp <- kpPlotBAMDensity(kp, data = bam_file, window.size = window_size, col = "darkorange")
    kp <- kpAxis(kp, ymax = kp$latest.plot$computed.values$max.density)

    dev.off()
}

# Plot and save ----------------------------------------------------------------
plot_BAM_density(opt$fasta_file, opt$bam_file, opt$output, opt$window_size, opt$min_contig_size, opt$plot_width, opt$plot_height)

################################################################################
