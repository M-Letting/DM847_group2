# Install package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pdfLaTeX")

BiocManager::install("msa")
library(PdfLaTeX)
library(msa)
library(Biostrings)
BiocManager::install("Biostrings")
?msa::msaConvert
# Read in sequences
human_seqs <- "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/human.fa"
mouse_seqs <- "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/mouse.fa"

# Read in BLOSSUM62 substitution matrix
data(BLOSUM62)

# Align sequences
human_out <- msa::msa("/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/human.fa",
                      method = c("Muscle"), maxiters = 10, 
                      substitutionMatrix = BLOSUM62, type = "protein",
                      verbose = TRUE)

mouse_out <- msa::msa(mouse_seqs, method = c("Muscle"), maxiters = 10,
                      substitutionMatrix = BLOSUM62, type = "protein",
                      verbose = TRUE)
# Saving the alignments as r objects
save(human_out, 
     file = "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/human_out.RData")
save(mouse_out, 
     file = "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/mouse_out.RData")

# Saving the alignments as fasta files
alignment1 <- as(human_out, "AAStringSet")
writeXStringSet(alignment1, 
                filepath = "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/human_fasta.fasta", 
                format = "fasta")
alignemnt2 <- as(mouse_out, "AAStringSet")
writeXStringSet(alignemnt2, 
                filepath = "/Users/malthe/Desktop/0Kandidat/DM847/DM847_project/project_sequences/mouse_fasta.fasta", 
                format = "fasta")


# Function to calculate consensus sequence
consensus_sequence <- function(x) {
  # Initialize an empty string to store the consensus sequence
  consensus <- ""
  
  # Loop over each column (position)
  for (col in 1:ncol(x)) {
    # Find the residue with the highest frequency in this column
    max_residue <- rownames(x)[which.max(x[, col] / sum(x[, col]))]
    max_freq <- max(x[, col] / sum(x[, col]))
    
    # Determine the consensus character based on the rules provided
    if (max_freq > 0.8) {
      # More than 80% of the values are the same residue, use uppercase
      consensus <- paste0(consensus, toupper(max_residue))
    } else if (max_freq > 0.2) {
      # More than 20% but less than 80%, use lowercase
      consensus <- paste0(consensus, tolower(max_residue))
    } else {
      # No residue makes up more than 20% of the values, use '.'
      consensus <- paste0(consensus, ".")
    }
  }
  
  # Return the consensus sequence as a single string
  return(consensus)
}

# Calculate the consensus sequence for the alignments
human_cons_matrix <- as.data.frame(consensusMatrix(human_out))
consesus_seq_human <- consensus_sequence(human_cons_matrix)

mouse_cons_matrix <- as.data.frame(consensusMatrix(mouse_out))
consesus_seq_mouse <- consensus_sequence(mouse_cons_matrix)

consesus_seq_human
