# Load the ggplot2 library if you haven't already
library(ggplot2)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the number of arguments is correct
if (length(args) != 3) {
  cat("Usage: Rscript density_plots.R <cell_name> <fill_color> <input_file>\n")
  q()
}

# Extract the arguments
cell_name <- args[1]
fill_color <- args[2]
input_file <- args[3]

# Read the data from the specified input file
data <- read.table(input_file, header = FALSE, sep = ":")

# Extract the second column of numbers
numbers <- data$V2

# Create a density plot with the specified fill color and cell name
plot_title <- cell_name
plot <- ggplot(data, aes(x = numbers)) +
  geom_density(fill = fill_color, color = "black") +
  labs(x = "Sequence Counts", y = "Density") +
  ggtitle(plot_title)

# Generate an output file name based on the cell name and fill color
output_file <- paste0(
  cell_name, "_density", ".png"
)

# Save the plot to the output file
ggsave(output_file, plot, width = 6, height = 4)

cat(paste("Density plot saved as:", output_file, "\n"))
