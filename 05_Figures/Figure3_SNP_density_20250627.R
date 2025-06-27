# SNP density

# Clean environment
rm(list = ls())

# Load necessary libraries
library(vcfR)
library(tidyverse)
#library(dplyr)
#library(ggplot2)
library(scales)
library(cowplot)



#---------------------------------------------------------------------------------------------------------
#   BASED ON MULTIHETSEP
#---------------------------------------------------------------------------------------------------------
# Set directory
#setwd("~/Documents/03_Master/Thesis/CRIP/multihetsep-files/")

# Step 1: Read the multihetsep file
# Replace 'your_file.tabular' with the path to your multihetsep file
snps1 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr1/multihetsep_MF1_MF2_MF3_MF4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                   comment.char = "#")
snps2 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr2/multihetsep_MF1_MF2_MF3_MF4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps3 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr3/multihetsep_MF1_MF2_MF3_MF4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps4 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr4/multihetsep_MF1_MF2_MF3_MF4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps5 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr1/multihetsep_MG2_MG3_MG4_MG5_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps6 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr2/multihetsep_MG2_MG3_MG4_MG5_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps7 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr3/multihetsep_MG2_MG3_MG4_MG5_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps8 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr4/multihetsep_MG2_MG3_MG4_MG5_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps9 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr1/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps10 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr2/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps11 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr3/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps12 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr4/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps13 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr1/multihetsep_SI1_SI2_SI3_SI4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps14 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr2/multihetsep_SI1_SI2_SI3_SI4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps15 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr3/multihetsep_SI1_SI2_SI3_SI4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps16 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr4/multihetsep_SI1_SI2_SI3_SI4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

snps17 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr1/multihetsep_SS1_SS2_SS3_SS4_Chr1.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps18 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr2/multihetsep_SS1_SS2_SS3_SS4_Chr2.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps19 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr3/multihetsep_SS1_SS2_SS3_SS4_Chr3.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")
snps20 <- read.table("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/01_Revison/Sorted_final_files/Final_Supplement/Zenodo/multihetsep-files/multihetsep-Chr4/multihetsep_SS1_SS2_SS3_SS4_Chr4.txt", sep = "\t", header = FALSE, blank.lines.skip = TRUE,
                    comment.char = "#")

# Put them together
mf <- rbind(snps1,snps2,snps3,snps4)
mg <- rbind(snps5, snps6,snps7,snps8)
nmf <- rbind(snps9,snps10,snps11,snps12)
si <- rbind(snps13, snps14,snps15,snps16)
ss <- rbind(snps17, snps18,snps19,snps20)

# Add population label
mf$pop <- "MF"
mg$pop <- "MG"
nmf$pop <- "NMF"
si$pop <- "SI"
ss$pop <- "SS"

# Bind them all
snps <- rbind(mf,mg,nmf,si,ss)

# Step 2: Assign column names
colnames(snps) <- c("chr", "start", "nosegsitesincelast", "refallele", "pop")

# Step 3: Filter out rows with invalid SNP data (if necessary)

# Step 4: Define the correct order of chromosomes and populations
goodChrOrder <- c(paste("Chr", c(1:4), sep = ""))
snps$chr <- factor(snps$chr, levels = goodChrOrder)

goodPopOrder <- c("MG", "NMF", "MF", "SI", "SS")
snps$pop <- factor(snps$pop, levels = goodPopOrder)

snps <- snps %>% 
  mutate(pop = case_when(
    pop == "MG" ~ "Hesse (MG)",
    pop == "NMF" ~ "Lorraine (NMF)",
    pop == "MF" ~ "Rhône-Alpes (MF)",
    pop == "SI" ~ "Piedmont (SI)",
    pop == "SS" ~ "Andalusia (SS)",
    TRUE ~ pop
  )
  )


goodPopOrder <- c("Hesse (MG)", "Lorraine (NMF)", "Rhône-Alpes (MF)", "Piedmont (SI)", "Andalusia (SS)")
snps$pop <- factor(snps$pop, levels = goodPopOrder)


# Step 5: Plot the SNP density using ggplot2
snpDensity <- ggplot(snps, aes(x = start/1e6,  fill = pop)) + 
  geom_histogram(aes(y = after_stat(count)/(50e3/1e6)), binwidth = 50e3/1e6) +  # Adjust bin width as needed
  facet_grid(pop ~ chr, scales = "free_x") +  # Separate plots for each chromosome
  #ggtitle("Chromosome-wise SNP distribution based on multihetsep files") + 
  xlab("Position (Mb)") + 
  ylab("SNP density (per Mb)") + 
  labs(fill = "Population") +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20,face = "bold", color = "black"),
        axis.title.x = element_text(size = 20,face = "bold", color = "black"),
        title = element_text(size = 20, color = "black"),
        text = element_text(size = 20, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("white"),
        panel.grid.minor = element_line("white"),
        panel.grid.major.x = element_line("lightgrey", linetype = 3),
        panel.grid.major.y = element_line("lightgrey", linetype = 3),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.position="bottom",
        legend.text = element_text(size = 20, color = "black"),
        #strip.text.y = element_text(angle = 0),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 20))


my_palette <- c("Rhône-Alpes (MF)" = "#4C83EB",   # Assign colors to populations
                "Lorraine (NMF)" = "#CD70C6",
                "Hesse (MG)" = "#3AB6D8",
                "Piedmont (SI)" = "#FF941A",
                "Andalusia (SS)" = "#AB3232")

scientific_10 <- function(x) { 
  ifelse(x == 0, "0", parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))) 
} #define scale

snpDensity <- snpDensity + scale_fill_manual(values = my_palette) + scale_x_continuous(label = comma) + scale_y_continuous(labels = scientific_10) 

snpDensity <- snpDensity + theme(axis.text.y = element_text(angle = 45, hjust = 1), panel.spacing = unit(1.5, "lines"))

snpDensity

cowplot::save_plot("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/02_Revision_2/G3_GenomeReport/Plots/SNP_density.svg", snpDensity, nrow = 5, ncol = 5, base_asp = 1.2, dpi = 600, bg = "white", scale = 0.7)
cowplot::save_plot("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/02_Revision_2/G3_GenomeReport/Plots/SNP_density.png", snpDensity, nrow = 5, ncol = 5, base_asp = 1.2, dpi = 600, bg = "white", scale = 0.7)


# Total SNP density
mf_density <- snps %>%
  filter(pop == "MF")

total_snps <- nrow(mf_density)
genome_size_mb <- 192 
average_density <- total_snps / genome_size_mb
average_density



# Step 6: Plot the SNP count using ggplot2
snpDensity <- ggplot(snps, aes(x = start/1e6,  fill = pop)) + 
  geom_histogram(binwidth = 50e3/1e6) +  # Adjust bin width as needed
  facet_grid(pop ~ chr, scales = "free_x") +  # Separate plots for each chromosome
  #ggtitle("Chromosome-wise SNP distribution based on multihetsep files") + 
  xlab("Position (Mb)") + 
  ylab("SNP count") + 
  labs(fill = "Population") +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20,face = "bold", color = "black"),
        axis.title.x = element_text(size = 20,face = "bold", color = "black"),
        title = element_text(size = 20, color = "black"),
        text = element_text(size = 20, color = "black"),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line("white"),
        panel.grid.minor = element_line("white"),
        panel.grid.major.x = element_line("lightgrey", linetype = 3),
        panel.grid.major.y = element_line("lightgrey", linetype = 3),
        panel.border = element_rect("black", fill = NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_rect(fill="white"),
        legend.position="bottom",
        #strip.text.y = element_text(angle = 0),
        strip.text.y = element_blank())


my_palette <- c("Rhône-Alpes (MF)" = "#4C83EB",   # Assign colors to populations
                "Lorraine (NMF)" = "#CD70C6",
                "Hesse (MG)" = "#3AB6D8",
                "Piedmont (SI)" = "#FF941A",
                "Andalusia (SS)" = "#AB3232")

scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", 
                                                 scales::scientific_format()(x))) } # define scale

snpDensity <- snpDensity + scale_fill_manual(values = my_palette) + scale_x_continuous(label = comma) #+ scale_y_continuous(labels = scientific_10) 

snpDensity

cowplot::save_plot("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/02_Revision_2/G3_GenomeReport/Plots/SNP_count.svg", snpDensity, nrow = 5, ncol = 5, base_asp = 1.2, dpi = 600, bg = "white", scale = 0.7)
cowplot::save_plot("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/02_Revision_2/G3_GenomeReport/Plots/SNP_count.png", snpDensity, nrow = 5, ncol = 5, base_asp = 1.2, dpi = 600, bg = "white", scale = 0.7)

