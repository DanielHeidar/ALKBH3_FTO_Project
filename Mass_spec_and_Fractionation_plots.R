setwd("C:/Users/danielheidar/OneDrive - Menntaský")
ALKBH3_Mass <- read.csv2("ALKBH3_Mass.csv", header = TRUE, stringsAsFactors = FALSE)

# Load necessary libraries
library(ggplot2)


# Example structure of your data
# ALKBH3_Mass <- data.frame(
#   Protein = c("HSP70", "FTO", "HSP90", "ALKBH5", "NPM1"),
#   Unique_Peptides = c(5, 13, 7, 9, 4)
# )

# Create the barplot
ggplot(ALKBH3_Mass, aes(x = reorder(Protein, -Unique.Peptides), y = Unique.Peptides, fill = Protein == "FTO")) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(
    aes(label = Unique.Peptides),
    vjust = -0.5,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#FF6F61", "FALSE" = "#AEC6CF"),
    name = "Protein Highlight"
  ) +
  labs(
    title = "Unique Peptides Binding to ALKBH3",
    x = "Protein Names",
    y = "Number of Unique Peptides"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 15, face = "bold"),  # Tilt more, bold for clarity
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed")
  ) +
  ylim(0, max(ALKBH3_Mass$Unique.Peptides) + 3) +
  annotate(
    "text",
    x = which(ALKBH3_Mass$Protein == "FTO"),
    y = max(ALKBH3_Mass$Unique.Peptides) + 2,
    label = "Highlighted Protein",
    color = "#FF6F61",
    fontface = "italic",
    size = 4
  )



FTO_Mass <- read.csv2("FTO_Mass.csv", header = TRUE, stringsAsFactors = FALSE)


# Example structure of your data
# FTO_Mass <- data.frame(
#   Protein = c("HSP70", "ALKBH3", "HSP90", "ALKBH5", "NPM1"),
#   Unique.Peptides = c(5, 13, 7, 9, 4)
# )

# Create the barplot
ggplot(FTO_Mass, aes(x = reorder(Protein, -Unique.Peptides), y = Unique.Peptides)) +
  geom_bar(stat = "identity", fill = "#AEC6CF", color = "black", size = 0.3) +  # Uniform bar color
  geom_text(
    aes(label = Unique.Peptides),
    vjust = -0.5,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +  # Add labels above bars
  labs(
    title = "Unique Peptides Binding to FTO",
    subtitle = "Overview of protein interactions based on unique peptides",
    x = "Protein Names",
    y = "Number of Unique Peptides"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 15, face = "bold"),  # ← Improved label readability
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed")
  ) +
  ylim(0, max(FTO_Mass$Unique.Peptides) + 3)  # Extend y-axis for labels


library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

df <- read_excel("Fractionation_plot.xlsx")

# Convert to long format for ggplot
df_long <- pivot_longer(df, cols = everything(), names_to = "Knockdown", values_to = "Value")

# Calculate mean and standard error (optional)
df_summary <- df_long %>%
  group_by(Knockdown) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE),
            SE = sd(Value, na.rm = TRUE) / sqrt(n()))

# Generate dynamic color palette based on the number of knockdowns
color_palette <- RColorBrewer::brewer.pal(length(unique(df_summary$Knockdown)), "Set3")

# Create a clean bar plot
ggplot(df_summary, aes(x = Knockdown, y = Mean_Value, fill = Knockdown)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +  # Black outline for clarity
  geom_errorbar(aes(ymin = Mean_Value - SE, ymax = Mean_Value + SE), 
                width = 0.2, color = "black") +  # Error bars
  theme_minimal(base_size = 14) +  # Cleaner theme with larger font
  labs(title = "Average Values for Each Knockdown",
       x = "Knockdown Condition",
       y = "Nuclear:Cytoplasmic Ratio") +
  theme(legend.position = "top", legend.title = element_blank()) +  # Place legend at top
  scale_fill_manual(values = color_palette) +  # Dynamically assign colors
  annotate("text", x = 1.5, y = max(df_summary$Mean_Value) * 1.1, 
           label = "n = 3", size = 5, fontface = "italic")  # Legend text for n=3