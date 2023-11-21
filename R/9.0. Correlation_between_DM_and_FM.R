# ------------------------------------------------------------------------------
# Script to calculate the correlation between fresh mass and dry mass of three isogenic lines
# This script was prepared to answer one of the reviewer's concern
# ------------------------------------------------------------------------------
# Data collected by Félix P Leiva (e-mail: felixpleiva@gmail.com)
# Script created by Félix P Leiva (e-mail: felixpleiva@gmail.com)
# ------------------------------------------------------------------------------
# Cleaning working space
rm(list = ls())
# ------------------------------------------------------------------------------
# check directory
getwd()
# ------------------------------------------------------------------------------
# libraries
library(dplyr)
library(readxl)
# ------------------------------------------------------------------------------
# load data
dat <- read_excel("../Data/fresh mass and dry mass correlation 4 lines 20191206.xlsx")
# ------------------------------------------------------------------------------
# Remove rows with "NA" in the "sex" column
dat <- dat %>%
  filter(sex != "NA")

# Exclude the stock # 25203
dat <- dat %>%
  filter(stock != "25203")

# Calculate correlation for each "stock"
correlation_data <- dat %>%
  group_by(stock) %>%
  summarize(correlation = cor(dm, fw))

# Create a scatter plot with different shapes for males and females, color-coded by "stock," and facet by "stock"
combined_plot <- ggplot(data = dat, aes(x = dm, y = fw, color = stock)) +
  geom_point() +
  labs(x = "Dry mass (mg)", y = "Fresh mass (mg)", size = 2) +
  geom_text(data = correlation_data, aes(label = paste("Correlation =", round(correlation, 2)),
                                         x = max(dat$dm), y = min(dat$fw) + 0.2), hjust = 1, vjust = 0) +
  facet_wrap(~stock)

# Remove the legend for the "stock" variable
combined_plot <- combined_plot + guides(color = FALSE)

# Export the plot as a PNG file
ggsave("../Outputs/9.2.1. Combined_scatter_plot.png", combined_plot, width = 10, height = 4, units = "in", dpi = 600)
ggsave("../Outputs/9.2.1. Combined_scatter_plot.pdf", combined_plot, width = 10, height = 4)
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("../Outputs/9.1.2. Correlation_between_DM_and_FM_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################