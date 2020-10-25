## Instalation of Required Libraries

# BiomaRt

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# devtools

install.packages("devtools")

# SetupSequences custom package

devtools::install_github("JaimeMLegaz/setup-sequences")

# Install Keras

# Install Caret

# Install Python environment

reticulate::install_miniconda()

# Install Tensorflow

install.packages("tensorflow")
library(tensorflow)
install_tensorflow(version = "gpu") # Version can be left empty
