This directory contains all code used to generate figures in the paper, with TSI 2002-2012 data serving as the primary demonstration of the full processing workflow.

For all datasets (including TSI 2016-2024 and FUI data), identical processing procedures were applied. 
To facilitate reproduction and further analysis, pre-processed results for the additional datasets are provided in their respective subdirectories.

Within each dataset-specific folder:
Core execution scripts are prefixed with "A" followed by a sequential number (e.g., "A1_data_preprocessing.m", "A2_analysis.m")
Supplementary files contain helper functions and utilities
To recreate the workflow, please execute files in strict numerical order following the prefix numbering (beginning with A1, then A2, etc.)
