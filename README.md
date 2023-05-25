# White-Lab
Code I've written for my undergraduate research position in the White Lab at MIT.

### Code use:
Each folder contains a different subproject.
Please see .py files inside folders for code and descriptions about what the code does.
The rest of the files in each folder are inputs or outputs

### Project summary:
Currently, there is a challenge with bulk sequencing in tissue analysis. Researchers and doctors want to identify the types of immune cells found in tumor tissue samples and in what proportions the cells occur. The overall aim of the project was to develop a computational method that can predict the proportions of distinct cell types from bulk phosphoproteomic data collected from tumor tissue samples. This information could help us understand the biology of complex patient samples and the roles of various cells in tumor growth.

We use phosphoproteomic and immune cell data from tumor tissue samples. First, I created code to perform PCA on immune cells vs. samples and peptides vs. samples. I also performed correlation analysis on the two sets of data to identify relationships between specific peptides and immune cells. I found an interesting negative correlation between immune receptor substrate 1 (IRS-1) and natural killer T (NKT) cells. We realized that the dataset was too noisy to form accurate predictions of immune cells based on peptide signatures. Thus, I focused on specific peptides of I performed correlations between specific pairs of peptides and between specific peptides and immune cells to identify relationships.

See "Cell Type Deconvolution by Phosphoproteomics.pdf" for overview and results.

This is only the code from the first half of my UROP (Fall 2020 and not Spring 2021).
Unfortunately, I lost the second half of my code when my laptop broke. However, the slide deck PDF contains results from the entire UROP.
