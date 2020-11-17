# coral_deoxygenation_RNASeq

Project Summary
We experimentally exposed two common reef-building corals, Acropora tenuis and Acropora selago, to 12 h continuous deoxygenation (<2 mg O2 L-1) followed by 12 h re-oxygenation (~6 mg O2 L-1) via their natural night-day light cycle and followed gene expression changes.

Manuscript publised @ Global Change Biology: https://onlinelibrary.wiley.com/doi/10.1111/gcb.15436

RNA-Seq data @ NCBI: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA635763


## Workflow
1. The script `coral_deoxygenation_RNASeq.sh` was used to produce transcript count tables
2. Differential expression analysis was done using the script `Coral_deoxigenation_DESeq.R`
3. PCA plots were produced using the script `Coral_deoxigenation_PCA_plots.R`
