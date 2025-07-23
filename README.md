# scRNAseq-checker
Shiny R application for seurat data check

This application was made as a temporary solution for biologists. Gives oportunity to check pathways, clusterization, cell distribution, gene expression with no programming skills. Uses seurat files as an imput. You can save pictures with plotly.

Note: if your seurat file is big - it might take time. Little optimisation is proposed in comments inside the script

# How to use
1) Upload seurat file via load("YourFile.RData")
2) Download pathways. You can use human pathways from GSEA. Change names in pathways <- gmtPathways("m2.all.v2023.1.Mm.symbols.gmt.txt")
3) Run application
