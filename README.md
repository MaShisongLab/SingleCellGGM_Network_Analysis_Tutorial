# SingleCellGGM_Network_Analysis_Tutorial

## A tutorial on how to conduct single-cell gene co-expression network analysis using SingleCellGGM.

## 1. Conduct network single-cell gene co-expression network analysis using the SingleCellGGM algorithm.

Please refer to the SingleCellGGM package for details.

SingleCellGGM takes a log-normalized gene expression matrix, the number of iterations, the names of the genes, and the name of dataset as inputs. The expression matrix should have samples in rows and genes in columns. The sample numbers should be large and the low-expression genes should be filtered out first.

we use a mouse single-cell gene expression matrix obtained from the MCA project ([Han *et al*, 2018](#References)) as an example to demonstrate how to conduct single-cell GGM gene co-expression network analysis via scGGM. The matrix file "MCA_Figure2-batch-removed.txt.tar.gz" can be downloaded from [Figshare](https://figshare.com/ndownloader/files/10351110?private_link=865e694ad06d5857db4b) as provided by MCA. Unzip and place the file "Figure2-batch-removed.txt" into the MATLAB working folder. We also obtained the Ensembl gene IDs for the genes within the matrix and saved it in a file "data/MCA.ensembl.gene.ids.txt". 

```matlab
% MATLAB code
% Read in the single-cell gene expression matrix
expression_matrix = readtable('Figure2-batch-removed.txt','ReadRowNames',true);
expression_matrix = table2array(expression_matrix);

% Log normalization
expression_matrix = log2( expression_matrix ./ sum( expression_matrix ) * 10000 + 1 );
expression_matrix = expression_matrix';

% Read in gene names
gene = readcell('data/MCA.ensembl.gene.ids.txt');

% Out of 25133 genes in the matrix, 24802 have Ensembl gene IDs. Only the 
% genes with Ensembl gene IDs will be used for network construction.
idx = contains(gene,'ENSMUSG');

expression_matrix = expression_matrix(:,idx);
gene = gene(idx);

% Conduct single-cell gene co-expression analysis via scGGM
ggm = scggm( expression_matrix, 20000, gene, 'mca')

% Examine the results
ggm
ggm.SigEdges(1:5,:)

% Save all gene pairs to a file for gene co-expression construction
writetable(ggm.SigEdges(:,1:3),'mca.ggm.network.txt','Delimiter','tab','WriteVariableNames',FALSE)
```
## 2. Clustering the network into gene modules using the MCL algorithm. 

The [MCL algorithm](https://www.micans.org/mcl/) should be installed.

```shell
# Shell Code
# Network clustersing using MCL, and the result is saved to a file 'MCLresult'. 
# The MCL algorithm should be installed first.
mcl mca.ggm.network.txt -I 1.7 -scheme 7 -o MCLresult -te 20 --abc
```

## 3. Convert the MCLresult into a table with module information in R.

We will use R to do the job.

```R
# R Code

# Read in the MCLresult file and convert it to a table
# We only retain gene modules / GEPs with 10 or more genes.
# One can also set the limit at another number.
MCLresult  <- readLines("MCLresult")
gep_size_limit = 10
GGM_Modules <- as.data.frame(matrix(0, nrow = 0, ncol = 2))
for (i in 1:length(MCLresult)) {
  Module <- cbind(as.character(unlist(strsplit(MCLresult[i],'\t'))),i)
  if (nrow(Module) >= gep_size_limit) { GGM_Modules <- rbind(GGM_Modules,Module); j = i }
}
colnames(GGM_Modules) <- c("Gene_ID", "Module_GEP_ID") 
fmt = c("%01d","%02d","%03d","%04d","%05d")[min((floor(log(j, 10)) + 1), 5)]
GGM_Modules$Module_GEP_ID <- paste("M",sprintf(fmt,as.numeric(GGM_Modules$Module_GEP_ID)), sep = "")

# Save the table to an output file
write.csv(GGM_Modules,"mca.ggm.network.modules.csv",row.names = F)
```

## 4. GO and MP enrichment analysis for the identified GEPs.

To be done.

