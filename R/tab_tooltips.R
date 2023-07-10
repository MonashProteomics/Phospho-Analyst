# tool tips for tab panels

tooltips_ui <- function(tab){
  title = switch (tab,
                  # common plots
                  "Volcano plot" = "For each pairwise comparison, volcano plot illustrates statistically significant (p-value) versus change magnitude (fold change)",
                  "Heatmap" = "The heatmap provides an overview of all differentially expressed proteins (rows) across all samples (columns) and clustered hierarchically based on similar expressions",
                  "Individual Plot" = "Individual intensities of a given phosphosite/protein are plotted across all replicates of a condition either as box plot, 
                                       violin plot, interaction plot or intensity plot",
                  "PCA Plot" = "A Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. 
                                In brief, the more similar 2 samples are, the closer they cluster together",
                  "Sample Correlation" = "A correlation matrix is plotted as a heatmap to visualize the Pearson correlation coefficient between the various samples",
                  "Sample CVs" = "A histogram plot showing the distribution of protein level coefficient of variation (CV) for each condition. 
                                  Each plot also contains a vertical line, which indicates the median CV percentage for that condition",
                  "Phosphosite Numbers" = "The number of identified and quantified phosphosites/proteins in each sample",
                  "Sample coverage" = "A summary of how many phosphosites/proteins have been quantified consistently in how many samples",
                  "Normalization" = "These two plots represent the effect of the variant stabilizing normalization (vsn) method on 
                                     the phosphosite/protein intensity distribution in each sample",
                  "Missing values - Heatmap" = "Indicates whether a value of a given protein (rows) in a given sample (columns) is missing (0; white) or not (1; black). 
                                                Only phosphosites/proteins with at least one missing value are visualized",
                  "Imputation" = "A density plot of phosphosite/protein intensity (log2) distribution for each condition after and before missing value imputation being performed",
                  "Gene Ontology/ Pathway" = "Enrichment analyses performed by using the protein IDs of significantly regulated phosphosites, 
                                              three Gene Ontology terms and two pathway databases can be selected",
                  "Kinase-Substrate enrichment" = "A summary bar plot of kinase active scores based on significantly regulated phosphosites",
                  
                  # special on Phosphosite page
                  "Abundance rank" = "A rank/abundance plot illustrates a dynamic distribution of phosphosites by plotting the mean abundance values against 
                                      corresponding rank",
                  "Abundance comparison" = "A scatter plot to represent the comparison of mean phosphosite abundance between each pairwise across all phosphosites 
                                            with a 1:1 dashed line as reference",
                
                  # special on ProteinGroup page
                  "Protein Numbers" = "The number of identified and quantified phosphosites/proteins in each sample",
                  "Gene Ontology" = "Perform enrichment analyses of significantly regulated proteins by a selection of three Gene Ontology terms",
                  "Pathway enrichment" = "Perform enrichment analyses of significantly regulated proteins by a selection of two pathway enrichment databases",
                  
                  # special on Comparison page
                  "Numbers" = "The number of identified and quantified phosphosites/proteins in each sample",
                  "Interaction plot" = "Combined data to show the interaction between different replicates. Enabling distinction of abundance at the protein and phosphosite level",
                  "Bubble plot" = "Displays the phosphosite intensity values of a selected protein group. Phosphosites are coloured by experiment conditions, 
                                    and the size of points represents the mean intensity values of a same phosphosite",
                  # special on Phosphosite(corrected) page
                  "Normalization (normal vs corrected)" = "These two plots represent the normal and corrected phosphosite intensity distribution in each sample",
                  "Imputation (normal vs corrected)" = "A density plot of phosphosite intensity (log2) distribution for each condition after and before protein abundances correction being performed"
                  
  )
  span(tab, title = title,`data-toggle`="tooltip")
}