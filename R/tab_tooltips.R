# tool tips for tab panels

tooltips_ui <- function(tab){
  title = switch (tab,
                  # common plots
                  "Volcano plot" = 'A volcano plot maps the “log2 fold changes” on the x-axis against the –log10 “p-values” on the y-axis.
                                    Potentially interesting candidate proteins are located in the left and right upper quadrant',
                  "Heatmap" = "The heatmap provides an overview of all differentially expressed phosphosites or proteins (rows) across all samples (columns) and clustered hierarchically based on similar expressions",
                  "Individual Plot" = "Individual intensities of a given phosphosite or protein are plotted across all replicates of a condition either as box plot, 
                                       violin plot, interaction plot or intensity plot",
                  "PCA Plot" = "A Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. 
                                In brief, the more similar 2 samples are, the closer they cluster together",
                  "Sample Correlation" = "A correlation matrix is plotted as a heatmap to visualize the Pearson correlation coefficient between the various samples",
                  "Sample CVs" = "A histogram plot showing the distribution of protein level coefficient of variation (CV) for each condition. 
                                  Each plot also contains a vertical line, which indicates the median CV percentage for that condition",
                  "Phosphosite Numbers" = "The number of identified and quantified phosphosites in each sample",
                  "Sample coverage" = "A summary of how many phosphosites or proteins have been quantified consistently in how many samples",
                  "Normalization" = "These two plots represent the effect of the variant stabilizing normalization (vsn) method on 
                                     the phosphosite or protein intensity distribution in each sample",
                  "Missing values - Heatmap" = "Indicates whether a value of a given protein (rows) in a given sample (columns) is missing (0; white) or not (1; black). 
                                                Only phosphosites or proteins with at least one missing value are visualized",
                  "Imputation" = "A density plot of phosphosite or protein intensity (log2) distribution for each condition after and before missing value imputation being performed",
                  "Gene Ontology/ Pathway" = "The protein IDs of significantly regulated phosphosites are used for the enrichment analysis. 
                                              Three Gene Ontology terms in addition to KEGG and Reactome  can be selected",
                  "Kinase-Substrate enrichment" = "A tool to estimate changes in kinase activities using bar plots",
                  
                  # special on Phosphosite page
                  "Abundance rank" = "A plot illustrating the dynamic distribution of all phosphosites",
                  "Abundance comparison" = "A scatter plot comparing mean phosphosite abundances across experimental conditions,  with a 1:1 dashed line as reference",
                
                  # special on ProteinGroup page
                  "Protein Numbers" = "The number of identified and quantified proteins in each sample",
                  "Gene Ontology" = "A tool to perform enrichment analyses of significantly regulated proteins using Gene Ontology terms",
                  "Pathway enrichment" = "A tool to erform enrichment analyses of significantly regulated proteins using the KEGG or Reactome database",
                  
                  # special on Comparison page
                  "Numbers" = "The number of identified and quantified phosphosites or proteins in each sample",
                  "Interaction plot" = "Side-by-side comparison of protein and associated phosphosite intensities plotted across all individual replicates of all conditions",
                  "Bubble plot" = "Displays phosphosite intensity values of selected proteins. Phosphosites are coloured by experimental conditions, 
                                  and the size of each bubble correlates with the phosphosite’s mean intensity",
                  # special on Phosphosite(corrected) page
                  "Normalization (normal vs corrected)" = "These two plots represent the normal and corrected phosphosite intensity distribution in each sample",
                  "Imputation (normal vs corrected)" = "Density plots of phosphosite intensity (log2) distribution for each condition before and after correction for protein abundances",
                  
                  # Absence/Presence page
                  "Venn Plot" = "Qualitative overlap of identified phosphosites between two or three selected conditions shown as overlapping circles",
                  "UpSet Plot" = "Qualitative overlap of identified phosphosites between any number of condition(s)"
                  
  )
  span(tab, title = title,`data-toggle`="tooltip")
}