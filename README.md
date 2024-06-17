# Heredity_pigment

Data files and R code underlying Tapanes and Rennison Paper "The genetic basis of divergent melanic pigmentation in benthic and limnetic threespine stickleback"
 
The aim of the study was to characterize the genetic archetecture of melanic pigment traits in benthic and limnetic threespine stickleback.

In this repository there is:
-BenLim-Pigment-QTL_Mapping_Final.Rmd all of the code for the heritability and QTL mapping analyses. 
-candidate_1.5LOD_8.txt an input file for "Funcitonal_Enrichment.R" containing all of the annotated genes falling withing a 1.5 LOD interval of the pigmentation QTL peak on chromosome 8. 
--candidate_1.5LOD_21.txt an input file for "Funcitonal_Enrichment.R" containing all of the annotated genes falling withing a 1.5 LOD interval of the pigmentation QTL peak on chromosome 21. 
-combined_maps.map input file of the linkage map used for QTL mapping in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- Covar.csv covariance input file with information on family and sex used for QTL mapping in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- finalped.csv pedigree input file formatted for the heredity analysis in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- finalped3.csv pedigree input file formatted for kinship matrix estimation in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- Functional_Enrichment.R code for functional enrichment analysis of candidate genes.
- Pheno-071822.csv Full phenotype input file used for QTL mapping analysis in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- Pheno.csv trimmed phenotype input file formatted for heretability analysis in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- pondsAll.fixedgeno_recode.sv geneotype input file for QTL mapping in "BenLim-Pigment-QTL_Mapping_Final.Rmd"
- stickleGO.txt GO term key from the BROAD stickleback genome annotation used as input for "Functional_Enrichment.R"


