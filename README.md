# Whole-tissue deconvolution and scRNAseq analysis identify altered endometrial cellular compositions and functionality associated with endometriosis
This repository is a companion to a study that carefully applies cell type deconvolution analysis to publicly available eutopic endometrial tissue transcriptional profiles. In addition to applying the xCell algorithm to this microarray data, we also describe and apply methods for assessing how well the standard signatures of the xCell algorithm, which were built on non-endometrial samples, might match with actual cognate cell type signatures of the target endometrial tissue.

If you use the repository, we ask that you cite the paper in any publications that follow:

Daniel G. Bunis*, Wanxin Wang*, Júlia Vallvé-Juanico, Sahar Houshdaran, Sushmita Sen, Isam Ben Soltane, Idit Kosti, Kim Chi Vo, Juan Irwin, Linda C. Giudice, Marina Sirota. "Whole-Tissue Deconvolution and scRNAseq Analysis Identify Altered Endometrial Cellular Compositions and Functionality Associated With Endometriosis", Front. Immunol., 05 January 2022. https://doi.org/10.3389/fimmu.2021.788315

## Additional Details:
GSE111976_ct_endo_10x.rds can be downloaded from GEO, GSE111976. All other \*.rds files used in the \*.R scripts can be found in the associated [figshare project](https://figshare.com/projects/Whole-tissue_deconvolution_and_scRNAseq_analysis_identify_altered_endometrial_cellular_compositions_and_functionality_associated_with_endometriosis/127208).
The R packages used for interpreting microarray probes in 'microarray_processing.Rmd' were downloaded from [this site](http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp), specifically the 'HGU133Plus2_Hs_ENTREZG' "C" and "A" source packages.
