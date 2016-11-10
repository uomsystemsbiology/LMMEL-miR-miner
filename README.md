# LMMEL-miR-miner
MATLAB scripts that use miR target prediction databases (TargetScan 7.1, DIANA-microT CDS) to enrich a statistical analysis 
(Pearson's correlation, mutual information) of miR and mRNA data from the Ludwig Melbourne melanoma (LM-MEL) cell line panel data


This script accompanies the manuscript:
#### MC Andrews/J Cursons, DG Hurley, M Anaka, JS Cebon, A Behren, EJ Crampin (2016). Systems analysis identifies miR-29b 
#### regulation of invasiveness in melanoma. *BMC Molecular Cancer*, (Accepted Nov 2016).
* doi: to-be-assigned


## Contacts
For further information, please contact:

#### Dr. Miles Andrews
* 
* *ex*: Cancer Immunobiology Laboratory, Olivia Newton-John Cancer Research Institute, Australia  
* miles.andrews (at) onjcri.org.au

#### Dr. Joe Cursons
* Bioinformatics Division, Walter and Eliza Hall Institute of Medical Research, Australia
* *ex*: Systems Biology Laboratory, University of Melbourne, Australia  
* cursons.j (at) wehi.edu.au

#### Dr. Daniel Hurley 
* 
* *ex*: Systems Biology Laboratory, University of Melbourne, Australia  
* daniel.hurley (at) unimelb.edu.au

#### Dr. Andreas Behren
* Cancer Immunobiology Laboratory, Olivia Newton-John Cancer Research Institute, Australia  
* andreas.behren (at) onjcri.org.au

#### Prof. Jonathan Cebon
* Cancer Immunobiology Laboratory, Olivia Newton-John Cancer Research Institute, Australia  
* jonathan.cebon (at) onjcri.org.au

#### Prof. Edmund Crampin
* Systems Biology Laboratory, University of Melbourne, Australia  
* edmund.crampin (at) unimelb.edu.au


### Virtual Reference Environment

For users unfamiliar with python, a Virtual Reference Environment will be available for this  
project, containing all scripts, data and documentation in an easily-deployed format.
* [GitHub page for VRE]( http://github.com/uomsystemsbiology/LMMEL-miR-miner_reference_environment )

For further information on Virtual Reference Environments, please refer to the
[Online Documentation]( http://uomsystemsbiology.github.io/research/reference-environments/ )


## Accompanying Data

These scripts use a number of data sources, including:

#### LM-MEL Cell Line Panel: miR and mRNA abundance data
* [Behren *et al.*, Pig. Cell Mel. Res., (2013)]( http://dx.doi.org/10.1111/pcmr.12097 )
* [mRNA abundance data (ArrayExpress)]( https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1496/ )
* [miR abundance data (GEO)]( http://to.be.fixed )

#### TCGA Skin and Cutaneous Melanoma (SKCM) data: miR and mRNA abundance data from clinical samples
* [Project Link]( http://cancergenome.nih.gov/ )
* [Data Link]( https://tcga-data.nci.nih.gov/tcga/tcgaCancerDetails.jsp?diseaseType=SKCM&diseaseName=Skin%20Cutaneous%20Melanoma )
* [Related Manuscript (Analytical)]( http://dx.doi.org/10.1111/pcmr.12374 )

#### TargetScan v7.1
* [Project Website]( http://http://www.targetscan.org/vert_71/ )
* [Data Link]( http://www.targetscan.org/cgi-bin/targetscan/data_download.vert71.cgi )
* [Agarwal *et al.*, eLife, (2015).]( http://dx.doi.org/10.7554/eLife.05005 )
* [Chiang *et al.*, Genes Dev., (2010).]( http://dx.doi.org/10.1101/gad.1884710 )
* [Friedman *et al.*, Genome Res., (2009).]( http://dx.doi.org/10.1101/gr.082701 )
* [Fromm *et al.*, Annu. Rev. Genet., (2015).]( http://dx.doi.org/10.1146/annurev-genet-120213-092023 )
* [Nam *et al.*, Mol. Cell, (2014).]( http://dx.doi.org/10.1016/j.molcel.2014.02.013 )
* [Garcia *et al.*, Nat. Struct. Mol. Biol., (2011).]( http://dx.doi.org/10.1038/nsmb.2115 )
* [Shin *et al.*, Mol. Cell, (2010).]( http://dx.doi.org/10.1016/j.molcel.2010.06.005 )
* [Grimson *et al.*, Mol. Cell, (2007).]( http://dx.doi.org/10.1016/j.molcel.2007.06.017 )
* [Lewis *et al.*, Cell, (2005).]( http://dx.doi.org/10.1016/j.cell.2004.12.035 )

#### DIANA-microT CDS (v5.0)
* [Project Website]( http://diana.imis.athena-innovation.gr/DianaTools/index.php )
* [Data Link]( http://to-be-fixed ) **NB**: a free account with DIANA Tools is required for download
* [DIANA Tools account creation]( http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=user/create )
* [Paraskevopoulou *et al.*, Nucleic Acids Res. (2013).]( http://dx.doi.org/10.1093/nar/gkt393 )
* [Alexiou *et al.*, Bioinformatics, (2009)]( http://dx.doi.org/10.1093/bioinformatics/btp565 )
* [Maragkakis *et al.*, BMC Bioinformatics, (2009).]( http://dx.doi.org/10.1186/1471-2105-10-295 )
* [Maragkakis *et al.*, Nucleic Acids Res., (2009).]( http://dx.doi.org/10.1093/nar/gkp292 )

#### miRTarBase
* [Project Website]( http://http://mirtarbase.mbc.nctu.edu.tw/ )
* [Data Link]( http://mirtarbase.mbc.nctu.edu.tw/php/download.php )
* [Chou *et al.*, Nucleic Acids Res., (2016)]( http://dx.doi.org/10.1093/nar/gkv1258 )
* [Hsu *et al.*, Nucleic Acids Res., (2014)]( http://dx.doi.org/10.1093/nar/gkt1266 )
