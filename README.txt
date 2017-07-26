shinyGISPA: Gene Integrated Set Profile Analysis with Shiny
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Authors
~~~~~~~
Bhakti Dwivedi & Jeanne Kowalski


Maintainer
~~~~~~~~~~
Bhakti Dwivedi


Organization
~~~~~~~~~~~~
Biostatistics & Bioinformatics Shared Resource
Winship Cancer Institute
Emory University
Atlanta, GA 30322


Purpose
~~~~~~~
shinyGISPA is a web-based tool intended for the researchers who are interested in defining gene sets with similar, a priori specified molecular profile. shinyGISPA is based on the GISPA method published in our previous paper (Kowalski et al., 2016; PMID: 26826710). The tool is developed using shiny, a web-application framework for R. 


Availability
~~~~~~~~~~~~
GISPA Bioconductor R package is available at https://www.bioconductor.org/packages/release/bioc/html/GISPA.html and shinyGISPA is available from the GitHub, https://github.com/BhaktiDwivedi/shinyGISPA.


Current Version
~~~~~~~~~~~~~~~
shinyGISPA version 1.0


Copyright
~~~~~~~~~
Copyright (c) 2016-1017, Bhakti Dwivedi


Bug Reports/Questions
~~~~~~~~~~~~~~~~~~~~~
Please contact Bhakti Dwivedi, bhakti.dwivedi@emory.edu for questions, comments, or feature requests.


Prerequisities
~~~~~~~~~~~~~~
1) Download and install R or RStudio (version 3.3. or later) from https://cran.r-project.organd
2) Open R and install the required R packages:

	install.packages( c("changepoint", “colourpicker”, "data.table”, “genefilter”, “ggplot2”, “graphics”,  “HH", “knitr”, ”latticeExtra", “plyr”, ”scatterplot3d", “stats”, “splitstackshape”, “openxlsx”) )  
  
3) Install the Bioconductor R packages required by shinyGISPA

	source("http://bioconductor.org/biocLite.R")
	biocLite(package.name)


Running shinyGISPA
~~~~~~~~~~~~~~~~~~
Users can run shinyGISPA locally by downloading the source code available from the GitHub: https://github.com/BhaktiDwivedi/shinyGISPA, by typing the below commands in R console:

    > library(shiny)
    > runApp("shinyGISPA")

Users can also download and run the app from GitHub directly as:

    > shiny::runGitHub('shinyGISPA', 'BhaktiDwivedi')
  
Please see shinyGISPA_manual.pdf for more details.  


Input File Requirements
~~~~~~~~~~~~~~~~~~~~~~~
1) ASCII formatted tab-delimited file only, where each row represents a gene (or related gene id) and each column a sample. 
2) First and second column must correspond to gene names and gene id’s (e.g., gene transcript, ensemble id or other) followed by the three sample classes data for each data type selected. Please see ‘File Format Requirements’ section of the manual for details.

Compatibility
~~~~~~~~~~~~~
Compatible with Firefox or Chrome browsers.


Funding
~~~~~~~
This work is funded by the Leukemia and Lymphoma Society Translational Research Program Award (to Jeanne Kowalski); Georgia Research Alliance Scientist Award (Jeanne Kowalski); a Team Science Seed Funding from the Winship Cancer Institute of Emory University (Lawrence H. Boise, Sagar Lonial, Michael R. Rossi); Biostatistics and Bioinformatics Shared Resource of Winship Cancer Institute of Emory University and NIH/NCI [Award number P30CA138292, in part]. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.


Citing the GISPA method:
~~~~~~~~~~~~~~~~~~~~~~~~
Kowalski J, Dwivedi B, Newman S, Switchenko JM, Pauly R, Gutman DA, Arora J, Gandhi K, Ainslie K, Doho G, Qin Z, Moreno CS, Rossi MR, Vertino PM, Lonial S, Bernal-Mizrachi L, Boise LH. Gene integrated set profile analysis: a context-based approach for inferring biological endpoints. Nucleic Acids Res. 2016 Apr 20;44(7):e69. doi: 10.1093/nar/gkv1503. Epub 2016 Jan 29. PubMed PMID: 26826710; PubMed Central PMCID: PMC4838358.
