![GLASS_logo_072219_08](https://user-images.githubusercontent.com/6731211/64618915-2ed59e80-d3af-11e9-8983-d41414379ad3.png)

## The GLASS consortium

### Overview
The Glioma Longitudinal AnalySiS (GLASS) consortium consists of clinical, bioinformaticians, and basic science researchers from leading institutions across the world striving to better understand glioma tumor evolution and to expose its therapeutic vulnerabilities. For more information, please read our position paper: **The Glioma Longitudinal Analysis Consortium. Glioma Through the Looking GLASS: Molecular Evolution of Diffuse Gliomas and the Glioma Longitudinal AnalySiS Consortium.** *Neuro Oncol. 2018* PMID: 29432615. 

The code in this respository was used to generate the figures and perform the analyses described in the [2021 bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2021.05.03.442486v1).

Previous projects and associated repositories are listed below:
Barthel, F.P., Johnson, K.C., Varn, F.S., Moskalik, A.D., Tanner, G., Kocakavuk, E., Anderson, K.J., Abiola, O., Aldape, K., Alfaro, K.D., et al. (2019). Longitudinal molecular trajectories of diffuse glioma in adults. Nature 576, 112-120. [Paper.](https://www.nature.com/articles/s41586-019-1775-1) [Code.] (https://github.com/TheJacksonLaboratory/GLASS)

### The GLASS Data Model

The GLASS data model was built to carefully take into the account the temporal nature of the dataset. The layers of the dataset is best explained using the diagram below. Here we see a single hypothetical patient (left-most column) from which we have obtained tumor samples at three different time points (second column), including a second tumor sub-sample (eg. multisector sample) from the primary tumor (third column). Moreover, a simultaneous DNA-RNA extraction was performed on 1st tumor recurrence and both were sequenced (fourth column). Lastly, some of the DNA aliquots were split off and separately sent off for whole genome and whole exome sequencing (fifth column).

The gray colored blocks in the diagram represent the abstraction layers we are using to represent this complexity. Firstly, the `case` level represents patient in the cohort. Secondly, the `sample` level represent bulk samples (tumor and control) in the cohort and finally the `aliquot` level is the aggregate of the sample portion used for differentiating between multi-sector samples, the analyte used for differentiating between DNA and RNA extractions and finally the analysis to which the aliquot was subjected, such as whole genome, whole exome or RNA sequencing.

![data-model](https://user-images.githubusercontent.com/9220167/48782460-3dca5380-ecac-11e8-8ac5-c3c2d71bb94a.png)

Internally, we are using the PostgreSQL database management system to manage the complex relationships between data entities. More information on the database scheme can be found [here](https://www.synapse.org/#!Synapse:syn17038081/wiki/585706).

#### GLASS Barcode

The GLASS barcodes are inspired by the TCGA barcodes and are constructed in a recognizable fashion:

![barcode](https://user-images.githubusercontent.com/9220167/48782475-491d7f00-ecac-11e8-9bba-c4ba02e8aa07.png)

Key parts of the data model are represented in the barcode.

### Data Release version spring 2022.
The data available here marks the third release of the GLASS project dataset, which is managed internally using the PostgreSQL database management system. The data available here represents the `spring-2022` snapshot of our dataset as described in the upcoming GLASS manuscript that will soon be published. These data are under active curation so future versions will include additional data as well as correct potential errors.

#### Data Download
The GLASS data can be downloaded from the `Tables` page [here](https://www.synapse.org/#!Synapse:syn17038081/tables/) and the `Files` page [here](https://www.synapse.org/#!Synapse:syn26465623). It is also possible to query the data directly using the the API by using queries. You can read more about that [here](https://docs.synapse.org/articles/tables.html).
