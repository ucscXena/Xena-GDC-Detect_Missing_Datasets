# Project Motivation

The UCSC Xena [1] browser is an online tool that allows users to analyze and visualize different types of genomic data. This data comes from multiple sources, most importantly the Genomic Data Commons (GDC) [2]. Data from the GDC is imported into UCSC Xena through Xena-GDC-ETL[6].  The GDC has tens of projects, hundreds of datasets, and thousands of samples, all of which take considerable time and effort to import into Xena. The GDC periodically adds and updates its data but only sometimes indicates in its release notes which datasets it has updated. Programmers who run the Xena-GDC-ETL code need to know which datasets need updating and which still need to change. While there is an automated testing suite that checks every data value in the GDC and compares it to Xena, it takes a long time to run for all projects and all datasets. The goal of this project is to provide users of the Xena-GDC-ETL code with information on which project(s) and datatype(s) require updating. This will allow the Xena browser to be updated quickly after the GDC releases data and to monitor the GDC for any dataset changes.  

## Background: GDC

The NCI's Genomic Data Commons (GDC) provides the cancer research community with a unified repository and cancer knowledge base that enables data sharing across cancer genomic studies in support of precision medicine. The GDC offers public access to data from important national and international cancer research projects such as The Cancer Genome Atlas (TCGA)[4], Therapeutically Applicable Research to Generate Effective Treatments (TARGET) [3], and others. This data, as well as its metadata, can be retrieved and downloaded using the GDC API’s various endpoints, such as “/cases”, “/files”, “/projects”, etc. These endpoints are used by the Xena-GDC-ETL to import data, transform the data into files compatible with the Xena browser, and transfer it to the Xena hub. The Xena hub is organized into Cohorts, each representing a project contained in the GDC. Individual cohorts also have datasets that correspond to the data types with samples imported from the GDC. 

## Background: Xena API

The Xena API [7] is a Python API that can be used to query the Xena hub. Users can use the Xena API to get metadata or data on specific samples, genes, datasets, and cohorts. This project uses the Xena API to query the sample submitter IDs that are contained in specific datasets. This way, they can be compared with those returned by the GDC.

## Overview

The objective of this project is to build a testing suite to independently verify the quantity of data between the data on the Xena data hub and the GDC portal. Specifically it will test 3 things:



* Whether any data sets on the GDC are missing from the Xena data hub.
* The number of samples in each data set on the Xena data hub matches the number of samples for that dataset on the GDC
* The individual samples contained in each dataset on the Xena data hub match those on the GDC.

All results including missing datasets, data sets with wrong sample numbers, and missing data types are saved to a TSV file. 

## Methods

### Testing for datasets missing from the Xena data hub.

The test first requests the relevant metadata of each file in a given project using the GDC API. This metadata includes information such as “data_type”, “workflow_type”, “platform” etc which is then used to match each file to a dataset. Next, a list of all possible datasets is made for each project, noting that some projects may not have a dataset for certain data types. Then individual files from the GDC API are matched to individual datasets. The parameters for matching a file to a dataset are derived from the xena_dataset script in the Xena-GDC-ETL code. Some files may not match any dataset and are saved as missing files with their data type presented as missing. Other files have a data type of  “slide images”, “Biospecimen Supplement”, and other clinical data; they are automatically removed because these data types are not imported to the Xena hub. This way, an accurate number of samples are placed into datasets, and we can remove data sets that don’t have any samples. 

Once there is a list of data sets compiled from the GDC API call, it is compared to the list of data sets returned by the Xena API. Data sets missing from the Xena API call are saved to the TSV under the column missing_datasets. 

### Comparing the number of samples in the Xena data hub to the GDC.

To check if the GDC has updated data for a specific project, the test checks for a change in the sample count in the GDC of a certain data type. This was done by comparing the number of submitter IDs in a data set compiled by a GDC API call and the number of submitter IDs in a dataset in the Xena hub compiled via a Xena API call. 

For mutation and copy number segment data, there are two submitter IDs associated with each file, one for the tumor and one for the normal. This is because there are normal mutations and copy number variation across the human population. These are subtracted out of the tumor data to ensure data quality and protect patient privacy. The Xena-GDC-ETL script maps these data files to the tumor submitter ID. To account for this the test checks the “tissue_type” metadata to ensure it only takes data from tumors.

Some projects, such as CPTAC-3 and BEATAML1.0-COHORT, do not use submitter IDs for various reasons. They instead use Case IDs and a truncated submitter ID, respectively. Thus, for these two projects, the comparison is not done with submitter IDs, but with Case IDs and a truncated submitter ID. 

### Testing for changes in version and redaction of samples in the GDC.

The names of each submitter ID are saved as part of the comparing samples test. When the comparison is conducted, extra submitter IDs are saved to a logging file and organized by project and data set. Submitter IDs that are missing in the Xena hub indicate that the GDC has added more samples to the dataset. Submitter IDs that are missing from the GDC indicate that the GDC has redacted certain samples from a dataset. 

The number of submitter IDs missing or redacted is printed on the screen. The test also saves submitter IDs that have been removed/updated to a logging file, which is saved in the same directory in which the test was run.


![alt_text](images/image1.png "image_tooltip")

![alt_text](images/image2.png "image_tooltip")


### Figure 1: Example results from running the script on the command line

Both images display what is printed on the command line when running the script. These images only show the projects TARGET-OS and REBC-THYR, but the code runs through all projects and prints with identical format. 

The command line prints out the results for each project one at a time. At the top, it displays the project ID as well as the number of samples that don’t belong into any data type that the Xena hub has. If that number is not 0, it may be due to an error in the code or a new data type that has never been imported to the Xena hub. In the middle, it prints out the results for each data set. First, the data type, and then the number of 

TARGET-OS is a project that has already been imported to the Xena hub but requires updating. The missing datasets and datasets that need updating are found at the bottom. The specific number of samples in Xena and the GDC for each dataset are shown at the top. 

REBC-THYR has not been imported to the Xena hub. However, the script still shows which data sets would need to be imported. 

## Results

The testing suite has been run on all 83 projects in the GDC with a total runtime of 366 seconds. It has already found many datasets that need updating such as 'somaticmutation_wxs', 'methylation27', 'methylation450', ‘methylation epic’, and 'somaticmutation_targeted'.  Imported projects that were found to need updating include BEATAML1.0-COHORT, CGCI-BLGSP, CGCI-HTMCP-CC, CPTAC-2, CPTAC-3, TARGET-ALL-P3, TARGET-WT, TCGA-COAD, TCGA-DLBC, TCGA-ESCA, TCGA-GBM, TCGA-HNSC, and others.

Some projects are also missing entire data sets which include these data sets: somaticmutation_wxs', 'methylation27', 'methylation450', and 'somaticmutation_targeted'. These projects were CGCI-HTMCP-CC, CPTAC-2, CPTAC-3, TARGET-ALL-P3, TARGET-CCSK, TARGET-WT, TCGA-BLCA, TCGA-DLBC, TCGA-ESCA, TCGA-LIHC, TCGA-READ, and TCGA-STAD.

Other projects are completely controlled and return no results. 

Lastly, some projects have not yet been imported to the Xena hub, so the TSV file displays the datasets that need to be imported in the missing datasets column for those projects. The specifics of which project is missing which dataset can be found in this file:

A more detailed look on which submitter IDs are missing from which data set can be found in the logging file created by the script: 

[Missing Submitter Ids](https://drive.google.com/file/d/1QBvQUb6bK0ijkiNwnSARtMMz000rBX7K/view)

## Usage instructions

To run the script on all projects in the GDC, run the file on the command line without any arguments

Example of running the script:

python3 DetectMissingDatasets.py

To run the script on an individual project:

[This file name] [project ID]

Example:

python3 DetectMissingDatasets.py TCGA-BRCA

Example results: 


![alt_text](images/image3.png "image_tooltip")

![alt_text](images/image4.png "image_tooltip")


## References

[1]Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0546-8

[2]Heath AP, Ferretti V, Agrawal S, An M, Angelakos JC, Arya R, Bajari R, Baqar B, Barnowski JHB, Burt J, Catton A, Chan BF, Chu F, Cullion K, Davidsen T, Do PM, Dompierre C, Ferguson ML, Fitzsimons MS, Ford M, Fukuma M, Gaheen S, Ganji GL, Garcia TI, George SS, Gerhard DS, Gerthoffert F, Gomez F, Han K, Hernandez KM, Issac B, Jackson R, Jensen MA, Joshi S, Kadam A, Khurana A, Kim KMJ, Kraft VE, Li S, Lichtenberg TM, Lodato J, Lolla L, Martinov P, Mazzone JA, Miller DP, Miller I, Miller JS, Miyauchi K, Murphy MW, Nullet T, Ogwara RO, Ortuño FM, Pedrosa J, Pham PL, Popov MY, Porter JJ, Powell R, Rademacher K, Reid CP, Rich S, Rogel B, Sahni H, Savage JH, Schmitt KA, Simmons TJ, Sislow J, Spring J, Stein L, Sullivan S, Tang Y, Thiagarajan M, Troyer HD, Wang C, Wang Z, West BL, Wilmer A, Wilson S, Wu K, Wysocki WP, Xiang L, Yamada JT, Yang L, Yu C, Yung CK, Zenklusen JC, Zhang J, Zhang Z, Zhao Y, Zubair A, Staudt LM, Grossman RL. The NCI Genomic Data Commons. Nat Genet. 2021 Mar;53(3):257-262. doi: 10.1038/s41588-021-00791-5. Erratum in: Nat Genet. 2021 Jun;53(6):935. PMID: 33619384.

[3]Ma X, Edmonson M, Yergeau D, Muzny DM, Hampton OA, Rusch M, Song G, Easton J, Harvey RC, Wheeler DA, Ma J, Doddapaneni H, Vadodaria B, Wu G, Nagahawatte P, Carroll WL, Chen IM, Gastier-Foster JM, Relling MV, Smith MA, Devidas M, Guidry Auvil JM, Downing JR, Loh ML, Willman CL, Gerhard DS, Mullighan CG, Hunger SP, Zhang J. Rise and fall of subclones from diagnosis to relapse in pediatric B-acute lymphoblastic leukaemia. Nat Commun. 2015 Mar 19;6:6604. doi: 10.1038/ncomms7604. PMID: 25790293; PMCID: PMC4377644.

[4]National Cancer Institute, 2023 The Cancer Genome Atlas Program (TCGA). www.cancer.gov/ccg/research/genome-sequencing/tcga Accessed 24 June 2024.

[5]Schmitz R, Wright GW, Huang DW, Johnson CA, Phelan JD, Wang JQ, Roulland S, Kasbekar M, Young RM, Shaffer AL, Hodson DJ, Xiao W, Yu X, Yang Y, Zhao H, Xu W, Liu X, Zhou B, Du W, Chan WC, Jaffe ES, Gascoyne RD, Connors JM, Campo E, Lopez-Guillermo A, Rosenwald A, Ott G, Delabie J, Rimsza LM, Tay Kuang Wei K, Zelenetz AD, Leonard JP, Bartlett NL, Tran B, Shetty J, Zhao Y, Soppet DR, Pittaluga S, Wilson WH, Staudt LM. Genetics and Pathogenesis of Diffuse Large B-Cell Lymphoma. N Engl J Med. 2018 Apr 12;378(15):1396-1407. doi: 10.1056/NEJMoa1801445. PMID: 29641966; PMCID: PMC6010183.

[6] Xena-GDC-ETL, https://github.com/ucscXena/xena-GDC-ETL Accessed 24 June 2024.

[7] Xena API, [https://github.com/ucscXena/xenaPython](https://github.com/ucscXena/xenaPython) Accessed 26 June 2024.
