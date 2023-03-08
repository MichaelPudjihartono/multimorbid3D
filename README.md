This repository contain scripts used for the PPIN and comorbidity analysis in the manuscript "Integrative Tissue-Specific Analysis of Melanoma Risk Loci Reveals Novel Biological Insights That Links to Driver Genes and Comorbidities". This study integrates data on multiple levels of tissue-specific biological information to gain novel insights into the mechanisms underlying melanoma risk.

The scripts contained in this repository is a version of the original multimorbid3D pipeline (https://github.com/Genome3d/multimorbid3D) that have been modified to fit the study design of the manuscript. The changes are listed below:
1. Remove the requirement to only query the top 5 STRING interaction partner of each query protein regardless of STRING score cut-off
2. Add the capability to handle LD expansion on chromosome 17q21.31 (--window-control) separately from the rest of the genome (--window-size) due to the known inversion at this locus (1). 


The base command used to generate level 0-4 PPIN and identify multimorbid conditions used in the manuscript is shown below:
python comorbid.py -g genes.txt --grn-dir filtered_GRN -o output_dir --ld --window-control 10000 --levels 4 --ppin string --string-score 0.9 --bootstrap --bootstraps 350 

The usage are otherwise identical to the original multimorbid3D pipeline. Please refer to (https://github.com/Genome3d/multimorbid3D) for a complete documentation of the pipeline.



References
1. de Jong S, Chepelev I, Janson E, Strengman E, van den Berg LH, Veldink JH, Ophoff RA. Common inversion polymorphism at 17q21.31 affects expression of multiple genes in tissue-specific manner. BMC Genomics. 2012 Sep 6;13:458. doi: 10.1186/1471-2164-13-458. PMID: 22950410; PMCID: PMC3582489.

