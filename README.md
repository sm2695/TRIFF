# **TRIFF**
Analysis code for the TRIFF consortium/projects

Collab between Sellgren, Khodosevich, Lyndby, and Meyer Labs. 

To run the analysis, follow these steps:  

### **1. Clone the Repository**  
```sh
git clone https://github.com/sm2695/TRIFF.git
cd TRIFF
```

### **2. Processing**  
Process in the following order: 
1. Run script to curate genes from different sources and standardise them.
```sh
Rscript code/Gene_preprocessing.r
```
2. Run script to test selected genes during brain development using bulk RNA-seq data from BrainSpan.
```sh
Rscript code/Gene_prioritization_brainspan.R
```
3. Run script to test prenatally-enriched genes in human fetal single-cell RNA seq data. The script returns both mean expression and specificity values across brain regions and cell types. 
```sh
Rscript code/Gene_prioritisation_sc.r
```
Helper functions are stored in R/utils.R. 

## **Contribute**  
To add/modify to the existing repository (if not main):  
1. Fork the repository.  
2. Create a feature branch (`git checkout -b feature-name`).  
3. Commit your changes (`git commit -m "Add feature"`).  
4. Push to the branch (`git push origin feature-name`).  
5. Open a pull request.  
 

