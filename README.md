# Feature Extraction: A Fourier and Numerical Mapping Approach ---- To extract features using 5 numeric mappings

Due to the high number of genomic sequencing projects, the number of RNA transcripts increased significantly, creating a huge volume of data. Thus, new computational methods are needed for the analysis and information extraction from these data. In particular, when parts of a genome are transcribed into RNA molecules, some specific classes of RNA are produced, such as mRNA and ncRNA with different functions. In this way, long non-coding RNAs have emerged as key regulators of many biological processes. Therefore, machine learning approaches are being used to identify this enigmatic RNA class. Considering this, we present a Fourier transform-based features extraction approach with 5 numerical mapping techniques (Voss, Integer, Real, EIIP and Z-curve), in order to classify lncRNAs from plants. We investigate four classification algorithms like Naive Bayes, Random Forest, Support Vector Machine and AdaBoost. Moreover, the proposed approach was compared with 4 competing methods available in the literature (CPC2, CNCI, PLEK, and RNAplonc). The experimental results demonstrated high efficiency for the classification of lncRNAs, providing competitive performance.


## Authors

* Robson Parmezan Bonidia, Fabrício Martins Lopes, Lucas Dias Hiera Sampaio and Danilo Sipoli Sanches.

* **Correspondence:** robservidor@gmail.com and danilosanches@utfpr.edu.br


## Publication

Bonidia R.P., Sampaio L.D.H., Lopes F.M., Sanches D.S. (2019) Feature Extraction of Long Non-coding RNAs: A Fourier and Numerical Mapping Approach. In: Nyström I., Hernández Heredia Y., Milián Núñez V. (eds) Progress in Pattern Recognition, Image Analysis, Computer Vision, and Applications. CIARP 2019. Lecture Notes in Computer Science, vol 11896. Springer, Cham.


## List of files

 - **Datasets:** Fasta and csv files used in our research;

 - **Examples:** Files of Example

 - **README:** Documentation;

 - **featureExtraction.py** Main File - Python.


## Dependencies

- Python (>=3.5)
- NumPy 
- SciPy
- Biopython
- Statistics


## Installing our tool

```sh
$ git clone https://github.com/Bonidia/FourierFeatureExtraction.git FeatureExtractionFourier
```

## Usange and Examples


```sh
Access folder: $ cd FeatureExtractionFourier
 
To run the tool (Example): $ python3.5 featureExtraction.py -i input -o output -l mRNA -r 2


Where:

-i = Input - Fasta format file, e.g., test.fasta

-o = output - CSV format file, e.g., test.csv

-l = Label - Dataset Label, e.g., lncRNA, mRNA, sncRNA

-r = representation/mappings, e.g., 1 = Binary, 2 = Z-curve, 3 = Real, 4 = Integer, 5 = EIIP.
```

## About

If you use this code in a scientific publication, we would appreciate citations to the following paper:

Bonidia R.P., Sampaio L.D.H., Lopes F.M., Sanches D.S. (2019) Feature Extraction of Long Non-coding RNAs: A Fourier and Numerical Mapping Approach. In: Nyström I., Hernández Heredia Y., Milián Núñez V. (eds) Progress in Pattern Recognition, Image Analysis, Computer Vision, and Applications. CIARP 2019. Lecture Notes in Computer Science, vol 11896. Springer, Cham.
