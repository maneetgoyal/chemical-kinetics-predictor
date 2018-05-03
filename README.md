# kineticsML
Machine Learning Experiments on the NIST Chemical Kinetics Database.

In this work, we examine how the presence of certain species and the formation/breaking of set of chemical bonds during gas phase homogeneous reactions correlates with the corresponding activation energy \& reaction order, using machine learning approaches. The reactions data was acquired from the NIST Kinetics database via web-scraping and HTML parsing to convert the online database to a more accessible dataframe (Pandas) format. We then engineered a feature vector that captures the presence of certain species and bonds, and utilized it to train our models to predict activation energy and reaction order. The Random Forest based classification model was chosen because this robust ensemble method has a very intuitive basis behind its working (decision trees), and is also resilient to over-fitting. Huge class imbalance was encountered in the favor of second order reactions and oversampling of the minority classes was done by random oversampling and synthetic minority over-sampling technique (SMOTE). A high prediction score of around 95\% was achieved. In regards to regression for predicting activation energy, the results weren't satisfactory and call for the engineering of better feature vectors.

### Packages you may need to install:
##### SciPy Stack | https://www.scipy.org/install.html
##### NumPy | https://www.scipy.org/install.html
##### Pandas | https://pandas.pydata.org/pandas-docs/stable/install.html
##### BeautifulSoup | https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup
##### Juypter | http://jupyter.org/install
##### ChemSpiPy | http://chemspipy.readthedocs.io/en/latest/guide/install.html
##### PubChemPy | https://pubchempy.readthedocs.io/en/latest/guide/install.html
##### Scikit-learn | http://scikit-learn.org/stable/install.html
##### imbalanced-learn | http://contrib.scikit-learn.org/imbalanced-learn/stable/install.html
##### seaborn | https://seaborn.pydata.org/installing.html#installing

### Project Package Skeleton
![alt text](https://github.gatech.edu/mgoyal35/kineticsML/blob/master/Report/PackageSkeleton.png)

### Project Workflow
![alt text](https://github.gatech.edu/mgoyal35/kineticsML/blob/master/Report/FlowChart.png)

### Raw Data Source
J. A. Manion, R. E. Huie, R. D. Levin, D. R. Burgess Jr., V. L. Orkin, W. Tsang, W. S. McGivern, J. W. Hudgens, V. D. Knyazev, D. B. Atkinson, E. Chai, A. M. Tereza, C.-Y. Lin, T. C. Allison, W. G. Mallard, F. Westley, J. T. Herron, R. F. Hampson, and D. H. Frizzell, NIST Chemical Kinetics Database, NIST Standard Reference Database 17, Version 7.0 (Web Version), Release 1.6.8, Data version 2015.09, National Institute of Standards and Technology, Gaithersburg, Maryland, 20899-8320.  Web address:  http://kinetics.nist.gov/
