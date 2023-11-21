# STAT846_FinalProject

# Introduction
The K-Nearest Neighbors (KNN) algorithm is a non-parametric classification technique that makes no assumptions about the shape of decision boundaries. It is easy to implement and often performs well despite its simplicity.  

Classification algorithms operate on data where each case belongs to a class (or category) along with having several other known categorical or quantitative attributes. We often know which class some cases belong to (these cases are called our training data) but do not know the classes of other cases (called our test data). KNN predicts classes in the test data by examining known variables and seeing which cases in the test set are the “closest” (using a defined distance metric) to cases with known classifications in the training data. The number of other cases considered in this comparison is called the “k value”, and the closest cases in the training set are called the “nearest neighbors”. Smaller k values generally produce more flexible models with easy-to-compute results.  

Many R packages offer functions that implement the KNN algorithm, presenting users with a wide range of choices to consider. This report examines seven available implementations and provides recommendations on the best overall function and specific situations where other functions should be considered. HIV treatment data is used to test the performance of each package so that direct comparisons can be made. 

## Description of Example Dataset 

The Data used is from “The HIV Drug Resistance Database” built by a Stanford team.
The sample size is 1246.
The outcome is Resistance information of 5 drugs in the NRTI class. 

In this data set, there are 228 mutations of each isolate/virus and five predictors that represent five different drugs. 
These drugs are:
  * Lamivudine (3TC)
  * Abacavir (ABC)
  * Zidovudine (AZT)
  * Stavudine (D4T)
  * Didanosine (DDI)
 
These response variables have been changed to binary variables using specified cut-off values suggested by biologists.

We used only the drug Abacavir (ABC) for our purpose here.

These data sets are accessible by package "MTPS". First, installing this package in R, as proposed in the R file, and then by using "data(HIV)" command this data will be uploaded in the R enviroment.


## Overview of Evaluation Criteria 

Each KNN function is evaluated based on the following criteria: 
* Prediction accuracy 
* Algorithm robustness 
* Runtime (speed)/computational efficiency  
* Ease of use 
* Other features of the functions (output options, built-in data processing, missing data handling, etc.) 

## Folders

There are 3 folders here containing our projects.

Under folder Final_R, the R code is stored.

Under plots folder, the plots for the reports are saved.

Under Saved_RData_Files, the outputs of each packages (Confusion Matrices) are saved.
