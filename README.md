# Master's of Statistics Thesis 2023

# Abstract

One of the key challenges hindering accurate statistical modeling of COVID-19 data is the dynamic nature of the evolving virus variants. Traditional statistical models typically assume key features of a phenomenon (such as associations between infection and mortality rates of a virus) remain constant, which is not the case with COVID-19. This thesis examines how local polynomial regression models, a powerful tool which accounts for dynamic relationships, can be applied to COVID-19 data. The investigated models also account for a lag in the relationship; that is, it takes several days before a person infected with COVID-19 is hospitalized or deceased. 

This thesis build upon the work of Liu et al. (2023) to further develop and refine these models. Several enhancements to previous work in this area are investigated, including assessing global bandwidth selection procedures, effectiveness of using different data types and the impact of different time series lengths.  An R implementation for these models is provided. 

# Description of Data used for Analysis

Data used for analysis is publicly available via the Our World in Data team (https://ourworldindata.org/coronavirus). This data set contains reported COVID-19 cases and outcomes (such as deaths, hospitalizations and ICU admissions) for individual countries since the start of the pandemic. The countries analysed below were chosen to include a range of country locations and demographics. Data availability also differs by country, resulting in some countries not being included for certain types of analysis. For example, South Korea reports daily ICU admission counts but not daily hospitalization counts. 

In January 2024, the main data source used by OWID changed their case count reporting frequency from daily to weekly and retroactively applied this change to historic data. Since the models in this thesis are built on daily data, the OWID data set as of December 31, 2023 was downloaded and used for the analysis. This data is available in this repository.

# Folders

The folder "Data" contains the data file used in the R script. This data is downloaded from the OWID COVID-19 data project as of December 31, 2023. 

The folder "Data Analysis Code" contains the R code used for data analysis in Chapter 6. 
