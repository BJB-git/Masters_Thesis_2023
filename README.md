# Master's of Statistics Thesis 2023

# Abstract

One of the key challenges hindering accurate statistical modeling of COVID-19 data is the dynamic nature
of the evolving virus variants. Traditional statistical models typically assume key features of a phenomenon
(such as associations between infection and mortality rates of a virus) remain constant, which is not the case
with COVID-19. This thesis examines how local polynomial regression models, a powerful tool which accounts
for dynamic relationships, can be applied to COVID-19 data. The investigated models also account for a lag
in the relationship; that is, it takes several days before a person infected with COVID-19 is hospitalized or
deceased.

# Description of Data used for Analysis

Data used for analysis is publicly available via the Our World in Data team (https://ourworldindata.org/coronavirus). This data set contains reported COVID-19 cases and outcomes (such as deaths, hospitalizations and ICU admissions) for individual countries since the start of the pandemic. The countries analysed below were chosen to include a range of country locations and demographics. Data availability also differs by country, resulting in some countries not being included for certain types of analysis. For example, South Korea reports daily ICU admission counts but not daily hospitalization counts. 

# Folders

The folder "Thesis" contains a final PDF and LaTeX source code for the thesis document.

The folder "Data" contains the data file used in the R script. This data is downloaded from the OWID COVID-19 data project as of December 31, 2023. 

The folder "Data Analysis Code" contains the R code used for data analysis in Section 3. 
