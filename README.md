# Dissertation
This repository is home to the R script I used in my MPH dissertation to analyze the ICEMR PRISM Cohort data set. The data are publicly available at https://clinepidb.org. There is a very great description of the data on the site which may be helpful for background.

# Code Flow
The script follows this general track, though hopefully it is clean enough to follow without consulting this overview! The four dashes (----) in R studio allow for the subsection of code between the dashes to be collapsed and ran all at once, and thus I used these dashes to signpost the different sections.

## Introduction
First, I load all the necessary libraries and read the .txt files which store the data. Then I factored all the categorical variables, created the necessary derived variables, and joined a few of the tables to make the master dataframe which was used in the analyses. I also created a few functions in this introduction section which helped limit repetition in the code later.

## Stage One
The shortest stage of the analysis was the population demographics, which mostly used the functions I created earlier. Then can both be run all at once and will results in two new data objects, which ended up being tables in my dissertation (a remove() call at the end of the code chunk cleans up unused objects). 

## Stages Two-Four
They each have an exploratory analysis chunk and a few chunks where models are fit. Exploratory analysis chunks can be run all at once and have remove() calls to clean up after themselves. The model chunks probably shouldn't be run all at once because 1) it would take a few minutes, and 2) you will want to look at each model individually by running the modOutputs() line after fitting each model
