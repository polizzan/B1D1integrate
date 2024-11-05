# An integrated formal demographic fertility framework

## Purpose of This Repository
The repository **B1D1integrate** was created to enable the replication of empirical findings reported in:

*An integrated formal demographic fertility framework*,

hereafter *our manuscript*.

## Repository Structure
The repository **B1D1integrate** contains two folders:

### 1. data
This folder is used to store the data files necessary to replicate our findings. We use female age-, cohort-, and country-specific fertility and death rates, as provided in the [Human Fertility Database](https://humanfertility.org) (HFD) and the [Human Mortality Database](https://mortality.org) (HMD). These data need to be downloaded from the HFD and HMD directly and stored in the **data** folder following our naming conventions. Cohort fertility rates can be downloaded from the [HFD](https://humanfertility.org/Data/ZippedDataFiles) as a *.zip* file called *asfr.zip*. Cohort death rates can be downloaded from the [HMD](https://www.mortality.org/Data/ZippedDataFiles) as a *.zip* file called *c_death_rates.zip*. Please note that our reported findings are based on HFD and HMD data downloaded on 04 November 2024 and that data distributed by the HFD and HMD may have been updated or revised in the meantime.

### 2. scripts
This folder contains the analysis file `00_analysis.R` necessary to replicate our empirical findings. This file loads the HFD and HMD data, calculates the density, survival, and hazard functions of fertility for the selected countries, and creates Figure 1 and Figure 2 reported in the main manuscript, as well as Figure 3 and Figure 4 reported in the Appendix. The file `00_analysis.R` saves these figures in *.pdf* format in an automatically generated folder called **output**.  

## How to Use This Repository
In order to run the `R` code provided in the repository **B1D1integrate**, please proceed in the following order:

1. Download the repository from `OSF` or `github`. If applicable, unzip the downloaded folder and place it in a location convenient for you. 
2. Store the necessary HFD and HMD data in the folder **data** following our naming conventions. 
3. Double click on the file `B1D1integrate.Rproj`. This should open `RStudio` on your machine.  
4. Within `RStudio`, click on `File`/`Open File...` and select the analysis file `00_analysis.R` located in the **scripts** folder.
5. You should now be able to run our code without adjusting any directories.

## License
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
