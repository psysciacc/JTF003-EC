# Data preparation notes

This sub-folder contains the code for wrangling and preparing the data.

There are several scripts, to be executed in order as documented in the main script `prepare_dataset.R`. 

1. `fetch_data_v4.R` this script load and join raw data  from all sources (including Norway data which was collected offline). This code is mainly the same as in the tracker app, and already apply some data cleaning steps. Running this require the Qualtrics API key, but is not necessary to run the rest of the analyses as the output of this script is in the file `all_withNorway_raw.csv`

2. `data_cleaning.R` More data clearning, and removal of `Test` responses

3. `fill_in_country.R` Fill in missing country values by cross-checking 

At the end of hte script the cleaned dataset is copied in the `analysis-code`  folder
