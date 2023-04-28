# baldur_workflow
The R process to generate the results presented in the manuscript.

Everything except the power analysis can be ran with the following:
```{r, eval = FALSE}
source("1_packages_and_function.R")
source("2_read_data.R")
source("3_lgmr.R")
source("4_significane_calling.R")
source("5_performance.R")
```

But, for the power analysis, one needs to first make `6_human_power.R` executable:
```
chmod +x 6_human_power.R
```

Then, you can run it as follows:
```
./6_human_power.R <number of columns> <number of replicates>
```
Note that the script will run the LGMR model with 20 parallel threads.
If your computer cannot handle that, please change that manually line 37.
each run will save an `RData` file that can (and will) later be read in.
After running all the subsets you can plot with:

```{r, eval = FALSE}
source("7_plotting.R")
```

Note that the plotting script is hardcoded for the following sequence of columns:
```{r, eval = FALSE}
3, 6, 9, 12, 15, 18, 21
```
