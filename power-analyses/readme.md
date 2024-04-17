# Notes on power analyses

We show here results for the analyses assuming Cohen's $d$ approximately equal to 0.2 for both effect of answer justification and feedback (see file `power-results.pdf`).

The power analyses script `power-script-cluster.R` is designed to be run via Rscript in command line and receive some parameters as input - see in particular at the top of the script:

```
cmd_args=commandArgs(TRUE)
N = as.numeric(cmd_args[1])
N_c = as.numeric(cmd_args[2])
N_rep = as.numeric(cmd_args[3])
sim_label = cmd_args[4]
```

- `N` is the number of participants per condition
- `N_c` is the number of countries
- `N_rep` is the number of iteration a simulation with a particular set of parameters is repeated
- `sim_label` is just a label appended to the results (useful is multiple set of simulations are runned in parallel with different settigns)
