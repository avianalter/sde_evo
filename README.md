# sde_evo
sde_evo.cpp contains the main class for running the SDE simulations and selection. Reference libraries are found in ./libs. 

$> g++ -I libs/ sde_evo.cpp -o sde_evo

$> ./sde_evo

This will run the software, but will produce a large amount of text data files. The largest files will be \*.x.tsv and \*.xis.tsv, which may be deleted, since a subset of that data is saved in the \*.x_trunc.tsv and \*.xis_trunc.tsv. 

load_in.R is the R script for loading the data files into an R environment, while pop_avgs.R is the R script for producing population averages using the objects loaded into the R workspace via load_in.R
