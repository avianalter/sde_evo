parm_name_1 <- "lr_lp_"
parm_name_2 <- ".parm"
parm_name_3 <- ".tsv"

mod_name <- "stoch_self_reg"

xi_name <- paste(".", mod_name, ".x_trunc.tsv", sep="")
x_name <- paste(".", mod_name, ".xis_trunc.tsv", sep="")

idx_run <- 0:99
idx_parm <- 0:3

type_list <- c("xis", "x")

#parm_0_0 <- read.csv("lr_lp_0.parm0.tsv", sep="\t", header=T)

type_name <- "parm"
for(j in idx_parm){
	for(k in idx_run){
		r_cmd <- paste(type_name, "_", j, "_", k, " <- read.csv('lr_lp_", k, ".parm", j, ".tsv', sep='\t', header=T)", sep="")

		print(r_cmd)
		eval(parse(text=r_cmd))
	}
}

#x_0 <- read.csv("lr_lp_0.stoch_self_reg.x_trunc.tsv", sep="")
for(type_name in type_list){
	for(k in idx_run){
		cmd_append <- "_trunc.tsv', sep='\t', header=T)"
		if(type_name == "x") cmd_append <- "_trunc.tsv', sep='\t', header=T, skip=2)"
	
		r_cmd <- paste(type_name, "_", k, " <- read.csv('lr_lp_", k, ".", mod_name, ".", type_name, cmd_append, sep="")

		print(r_cmd)
		eval(parse(text=r_cmd))
	}
}



