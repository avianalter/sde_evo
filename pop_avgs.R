idx_run <- 0:99
idx_parm <- 0:3
type_list <- c("xis", "x")
mod_name <- "stoch_self_reg"
generations <- 1:nrow(x_0)
rep_pops <- 1:100

times <- parm_0_0[,1]

#parm_0_0_means <- rowMeans(parm_0_0[,-1])
for(i in idx_parm){
	this_out_list <- "times"

	for(j in idx_run){
		this_idx_name <- paste(i, "_", j, sep="")
		this_idx_out <- paste("parm_", this_idx_name, "_means", sep="")
		r_cmd <- paste(this_idx_out, " <- rowMeans(parm_", this_idx_name, "[,-1])", sep="")
		print(r_cmd)
		eval(parse(text=r_cmd))
		
		this_out_list <- paste(this_out_list, ", ", this_idx_out, sep="")
	}
	
	#parm_0_collated <- do.call("rbind", list(times, parm_0_0_means, ....))
	r_cmd <- paste("parm_", i, "_collated <- do.call('rbind', list(", this_out_list, "))", sep="")
	#print(r_cmd)
	eval(parse(text=r_cmd))
}

#x_0_means <- rowMeans(x_0[,-1])
for(i in type_list){
	this_out_list <- "times"

	for(j in idx_run){
		this_idx_name <- paste(i, "_", j, sep="")
		r_cmd <- paste(this_idx_name, "_means <- rowMeans(", this_idx_name, "[,-1])", sep="")
		print(r_cmd)
		eval(parse(text=r_cmd))
		
		this_out_list <- paste(this_out_list, ",", this_idx_name, "_means", sep="")
	}
	
	#x_collated <- do.call("rbind", list(times, x_0_means, ...))
	r_cmd <- paste(i, "_collated <- do.call('rbind', list(", this_out_list, "))", sep="")
	#print(r_cmd)
	eval(parse(text=r_cmd))
}

png("tau.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(parm_0_collated[-1,], na.rm=T)
upper_bnd <- max(parm_0_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="Tau")
for(i in rep_pops+1) lines(generations, parm_0_collated[i,], col=i)
dev.off()

png("D.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(parm_1_collated[-1,], na.rm=T)
upper_bnd <- max(parm_1_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="D")
for(i in rep_pops+1) lines(generations, parm_1_collated[i,], col=i)
dev.off()

png("gamma.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(parm_2_collated[-1,], na.rm=T)
upper_bnd <- max(parm_2_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="Gamma")
for(i in rep_pops+1) lines(generations, parm_2_collated[i,], col=i)
dev.off()

png("epsilon.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(parm_3_collated[-1,], na.rm=T)
upper_bnd <- max(parm_3_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="Epsilon")
for(i in rep_pops+1) lines(generations, parm_3_collated[i,], col=i)
dev.off()

png("x.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(x_collated[-1,], na.rm=T)
upper_bnd <- max(x_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="X")
for(i in rep_pops+1) lines(generations, x_collated[i,], col=i)
dev.off()

png("xi.png", width=300/72*480, height=300/72*480, res=300)
lower_bnd <- min(xis_collated[-1,], na.rm=T)
upper_bnd <- max(xis_collated[-1,], na.rm=T)
plot(generations, generations, col="white", ylim=c(0.95*lower_bnd, 1.05*upper_bnd), xlab="Generations", ylab="Xi")
for(i in rep_pops+1) lines(generations, xis_collated[i,], col=i)
dev.off()
