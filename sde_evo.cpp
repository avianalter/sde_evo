//============================================================================
// Name        : selection.cpp
// Author      : UnJin Lee
// Description : Implementataion of selection over models
// Version     : 0.0.1
//============================================================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include <func.h>
#include <solver.h>
#include <population.h>
#include <individual.h>

#define NSAMPLES 50000
#define NEQNS 3
#define NPARM 7
//parms_flattened(i,j), where i on [0, NEQNS] and j on [0, NPARM],
//follows parms_flattened[i*NPARM + j]
//similarly, x(i,j), xis(i,j), and x_init(i,j), where i on [0, NEQNS] and j on [0, NSAMPLES],
//follows x[i*NSAMPLES + j]
#define BC_LEFT 0
#define TIME_FACTOR 0.01
#define VERBOSE true

using namespace std;

bool coin_flip(double p){
	double rv = 1.*rand()/RAND_MAX;
	return p > rv;
}

double * gen_ecdf(Func * cdf, double * fit_parm, double x_left, double x_right, double time, int n_bins){
	double * ecdf = new double[n_bins+3];
	double step_size = (x_right - x_left)/n_bins;
	
	ecdf[0] = 1.*n_bins;
	ecdf[1] = step_size;
	ecdf[2] = x_left;
	
	for(int i=0; i<n_bins; i++){
		ecdf[i+3] = (x_left+i*step_size) * (*cdf).apply(x_left + i*step_size, time, fit_parm);
	}
	
	//normalize
	for(int i=0; i<n_bins;i++){
		ecdf[i+3] = ecdf[i+3]/ecdf[n_bins+3-1];
	}
	
	return ecdf;
}

int runsim_plasticity(string pop_name, string func, string transf, Func * fitness, double * fit_parms, int solver_order, int solver_stoch, int solver_colored, double * parms_model, int * parms_seld, int n_par, int size, Func * parameter_dist, double ** mut_parm, double x_init, double bound_lower, double bound_upper, double p_sel, double p_mut, bool sel_flip, bool select_ind=false, bool select_all=true, bool recomb = true, bool pleiotrop = true,int ecdf_size=10000, double t_select = 0.5, double finaltime = 10., double timeStep = 0.0001, double gen_time = 0.01){

	srand (time(NULL));

	Func model(func);
	Func transform(transf);
	
	Solver solv(solver_order, solver_stoch, solver_colored);

	Population pop(pop_name, size, &model, &transform, &solv, parms_model, n_par, x_init, parameter_dist, parms_seld, recomb, mut_parm, timeStep, (solver_colored==1), sel_flip, pleiotrop);
	
	int step_interval = round(gen_time/timeStep);
	int clicker = 0;
	
	while(pop.get_time() <= finaltime){
		if(VERBOSE) cout << "Time: " << pop.get_time() << "\n";
		
		pop.bound_cond(bound_lower, bound_upper);

		
/*		if(coin_flip(p_sel) && pop.get_time() > t_select){
			if(select_ind) pop.selection_single(*fitness, fit_parms); 
			if(select_all) pop.selection_sweep(*fitness, fit_parms);
		}
		
		pop.repop(coin_flip(p_mut) && pop.get_time() > t_select);  */

		if(clicker == 0 && pop.get_time() > t_select) {
			pop.repop_wf_sel(*fitness, fit_parms, p_mut);
			if(VERBOSE) cout << "Repopulating" << endl;
		}
	
		pop.step();
		
		pop.write_out();
		
		clicker++;
		if(clicker > step_interval) clicker = 0;
	}
	
	pop.write_close();
	
}

/*
int runsim_moving_target(string pop_name, string func, string transf, Func * fitness, double * fit_parms, Func * fit_shift, double * shift_parms, int solver_order, int solver_stoch, int solver_colored, double * parms_model, int * parms_seld, int n_par, int size, Func * parameter_dist, double ** mut_parm, double x_init, double bound_lower, double bound_upper, double p_sel, double p_mut, bool sel_flip, bool select_ind=false, bool select_all=true, bool recomb = true, bool pleiotrop = true,int ecdf_size=10000, double t_select = 0.5, double finaltime = 5., double timeStep = 0.0001, double gen_time = 0.01){

	srand (time(NULL));

	Func model(func);
	Func transform(transf);
	
	Solver solv(solver_order, solver_stoch, solver_colored);

	Population pop(pop_name, size, &model, &transform, &solv, parms_model, n_par, x_init, parameter_dist, parms_seld, recomb, mut_parm, timeStep, (solver_colored==1), sel_flip, pleiotrop);
	
	int step_interval = round(gen_time/timeStep);
	int clicker = 0;
	
	while(pop.get_time() <= finaltime){
		if(VERBOSE) cout << "Time: " << pop.get_time() << "\n";
		
		pop.bound_cond(bound_lower, bound_upper);

		if(clicker == 0 && pop.get_time() > t_select) {
			pop.repop_wf_sel(*fitness, fit_parms, p_mut);
			if(VERBOSE) cout << "Repopulating" << endl;
		}
	
		pop.step();
		
		pop.write_out();
		
		clicker++;
		if(clicker > step_interval) clicker = 0;
	}
	
	pop.write_close();
	
}
*/

int main(){
	string mufasa = "direct_lowv_higher";
	string simba = "stoch_self_reg";
	string transf = "trans_minv";
	Func fit("gaussian");
	Func mutation("neighbor_disc");
	
	double left_bound = 0.0001;
	double right_bound = 100.01;
	
	double stoch_parms[4] = {0.01, 100., 1., 67.};
	
	bool select_flip = false;
	
	int steps_cdf = 10000;
	//this is in x-space
	//double * fit_p = gen_ecdf(&stoch_cdf, stoch_parms, left_bound, right_bound, 999., steps_cdf);
	double fit_p[2] = {1./90., 0.0001};
	
	int n_parms = 4;
	int n_pop = 100;
	int sel_parms[4] = {0,0,0,1};
	
	double mut_tau[2] = {1./2., 0.0005};
	double mut_D[2] = {1./2., 1.};
	double mut_gamma[1] = {0.05};
	double mut_epsilon[2] = {1./2., 1.};
	
	double ** mutation_parms = new double*[n_parms];
	mutation_parms[0] = mut_tau;
	mutation_parms[1] = mut_D;
	mutation_parms[2] = mut_gamma;
	mutation_parms[3] = mut_epsilon;
	
	double prob_selection = 0.01;
	double prob_mutation = 0.1;

	stringstream s;
	
	for(double x=2.; x < 2.1; x += 0.25){
		for(int i=0; i < 100; i++){
			s << "lr_lp_" << i ;
			mufasa = s.str();

			stoch_parms[3] = 67.;
		
			runsim_plasticity(mufasa, simba, transf, &fit, fit_p, 2, 1, 1, stoch_parms, sel_parms, n_parms, n_pop, &mutation, mutation_parms, 1./67., left_bound, right_bound, prob_selection, prob_mutation, select_flip, false, true, true, true, 10000, 3.5, 4.2+3.5, 0.0001, x*stoch_parms[0]);

			s.str("");
			s.clear();
		}
	}
	
	
/*	
	double x_vals[50000];
	double gaus_p[50000];
	double x_left = 0.005;
	double x_right = 0.025;
	double x_step = (x_right - x_left)/(1.*50000);
	
	ofstream gaus_control;
	gaus_control.open("gaus_control.csv");
	
	for(int i =0;i< 50000; i++){
		x_vals[i]=x_left+1.*i*x_step;
		gaus_p[i]= fit.apply(x_vals[i],  1., fit_p);
		gaus_control << gaus_p[i] << "\t";
	}
	
	gaus_control << "\r\n";
	gaus_control.close();

	gaus_control.open("neighbor_control.csv");

	double parm_val[500];
	Func nei("neighbor_disc");
	for(int i=0; i < 500; i++){
		parm_val[i] = 60.;
		parm_val[i] = nei.apply(parm_val[i], 1., mut_epsilon, 1);
		gaus_control << parm_val[i] << "\t";
	}
	gaus_control << "\r\n";
	gaus_control.close();
*/	
	return 0;
	}













