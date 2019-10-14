//============================================================================
// Name        : individual.h
// Author      : UnJin Lee
// Description : Individuals class and helpers (parameter set)
//============================================================================

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <math.h>
#include <stdlib.h>

#include <func.h>
#include <solver.h>

using namespace std;

class Individual{
	Func * model_func;
	Solver * solv;
	double * parameters;
	double x;
	double xi;
	int n_parms;
	bool flip;
	
	public:
		Individual(Func * model, Solver * solver, double * parm, int npar, double x_0, double xi_0, bool dir){
			model_func = model;
			solv = solver;
			parameters = parm;
			n_parms = npar;
			x = x_0;
			xi = xi_0;
			flip = dir;
		};
		
		void copy(Individual * other_guy);
		
		double get_x();
		double get_xi();
		
		Func * get_func(){return model_func;};
		Solver * get_solv(){return solv;};

		int get_n_parms(){return n_parms;};
		
		bool get_flip(){return flip;};
	//functions for changing specific parameters
	//functions for extracting specific parameters
		double get_parm(int idx);
	//functions for the number of parameters 
	//functions for returning parameters
		double * get_parms();
		
	//functions for setting parms
		void set_x(double x_in);
		void set_xi(double xi_in);
		
	//functions for returning model
	//functions for returning det
	//function for selection
		bool selection(Func fitness, double * selection_parms, double time);
	//function for mutation
		void mutate(Func * dists, double ** parms_array, double time, int * p_seld);
		void mutate_oneparm(Func *dists, double ** parms_array, double time, int * p_seld);
	//function for recombination
		void recombine(Individual momz, Individual popz);	
	
	private:
	
	
};

void Individual::copy(Individual * other_guy){
	model_func->set_type((other_guy->get_func())->get_type());
	solv->set_order((other_guy->get_solv())->get_order());
	solv->set_stochastic((other_guy->get_solv())->get_stochastic());
	solv->set_colored_noise((other_guy->get_solv())->get_colored_noise());
	std::copy(other_guy->get_parms(), other_guy->get_parms()+(other_guy->get_n_parms()), parameters);
	x = other_guy->get_x();
	xi = other_guy->get_xi();
	n_parms = other_guy->get_n_parms();
	flip = other_guy->get_flip();
}

double Individual::get_x(){return x;}

double Individual::get_xi(){return xi;}

double Individual::get_parm(int idx) {return parameters[idx];}

double * Individual::get_parms(){return parameters;}

void Individual::set_x(double x_in){x = x_in;};

void Individual::set_xi(double xi_in){xi = xi_in;};

bool Individual::selection(Func fitness, double * selection_parms, double time){
	double rv = rand()/RAND_MAX;
	
	if(flip) return rv < fitness.apply(x, time, selection_parms);
	else return rv > fitness.apply(x, time, selection_parms);
}

void Individual::mutate(Func * dists, double ** parms_array, double time, int * p_seld){
	for(int i=0; i < n_parms; i++){
		if(p_seld[i] == 1) parameters[i] = dists->apply(parameters[i], time, parms_array[i], 1);
	}
}

void Individual::mutate_oneparm(Func *dists, double ** parms_array, double time, int * p_seld){
	int n_seld_p = 0;
	for(int i=0; i < n_parms; i++){
		n_seld_p += p_seld[i];
	}

	double rv = 1.*rand()/RAND_MAX;
	double weighted = 0.;
	
	for(int i=0; i < n_parms; i++){
		weighted += 1.*p_seld[i];
		if(weighted/n_seld_p*1. > rv){
			parameters[i] = dists->apply(parameters[i], time, parms_array[i], 1);
			i = n_parms;
		}
	}
}

void Individual::recombine(Individual momz, Individual popz){
	for(int i=0; i < n_parms; i++){
		parameters[i] = momz.get_parm(i);
		if(rand() < RAND_MAX/2.) parameters[i] = popz.get_parm(i);
	}
	
	x = momz.get_x();
	xi = momz.get_xi();
	
	if(rand() < RAND_MAX/2.) x = popz.get_x();
	if(rand() < RAND_MAX/2.) xi = popz.get_xi();	
}


#endif //INDIVIDUAL_H
