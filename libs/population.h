//============================================================================
// Name        : population.h
// Author      : UnJin Lee
// Description : Population class and helpers
//============================================================================

#ifndef POPULATION_H
#define POPULATION_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include <individual.h>
#include <func.h>

#define TIME_FACTOR 0.01

using namespace std;

class Population{
	string pop_name;
	Individual ** popln;
	int pop_size;
	Func * model;
	Func * transform;
	double * pop_xs;
	double * pop_xis;
	bool * selected;
	Func * parm_dists;
	bool recomb;
	double * parm_mod;
	int nparm;
	double ** mut_parms;
	double t;
	double time_step;
	double time_fac;
	bool col;
	Solver * sol;
	int counter;
	int time_chunk;
	bool flip;
	double * sel_coeff;
	int * parent;
	int * p_under_sel;
	bool pleiotropy;
	int n_seld;
	bool io;

	ofstream * outfile_x;
	ofstream * outfile_x_trunc;
	ofstream * outfile_xis;
	ofstream * outfile_xis_trunc;
	ofstream * outfile_parm;
	ofstream * outfile_parent;
	ofstream * outfile_fitness;
	
	stringstream * s;
	
	public:
		Population(string name, int size, Func * mod, Func * trans, Solver * solver, double * parm_model, int n_parm, double x_0, Func * param_dist, int * parm_seld, bool recombine, double ** mutation_parm, double timeStep, bool colored, bool flippy, bool pleio, bool inout=true){
			pop_name = name;
			pop_size = size;
			
			model = mod;
			transform=trans;

			stringstream p;
			s=&p;

			parm_mod = parm_model;
			nparm = n_parm;
			
			pop_xs = new double[pop_size];
			pop_xis = new double[pop_size];
			
			selected = new bool[pop_size];
			
			sel_coeff = new double[pop_size];
			parent = new int[pop_size];
			
			Individual ** population;
			population = new Individual*[pop_size];
			popln = population;
			
			flip = flippy;
			pleiotropy = pleio;
			
			for(int i=0; i < pop_size; i++){
				double * my_model_parm = new double[n_parm];
				std::copy(parm_model, parm_model + n_parm, my_model_parm);
				popln[i] = new Individual(mod, solver, my_model_parm, n_parm, x_0, 0., flippy);
				parent[i] = 0;
			}
			
			parm_dists = param_dist;
			mut_parms = mutation_parm;
			p_under_sel = parm_seld;
			
			recomb = recombine;
			col = colored;
			
			sol = solver;
			
			t = 0.;
			time_step = timeStep;
			time_fac = TIME_FACTOR;
			if(TIME_FACTOR > timeStep) time_fac = timeStep;
			time_chunk = round(time_fac / time_step);
			n_seld = 0;
			
			update_xs();
			update_xis();
			

			//file output stuffs
			counter=0;

			outfile_x = new ofstream;
			outfile_x_trunc = new ofstream;
			outfile_xis = new ofstream;
			outfile_xis_trunc = new ofstream;
			outfile_parent = new ofstream;
			outfile_parm = new ofstream[nparm];
			outfile_fitness = new ofstream;
			
			io = inout;
			if(inout==true) write_init();			
		};
	
	//functions for extracting a single parameter for the whole population (array or summary output)
		double* get_xs();
		double* get_xis();
		double get_time(){return t;};
		double get_parm_ind(int indiv, int p_idx){return popln[indiv]->get_parm(p_idx);};
		double * get_parms_ind(int indiv){return popln[indiv]->get_parms();};
		double get_sel_coeff(int indiv){return sel_coeff[indiv];};
		
		int get_n_seld(){return n_seld;};
	//functions for setting a single parameter for the whole population (input is array or single value)
		void set_xs(double * xs);
		void set_xis(double * xis);
		void update_xs();
		void update_xis();
	
	//functions for writing outputs
		void write_init();
		void write_out();
		void write_out_repop();
		void write_close();
	
	//functions for returning summary values of the population (like size, etc)
	
	
	//function for selection
		void selection_single(Func fitness, double * parm_fit);
		void selection_sweep(Func fitness, double * parm_fit);
	//function for checking bounds
		void bound_cond(double lower, double upper);
	//function for repopulation
		void repop(bool mutation);
	//function for drift-selection based repopulation
		void repop_wf_sel(Func fitness, double * parm_fit, double mutation);
	//function for drift-selection based repopulation for genetic distances	
		int repop_wf_sel_gen_dist(Func fitness, Func * pdf, double * parm_fit, double mutation, double target, double target_margin, double right_bound);		
	//function for stepping simulation, including selection, checking bounds, and repop
		void step();		

		bool isbetween(double a, double b, double epsilon);	

//	private:
	private:

};

	//functions for extracting a single parameter for the whole population (array or summary output)
	//warning: these functions return a pointer, so watch for strange behavior
double* Population::get_xs(){
	return pop_xs;
}
double* Population::get_xis(){
	return pop_xis;
}
	//functions for setting a single parameter for the whole population (input is array or single value)
void Population::set_xs(double * xs){
	for(int i=0; i<pop_size; i++){
		popln[i]->set_x(xs[i]);
		pop_xs[i] = xs[i];
	}
}
void Population::set_xis(double * xis){
	for(int i=0; i<pop_size; i++){
		popln[i]->set_xi(xis[i]);
		pop_xis[i] = xis[i];
	}
}
void Population::update_xs(){
	for(int i=0; i<pop_size; i++){
		pop_xs[i] = popln[i]->get_x();
	}
}
void Population::update_xis(){
	for(int i=0; i<pop_size; i++){
		pop_xis[i] = popln[i]->get_xi();
	}
}
	
	//functions for returning summary values of the population (like size, etc)
	
void Population::write_init(){
	for(int i=0; i<nparm; i++){
		*s << pop_name << ".parm" << i << ".tsv";
		outfile_parm[i].open(s->str().c_str());
		s->str("");
		s->clear();
	}

	*s << pop_name << "." << model->get_type() << ".x.tsv";
	outfile_x->open(s->str().c_str());
	s->str("");
	s->clear();

	*s << pop_name << "." << model->get_type() << ".x_trunc.tsv";
	outfile_x_trunc->open(s->str().c_str());
	s->str("");
	s->clear();

	*s << pop_name << "." << model->get_type() << ".xis.tsv";
	outfile_xis->open(s->str().c_str());
	s->str("");
	s->clear();

	*s << pop_name << "." << model->get_type() << ".xis_trunc.tsv";
	outfile_xis_trunc->open(s->str().c_str());
	s->str("");
	s->clear();

	*s << pop_name << "." <<model->get_type() << ".parent.tsv";
	outfile_parent->open(s->str().c_str());
	s->str("");
	s->clear();
	
	*s << pop_name << "." << model->get_type() << ".fitness.tsv";
	outfile_fitness->open(s->str().c_str());
	s->str("");
	s->clear();
	
	*outfile_x << model->get_type() << " // s_ord // " << sol->get_order() << " // s_stoch // " << sol->get_stochastic() << " // s_col // "  << sol->get_colored_noise() << " // t_step // " << time_step << "\n";
	
	*outfile_x_trunc << model->get_type() << " // s_ord // " << sol->get_order() << " // s_stoch // " << sol->get_stochastic() << " // s_col // "  << sol->get_colored_noise() << " // t_step // " << time_step << "\n";
				
	
	*outfile_x << "\ntime";
	*outfile_x_trunc << "\ntime";
	*outfile_xis << "\ntime";
	*outfile_xis_trunc << "\ntime";
	*outfile_parent << "\ntime";
	*outfile_fitness << "\ntime";

	for(int i=0; i<nparm; i++){
		outfile_parm[i] << "\ntime";
	}

	for(int i=0; i<pop_size; i++){
		*outfile_x << "\tx_" << i;
		*outfile_x_trunc << "\tx_" << i;
		*outfile_xis << "\txi_" << i;
		*outfile_xis_trunc << "\txi_" << i;
		*outfile_parent << "\tsol_" << i;
		*outfile_fitness << "\tx_" << i;
		for(int j=0; j<nparm; j++){
			outfile_parm[j] << "\tx_" << i;
		}
	}
}
	
	//function for writing out
void Population::write_out(){
	if(counter % 10 == 0){
		*outfile_x << "\n" << t;
		*outfile_xis << "\n" << t;
		
		
		for(int i=0; i<pop_size; i++){
			*outfile_x << "\t" << popln[i]->get_x();
			*outfile_xis << "\t" << popln[i]->get_xi();
		}
	}
}

void Population::write_out_repop(){
	*outfile_parent << "\n" << t;
	*outfile_fitness << "\n" << t;
	*outfile_x_trunc << "\n" << t;
	*outfile_xis_trunc << "\n" << t;
	
	for(int i=0; i<nparm; i++){
		outfile_parm[i] << "\n" << t;		
	}
		
	for(int i=0; i<pop_size; i++){
		*outfile_parent << "\t" << parent[i];
		*outfile_x_trunc << "\t" << popln[i]->get_x();
		*outfile_xis_trunc << "\t" << popln[i]->get_xi();
		*outfile_fitness << "\t" << sel_coeff[i];
		for(int j=0; j<nparm; j++){
			outfile_parm[j] << "\t" << popln[i]->get_parm(j);
		}
	}
}

void Population::write_close(){
	*outfile_x << "\r\n";
	*outfile_xis << "\r\n";

	*outfile_x_trunc << "\r\n";
	*outfile_xis_trunc << "\r\n";
	
	*outfile_fitness << "\r\n";
	
	for(int i=0; i<nparm; i++){
		outfile_parm[i] << "\r\n";
	}

	outfile_x->close();
	outfile_xis->close();
	outfile_x_trunc->close();
	outfile_xis_trunc->close();
	
	outfile_parent->close();
	
	outfile_fitness->close();
	
	for(int i=0; i<nparm; i++){
		outfile_parm[i].close();
	}
}	
	
	//functions for selection
void Population::selection_single(Func fitness, double * parm_fit){
	int selection_ind = rand() % pop_size;
	selected[selection_ind] = popln[selection_ind]->selection(fitness, parm_fit, t);
}
void Population::selection_sweep(Func fitness, double * parm_fit){
	for(int i=0; i<pop_size; i++){
		selected[i] = popln[i]->selection(fitness, parm_fit, t);
	}
}
	//function for checking bounds
void Population::bound_cond(double lower, double upper){
	for(int i=0; i<pop_size; i++){
		selected[i] = (popln[i]->get_x() <= lower) || (popln[i]->get_x() >= upper);
	}
}
	//function for mutation, recombination, and repopulation
void Population::repop(bool mutation){
	for(int i=0; i<pop_size; i++){
		if(selected[i]){
			int repop_idx0 = 0;
			int repop_idx1 = -1;
			
			//draw until a surviving individual is found or max attempts tried
			try{
				for(int count=0; count <= pop_size && selected[repop_idx0]; count++){
					repop_idx0 = rand() % pop_size;
					
					//if recombining, need to select 2nd individual
					if(!selected[repop_idx0] && recomb){
						for(int count1=0; count1 <= pop_size && selected[repop_idx1]; count1++){
							repop_idx1 = rand() % pop_size;
						}
						
						if(selected[repop_idx1]) throw 1;
					}
				}
				
				if(selected[repop_idx0]) throw 0;
			}
			
			
			//no successful draws at all, checking for extinction
			catch(int e){
				//cant find a single solution
				if(e == 0){
					cerr << "\nMax iterations hit, no solution, checking if survival";
					repop_idx0 = -1;
					for(int count=0; count <= pop_size && repop_idx0 == -1; count++){
						if(!selected[count]) repop_idx0 = count;
					}
			
					if(repop_idx0 == -1){
						cerr << "\nExtinction!!";
						terminate();
					}
					else{
						if(!recomb)	cerr << "\nUsing lowest indexed survivor";
						else{
							for(int count=repop_idx0; count <= pop_size && repop_idx1==-1; count++){
								if(!selected[count]) repop_idx1 = count;
							}
							
							if(!repop_idx1 == -1) cerr << "\nUsing lowest two indexed survivors";
							else{
								cerr << "\nRepopulating only using single survivor";
								repop_idx1 = repop_idx0;
							}
						}
					}
				}
			}
			
			//recombine or replicate
			if(recomb) popln[i]->recombine(*popln[repop_idx0], *popln[repop_idx1]);
			else popln[i]->copy(popln[repop_idx0]);
		}
		
		//mutate, regardless of selection status
		if(mutation){
			for(int i=0; i < pop_size; i++){
				popln[i]->mutate(parm_dists, mut_parms, t, p_under_sel);
			}
		}
		
		//reset selection mask
		selected[i] = false;
	}
}

void Population::repop_wf_sel(Func fitness, double * parm_fit, double mutation){
	double sum = 0.;
	int idx = 0;
	
	for(int i=0; i<pop_size; i++){
		sel_coeff[i] = fitness.apply(popln[i]->get_x(), t, parm_fit);
		if(selected[i]) sel_coeff[i] = 0.;
		sum += sel_coeff[i];
		//cout << sel_coeff[i] << "\t" << popln[i]->get_x() << endl;
	}
	
	if(sum == 0.){
		cerr << "\nExtinction!" << endl;
		for(int i=0; i<pop_size; i++) sel_coeff[i] = 1.;
		sum = pop_size*1.;
		//terminate();
		//write_close();
		//terminate();
	}
	
	for(int i=0; i<pop_size; i++){
		double rv = 1.*rand()/RAND_MAX;
		double sel_prop = 0.;
		int idx = 0;
		for(int j=0;sel_prop < rv; j++){
			sel_prop += sel_coeff[j]/sum;
			idx = j;
		}
		
		popln[i]->copy(popln[idx]);
		parent[i] = idx;
		
		selected[i] = false;
	}
	
	
	for(int i=0; i<pop_size; i++){
		double rv = 1.*rand()/RAND_MAX;
		if(rv < mutation){
			if(pleiotropy) popln[i]->mutate(parm_dists, mut_parms,t, p_under_sel);
			else popln[i]->mutate_oneparm(parm_dists, mut_parms,t, p_under_sel);
		}
	}
	
	n_seld++;
	
	if(io) write_out_repop();
}

int Population::repop_wf_sel_gen_dist(Func fitness, Func * pdf, double * parm_fit, double mutation, double target, double target_margin, double right_bound){
	double sum = 0.;
	double avg = 0.;
	int idx = 0;
	double best_soln = 0.;
	double parms_to_send[12];
	parms_to_send[0] = double(4);
	parms_to_send[1] = double(50000);
	parms_to_send[2] = 0.0001;
	parms_to_send[3] = right_bound;
	
	int result = -1;
	
	for(int i=0; i < nparm; i++){
		parms_to_send[i+nparm] = parm_fit[i];
	}
	
	for(int i=0; i<pop_size; i++){
		for(int j=0; j<nparm; j++){
			parms_to_send[j+nparm+4] = popln[i]->get_parm(j);
		}
	
		sel_coeff[i] = fitness.apply(popln[i]->get_x(), 999., parms_to_send, pdf);
		avg += sel_coeff[i];
		cout << "\t" << sel_coeff[i];
		//sel_coeff[i] = target-sel_coeff[i];
		//sel_coeff[i] = sqrt(sel_coeff[i]*sel_coeff[i]);
//		cout << "\t" << sel_coeff[i];
		//if(sel_coeff[i] < target_margin) result = i;
		//sel_coeff[i] = 1/sel_coeff[i];
		
//		cout << "\t" << sel_coeff[i];
		
		if(sel_coeff[i] > best_soln) {
			best_soln = sel_coeff[i];
			result = i;	
		}
		
		if(selected[i]) sel_coeff[i] = 0.;
		sum += sel_coeff[i];
//		cout << sel_coeff[i] << "\t" << sum << endl;
	}
	
	cout << "\tavg:\t" << avg/pop_size << endl;
	
/*	if(sum == 0.){
		cerr << "\nExtinction!" << endl;
		terminate();
	}*/
	
	if(result >= -1){
		for(int i=0; i<pop_size; i++){
			double rv = 1.*rand()/RAND_MAX;
			double sel_prop = 0.;
			int idx = 0;
			for(int j=0;sel_prop < rv; j++){
				sel_prop += sel_coeff[j]/sum;
				idx = j;
			}
		
			popln[i]->copy(popln[idx]);
			parent[i] = idx;
		
			selected[i] = false;
		}
	
	
		for(int i=0; i<pop_size; i++){
			double rv = 1.*rand()/RAND_MAX;
			if(rv < mutation){
				if(pleiotropy) popln[i]->mutate(parm_dists, mut_parms,t, p_under_sel);
				else popln[i]->mutate_oneparm(parm_dists, mut_parms,t, p_under_sel);
			}
		}
	}
	
	n_seld++;
	
	return result;
}

	//function for stepping simulation, including selection, mutation, checking bounds, and repop
void Population::step(){
	counter++;
	t += time_step;
	double forward[2] = {1., 0.};
	double reverse[2] = {0., 0.};

	for(int i=0; i < pop_size; i++){

		//colored noise section	
		if(col){
			double rv = model->rand_normal(0, 1);

			double c_now = popln[i]->get_x();
			double c_now_trans = transform->apply(c_now, t, forward);

			double c_next = sol->step(c_now_trans, popln[i]->get_xi(), t, time_step, *model, popln[i]->get_parms(), rv, 1);
			double gamma = sol->step(c_now_trans, popln[i]->get_xi(), t, time_step, *model, popln[i]->get_parms(), rv, 0);
			double c_next_revtrans = transform->apply(c_next, t, reverse);

			popln[i]->set_x(c_next_revtrans);
			popln[i]->set_xi(gamma);
		}
		
		//non-colored noise section
		else{
			double rv = model->rand_normal(0, sqrt(time_step));
			double c_now = popln[i]->get_x();
			double c_now_trans = transform->apply(c_now,t, forward);
			double c_next = sol->step(c_now_trans, t, time_step, *model, popln[i]->get_parms(), rv);
			double c_next_revtrans = transform->apply(c_next, t, reverse);
			popln[i]->set_x(c_next_revtrans);
			popln[i]->set_xi(rv);
		}
	}
}

bool Population::isbetween(double a, double b, double epsilon){
	return((a > b-epsilon)&&(a < b+epsilon));
}


#endif //POPULATION_H


