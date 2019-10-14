//============================================================================
// Name        : func.h
// Author      : UnJin Lee
// Description : Numerical function class and methods
//============================================================================

#ifndef FUNC_H
#define FUNC_H
//#include <individual.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <math.h>

using namespace std;

class Func{
	string type;
	
	public:
		Func() {type = "";};
		Func(string s) {type = s;};
		void set_type(string s) {type = s;};
		string get_type() {return type;};
		double apply(double y, double time, double parm[]);
		double apply(double y, double time, double parm[], int det);
		double apply(double y, double time, double parm[], Func * funky_town);
//		double apply(double y, double time, Individual indiv);
//		double apply(double y, double time, Individual indiv, int det);
		double expCorFcn(double D, double time, double time_prime, double cor_time);
		double rand_normal(double mean, double stddev);
		double genCorNoise(double x, double sd, double cor);
	private:
		//declare all function types
		double testcase(double y, double time, double parm[]);
		double testcase1(double y, double time, double parm[]);
		double onlynoise(double y, double time, double parm[]);
		double cg_exp(double y, double time, double parm[]);
		double cg_fullmodel(double y, double time, double parm[]);
		double cg_tt(double y, double time, double parm[]);
		double cg_tt_nt(double y, double time, double parm[]);
		double cg_tt_noise(double y, double time, double parm[]);
		double sine(double y, double time, double parm[]);
		
		double nnet(double y, double time, double parm[]);	
		
		double stoch_self_reg_pdf(double y, double time, double parm[]);
		double gaussian(double y, double time, double parm[]);
		double comp_to_cdf(double y, double time, double parm[]);

		double trans_none(double y, double time, double par[]);
		double trans_minv(double y, double time, double par[]);
		
		double line(double y, double time, double parm[]);

		double selection_control(double y, double time, double parm[]);
		
		double cle(double y, double time, double parm[], int det);	
		
		double noise(double y, double time, double parm[], int det);
		double linfric(double y, double time, double parm[], int det);
		double multnoise(double y, double time, double parm[], int det);
		double multnoise_col_kim(double y, double time, double parm[], int det);
		double stoch_self_reg(double y, double time, double parm[], int det);
		double bistable_potential(double y, double time, double parm[], int det);
		double twohit_alpha(double y, double time, double parm[], int det);
		double twohit_beta(double y, double time, double parm[], int det);
		double twohit_mu(double y, double time, double parm[], int det);
		double twohit_repsca(double y, double time, double parm[], int det);
		double twohit_repind(double y, double time, double parm[], int det);
		
		double neighbor(double y, double time, double parm[], int det);
		double neighbor_unif(double y, double time, double parm[], int det);
		double neighbor_disc(double y, double time, double parm[], int det);


		
		double gen_dist(double y, double time, double parm[], Func * pdf);
};

//correlation function as defined by equation 2 of ./reference_documents/UNJIN.pdf
double Func::expCorFcn(double D, double time, double time_prime, double cor_time){
	return(D*exp(-abs(time - time_prime)/cor_time));
}

//update all private functions into this wrapper
double Func::apply(double y, double time, double parm[]){
	if(type == "testcase") return(testcase(y, time, parm));
	if(type == "testcase1") return(testcase1(y, time, parm));
	if(type == "cg_exp") return(cg_exp(y, time, parm));
	if(type == "cg_fullmodel") return(cg_fullmodel(y, time, parm));
	if(type == "cg_tt") return(cg_tt(y, time, parm));
	if(type == "cg_tt_nt") return(cg_tt_nt(y, time, parm));
	if(type == "cg_tt_noise") return(cg_tt_noise(y, time, parm));
	if(type == "onlynoise") return(onlynoise(y, time, parm));
	if(type == "sin") return(sine(y, time, parm));
	if(type == "nnet") return(nnet(y, time, parm));
	if(type == "stoch_self_reg_pdf") return(stoch_self_reg_pdf(y, time, parm));
	if(type == "gaussian") return(gaussian(y, time, parm));
	if(type == "comp_to_cdf") return(comp_to_cdf(y, time, parm));
	if(type == "trans_none") return(trans_none(y, time, parm));
	if(type == "trans_minv") return(trans_minv(y, time, parm));
	if(type == "line") return(line(y, time, parm));
	if(type == "selection_control") return(selection_control(y, time, parm));
	else return(0.);
}

//this wrapper is for functions with stochastic and deterministic parts
double Func::apply(double y, double time, double parm[], int det){
	if(type == "noise") return(noise(y, time, parm, det));
	if(type == "linfric") return(linfric(y, time, parm, det));
	if(type == "multnoise") return(multnoise(y, time, parm, det));
	if(type == "multinoise_col_kim") return(multnoise_col_kim(y, time, parm, det));
	if(type == "stoch_self_reg") return(stoch_self_reg(y, time, parm, det));
	if(type == "bistable_potential") return(bistable_potential(y, time, parm, det));
	if(type == "twohit_alpha") return(twohit_alpha(y, time, parm, det));
	if(type == "twohit_beta") return(twohit_beta(y, time, parm, det));
	if(type == "twohit_mu") return(twohit_mu(y, time, parm, det));
	if(type == "twohit_repsca") return(twohit_repsca(y, time, parm, det));
	if(type == "twohit_repind") return(twohit_repind(y, time, parm, det));
	if(type == "neighbor") return(neighbor(y, time, parm, det));
	if(type == "neighbor_unif") return(neighbor_unif(y, time, parm, det));
	if(type == "neighbor_disc") return(neighbor_disc(y, time, parm, det));
	else return(0.);
}

//this wrapper is for meta-functions that call other Func objects
double Func::apply(double y, double time, double parm[], Func * funky_town){
	if(type == "gen_dist") return(gen_dist(y, time, parm, funky_town));
	else return(0.);
}

//this wrapper is for using the Individual class
//double Func::apply(double y, double time, Individual ind){
	

//}

//code adapted from http://en.literateprograms.org/Box-Muller_transform_%28C%29
double Func::rand_normal(double mean, double stddev) {
	//choose a point x,y in the unit circle uniformly at random
	double x, y, r;
	do {
		//scale two random integers to doubles between -1 and 1
		x = 2.0*rand()/RAND_MAX - 1;
		y = 2.0*rand()/RAND_MAX - 1;
		r = x*x + y*y;
	} while (r == 0.0 || r > 1.0);
        
	//Apply Box-Muller transform on x, y
	double d = sqrt(-2.0*log(r)/r);
	double n1 = x*d;
	double n2 = y*d;
	//scale and translate to get desired mean and standard deviation
	double result = n1*stddev + mean;
	
	return result;
}

//generates a single random variable gamma(time_prime) with non-zero time correlation 
//w.r.t. gaussian noise (~N(0,1)) 
double Func::genCorNoise(double x, double sd, double cor){
	double y = rand_normal(0., sd);
	double z = cor*x + sqrt(1-cor*cor)*y;
	
	return(z);
}

//////////////////////////////////////////////////////////////////
//Deterministic models below
//////////////////////////////////////////////////////////////////

//test function with no external parameters
double Func::testcase(double y, double time, double parm[]){
	return(y + 2*time);
}

//test function just returning first parameter
double Func::onlynoise(double y, double time, double parm[]){
	return(parm[0]);
}

//test function with external parameters, parm[] = {1., 2.} is identical to testcase(...)
double Func::testcase1(double y, double time, double parm[]){
	return(parm[0]*y + parm[1]*time);
}

//minimal model for growth (exponential), parm[0] is r
double Func::cg_exp(double y, double time, double parm[]){
	return(parm[0]*y);
}

//full model for cancer growth, parm[] = {r, N, p}
double Func::cg_fullmodel(double y, double time, double parm[]){
	return(parm[0]*y - (y*y + parm[1]*y) + parm[2]*parm[1]);
}

//model with only tumor-tumor interaction, parm[0] is r
double Func::cg_tt(double y, double time, double parm[]){
	return(parm[0]*y - y*y);
}

//model with tumor-tumor and normal-tumor interactions, parm[] = {r, N}
double Func::cg_tt_nt(double y, double time, double parm[]){
	return(parm[0]*y - (y*y + parm[1]*y));
}

//model with tumor-tumor and random mutations, parm[] = {r, N, p}
double Func::cg_tt_noise(double y, double time, double parm[]){
	return(parm[0]*y - (y*y) + parm[1]*parm[2]);
}

double Func::sine(double y, double time, double parm[]){
	return(sin(time));
}

//computes pdf for a given y @ t, parm[] = {tau, D, gamma, epsilon}
double Func::stoch_self_reg_pdf(double y, double time, double parm[]){
	double expt_z_sq1 = 1/(parm[2] + 1/(parm[0]));
	double expt_z_sq2 = (2*parm[2])/(parm[2]*parm[2] - 1/(parm[0]*parm[0]));
	double expt_z_sq2_1 = exp(-(parm[2] + 1/parm[0])*time);
	double expt_z_sq3 = 1/(parm[2]-1/parm[0])*exp(-2*parm[2]*time);
	double expt_z_sq = (parm[1]/parm[2])*(expt_z_sq1- expt_z_sq2*expt_z_sq2_1 + expt_z_sq3);
	double beta = 1/(2*expt_z_sq);
	
	double p_c_ins = 1/y - exp(-parm[2]*time) - parm[3]/parm[2]*(1-exp(-parm[2]*time));
	double p_c = 1/(y*y) * exp(-beta*p_c_ins*p_c_ins);
	
	return p_c;
}

//computes pdf for a given y @ t, parm[] = {mu, sigma}
double Func::gaussian(double y, double time, double parm[]){
	double mu = parm[0];
	double sigma = parm[1];
	double norm_fac = 1/(sigma*sqrt(2*3.14));
	double ins_exp = (y - mu)/sigma;
	ins_exp = -0.5 * ins_exp * ins_exp;
	
	return norm_fac * exp(ins_exp);
}

//parm[0] is length of integration, parm[1] is x step, parm[2] is x_0, parm length is parm[0]+3
double Func::comp_to_cdf(double y, double time, double parm[]){
	int len = int(parm[0]);
	double x_step = parm[1];
	double x_0 = parm[2];
	
	int first = 0;
	
	for(int i=3; i<(len+3) && first == 0; i++){
		if(parm[i] > y) first = i-3;
	}
	
	return x_step*first + x_0;
	
}

double Func::trans_none(double y, double time, double parm[]){return y;};

double Func::trans_minv(double y, double time, double parm[]){
	if(parm[0] == 1.) return 1./y;
	else return 1./y;
}

//parm is {m, b}
double Func::line(double y, double time, double parm[]){
	return parm[0]*time + parm[1];
}

//returns gamma/epsilon, or mean of the distribution
double Func::selection_control(double y, double time, double parm[]){
	return parm[2]/parm[3];
}

//////////////////////////////////////////////////////////////////
//Neural network model below
//////////////////////////////////////////////////////////////////

double Func::nnet(double y, double time, double parm[]){
	//parm[] as follows: parm[0] is number of dimensions, 
	//parm[1] to parm[dim] is R
	//parm[dim+1] to parm[2*dim] is c
	//parm[2*dim+1] to parm[3*dim] is decay,
	//parm[3*dim+1] to parm[4*dim] is 1st row of T, 
	//parm[4*dim+1] to parm[5*dim+1] is 2nd row of T, ...
	//parm[(dim+2)*dim+1] to parm[(dim+3)*dim] is dimth row of T
	int dim = parm[0];
	double R = parm[1];
	double c = parm[2];
	double decay = parm[3];
	double vminus[dim];
	double T_y[dim];
	double rawinput = 0.;
	
	for(int i=0; i < dim; i++){
		vminus[i] = parm[i+4];
		T_y[i] = parm[i+dim+4];
	}

	for(int i=0; i < dim; i++){
		rawinput += vminus[i]*T_y[i];
	}
	
	double y_next = R*(1+tanh(rawinput+c))/2 - decay*y;
	
	return(y_next);
}

//////////////////////////////////////////////////////////////////
//Chemical Langevin model below
//////////////////////////////////////////////////////////////////

double Func::cle(double y, double time, double parm[], int det){
	return(0.);
}



//////////////////////////////////////////////////////////////////
//Stochastic models below
//////////////////////////////////////////////////////////////////

double Func::noise(double y, double time, double parm[], int det){
	if(det == 1) return(0.);
	if(det == 0) return(1.);
}

//parm[0] is tau_c as in eqn 4 of Newton, Kim, Liu (2013)
double Func::linfric(double y, double time, double parm[], int det){
	if(det == 1) return((-1/parm[0])*y);
	if(det == 0) return(1);
}

double Func::multnoise(double y, double time, double parm[], int det){
	if(det == 1) return(0.);
	if(det == 0) return(y);
}

double Func::multnoise_col_kim(double y, double time, double parm[], int det){
	if(det == 1) return(0.);
	if(det == 0) return(1.);
}

//parm[2] is gamma, parm[3] is epsilon (parm[0] and [1] are used by solver)s
double Func::stoch_self_reg(double y, double time, double parm[], int det){
	if(det == 1) return(-parm[2]*y + parm[3]);
	if(det == 0) return(1.);
}

double Func::bistable_potential(double y, double time, double parm[], int det){
	if(det == 1) return(y - y*y*y);
	if(det == 0) return(1.);
}

//twohit_* functions are used in non-colored noise models, therefore all parm[] will
//be directly forwarded to functions

//parm is {r_1, u_1, <d_1>}
double Func::twohit_alpha(double y, double time, double parm[], int det){
	if(det == 1) return(-parm[0]*(1-parm[1]) + parm[2]);
	if(det == 0) return(1.);
}

//parm is {r_1, u_1, alpha, r_2, u_2, <d_2>}
double Func::twohit_beta(double y, double time, double parm[], int det){
	if(det == 1) return(-parm[0]*parm[1]*exp(-parm[2])*exp(y) - parm[3]*(1-parm[4]) + parm[5]);
	if(det == 0) return(1.);
}

//parm is {r_2, u_2, beta, r_3, <d_3>}
double Func::twohit_mu(double y, double time, double parm[], int det){
	if(det == 1) return(-parm[0]*parm[1]*exp(-parm[2])*exp(y) - parm[3] + parm[4]);
	if(det == 0) return(1.);
}

//parm is {x_0, u_0, r_0, u_1, r_1, d_1, C} (model is log transformed)
double Func::twohit_repsca(double y, double time, double parm[], int det){
	if(det == 1) return(parm[4] - parm[5]);
	if(det == 0){
		double alpha_0 = parm[1] * parm[2];
		double alpha_1 = parm[3] * parm[4];
		return(sqrt( pow((alpha_0  * parm[0] * exp(-y)),2.) + parm[6] * pow(alpha_1,2.)));
	}
}

//parm is {u_0, r_0, x_0, u_1, r_1, d_1, C} (model is log transformed)
double Func::twohit_repind(double y, double time, double parm[], int det){
	if(det == 1) return(parm[4] - parm[5]);
	if(det == 0){
		double alpha_0 = parm[1];
		double alpha_1 = parm[3];
		return(sqrt( pow((alpha_0  * parm[0] * exp(-y)),2.) + parm[6] *  pow(alpha_1,2.)));
	} 
}

//parm is {perc}, (model randomly selects on [(1-perc)*y, (1+perc)*y] using uniform dev
double Func::neighbor(double y, double time, double parm[], int det){
	double rv = 1.*rand()/RAND_MAX;	
	rv = rv * 2 * parm[0] + (1-parm[0]);
	return rv*y;
}

//parm is {perc}, (model randomly selects on [y-perc, y+perc] using uniform dev
double Func::neighbor_unif(double y, double time, double parm[], int det){
	double rv = 1.*rand()/RAND_MAX;	
	rv = (rv-0.5)* 2. * parm[0];
	return rv + y;
}

//parm is {p, step}, (model randomly selects on {-n, -n+1, ..., -1, 1, ..., n-1, n} using a symmetric geometricdist)
double Func::neighbor_disc(double y, double time, double parm[], int det){
	double p = parm[0];
	double rv = 1.*rand()/RAND_MAX;
	double n_step = ceil(log(1.-rv)/log(1.-p));
	
	if(0.5 < 1.*rand()/RAND_MAX) n_step = -n_step;
	
	return y+ 1.*n_step*parm[1];
}



//////////////////////////////////////////////////////////////////
//Meta-models below
//////////////////////////////////////////////////////////////////

//parm[0] is number of parms to send to pdf, [1] is n_bins, [2] is x_left, [3] is x_right, 
//parm[4:(parm[0]+3)] is one parm set, parm[(parm[0]+4):(2*parm[0]+3)] is the other
double Func::gen_dist(double y, double time, double parm[], Func * pdf){
	int n_parms = int(parm[0]);
	int n_bins = int(parm[1]);
	double x_left = parm[2];
	double x_right = parm[3];
	
	double * parm_left = new double[n_parms];
	double * parm_right = new double[n_parms];
	
	for(int i=0; i < n_parms; i++){
		*(parm_left+i) = parm[i+4];
		*(parm_right+i) = parm[i+4+n_parms];
	}

	double distance = 0.;
	double step_size = (x_right - x_left)/(n_bins*1.);
	double x= x_left;
	
	double funky_left[n_bins];
	double funky_right[n_bins];
	double funky_left_norm = 0.;
	double funky_right_norm = 0.;
	
	for(int i=0; i < n_bins; i++){
		funky_left[i] = pdf[0].apply(x, time, parm_left);
		funky_right[i] = pdf[1].apply(x, time, parm_right);
				
		funky_left_norm += funky_left[i]*step_size;
		funky_right_norm += funky_right[i]*step_size;	

		x += step_size;
	}
	
	for(int i=0; i< n_bins; i++){
		funky_left[i] = funky_left[i]/funky_left_norm;
		funky_right[i] = funky_right[i]/funky_right_norm;
	}
	
	for(int i=0; i < n_bins; i++){
		distance += sqrt(funky_left[i])*sqrt(funky_right[i])*step_size;
	}
	
	return distance;	

}

#endif //FUNC_H



