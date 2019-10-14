//============================================================================
// Name        : solver.h
// Author      : UnJin Lee
// Description : ODE solver class and methods
//============================================================================

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <math.h>

using namespace std;

class Solver{
	int order;
	int stochastic;
	int colored_noise;
	
	public:
		Solver() {
			order = 1;
			stochastic = 0;			
			colored_noise = 0;			
			};
		Solver(int o, int s, int c = 0) {
			order = o;
			stochastic = s;
			colored_noise = c;
			};
		void set_order(int o) {order = o;};	
		void set_stochastic(int s) {stochastic = s;};
		void set_colored_noise(int cn) {colored_noise = cn;};
		
		int get_order(){return order;};
		int get_stochastic(){return stochastic;};
		int get_colored_noise(){return colored_noise;};
		
		double step(double y, double time, double timeStep, Func f, double parm[]);
		double step(double y, double time, double timeStep, Func f, double parm[], double deltaW);
		double step(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det);
		double genCorNoise(double x, double cor);
		double rand_normal(double mean, double stddev);
	private:
		double forwardEuler(double y, double time, double timeStep, Func f, double parm[]);
		double RK2(double y, double time, double timeStep, Func f, double parm[]);
		double RK3(double y, double time, double timeStep, Func f, double parm[]);
		double RK4(double y, double time, double timeStep, Func f, double parm[]);
		double sRK(double y, double time, double timeStep, Func f, double parm[], double deltaW);

		double sRK2_col(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det);			
		
		double no_integration(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det);
};

//if order is not stated, will automatically assume order is 1
double Solver::step(double y, double time, double timeStep, Func f, double parm[]){
	if(stochastic == 0){
		if(order == 2) return(RK2(y, time, timeStep, f, parm));
		else if(order == 3) return(RK3(y, time, timeStep, f, parm));
		else if(order == 4) return(RK4(y, time, timeStep, f, parm));
		else return(forwardEuler(y, time, timeStep, f, parm));
	}
	else{
		cout << "\nIncorrect solver";
		return(0.);
	}
}

//step for stochastic solvers w/ no colored noise
double Solver::step(double y, double time, double timeStep, Func f, double parm[], double deltaW){
	if(stochastic == 1){
		if(order == 1)  return(sRK(y, time, timeStep, f, parm, deltaW));
		else{
			cout << "\nIncorrect solver";
			return(0.);
		}
	}
	else{
		cout << "\nIncorrect solver";
		return(step(y, time, timeStep, f, parm));
	}
}

//step for stochastic solvers w/ colored noise
double Solver::step(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det){
	if(colored_noise == 1 && stochastic == 1){
		if(order == 2) return(sRK2_col(y, gamma, time, timeStep, f, parm, deltaW, det));
		if(order == 0) return(no_integration(y, gamma, time, timeStep, f, parm, deltaW, det));
		else{
			cout << "\nIncorrect solver";
			return(0.);
		}
	}
	else{
		cout << "\nIncorrect solver";
		return(0.);
	}
}

double Solver::forwardEuler(double y, double time, double timeStep, Func f, double parm[]){
	double k_1 = timeStep * f.apply(y, time, parm);
	
	double y_next = y + k_1;
	return(y_next);
}

double Solver::RK2(double y, double time, double timeStep, Func f, double parm[]){
	double k_1 = timeStep * f.apply(y, time, parm);
	
	double k_2 = timeStep * f.apply(y + k_1/2, time + timeStep/2, parm);
	
	double y_next = y + k_2;
	return(y_next);
}

double Solver::RK3(double y, double time, double timeStep, Func f, double parm[]){
	double k_1 = timeStep * f.apply(y, time, parm);
	double k_2 = timeStep * f.apply(y + k_1, time + timeStep, parm);
	double k_3 = timeStep * f.apply(y + (k_1 + k_2)/4, time + timeStep/2, parm);
	
	double y_next = y + (k_1 + k_2 + 4*k_3)/6;
	return(y_next);
}

double Solver::RK4(double y, double time, double timeStep, Func f, double parm[]){
	double k_1 = timeStep * f.apply(y, time, parm);
	double k_2 = timeStep * f.apply(y + k_1/2, time + timeStep/2, parm);
	double k_3 = timeStep * f.apply(y + k_2/2, time + timeStep/2, parm);
	double k_4 = timeStep * f.apply(y + k_3, time + timeStep, parm);
	
	double y_next = y + k_1/6 + k_2/3 + k_3/3 + k_4/6;
	return(y_next);
}

//here, deltaW should be deltaW ~ N(0, sqrt(timeStep))

double Solver::sRK(double y, double time, double timeStep, Func f, double parm[], double deltaW){
	double a = f.apply(y, time, parm, 1);
	double b = f.apply(y, time, parm, 0);
	double gammaHat = y + a*timeStep + b*sqrt(timeStep);
	double b_gammaHat = f.apply(gammaHat, time, parm, 0);
	
	double y_next = y + a*timeStep + b*deltaW + 0.5*(b_gammaHat - b)*(deltaW*deltaW - timeStep)/sqrt(timeStep);
	return(y_next);
}

//note, this is a special case of the sRK2 algorithm for colored Gaussian noise (see Honeycutt 1992)
//parm[0] is defined as tau, where <e(t), e(t')> = D'*lambda*exp(-lambda*|t-t'|)
//to assure compatibility between Kim and Honeycutt, lambda is 1/tau 
//parm[1] is D from Kim's model - to assure compatibility between Kim's model and Honeycutt's model, D' is (Kim's D)/lambda
//deltaW ~ N(0, 1), det == 1 gives the determinate portion, det == 0 gives our correlated gammas
double Solver::sRK2_col(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det){
	double lambda = 1/parm[0];
	double D_prime = parm[1]/lambda;
	double H1 = -lambda*gamma;
	double H2 = -lambda*(gamma + timeStep*H1 + sqrt(2*D_prime*lambda*lambda*timeStep)*deltaW);

	double F1 = f.apply(y, time, parm, 1) + gamma;
	double F2 = f.apply(y + timeStep*F1, time, parm, 1) + -H2/lambda;
	
	double y_next = y + 0.5*timeStep*(F1 + F2);
	double gamma_next = gamma + 0.5*timeStep*(H1 + H2) + sqrt(2*D_prime*lambda*lambda*timeStep)*deltaW;
	
	if(det == 1) return(y_next);
	if(det == 0) return(gamma_next);
}

double Solver::no_integration(double y, double gamma, double time, double timeStep, Func f, double parm[], double deltaW, int det){
	if(det == 1) f.apply(y, time, parm);
	if(det == 0) return(0.);
}

#endif //SOLVER_H

