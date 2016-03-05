#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <math.h>
#include <vector>

// Author Jos Gibbons, 2016

/* The int.main() below is divided into six "parts"
When a PART begins, its title appears alone on a comment line
Some descriptions of sub-PARTs, which would print to console for the user, are commented out (reinstate these in testing if processes are slow)
The comment lines below a PART's title explain it further, before a blank line: this includes describing possible future extensions
PART ONE imports historical financial data
PART TWO computes statistics for distributions used herein (see Documentation.pdf)
PART THREE estimates parameters, with errors
PART FOUR tests CDF fits
Before int.main() begins, several other functions must be defined for use in PART FOUR; these are CDFs
PART FIVE reports findings of data-fitting to an external file
PART SIX outputs 100,000 samples for each distribution considered
Before int main() begins, several other functions must be defined for use in PART SIX; these are inverse CDFs */

// CDFs for PART FOUR:

double NormalCDF(double value)
{
	return 0.5 * erfc(-value * sqrt(0.5));
}

double flambda;// fs at the beginning of variable names are simply field flags

double ExponentialCDF(double value)
{
	return 1 - exp(-flambda*value);
}

double fk;// f is a field flag
double ftheta;

double GammaCDF(double value)
{
	double factor = 0;
	double addend = 1.0 / fk; // fk's value will be set in PART THREE: this function will be used in PART FOUR
	int a = 1;// Each function herein requires at most one iterated int; these differ for modularity
	double precision = 0.000001;	
	while(addend > precision)
	{
		factor += addend;
		addend *= value / (fk + a);
		++a;
	}
	return pow(value, fk)*exp(-value) * factor / tgamma(fk);
}

double finverse_k;
double inverse_ftheta;

double InversegammaCDF(double value)
{
	double inverse_factor = 0;
	double inverse_addend = 1.0 / finverse_k;
	int b = 1;
	double inverse_precision = 0.000001;
	while(inverse_addend > inverse_precision)
	{
		inverse_factor += inverse_addend;
		inverse_addend *= value / (finverse_k + b);
		++b;
	}
	return pow(value, finverse_k)*exp((-1.0)*value) * inverse_factor / tgamma(finverse_k);
}

double pareto_minimum;
double pareto_alpha;

double ParetoCDF(double value)
{
	return 1.0 - exp(pareto_alpha * log(pareto_minimum / value));
}

// Inverse CDFs, obtained by Newton-Raphson, for PART SIX

// Normal and Lognormal can both be done with this code

double Normal_Inverse(double value)
{
	double x;
	x = 0;
	for(int c = 0; c < 10; ++c)
	{
		x += 2.50662827463 * exp(0.5 * x * x) * (value - NormalCDF(x));
	}
	return x;
}

// Gamma (but not inverse-Gamma, because of a parameter dependence)

double Gamma_Inverse(double value)
{
	double y;
	if(value <= 0.5)
	{
		y = exp(log(fk * value * tgamma(fk)) * 1.0 / fk);
	}
	else
	{
		y = -log(tgamma(fk) * (1.0 - value));
	}
	for(int d = 0; d < 10; ++d)
	{
		y += exp((1.0 - fk) * log(y)) * exp(y) * value * tgamma(fk) * (1.0 - GammaCDF(y));
	}
	return y;
}

// Inverse-gamma

double Inverse_gamma_Inverse(double value)
{
	double z;
	if(value > 0.5)
	{ 
		z = exp(log(finverse_k * (1.0 - value) * tgamma(finverse_k)) * 1.0 / finverse_k);
	}
	else
	{ 
		z = -log(tgamma(finverse_k) * value);
	}
	for(int e = 0; e < 10; ++e)
	{
		z += exp((1.0 - finverse_k) * log(z)) * exp(z) * value * tgamma(finverse_k) * (1.0 - InversegammaCDF(z));
	}
	return z;
}

// Analytic inverse CDFs for Exponential and Pareto are obtained later

int main()
{ // int main() opens
    // std::cout << "Collecting data...\n";
	/* PART ONE
	This "part" imports data, which for convenience is also sorted
	A future extension could ask the user to specify a custom CSV's file location, and/or explain how to create one from an XLS or XLSX
	First we need to import data; I deleted cells from Data.xlsx not containing dollar values, and saved as CSV.csv */

	std::ifstream file("CSV.csv"); // declare file stream
	double value;
	std::string valuetmp;
	std::vector<double> dataset;
	while(getline(file, valuetmp))
	{
		value = ::atof(valuetmp.c_str());
		dataset.push_back(value);
	}
	sort(dataset.begin(), dataset.end());// Sorted data will be needed to compute Pareto's x_m, and for PARTs FOUR and SIX's CDF analysis
	
	/* PART TWO
	This "part" computes statistics from the data set of the form E(f(x))
	The functions f(x) for which such means are computed are those needed for parameter estimation in PART THREE
	Our models and the E(f(x)) they need are listed in Documentation.pdf
	A future extension could consider more models
	Such additional models would require more code in PART THREE, but may also require more code in PART TWO */
	
	// std::cout << "Counting data points...\n";
	double n = (double) dataset.size();// Casting prevents quotient problems later
	
	// std::cout << "Computing the means of several functions...\n";
	double sigma_xi = 0;
	double sigma_xi_squared = 0;
	double sigma1_over_xi = 0;
	double sigmaln_xi = 0;
	double sigmaln_xi_squared = 0;
	int i;
	for(i = 0; i < n; ++i)
	{
		sigma_xi += dataset[i];
		// For efficiency, I use += so a temporary variable isn't called
		sigma_xi_squared += dataset[i] * dataset[i];
		sigma1_over_xi += 1.0 / dataset[i];
		sigmaln_xi += log(dataset[i]);
		sigmaln_xi_squared += log(dataset[i]) * log(dataset[i]);
	}
	double mean_xi = sigma_xi / n;
	double mean_xi_squared = sigma_xi_squared / n;
	double mean_1_over_xi = sigma1_over_xi / n;
	double mean_ln_xi = sigmaln_xi / n;
	double mean_ln_xi_squared = sigmaln_xi_squared / n;

	/* PART THREE
	This "part" computes estimators and errors thereof; see Documentation.pdf
	I'll group by model; the models' order here is reproduced in later parts
	I've grouped Lognormal with Normal and Inverse-gamma with Gamma
	A future extension could also add more models */

	// std::cout << "Parameter estimation for the Normal model...\n";
	double mu = mean_xi;
	double sigma = sqrt(mean_xi_squared - mean_xi * mean_xi);
	double error_in_mu = sigma;
	double error_in_sigma = sigma * sqrt(1.0 / 6.0);

	// std::cout << "Parameter estimation for the Lognormal model...\n"; 
	double m = mean_ln_xi;
	double s = sqrt(mean_ln_xi_squared - mean_ln_xi * mean_ln_xi);
	double error_in_m = s;
	double error_in_s = s / sqrt(6.0);

	// std::cout << "Parameter estimation for the Exponential model...\n";
	flambda = (1.0 + 1.0 / n) / mean_xi;
	double error_in_flambda = sqrt(1.0 + 1.0 / n) / mean_xi;

	// std::cout << "Parameter estimation for the Gamma model...\n";
	double S = mean_ln_xi - log(mean_xi);
	fk = -0.5 / S;
	double numerator;
	double denominator;
	double step = 0.000001;
	for(i = 0; i < 10; ++i)
	{
		numerator = S + log(fk) - (lgamma(fk + step) - lgamma(fk)) / step;
		denominator = (lgamma(fk + 2 * step) + lgamma(fk) - 2 * lgamma(fk + step)) / (step * step) - 1.0 / fk;
		fk += numerator / denominator;
	}
	double error_in_k = fk / sqrt(n);
	ftheta = mean_xi / fk;
	double error_in_ftheta = ftheta / sqrt(n);

	// Inverse-gamma is similar to Gamma, respectively replacing mean_xi, mean_ln_xi with mean_1_over_xi, -mean_ln_xi
	
	// std::cout << "Parameter estimation for the Inverse-gamma model...\n";
	double inverse_S = -mean_ln_xi - log(mean_1_over_xi);
	finverse_k = -0.5 / inverse_S;
	for(i = 0; i < 10; ++i)
	{
		numerator = inverse_S + log(finverse_k) - (lgamma(finverse_k + step) - lgamma(finverse_k)) / step;
		denominator = (lgamma(finverse_k + 2 * step) + lgamma(finverse_k) - 2 * lgamma(finverse_k + step)) / (step * step) - 1.0 / finverse_k;		
		finverse_k += numerator / denominator;
	}
	double error_in_finverse_k = finverse_k / sqrt(n);
	inverse_ftheta = mean_1_over_xi / finverse_k;
	double error_in_inverse_ftheta = inverse_ftheta / sqrt(n);
	
	// std::cout << "Parameter estimation for the Pareto model...\n";
	pareto_minimum = dataset[0];
	pareto_alpha = 1.0 / (mean_ln_xi - log(pareto_minimum));
	double error_in_pareto_minimum = pareto_minimum / pareto_alpha;
	double error_in_pareto_alpha = pareto_alpha / sqrt(n);

	/* PART FOUR
	This part computes sums of squared vertical errors in CDFs, so it's a good thing the data has already been sorted
	A future extension of this code would cover any models added in PART THREE */

	// std::cout << "Testing quality of data fit to the Normal model...\n";
	double NormalCDF_sum_of_squared_vertical_errors = 0.0;
	double NormalCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		NormalCDF_latest_vertical_error = (i  * 1.0 + 1.0) / n - NormalCDF((dataset[i] - mu) / sigma);
		// The + 1.0 is needed because i = 0 corresponds to the 1st data point etc.
		NormalCDF_sum_of_squared_vertical_errors += NormalCDF_latest_vertical_error * NormalCDF_latest_vertical_error;
	}

	// Lognormal doesn't need a new CDF
	
	// std::cout << "Testing quality of data fit to the Lognormal model...\n";
	double LognormalCDF_sum_of_squared_vertical_errors = 0.0;
	double LognormalCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		LognormalCDF_latest_vertical_error = (i  * 1.0 + 1.0) / n - NormalCDF((log(dataset[i]) - m) / s);
		LognormalCDF_sum_of_squared_vertical_errors += LognormalCDF_latest_vertical_error * LognormalCDF_latest_vertical_error;
	}
	
	// std::cout << "Testing quality of data fit to the Exponential model...\n";
	double ExponentialCDF_sum_of_squared_vertical_errors = 0.0;
	double ExponentialCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		ExponentialCDF_latest_vertical_error = (i  * 1.0 + 1.0) / n - ExponentialCDF(dataset[i]);
		ExponentialCDF_sum_of_squared_vertical_errors += ExponentialCDF_latest_vertical_error * ExponentialCDF_latest_vertical_error;
	}

	// std::cout << "Testing quality of data fit to the Gamma model...\n";
	double GammaCDF_sum_of_squared_vertical_errors = 0.0;
	double GammaCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		GammaCDF_latest_vertical_error = (i  * 1.0 + 1.0) / n - GammaCDF(dataset[i] / ftheta);
		GammaCDF_sum_of_squared_vertical_errors += GammaCDF_latest_vertical_error * GammaCDF_latest_vertical_error;
	}

	// Inverse-gamma has a subtlety regarding the CDF
	
	// std::cout << "Testing quality of data fit to the Inverse-gamma model...\n";
	double InversegammaCDF_sum_of_squared_vertical_errors = 0.0;
	double InversegammaCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		InversegammaCDF_latest_vertical_error = (n - 1.0 * i) / n - InversegammaCDF(1.0 / (inverse_ftheta * dataset[i]));
		// i = 0 corresponds to the nth data point etc.
		InversegammaCDF_sum_of_squared_vertical_errors += InversegammaCDF_latest_vertical_error * InversegammaCDF_latest_vertical_error;
	}

	// std::cout << "Testing quality of data fit to the Pareto model...\n";
	double ParetoCDF_sum_of_squared_vertical_errors = 0.0;
	double ParetoCDF_latest_vertical_error;
	for(i = 0; i < n; ++i)
	{
		ParetoCDF_latest_vertical_error = (i  * 1.0 + 1.0) / n - ParetoCDF(dataset[i]);
		ParetoCDF_sum_of_squared_vertical_errors += ParetoCDF_latest_vertical_error * ParetoCDF_latest_vertical_error;
	}

	/* PART FIVE
	This part prints the findings of the analysis to Report.txt
	A future extension of this code would cover any models added in PART THREE */

	// First I explain the output the user will see
	std::cout << "This program has tried to fit several parametrised models to the provided data using a combination of Bayesian and maximum-likelihood methods. ";
	// I've written most of this one line at a time to keep code lines short enough for editors to read. However, I use line breaks in the output sparingly.
	std::cout << "The models used are Normal, Lognormal, Exponential, Gamma, Inverse-gamma and Pareto. ";
	// Note spaces before closed quotation marks (unless new line is being started) to prevent word.word strings
	std::cout << "For each model Report.txt presents estimates, with absolute errors, of the parameters. ";
	std::cout << "We can consistently compare all models' data fit by summing the vertical squared errors when their CDFs are compared with the empirical CDF. ";
	std::cout << "Each model's sum of vertical squared errors, hereafter Sigma Epsilon Squared, is also presented therein. ";
	std::cout << "For more information on the calculations this program used, consult Documentation.pdf.\n";// Line break needed prior to cout in PART SIX

	// Output to file

	std::ofstream Report;
	Report.open("Report.txt");
	Report << "The mu and error thereof, sigma and error thereof, and Sigma Epsilon Squared of the Normal fit are respectively\n";
	Report << mu << " " << error_in_mu << " " << sigma << " " << error_in_sigma << " " << NormalCDF_sum_of_squared_vertical_errors << "\n";
	Report << "The m and error thereof, s and error thereof, and Sigma Epsilon Squared of the Lognormal fit are respectively\n";
	Report << m << " " << error_in_m << " " << s << " " << error_in_s << " " << LognormalCDF_sum_of_squared_vertical_errors << "\n";
	Report << "The flambda, error thereof, and Sigma Epsilon Squared of the Exponential fit are respectively\n";
	Report << flambda << " " << error_in_flambda << " " << ExponentialCDF_sum_of_squared_vertical_errors << "\n";
	Report << "The shape parameter and error thereof, scale parameter and error thereof, and Sigma Epsilon Squared of the Gamma fit are respectively\n";
	Report << fk << " " << error_in_k << " " << ftheta << " " << error_in_ftheta << " " << GammaCDF_sum_of_squared_vertical_errors << "\n";
	Report << "The shape parameter and error thereof, scale parameter and error thereof, and Sigma Epsilon Squared of the Inverse-gamma fit are respectively\n";
	Report << finverse_k << " " << error_in_finverse_k << " " << inverse_ftheta << " " << error_in_inverse_ftheta << " " << InversegammaCDF_sum_of_squared_vertical_errors << "\n";
	Report << "The shape parameter and error thereof, scale parameter and error thereof in dollars, and Sigma Epsilon Squared of the Pareto fit are respectively\n";
	Report << pareto_minimum << " " << error_in_pareto_minimum << " " << pareto_alpha << " " << error_in_pareto_alpha << " " << ParetoCDF_sum_of_squared_vertical_errors << "\n";
	Report.close();

	/* PART SIX
	This part outputs to files
	A future extension of this code would cover any models added in PART THREE
	The number of samples per distribution may be changed too */

	std::cout << "Please wait while these distributions are randomly sampled, with output to Random_numbers.csv (we'll let you know when it's done).\n";
	std::ofstream Random_numbers;
	Random_numbers.open("Random_numbers.csv");
	Random_numbers << "Normal,Lognormal,Exponential,Gamma,Inverse Gamma,Pareto,\n";
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	unsigned int N = 100000; // Number of random samples per distribution
	for(i = 0; i < N; ++i)
	{
		double Normal_Output = mu + sigma * Normal_Inverse(dis(gen));
		double Lognormal_Output = exp(m + s * Normal_Inverse(dis(gen)));
		double Exponential_Output = log(1.0 / dis(gen)) / flambda;
		double Gamma_Output;
		Gamma_Output = ftheta * Gamma_Inverse(dis(gen));
		while(isinf(Gamma_Output) || isnan(Gamma_Output))
		{// Without this loop, a minority of outcomes diverge
			Gamma_Output = ftheta * Gamma_Inverse(dis(gen));// Testing showed which distributions this problem affects
		}
		double Inversegamma_Output;
		Inversegamma_Output = 1.0 / (inverse_ftheta * Inverse_gamma_Inverse(dis(gen)));
		while(isinf(Inversegamma_Output) || isnan(Inversegamma_Output))
		{
			Inversegamma_Output = 1.0 / (inverse_ftheta * Inverse_gamma_Inverse(dis(gen)));
		}
		double Pareto_Output = pareto_minimum * exp(-log(dis(gen)) / pareto_alpha);
		Random_numbers << Normal_Output << "," << Lognormal_Output << "," << Exponential_Output << "," << Gamma_Output << "," << Inversegamma_Output << "," << Pareto_Output << ",\n";
		// The next few lines are helpful if, during testing, a slow output of random numbers is suspected
		// if(i % 1000 == 999)
		// {
			// std::cout << i + 1 << " random values for each model have been output so far.\n";
		// }
	}
	Random_numbers.close();
	std::cout << "The file Random_numbers.csv now contains ";
	std::cout << N << " random values from each model's distribution with the above parameters.\n";

	// These three lines keep the program from closing, so the reader can read about the output files (it also helps with debugging)
	std::cout << "Press ENTER to close.";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return 0;

}// int main() closes