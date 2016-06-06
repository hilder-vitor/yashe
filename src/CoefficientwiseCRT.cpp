/*
This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

#include "CoefficientwiseCRT.h"

using namespace flint;
using namespace std;

flint::fmpzxx CoefficientwiseCRT::integerCRT(vector<flint::fmpzxx> remainders) const{
	fmpzxx x("0");
    fmpzxx t = modulus;
	fmpzxx t_;
	unsigned int N = coprimes.size();
	for (unsigned int i = 0; i < N; i++){
		t_ = t / coprimes[i];
		x += remainders[i] * t_ * invmod(t_, coprimes[i]);
	}
	x = x % t;
	return x;
}

vector<flint::fmpzxx> CoefficientwiseCRT::invIntegerCRT(flint::fmpzxx x) const{
	unsigned int N = coprimes.size();
	vector<fmpzxx> vec(coprimes.size());
	for (unsigned int i = 0; i < N; i++){
		vec[i] = x % coprimes[i];	
	}
	return vec;
}

fmpzxx product(const vector<fmpzxx>& numbers){
	fmpzxx x("1");
	unsigned int N = numbers.size();
	for (unsigned int i = 0; i < N; i++){
		x *= numbers[i];
	}
	return x;
}


CoefficientwiseCRT::CoefficientwiseCRT(const vector<flint::fmpzxx>& _coprimes, const flint::fmpz_mod_polyxx& cyclo) : cyclotomic(cyclo), coprimes(_coprimes), modulus(product(_coprimes)) {
}

const std::vector<flint::fmpzxx>& CoefficientwiseCRT::get_coprimes() const{
	return coprimes;
}

fmpzxx CoefficientwiseCRT::get_modulus() const {
	return modulus;
}


unsigned int maximum_degree(const vector<Plaintext>& polys){
	unsigned int N = polys.size();
	unsigned int max = 0;
	for (unsigned int i = 1; i < N; i++){
		if (!polys[i].is_zero() && polys[i].degree() > max)
			max = polys[i].degree();
	}
	return max;
}


Plaintext CoefficientwiseCRT::pack(vector<Plaintext> polys) const{
	unsigned int N = polys.size();
	fmpzxx t = modulus;
	vector<fmpzxx> coefficients(N);
	Plaintext resp(t, cyclotomic);
	unsigned int max_degree = maximum_degree(polys);
	for (unsigned int i = 0; i <= max_degree; i++){
		// pick the i-th coefficient of each polynomial
		for (unsigned int j = 0; j < N; j++)
			coefficients[j] = polys[j].get(i);
		fmpzxx coeff = integerCRT(coefficients);
		resp.set(i, coeff);	
	}
	return resp;
}

vector<Plaintext> CoefficientwiseCRT::unpack(Plaintext packaged) const{
	unsigned int N = coprimes.size();
	vector<Plaintext> polys;
	
	for (unsigned int i = 0; i < N; i++)
		polys.push_back(Plaintext(coprimes[i], cyclotomic));
	
	unsigned int deg = (packaged.is_zero() ? 1 : packaged.degree());
	for (unsigned int i = 0; i <= deg; i++){
		fmpzxx p_coefficient(packaged.get(i));
		vector<fmpzxx> coefficients = invIntegerCRT(p_coefficient);
		
		// set the j-th coefficient of each polynomial
		for (unsigned int j = 0; j < N; j++){
			polys[j].set(i, coefficients[j]);
		}
	}

	return polys;
}
