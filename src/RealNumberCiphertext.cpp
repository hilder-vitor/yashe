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

#include "RealNumberCiphertext.h"
#include "Yashe.h"
#include <iostream>
#include <assert.h>
#include <cmath>

int is_power_of_two(double val){
	double log2Val = log2(val);
	if (-0.000001 < (log2Val - std::round(log2Val)) && (log2Val - std::round(log2Val)) < 0.000001){
		return std::round(log2Val);	
	}
	return 0;
}

/* Constructors */
RealNumberCiphertext::RealNumberCiphertext(Yashe& pk, const fmpz_mod_polyxx& poly, long int _exp_shift) : 
    Ciphertext(pk,poly), exp_shift(_exp_shift){
}

RealNumberCiphertext::RealNumberCiphertext(Yashe& pk, const std::string& filename) : Ciphertext(pk, filename), exp_shift(0){
	std::string line;
	std::ifstream inFile(filename);
	while (inFile.good()){
		getline(inFile, line, ':');
		if ("exp_shift" == line){
			inFile >> exp_shift;
			getline(inFile, line); // ignores end of line
		}else{
			getline(inFile, line); // ignores the line
		}
	}
}

RealNumberCiphertext::RealNumberCiphertext(const RealNumberCiphertext& c) : Ciphertext(c), exp_shift(c.expoent_shift()){
}

RealNumberCiphertext::RealNumberCiphertext(const Ciphertext& c, long int _exp_shift) : Ciphertext(c), exp_shift(_exp_shift){
}

RealNumberCiphertext& RealNumberCiphertext::operator=(const RealNumberCiphertext& c) {
	Ciphertext::operator=(c);
	exp_shift = c.expoent_shift();
    return *this;
}

/* Operators */
std::ostream& operator<<(std::ostream& os, const RealNumberCiphertext& c) {
	os << "<RealNumberCiphertext: ";
	os << static_cast<Ciphertext>(c);
	os << " with expoent_shift =" << c.expoent_shift();
    os << ">";
    return os;
}

void RealNumberCiphertext::serialize(const std::string& filename) const{
	Ciphertext::serialize(filename);
	std::ofstream outFile(filename, std::ofstream::app);
	outFile << "exp_shift:" << exp_shift << endl;
}

/* Addition */
RealNumberCiphertext& RealNumberCiphertext::operator+=(const RealNumberCiphertext& c) {
	unsigned int complementar_precision;
	if (c.expoent_shift() == exp_shift){
		Ciphertext::operator+=(c); // add as simple ciphertexts
	}else if (c.expoent_shift() > exp_shift){
		complementar_precision = c.expoent_shift() - exp_shift;
		RealNumberPlaintext p = yashe.encode(1, complementar_precision);
		Ciphertext::operator*=(p);
		exp_shift = c.expoent_shift();
		Ciphertext::operator+=(c);
	}else{
		complementar_precision =  exp_shift - c.expoent_shift();
		RealNumberPlaintext p = yashe.encode(1, complementar_precision);
		Plaintext plain(p);
		Ciphertext::operator+=(c.Ciphertext::operator*(plain));
	}
    return *this;
}

RealNumberCiphertext RealNumberCiphertext::operator+(const RealNumberCiphertext& c) const {
    RealNumberCiphertext c2(*this);
    return (c2 += c);
}

RealNumberCiphertext& RealNumberCiphertext::operator-=(const RealNumberCiphertext& c) {
	(*this) += c * (-1.0);
	return *this;	
}

RealNumberCiphertext RealNumberCiphertext::operator-(const RealNumberCiphertext& c) const{
	RealNumberCiphertext c2(*this);
	return (c2 -= c);
}


/* Multiply  */
RealNumberCiphertext& RealNumberCiphertext::operator*=(const RealNumberCiphertext& c) {
	Ciphertext::operator*=(c);
	exp_shift += c.expoent_shift();
    return *this;
}

RealNumberCiphertext& RealNumberCiphertext::operator*=(const RealNumberPlaintext& p) {
	Ciphertext::operator*=(p);
	exp_shift += p.expoent_shift();
    return *this;
}

RealNumberCiphertext RealNumberCiphertext::operator*(const RealNumberCiphertext& c) const {
	RealNumberCiphertext c2(*this);
    c2 *= c;
    return c2;
}

RealNumberCiphertext RealNumberCiphertext::operator*(const RealNumberPlaintext& p) const {
	RealNumberCiphertext c2(*this);
    c2 *= p;
    return c2;
}


/* -------- Plaintext (double) operations ---- */

RealNumberCiphertext& RealNumberCiphertext::operator+=(const double& d) {
	RealNumberPlaintext p = yashe.encode(d, exp_shift);
	Ciphertext::operator+=(p);
    return *this;
}

RealNumberCiphertext RealNumberCiphertext::operator+(const double& d) const{
    RealNumberCiphertext c2(*this);
    c2 += d;
    return c2;
}

RealNumberCiphertext& RealNumberCiphertext::operator-=(const double& d) {
	(*this) += (-1.0 * d);
	return *this;
}

RealNumberCiphertext RealNumberCiphertext::operator-(const double& d) const{
    RealNumberCiphertext c2(*this);
    c2 -= d;
    return c2;
}


RealNumberCiphertext& RealNumberCiphertext::operator*=(const double& d){
	int prec;
	if (-1 == d || 1 == d)
		prec = 0;
	else
		prec = 64;
	RealNumberPlaintext p = yashe.encode(d, prec);
	Ciphertext::operator*=(p);
	exp_shift += prec;
    return *this;
}

RealNumberCiphertext RealNumberCiphertext::operator*(const double& d) const{
	RealNumberCiphertext c2(*this);
	c2 *= d;
	return c2;
}



/* Divide */
RealNumberCiphertext& RealNumberCiphertext::operator/=(const double& N){
	int expoent_power_of_two;
	int signal = (N < 0 ? -1 : 1); 
	if ((expoent_power_of_two = is_power_of_two(signal * N))){
		divide_by_power_of_two(expoent_power_of_two);
		if (signal == -1){
			Ciphertext::operator*=(yashe.encode(-1.0, 0));
		}
	}else{
		(*this) *= (1.0/N);
	}
	return *this;
}


RealNumberCiphertext RealNumberCiphertext::operator/(const double& N) const{
	RealNumberCiphertext c2(*this);
	c2 /= N;
	return c2;
}



void RealNumberCiphertext::divide_by_power_of_two(int k){ 
	exp_shift += k;
}

