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

#include "RealNumberPlaintext.h"
#include "Yashe.h"
#include <iostream>
#include <assert.h>

/* Constructors */
RealNumberPlaintext::RealNumberPlaintext(const fmpzxx& coeff_modulus, const fmpz_mod_polyxx& _cyclotomic, int _exp_shift) : Plaintext(coeff_modulus, _cyclotomic), exp_shift(_exp_shift){
}

RealNumberPlaintext::RealNumberPlaintext(const fmpz_mod_polyxx& poly, const fmpz_mod_polyxx& _cyclotomic, int exp_shift) : Plaintext(poly, _cyclotomic), exp_shift(exp_shift) {
}

RealNumberPlaintext::RealNumberPlaintext(const Plaintext& p, long int _exp_shift) : Plaintext(p), exp_shift(_exp_shift) {
}

RealNumberPlaintext::RealNumberPlaintext(const RealNumberPlaintext& p) : Plaintext(p), exp_shift(p.expoent_shift()){
}

bool RealNumberPlaintext::operator==(const RealNumberPlaintext& p) const{
	return p.expoent_shift() == exp_shift && Plaintext::operator==(p);
}

RealNumberPlaintext& RealNumberPlaintext::operator=(const RealNumberPlaintext& p){
	Plaintext::operator=(p);
	exp_shift = p.expoent_shift();
    return *this;
}

/* Addition */
RealNumberPlaintext& RealNumberPlaintext::operator+=(const RealNumberPlaintext& p) {
	if (p.expoent_shift() != exp_shift)
		throw std::invalid_argument("Addition of plaintexts with different expoent shift not implemented yet.");
	Plaintext::operator+=(p);
    return *this;
}

RealNumberPlaintext RealNumberPlaintext::operator+(const RealNumberPlaintext& p) const{
	RealNumberPlaintext p2(*this);
	return (p2 += p);
}

/* Multiplication */
RealNumberPlaintext& RealNumberPlaintext::operator*=(const RealNumberPlaintext& p) {
	Plaintext::operator*=(p);
	exp_shift += p.expoent_shift();
    return *this;
}

RealNumberPlaintext& RealNumberPlaintext::operator*=(const unsigned int& y){
	Plaintext::operator*=(y);
	return *this;
}

RealNumberPlaintext& RealNumberPlaintext::operator*=(const fmpzxx& y){
	Plaintext::operator*=(y);
	return *this;
}

RealNumberPlaintext RealNumberPlaintext::operator*(const RealNumberPlaintext& p) const {
	RealNumberPlaintext p2(*this);
	p2 *= p;
	return p2;
}

RealNumberPlaintext RealNumberPlaintext::operator*(const unsigned int& y) const{
	RealNumberPlaintext prod(*this);
	prod *= y;
	return prod;
}

RealNumberPlaintext RealNumberPlaintext::operator*(const fmpzxx& y) const{
	RealNumberPlaintext prod(*this);
	prod *= y;
	return prod;
}

int RealNumberPlaintext::expoent_shift() const{
	return exp_shift;
}


double RealNumberPlaintext::double_value() const{
	mpf_set_default_prec(2048);
	mpf_class value = 0;
	mpf_class tmp = 0;

	mpf_class double_power = 0.5;
	fmpzxx coeff;
	fmpzxx t = modulus();
	
	// do not consider expoents smaller than -128
	//for (int i = exp_shift - 1; i >= exp_shift - 128 && i >= 0; i--){
	for (int i = exp_shift - 1; i >= exp_shift - 64 && i >= 0; i--){
		coeff = get(i);
		if (2*coeff > t)
			coeff -= t;
		if (!coeff.is_zero()){
			tmp = coeff.to_string();
			value += tmp * double_power; 
		}
		double_power /= 2;
	}

	// decode the integer part
	unsigned int deg = (cyclo.degree() > degree() ? degree() : cyclo.degree());
	double_power = 1;
	for (unsigned int i = exp_shift; i <= deg; i++){
		coeff = get(i);
		if (!coeff.is_zero()){
			if (2*coeff > t)
				coeff -= t;
			tmp = coeff.to_string();
			value += tmp * double_power; 
		}
		double_power *= 2;
	}
	value.set_prec(64);
	return value.get_d();
}

std::ostream& operator<<(std::ostream& os, const RealNumberPlaintext& p) {
	os << "<RealNumberPlaintext: ";
	os << static_cast<Plaintext>(p);
	os << " with expoent_shift =" << p.expoent_shift();
    os << ">";
    return os;
}
