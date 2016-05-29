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

#include "Plaintext.h"
#include "Yashe.h"
#include <iostream>
#include <assert.h>

/* Constructors */
Plaintext::Plaintext(const fmpzxx& coeff_modulus, const fmpz_mod_polyxx& _cyclotomic) : pval(coeff_modulus), cyclo(_cyclotomic){
}

Plaintext::Plaintext(const fmpz_mod_polyxx& poly, const fmpz_mod_polyxx& _cyclotomic) :pval(poly), cyclo(_cyclotomic) {
}
	
Plaintext::Plaintext(const Plaintext& p) : pval(p.modulus()), cyclo(p.cyclotomic())
{
	pval = p.polynomial();
}

bool Plaintext::operator==(const Plaintext& p) const{
	return p.polynomial() == pval && p.modulus() == modulus();
}

Plaintext& Plaintext::operator=(const Plaintext& p){
    pval = fmpz_mod_polyxx(p.polynomial());
    return *this;
}


std::ostream& operator<<(std::ostream& os, const Plaintext& p) {
	unsigned int deg = p.cyclo.degree();
	fmpzxx zero("0");
	bool is_null_poly = true;
	for (unsigned int i = deg; i >= 1; i--){
		if (p.get(i) != zero){
			if (p.get(i) * 2 >= p.modulus())
				os << p.get(i) - p.modulus() << "x^" << i << " ";
			else
				os << "+" << p.get(i) << "x^" << i << " ";
			is_null_poly = false;
		}
	}
	
	if (p.get(0) != zero){
		if (p.get(0) * 2 >= p.modulus())
			os << p.get(0) - p.modulus() << " ";
		else
			os << "+" << p.get(0) << " ";

		is_null_poly = false;
	}
	if (is_null_poly){
		os << "0";
	}
	return os; 
}

/* Addition */
Plaintext& Plaintext::operator+=(const Plaintext& p) {
    pval += (p.polynomial() % cyclo);
    return *this;
}

Plaintext Plaintext::operator+(const Plaintext& p) const{
	Plaintext p2(*this);
	return (p2 += p);
}

/* Multiplication */
Plaintext& Plaintext::operator*=(const Plaintext& p) {
	pval *= (p.polynomial() % cyclo);
    return *this;
}

Plaintext& Plaintext::operator*=(const unsigned int& y){
	pval *= fmpzxx(y);
	return *this;
}

Plaintext& Plaintext::operator*=(const fmpzxx& y){
	pval *= y;
	return *this;
}


Plaintext Plaintext::operator*(const Plaintext& p) const {
	Plaintext p2(*this);
	p2 *= p;
	return p2;
}

Plaintext Plaintext::operator*(const unsigned int& y) const{
	Plaintext prod(*this);
	prod *= y;
	return prod;
}

Plaintext Plaintext::operator*(const fmpzxx& y) const{
	Plaintext prod(*this);
	prod *= y;
	return prod;
}


unsigned int Plaintext::degree() const{
	return pval.degree();
}

fmpzxx Plaintext::get(unsigned int i) const{
	fmpzxx coef(pval.get_coeff(i));
	return coef;
}

void Plaintext::set(unsigned int i, fmpzxx coefficient){
	pval.set_coeff(i, coefficient);
}

void Plaintext::set(unsigned int i, int coefficient){
	pval.set_coeff(i, fmpzxx(coefficient));
}

fmpzxx Plaintext::modulus() const{
	fmpzxx _mod(pval.modulus());
	return _mod;
}

const fmpz_mod_polyxx& Plaintext::cyclotomic() const{
	return cyclo;
}
fmpz_mod_polyxx Plaintext::polynomial() const{
	fmpz_mod_polyxx p(pval);
	return p;
}

