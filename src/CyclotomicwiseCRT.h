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

#ifndef __YASHE__CYCLOTOMIC_CRT
#define __YASHE__CYCLOTOMIC_CRT

#include "flint/fmpz_mod_polyxx.h"
#include "flint/fmpqxx.h"
#include "flint/arith.h"
#include <vector>
#include "Plaintext.h"

class CyclotomicwiseCRT {
    private:
	// all the coefficients of all the polynomials here belong to Z_t
	flint::fmpzxx t;
	// vector with the cyclotomic polynomials
	std::vector<flint::fmpz_mod_polyxx> coprimes;
	// product of the cyclotomic polynomials
	flint::fmpz_mod_polyxx poly_modulus;

	// each position k has the poly_modulus divided by coprimes[k]
	std::vector<flint::fmpz_mod_polyxx> partial_modulus;
	// each position k has the inverse of partial_modulus[k] mod coprimes[k]	
	std::vector<flint::fmpz_mod_polyxx> inverses;

    public:
	CyclotomicwiseCRT(const std::vector<flint::fmpz_mod_polyxx>& coprimes);

	/* This function returns the n-th cyclotomic polynomial with coefficients reduced mod t
	 * (which means it is represented as an element of Z_t[x]) */
	static flint::fmpz_mod_polyxx n_th_cyclotomic_mod_t(unsigned int n, const fmpzxx t);
	/* This method receives a cyclotomic polynomial with coefficients in Z_t and factors it
	 * in polynomials in Z_t[x].                                                              */
	static std::vector<flint::fmpz_mod_polyxx> factorize_cyclotomic(const fmpz_mod_polyxx& cyclotomic);

	const flint::fmpz_mod_polyxx& get_modulus() const;

	Plaintext pack(const std::vector<Plaintext>& remainders);
	
	std::vector<Plaintext> unpack(const Plaintext& remainders);
	
};	

#endif
