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

#ifndef __YASHE__PLAINTEXT
#define __YASHE__PLAINTEXT

#include <iostream>
#include <sstream>
#include <string>
#include "flint/fmpz_mod_polyxx.h"
#include "flint/fmpqxx.h"
using namespace flint;

class Yashe;

class Plaintext {
    protected:
    fmpz_mod_polyxx pval;
    const fmpz_mod_polyxx& cyclo;
    
    public:
	Plaintext(const fmpzxx& coeff_modulus, const fmpz_mod_polyxx& cyclotomic);
	Plaintext(const fmpz_mod_polyxx& poly, const fmpz_mod_polyxx& cyclotomic);
	Plaintext(const Plaintext& p);
	
	friend std::ostream& operator<<(std::ostream&, const Plaintext&);

	bool operator==(const Plaintext&) const;

	Plaintext& operator=(const Plaintext&);
	Plaintext& operator+=(const Plaintext&);
	Plaintext& operator*=(const Plaintext&);
	Plaintext& operator*=(const unsigned int&);
	Plaintext& operator*=(const fmpzxx&);

	Plaintext operator+(const Plaintext& p) const;
	Plaintext operator*(const Plaintext& p) const;
	Plaintext operator*(const unsigned int& y) const;
	Plaintext operator*(const fmpzxx& y) const;

	unsigned int degree() const;
	fmpzxx get(unsigned int i) const;
	void set(unsigned int i, fmpzxx coefficient);
	void set(unsigned int i, int coefficient);
	fmpzxx modulus() const;
	const fmpz_mod_polyxx& cyclotomic() const;
	fmpz_mod_polyxx polynomial() const;
	bool is_zero() const;
};

#endif
