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

#ifndef __YASHE__REAL_NUMBER_PLAINTEXT
#define __YASHE__REAL_NUMBER_PLAINTEXT

#include "Plaintext.h"
using namespace flint;


class RealNumberPlaintext: public Plaintext {
    private:
	int exp_shift = 64;
    
    public:
	RealNumberPlaintext(const fmpzxx& coeff_modulus, const fmpz_mod_polyxx& cyclotomic, int exp_shift = 64);
	RealNumberPlaintext(const fmpz_mod_polyxx& poly, const fmpz_mod_polyxx& cyclotomic, int exp_shift = 64);
	RealNumberPlaintext(const Plaintext& p, long int _exp_shif = 64);
	RealNumberPlaintext(const RealNumberPlaintext& p);
	
	bool operator==(const RealNumberPlaintext&) const;

	RealNumberPlaintext& operator=(const RealNumberPlaintext&);
	RealNumberPlaintext& operator+=(const RealNumberPlaintext&);
	RealNumberPlaintext& operator*=(const RealNumberPlaintext&);
	RealNumberPlaintext& operator*=(const unsigned int&);
	RealNumberPlaintext& operator*=(const fmpzxx&);

	RealNumberPlaintext operator+(const RealNumberPlaintext& p) const;
	RealNumberPlaintext operator*(const RealNumberPlaintext& p) const;
	RealNumberPlaintext operator*(const unsigned int& y) const;
	RealNumberPlaintext operator*(const fmpzxx& y) const;

	int expoent_shift() const;

	double double_value() const;
};

std::ostream& operator<<(std::ostream& os, const RealNumberPlaintext& p);
#endif
