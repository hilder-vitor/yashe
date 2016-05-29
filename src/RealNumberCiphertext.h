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

#ifndef __YASHE__REAL_NUMBER_CIPHERTEXT
#define __YASHE__REAL_NUMBER_CIPHERTEXT

#include "Ciphertext.h"
#include "RealNumberPlaintext.h"

using namespace flint;

class RealNumberCiphertext : public Ciphertext {
    private:
		long int exp_shift;
    public:
	RealNumberCiphertext(Yashe&, const fmpz_mod_polyxx& poly, long int expoent_shift);
	RealNumberCiphertext(Yashe&, const std::string& filename);
	RealNumberCiphertext(const RealNumberCiphertext& c);
	RealNumberCiphertext(const Ciphertext& c, long int _exp_shift = 64);

	void serialize(const std::string& filename) const;

	RealNumberCiphertext& operator=(const RealNumberCiphertext&);
	RealNumberCiphertext& operator+=(const RealNumberCiphertext&);
	RealNumberCiphertext& operator-=(const RealNumberCiphertext&);
	RealNumberCiphertext& operator*=(const RealNumberCiphertext&);
	RealNumberCiphertext& operator*=(const RealNumberPlaintext&);

	RealNumberCiphertext& operator+=(const double&);
	RealNumberCiphertext& operator-=(const double&);
	RealNumberCiphertext& operator*=(const double&);

	RealNumberCiphertext& operator/=(const double& N);
	

	inline long int expoent_shift() const { return exp_shift; }

	RealNumberCiphertext operator+(const RealNumberCiphertext& c) const ;
	RealNumberCiphertext operator-(const RealNumberCiphertext& c) const;
	RealNumberCiphertext operator*(const RealNumberCiphertext& c) const;
	RealNumberCiphertext operator*(const RealNumberPlaintext& c) const;
	RealNumberCiphertext operator+(const double& d) const;
	RealNumberCiphertext operator-(const double& d) const;
	RealNumberCiphertext operator*(const double& d) const;
	RealNumberCiphertext operator/(const double& N) const;

	void divide_by_power_of_two(int k);
};	


std::ostream& operator<<(std::ostream& os, const RealNumberCiphertext& c);

#endif
