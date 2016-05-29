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

#ifndef __YASHE__COEFFICIENT_CRT
#define __YASHE__COEFFICIENT_CRT

#include "flint/fmpz_mod_polyxx.h"
#include "flint/fmpqxx.h"
#include <vector>
#include "Plaintext.h"

class CoefficientwiseCRT {
    private:
	const flint::fmpz_mod_polyxx& cyclotomic;
	std::vector<flint::fmpzxx> coprimes;
	flint::fmpzxx modulus; // product of the coprimes

	flint::fmpzxx integerCRT(std::vector<flint::fmpzxx> remainders);
	std::vector<flint::fmpzxx> invIntegerCRT(flint::fmpzxx x);
    
    public:
	CoefficientwiseCRT(const std::vector<flint::fmpzxx>& coprimes, const flint::fmpz_mod_polyxx& cyclo); 

	fmpzxx get_modulus() const;

	Plaintext pack(std::vector<Plaintext> remainders);
	
	std::vector<Plaintext> unpack(Plaintext remainders);
	
};	

#endif
