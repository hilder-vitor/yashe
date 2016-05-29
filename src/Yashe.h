/*

Copyright or © or Copr. Tancrede Lepoint.

Tancrede.Lepoint@cryptoexperts.com

This software is a computer program whose purpose is to provide to the 
research community a proof-of-concept implementation of the homomorphic 
evaluation of the lightweight block cipher SIMON, describe in the paper
"A Comparison of the Homomorphic Encryption Schemes FV and YASHE", of
Tancrede Lepoint and Michael Naehrig, available at
http://eprint.iacr.org/2014.

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

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

#ifndef __YASHE__YASHEKEY
#define __YASHE__YASHEKEY

#include <ostream>
#include <vector>
#include <map>
#include "Ciphertext.h"
#include "Plaintext.h"
#include "RealNumberPlaintext.h"
#include "RealNumberCiphertext.h"
#include "Sampler.h"
#include "macros.h"
#include "matrixutils.h"
#include "timing.h"

#include "gmpxx.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include "flint/fmpz_mod_polyxx.h"
#include "flint/fmpz_mod_poly_factorxx.h"

using namespace flint;

/**********************************************************************
 * YASHE Parameters.
 **********************************************************************/
struct YASHEParams {
    unsigned long n, sigmakey, sigmaerr;
    fmpzxx q;
    fmpzxx t;
};
#define WORDLENGTH 72

class Ciphertext;

typedef std::vector<double> DoubleVector;
typedef std::vector<Plaintext> PlaintextVector;
typedef std::vector<RealNumberPlaintext> RealNumberPlaintextVector;
typedef std::vector<Ciphertext> CiphertextVector;
typedef std::vector<RealNumberCiphertext> RealNumberCiphertextVector;

typedef std::vector<std::vector<double> > DoubleMatrix;
typedef std::vector<std::vector<Plaintext> > PlaintextMatrix;
typedef std::vector<std::vector<RealNumberPlaintext> > RealNumberPlaintextMatrix;
typedef std::vector<std::vector<Ciphertext> > CiphertextMatrix;
typedef std::vector<std::vector<RealNumberCiphertext> > RealNumberCiphertextMatrix;


class Yashe {
    public:
	unsigned long	n, ell, sigmakey, sigmaerr, logwq;
	fmpzxx q, t, qdivt, tmodq, qdiv2t;

	fmpz_polyxx poly;
	fmpz_mod_polyxx * phi;
	fmpz_mod_polyxx * h, * f;
	std::vector<fmpz_mod_polyxx> gamma;

	Entropy random;
	Sampler *sampler;

	Yashe(const struct YASHEParams& params);
	Yashe(const Yashe& y);
	Yashe(const std::string& filename);
	void serialize(const std::string& filename);
	void unserialize(const std::string& filename);

	inline size_t	get_ell() const { return ell; }
	inline fmpz_mod_polyxx	get_phi() const { return *phi; }
	inline unsigned long get_logwq() const { return logwq; }
	inline fmpzxx get_q() const { return q; }
	inline fmpzxx get_t() const { return t; }
	inline fmpz_mod_polyxx get_gamma(unsigned i) const { return gamma[i]; }

	friend std::ostream& operator<<(std::ostream&, const Yashe&);

	/* function to encrypt */
	Ciphertext encrypt(const Plaintext&);
	RealNumberCiphertext encrypt(const RealNumberPlaintext&);

	define_function( encrypt, PlaintextVector, CiphertextVector);
	define_function( encrypt, PlaintextMatrix, CiphertextMatrix);
	define_function( encrypt, RealNumberPlaintextVector, RealNumberCiphertextVector);
	define_function( encrypt, RealNumberPlaintextMatrix, RealNumberCiphertextMatrix);

	/* function to decrypt */
	Plaintext decrypt(const Ciphertext&);
	RealNumberPlaintext decrypt(const RealNumberCiphertext&);
	define_function(decrypt, CiphertextVector, PlaintextVector);
	define_function(decrypt, CiphertextMatrix, PlaintextMatrix);
	define_function(decrypt, RealNumberCiphertextVector, RealNumberPlaintextVector);
	define_function(decrypt, RealNumberCiphertextMatrix, RealNumberPlaintextMatrix);

	/* function to encode */
	RealNumberPlaintext encode(const double& value, int exp_shift = 64);
	define_function(encode, DoubleVector, RealNumberPlaintextVector);
	define_function(encode, DoubleMatrix, RealNumberPlaintextMatrix);

	/* function to decode */
	double decode(const RealNumberPlaintext& value);
	define_function(decode, RealNumberPlaintextVector, DoubleVector);
	define_function(decode, RealNumberPlaintextMatrix, DoubleMatrix);

	fmpz_mod_polyxx	convert(const fmpz_mod_polyxx&);
	unsigned 		noise(const Ciphertext&);

	/* decrypt, decode, encode and encrypt.
	 * This function can receive RealNumberCiphertext,
	 * vectors and matrices of this type, etc... */
	template<typename Ciphertext_like>
	void recrypt(Ciphertext_like& c){
		timing timing;
		timing.start();
		c = encrypt(encode(decode(decrypt(c))));
		timing.stop("recrypt: ");
	}



	SymmetricMatrix<RealNumberCiphertext> encrypt(const SymmetricMatrix<RealNumberPlaintext>& plain);

	SymmetricMatrix<RealNumberPlaintext> decrypt(const SymmetricMatrix<RealNumberCiphertext>& c);

	SymmetricMatrix<RealNumberPlaintext> encode(SymmetricMatrix<double> messages, int precision = 64);

	SymmetricMatrix<double> decode(SymmetricMatrix<RealNumberPlaintext> plain);

};	
#endif

