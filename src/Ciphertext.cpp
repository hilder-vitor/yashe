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

#include "Ciphertext.h"
#include <iostream>
#include <assert.h>
#include "Yashe.h"

/* Constructors */
Ciphertext::Ciphertext(const Yashe& pk, const fmpz_mod_polyxx& poly) : 
    yashe(pk), level(0), cval(poly), aftermult(false)
{
}

Ciphertext::Ciphertext(const Yashe& scheme, const std::string& filename) : yashe(scheme), level(0), cval(scheme.q), aftermult(false) {
	unserialize(filename);
}

Ciphertext::Ciphertext(const Ciphertext& c) : 
    yashe(c.yashe), level(c.level),
    cval(c.cval), aftermult(c.aftermult) {
}

Ciphertext& Ciphertext::operator=(const Ciphertext& c) {
    assert(&yashe == &c.yashe);
    level = c.level;
    cval = c.cval;
    aftermult = c.aftermult;
    return *this;
}

/* Operators */
std::ostream& operator<<(std::ostream& os, const Ciphertext& c) {
    os << "<Ciphertext of level "
        << std::dec << c.level;
    
    if (c.aftermult)
        os << " of type II";
	
	os << " with max(noise)="  << std::dec << c.real_noise() << " bits";

    os << " of m=" << c.yashe.decrypt(c);
    os << ">";
    return os;
}

/* Addition */
Ciphertext& Ciphertext::operator+=(const Ciphertext& c) {
    if(aftermult != c.aftermult) {
    if(aftermult)
    {
        convert_self();
    }
    else
    {
        c.convert_self();
    }
    }

    cval += c.cval;
    level = (c.level>level) ? c.level : level;
    
    return *this;
}

Ciphertext Ciphertext::operator+(const Ciphertext& c) const {
	Ciphertext c2(*this);
	return (c2 += c);
}

/* Subtraction */ 
Ciphertext& Ciphertext::operator-=(const Ciphertext& c){
    if(aftermult != c.aftermult){
	    if(aftermult)
        	convert_self();
    	else
    	    c.convert_self();
    }

	cval -= c.cval;
    level = (c.level>level) ? c.level : level;
    return *this;
}
Ciphertext Ciphertext::operator-(const Ciphertext& c) const {
	Ciphertext c2(*this);
	return (c2 -= c);
}

/* Multiply 

cval = [near(t/q * cval * c.cval)]_q;
*/
Ciphertext& Ciphertext::operator*=(const Ciphertext& c) {
    if (aftermult)
        convert_self();
    if (c.aftermult)
        c.convert_self();
    fmpzxx one = fmpzxx(1), q=yashe.get_q();
    fmpzxx coeff, diff;
    fmpz_polyxx g;

    g = yashe.get_t()*(cval.to<fmpz_polyxx>()*c.cval.to<fmpz_polyxx>())%(yashe.poly);

    for (unsigned i=0; i<=g.degree(); i++)
    {
        ltupleref(coeff, diff) = fdiv_qr(g.get_coeff(i), q);
        if ( 2*diff > q )
            cval.set_coeff(i, coeff+one);
        else
            cval.set_coeff(i, coeff);
    }

    level = (c.level>level) ? c.level+1 : level+1;
    aftermult = true;

    return *this;
}

Ciphertext Ciphertext::operator*(const Ciphertext& c) const {
	Ciphertext c2(*this);
	c2 *= c;
	return c2;
}

/* -- Plaintexts operations -- */
Ciphertext& Ciphertext::operator+=(const Plaintext& p){
	if(aftermult)
		convert_self();

	// operations has to be done modulo q (not t)
	fmpz_mod_polyxx _poly_p(cval.modulus());
	_poly_p = p.polynomial();
	cval += yashe.qdivt * _poly_p;
	return *this;
}
	
Ciphertext Ciphertext::operator+(const Plaintext& p) const{
	Ciphertext c2(*this);
	return (c2 += p);
}

Ciphertext& Ciphertext::operator-=(const Plaintext& p){
	if(aftermult)
		convert_self();

	fmpz_mod_polyxx p_mod_q(cval.modulus());
	p_mod_q = p.polynomial();
	cval -= yashe.qdivt * p_mod_q;
	return *this;
}

Ciphertext Ciphertext::operator-(const Plaintext& p) const{
	Ciphertext c2(*this);
	return (c2 -= p);
}

Ciphertext& Ciphertext::operator*=(const Plaintext& p){
	if(aftermult)
		convert_self();

	fmpz_mod_polyxx p_mod_q(cval.modulus());
	p_mod_q = p.polynomial();
    cval = (cval.to<fmpz_polyxx>()*p_mod_q.to<fmpz_polyxx>())%(yashe.poly);
	return *this;
}

Ciphertext Ciphertext::operator*(const Plaintext& p) const{
	Ciphertext c2(*this);
	return (c2 *= p);
}

/* Convert */
void Ciphertext::convert_self(void) const {
    if(!aftermult)
    return;

    Ciphertext *that = const_cast<Ciphertext*>(this);
    that->cval = yashe.convert(cval);
    that->aftermult = false;
}

/* Get real noise */
unsigned Ciphertext::real_noise(void) const {
    if (!aftermult) return yashe.noise(*this);
    else return 0;
}


void Ciphertext::serialize(const std::string& filename) const{
	std::ofstream outFile(filename);
	outFile << "aftermult:" << (aftermult ? 1 : 0) << endl;
	outFile << "level:" << level << endl;
	outFile << "cval:" << cval.to<fmpz_polyxx>().to_string() << endl;
}

void Ciphertext::unserialize(const std::string& filename){
	std::string line;
	std::ifstream inFile(filename);
	while (inFile.good()){
		getline(inFile, line, ':');
		if ("aftermult" == line){
			inFile >> aftermult;
			getline(inFile, line); // ignores end of line
		}else if ("level" == line){
			inFile >> level;
			getline(inFile, line); // ignores end of line
		}else if ("cval" == line){
			getline(inFile, line);
			fmpz_polyxx tmp(line.c_str());
		    cval = fmpz_mod_polyxx(yashe.q);
		    cval = tmp;
		}else{
			getline(inFile, line); // ignores the line
		}
	}
}
