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

#include "Yashe.h"
#include "Sampler.h"
#include <iostream>
#include "flint/arith.h"
#include "timing.h"
#include <fstream>

/* Static values */
fmpzxx W((fmpzxx(1) << WORDLENGTH));
fmpzxx MASKING((fmpzxx(1) << WORDLENGTH)-fmpzxx(1));

/* Print Key */
std::ostream& operator<<(std::ostream& os, const Yashe& k) {
    os << "<Yashe with n=" << k.n << " cyclotomic_degree=" << k.ell << " q=" << k.q 
	<< " t=" << k.t << " sigma_key=" << k.sigmakey << " sigma_err=" << k.sigmaerr 
	<< " gamma_size=" << k.gamma.size()
	<< ">";
    return os;
}

/* Small useful functions */
bool isPowerOfTwo(int n)
{
    return (n) && !(n & (n - 1)); //this checks if the integer n is a power of two or not
}

void binaryGen(fmpz_mod_polyxx& f, unsigned degree)
{
    for (unsigned i=0; i<=degree; i++)
        f.set_coeff(i, fmpzxx((rand()%3)-1));
}

unsigned noise_from_poly(const fmpz_mod_polyxx& cval, const fmpzxx &q, unsigned ell)
{
    unsigned bitnoise = 0;
    fmpzxx coeff;
    for (unsigned i=0; i<ell; i++)
    {
        coeff = (cval.get_coeff(i).to<fmpzxx>());
        if (2*coeff > q)
            coeff = coeff - q;
        if (coeff.sizeinbase(2)>bitnoise)
            bitnoise = coeff.sizeinbase(2);
    }
    return bitnoise;
}

void Yashe::serialize(const std::string& filename){
	std::ofstream outFile(filename);
	outFile << "n:" << n << endl;
	outFile << "ell:" << ell << endl;
	outFile << "sigmakey:" << sigmakey << endl;
	outFile << "sigmaerr:" << sigmaerr << endl;
	outFile << "q:" << q << endl;
	outFile << "t:" << t << endl;

	outFile << "poly:" << poly.to_string() << endl;
	outFile << "h:" << h->to<fmpz_polyxx>().to_string() << endl;
	outFile << "f:" << f->to<fmpz_polyxx>().to_string() << endl;
	unsigned int _N = gamma.size();
	for (unsigned int i = 0; i < _N; i++)
		outFile << "gamma:" << gamma[i].to<fmpz_polyxx>().to_string() << endl;
}


Yashe::Yashe(const Yashe& y) :	
	n(y.n), ell(y.phi->degree()), sigmakey(y.sigmakey), sigmaerr(y.sigmaerr), logwq(y.q.sizeinbase(2)/WORDLENGTH+1),
	q(y.q), t(y.t),  qdivt(q/t), tmodq(y.tmodq), qdiv2t(y.q/(2*(y.t))), gamma(y.gamma) 
{
	phi = new fmpz_mod_polyxx(*(y.phi));
	h = new fmpz_mod_polyxx(*(y.h));
	f = new fmpz_mod_polyxx(*(y.f));
	sampler = new Sampler(sigmaerr*0.4, 1., &random); // 1/sqrt(2*pi) ~ 0.4
	poly = y.poly;
}

void Yashe::unserialize(const std::string& filename){
	std::string line;
	std::ifstream inFile(filename);
	while (inFile.good()){
		getline(inFile, line, ':');
		if ("n" == line){
			inFile >> n;
			getline(inFile, line); // ignores end of line
		}else if ("t" == line){
			getline(inFile, line);
			t = fmpzxx(line.c_str());
		}else if ("q" == line){
			getline(inFile, line);
			q = fmpzxx(line.c_str());
		}else if ("sigmakey" == line){
			inFile >> sigmakey;
			getline(inFile, line);
		}else if ("sigmaerr" == line){
			inFile >> sigmaerr;
			getline(inFile, line);
		}else if ("poly" == line){
			getline(inFile, line);
			fmpz_polyxx tmp(line.c_str());
			poly = tmp;
		    phi = new fmpz_mod_polyxx(q);
		    *phi = poly;
		}else if ("h" == line){
			getline(inFile, line);
    		h = new fmpz_mod_polyxx(q);
			fmpz_polyxx tmp(line.c_str());
			*h = tmp;
		}else if ("f" == line){
			getline(inFile, line);
    		f = new fmpz_mod_polyxx(q);
			fmpz_polyxx tmp(line.c_str());
			*f = tmp;
		}else if ("gamma" == line){
			getline(inFile, line);
			fmpz_polyxx tmp(line.c_str());
    		fmpz_mod_polyxx gamma_i(q);
			gamma_i = tmp;
			gamma.push_back(gamma_i);
		}else{
			getline(inFile, line); // ignores the line
		}
	}
	logwq = q.sizeinbase(2)/WORDLENGTH+1;
	qdivt = q/t;
	qdiv2t = q/(2*t);

	ell = phi->degree();
	sampler = new Sampler(sigmaerr*0.4, 1., &random); // 1/sqrt(2*pi) ~ 0.4
}

/* Constructor */
Yashe::Yashe(const std::string& filename){
	unserialize(filename);
}

/* Constructor */
Yashe::Yashe(const struct YASHEParams& params){
    n = params.n;
    sigmakey = params.sigmakey;
    sigmaerr = params.sigmaerr;
    q = params.q;
    t = params.t;

    logwq = q.sizeinbase(2)/WORDLENGTH+1;

    assert ( logwq*WORDLENGTH >= q.sizeinbase(2) );

    fmpz_mod_polyxx one(q);
    one.set_coeff(0, 1);

    // Define polynomial modulus
    arith_cyclotomic_polynomial(poly._data().inner, n);
    phi = new fmpz_mod_polyxx(q);
    *phi = poly;
    ell = phi->degree();

    fmpz_mod_polyxx finv(q);
    qdivt = q/t;
    qdiv2t = q/(2*t);

    // Creating sk/pk

    f = new fmpz_mod_polyxx(q);
    h = new fmpz_mod_polyxx(q);

    sampler = new Sampler(sigmaerr*0.4, 1., &random); // 1/sqrt(2*pi) ~ 0.4

    if (sigmakey == 1)
    {
        // Sample g, f1 with coefficients in {-1,0,1} 
        fmpz_mod_polyxx g(q);
        fmpz_mod_polyxx f1(q);
        binaryGen(g, ell-1);
        do
        {
            binaryGen(f1, ell-1);
            *f = t*f1+one;
            finv = (*f).invmod(*phi);
        } while (((*f)*finv)%(*phi) != one);
        *h = (t*g*finv)%(*phi);
    }
    else
    {
        exit(0);
    }

    // Create evaluation key
    fmpz_mod_polyxx pe(q);
    fmpz_mod_polyxx ps(q);
    gamma.resize(logwq, fmpz_mod_polyxx(q));
	for (unsigned k=0; k<logwq; k++){
        gamma[k] = *f;
        for (unsigned j=0; j<k; j++)
            gamma[k] = gamma[k]*W;

        for (unsigned i=0; i<ell; i++)
        {
            long value;

            value = sampler->SamplerGaussian();
            if (value>=0)   pe.set_coeff(i, fmpzxx(value));
            else            pe.set_coeff(i, q-fmpzxx(-value));
            value = sampler->SamplerGaussian();
            if (value>=0)   ps.set_coeff(i, fmpzxx(value));
            else            ps.set_coeff(i, q-fmpzxx(-value));
        }

        gamma[k] += pe+((*h)*ps)%(*phi);
    }
}
/* Encrypt a Plaintext */
Ciphertext Yashe::encrypt(const Plaintext& m) {
    fmpz_mod_polyxx cval(q);
    fmpz_mod_polyxx ps(q);
    fmpz_mod_polyxx delta_times_message(q);
    fmpzxx coeff;
    for (unsigned i=0; i<ell; i++)
    {
        cval.set_coeff(i,q+fmpzxx(sampler->SamplerGaussian())) ;
        ps.set_coeff(i,q+fmpzxx(sampler->SamplerGaussian()));
    }
	delta_times_message = (qdivt*(m.polynomial().to<fmpz_polyxx>()));
    cval += (((*h)*ps) + delta_times_message) % (*phi);
    return Ciphertext(*this, cval);
}


RealNumberCiphertext Yashe::encrypt(const RealNumberPlaintext& m) {
	RealNumberCiphertext c(encrypt(static_cast<Plaintext>(m)), m.expoent_shift());
	return c;
}

/* Decrypt */
Plaintext Yashe::decrypt(const Ciphertext& c) {
    fmpzxx coeff, diff;
    fmpz_polyxx g;
    if (c.aftermult)
        g = t*((((*f)*(*f)*c.cval)%(*phi)).to<fmpz_polyxx>());
    else
        g = t*((((*f)*c.cval)%(*phi)).to<fmpz_polyxx>());

	Plaintext m(t, *phi);

	unsigned int max_degree = phi->degree();
	for (unsigned i=0; i < max_degree; i++){
		// getting [g.coefficient(i) * q/t] (where [.] means closest integer)
		ltupleref(coeff, diff) = fdiv_qr(g.get_coeff(i), q);
		if (2*diff > q){
			if (coeff == t-fmpzxx(1))
				m.set(i, 0);
			else{
				coeff = coeff+fmpzxx(1);
				m.set(i, coeff);
			}
		}else
			m.set(i, coeff);
	}
    return m;
}

RealNumberPlaintext Yashe::decrypt(const RealNumberCiphertext& c){
	RealNumberPlaintext p(decrypt(static_cast<Ciphertext>(c)), c.expoent_shift());
	return p;
}

/* Convert */
void WordDecomp( std::vector<fmpz_mod_polyxx> &P, const fmpz_mod_polyxx &x )
{
    fmpzxx c;
    unsigned i,j;
    for (i=0; i <= x.degree(); i++)
    {
        c = x.get_coeff(i);
        j=0;
        while ( c > 0 )
        {
            P[j].set_coeff(i, c&MASKING);
            c = (c>>WORDLENGTH);
            j++;
        }
    }
}

fmpz_mod_polyxx Yashe::convert(const fmpz_mod_polyxx& cval) {
    fmpz_mod_polyxx result(q);

    std::vector<fmpz_mod_polyxx> P(logwq, fmpz_mod_polyxx(q));
    WordDecomp(P, cval);
    
    result = (P[0]*gamma[0]);
    for (unsigned i=1; i<logwq; i++)
    {
        result = result + (P[i]*gamma[i]);
    }
    result = result%(*phi);

    return result;
}

/* Get real noise */
unsigned Yashe::noise(const Ciphertext& _c){
	Ciphertext c = _c;
	c.convert_self();
	unsigned bitnoise = 0;
	fmpzxx coeff;
	fmpz_mod_polyxx g(q);
	fmpz_mod_polyxx m_bar(q);
	Plaintext m = decrypt(c);
	m_bar = m.polynomial();
	g = ((*f)*c.get_cval() - qdivt*m_bar) % (*phi);
	
	for (unsigned i=0; i<ell; i++){
		coeff = (g.get_coeff(i).to<fmpzxx>());
		if (2*coeff > q)
			coeff = coeff - q;

		if (coeff.sizeinbase(2)>bitnoise)
			bitnoise = coeff.sizeinbase(2);
	}
	return bitnoise;	
}


RealNumberPlaintext Yashe::encode(const double& _message, int exp_shift){
	double message = _message;
	int signal = (message > 0 ? 1 : -1);
	message *= signal;
	RealNumberPlaintext plain(t, *phi, exp_shift);
	// using mpz_class because values stored in a double may be greater than an int supports
	mpz_class integer_part = message;
	mpf_class __decimal_part = message;
	__decimal_part -= integer_part;
	double decimal_part = __decimal_part.get_d();//message - integer_part;
	// encoding the decimal part
	double cmp = 0.5;
	for (int i = plain.expoent_shift() - 1; i >= 0; i--){
		if (decimal_part >= cmp){ 
			plain.set(i, signal);
			decimal_part -= cmp;
		}
		cmp = cmp / 2.0;
	}
	// encode the integer part
	for (unsigned int i = plain.expoent_shift(); integer_part != 0; i++){
		mpz_class bit = integer_part % 2;
		plain.set(i, signal * bit.get_si());
		integer_part = integer_part / 2;
	}
	return plain;
}

double Yashe::decode(const RealNumberPlaintext& value){
	return value.double_value();
}

SymmetricMatrix<RealNumberCiphertext> Yashe::encrypt(const SymmetricMatrix<RealNumberPlaintext>& plain){
	unsigned int N = plain.size();
	RealNumberCiphertext c = encrypt(plain.get(0, 0));
	SymmetricMatrix<RealNumberCiphertext> C(N, c);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			c = encrypt(plain.get(i, j));	
			C.set(c, i, j);
		}
	}
	return C;
}

SymmetricMatrix<RealNumberPlaintext> Yashe::decrypt(const SymmetricMatrix<RealNumberCiphertext>& c){
	unsigned int N = c.size();
	RealNumberPlaintext plain = encode(1.0, 0);
	SymmetricMatrix<RealNumberPlaintext> P(N, plain);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			plain = decrypt(c.get(i, j));
			P.set(plain, i, j);	
		}
	}
	return P;
}

SymmetricMatrix<RealNumberPlaintext> Yashe::encode(SymmetricMatrix<double> messages, int precision){
	unsigned int N = messages.size();
	RealNumberPlaintext plain = encode(1.0, 0);
	SymmetricMatrix<RealNumberPlaintext> P(N, plain);
	
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			plain = encode(messages.get(i, j), precision);
			P.set(plain, i, j);	
		}
	}
	return P;
}

SymmetricMatrix<double> Yashe::decode(SymmetricMatrix<RealNumberPlaintext> plain){
	unsigned int N = plain.size();
	SymmetricMatrix<double> messages(N, 0.0);
	double value;
	
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			value = decode(plain.get(i, j));
			messages.set(value, i, j);
		}
	}
	return messages;
}

