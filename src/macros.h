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

/* ---------------------------------------------------------------------------
 *      Macro to define the encoding, encryption and decryption
 *  of vectors and matrices in the YASHE class.
 *
 *  	It is used with name equals to encrypt and decrypt
 *  and from_type and to_type equals to
 *  vector<Plaintext>, vector<Ciphertext>, vector<vector<Plaintext> >, etc
 *----------------------------------------------------------------------------*/
#define define_function(name, from_type, to_type)\
	to_type name ( const from_type & _vec_or_matrix_ ) const{\
		unsigned int N = _vec_or_matrix_.size();\
		to_type dest;\
		for (unsigned int i =0; i < N; i++)\
			dest.push_back( name (_vec_or_matrix_[i]));\
		return dest;\
	}

