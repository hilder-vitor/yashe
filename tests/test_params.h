#ifndef __YASHE___TEST_PARAMS__
#define __YASHE___TEST_PARAMS__


// -------------------------------------------------------------------------------
// --------- PARAMS FOR L = 4 WITH t = 2**32 (4294967296) --------------------------------
// -------------------------------------------------------------------------------
const struct YASHEParams params_L4_t2to32 = { 
	16384, // n = 2**14  (totient(n) = 2**13 is the degree of the cyclotomic polynomial)
	1,  // Key distribution's standard deviation
	8,  // Error distribution's standard deviation
    fmpzxx("2839213766779714416208296124562517712318911565184836172974571090549372219192960637992933791850638927971728600024477257552869537611963"), // q  (440-bit prime)   Ciphertext: ring R/qR 
    fmpzxx("4294967296") // t: Cleartext:  ring R/(moduli)R  (for instance, if t = 7, R/7R) 
};

// -------------------------------------------------------------------------------
// --------- PARAMS FOR L = 5 WITH t = 2**32 (4294967296) --------------------------------
// -------------------------------------------------------------------------------
const struct YASHEParams params_L5_t2to32 = { 
	16384, // n = 2**14  (totient(n) = 2**13 is the degree of the cyclotomic polynomial)
	1,  // Key distribution's standard deviation
	8,  // Error distribution's standard deviation
    fmpzxx("3773962424821541352241554580988268890916921220416440428376206300245624162392148852086126725177658767541468375030763844899770584629924792632561434251432696043649395327187"), // q  (560-bit prime)   Ciphertext: ring R/qR 
    fmpzxx("4294967296") // t: Cleartext:  ring R/(moduli)R  (for instance, if t = 7, R/7R) 
};


//-------------------------------------------------------------------------------
// --------- PARAMS FOR L = 5 WITH t = 2**64 (18446744073709551616) -------------
// ------------------------------------------------------------------------------
const struct YASHEParams params_L5_t2to64 = { 
	24593, // n 
	1,  // Key distribution's standard deviation
	8,  // Error distribution's standard deviation
	fmpzxx("10715086071862673209484250490600018105614048117055336074437503883703510511249361224931983788156958581275946729175531468251871452856923140435984577574698574803934567774824230985421074605062371141877954182153046474983581941267398767559165543946077062914571196477686542167660429831652624386837205668069673"), // q (1000-bit prime) Ciphertext: ring R/qR 
	fmpzxx("18446744073709551616") // t = 2**64
};

#endif
