#include "yashe_vectorial_tests.h"
#include <random>
#include <stdio.h>

#define EPSILON 0.00001

using namespace std;

Yashe* yashe = NULL;

vector<vector<double> > createRandomFloatMatrix(unsigned int N, unsigned int M){
	static std::default_random_engine e{((unsigned int)time(NULL))};
	static std::uniform_real_distribution<float> d{0, 10};

	vector<vector<double> > A(N);
	for (unsigned int i = 0; i < N; i++){
		A[i] = vector<double>(M);
		for (unsigned int j = 0; j < M; j++){
			A[i][j] = d(e) - d(e);
		}
	}
	return A;
}

vector<double> createRandomFloatVector(unsigned int N){
	static std::default_random_engine e{((unsigned int)time(NULL))};
	static std::uniform_real_distribution<float> d{0, 10};
	vector<double> v(N);
	for (unsigned int i = 0; i < N; i++){
		v[i] = d(e) - d(e);
	}
	return v;
}

void EXPECT_VECTOR_DOUBLE_EQ(vector<double> a, vector<double> b, string message = ""){
	EXPECT_EQ(a.size(), b.size()) << "Vectors' size are different. " << message;
	unsigned int n = a.size();
	for (unsigned int i = 0; i < n; i++){
		EXPECT_FLOAT_EQ(a[i], b[i]) << message;
	}
}

TEST (YASHE_encoding_decoding_vector, random_and_fixed_values) {
	vector<double> vec = {0.0, 1.2, -1, -1.232, -19221.28314, 827182.17261999997, -1.3864014009099695e-05, 1.0/(-72129.182742000005)};
	RealNumberPlaintextVector pVec = yashe->encode(vec);
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(pVec));

	vec = createRandomFloatVector(10);
	pVec = yashe->encode(vec);
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(pVec));
}

TEST (YASHE_encode_encrypt_decrypt_vector, random_and_fixed_values) {
	vector<double> vec = {0.0, 1.2, -1, 827182.17261999997, -1.3864014009099695e-05, 1.0/(-72129.182742000005), 9.09598e+10, -1.27908e+11, -5.17567e+10, -3.54928e+11, 6.97282e+10};
	RealNumberPlaintextVector pVec = yashe->encode(vec);
	RealNumberCiphertextVector cVec = yashe->encrypt(pVec);
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(cVec)));

	vec = createRandomFloatVector(10);
	pVec = yashe->encode(vec);
	cVec = yashe->encrypt(pVec);
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(cVec)));
}

TEST (YASHE_recrypt_vector, fixed_values) {
	vector<double> u = {0.0, 1.2, -1, -1.3864014009099695e-05, 1.0/(-72129.182742000005)};
	vector<double> v = {9.09598e+10, -1.27908e+11, -5.17567e+10, -3.54928e+11, 6.97282e+10};
	RealNumberCiphertextVector c_u = yashe->encrypt(yashe->encode(u));
	yashe->recrypt(c_u);
	EXPECT_VECTOR_DOUBLE_EQ(u, yashe->decode(yashe->decrypt(c_u)));

	RealNumberCiphertextVector c_v = yashe->encrypt(yashe->encode(v));
	yashe->recrypt(c_v);
	EXPECT_VECTOR_DOUBLE_EQ(v, yashe->decode(yashe->decrypt(c_v)));

	RealNumberCiphertextVector c_vec = c_u[3]*c_v[2]*(c_u + c_v);
	yashe->recrypt(c_vec);
	EXPECT_VECTOR_DOUBLE_EQ(u[3]*v[2]*(u + v), yashe->decode(yashe->decrypt(c_vec)));
}

TEST (YASHE_addtion_vector, vector_and_scalar) {
	vector<double> vec = {0.0, -1.22, 92.381, 82.177, -231.3, 9.1};
	vector<double> vecPlusPoint98 = {0.98, -0.24, 93.361, 83.157, -230.32, 10.08};
	vector<double> vecPlus23Point182 = {23.182, 21.962, 115.563, 105.359, -208.118, 32.282};
	vector<double> vecPlus5Plus1Point2 = {6.2000, 4.9800, 98.5810, 88.3770, -225.1000, 15.3000};

	RealNumberPlaintextVector pVec = yashe->encode(vec);
	RealNumberCiphertextVector cVec = yashe->encrypt(pVec);
	RealNumberCiphertext c = yashe->encrypt(yashe->encode(0.98));
	RealNumberCiphertextVector vecSum = cVec + c;
	EXPECT_VECTOR_DOUBLE_EQ(vecPlusPoint98, yashe->decode(yashe->decrypt(vecSum)));

	c = yashe->encrypt(yashe->encode(23.182, 64));
	vecSum = cVec + c;
	EXPECT_VECTOR_DOUBLE_EQ(vecPlus23Point182, yashe->decode(yashe->decrypt(vecSum)));

	vecSum = c + cVec;
	EXPECT_VECTOR_DOUBLE_EQ(vecPlus23Point182, yashe->decode(yashe->decrypt(vecSum)));

	RealNumberCiphertext c5 = yashe->encrypt(yashe->encode(5, 64));
	RealNumberCiphertext c1point2 = yashe->encrypt(yashe->encode(1.2, 64));
	vecSum = cVec + c5;
	vecSum = vecSum + c1point2;
	EXPECT_VECTOR_DOUBLE_EQ(vecPlus5Plus1Point2, yashe->decode(yashe->decrypt(vecSum)));

	// add zero to the vector
	c = yashe->encrypt(yashe->encode(0, 64));
	vecSum = cVec + c;
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(vecSum)));
}


TEST (YASHE_addtion_vector, vector_and_vector) {
	vector<double> u = {0.0,   -1.22,     92.381, 82.177,    -231.3,     9.1};
	vector<double> v = {-2.2,   183.29,   83.12, -193,       1,          0.1928};
	vector<double> w = {8.2912, 918231.1, 91.2,  183123.124, 18371.2948, 11.23};
	vector<double> uPlusV = {-2.2, 182.07, 175.501, -110.823, -230.3, 9.2928};
	vector<double> uPlusW = {8.2912, 918229.88, 183.581, 183205.301, 18139.9948, 20.33};
	vector<double> wPlusV = {6.0912, 918414.39, 174.32,  182930.124, 18372.2948, 11.4228};
	vector<double> wPlusUPlusV = {6.0912, 918413.17, 266.701, 183012.301, 18140.9948, 20.5228};

	RealNumberPlaintextVector plainU = yashe->encode(u);
	RealNumberPlaintextVector plainV = yashe->encode(v);
	RealNumberPlaintextVector plainW = yashe->encode(w);
	RealNumberCiphertextVector cipU = yashe->encrypt(plainU);
	RealNumberCiphertextVector cipV = yashe->encrypt(plainV);
	RealNumberCiphertextVector cipW = yashe->encrypt(plainW);

	RealNumberCiphertextVector vecSum = cipU + cipV;
	EXPECT_VECTOR_DOUBLE_EQ(uPlusV, yashe->decode(yashe->decrypt(vecSum)));

	vecSum = cipV + cipU;
	EXPECT_VECTOR_DOUBLE_EQ(uPlusV, yashe->decode(yashe->decrypt(vecSum)));

	vecSum = cipW + cipU;
	EXPECT_VECTOR_DOUBLE_EQ(uPlusW, yashe->decode(yashe->decrypt(vecSum)));

	vecSum = cipW + cipV;
	EXPECT_VECTOR_DOUBLE_EQ(wPlusV, yashe->decode(yashe->decrypt(vecSum)));

	vecSum = cipW + cipU;
	vecSum = vecSum + cipV;
	EXPECT_VECTOR_DOUBLE_EQ(wPlusUPlusV, yashe->decode(yashe->decrypt(vecSum)));
}
TEST (YASHE_subtraction_vector, vector_and_scalar) {
	vector<double> vec = {0.0, -1.22, 92.381, 82.177, -231.3, 9.1};
	vector<double> vecMinusPoint98 = {-0.98, -2.2, 91.401, 81.197, -232.28, 8.12000};
	vector<double> vecMinus23Point182 = {-23.182, -24.402, 69.199, 58.995, -254.482, -14.082};
	vector<double> vecMinus6Minus3Point14 = {-9.14, -10.36, 83.241, 73.037, -240.44, -0.04};

	RealNumberPlaintextVector pVec = yashe->encode(vec);
	RealNumberCiphertextVector cVec = yashe->encrypt(pVec);
	RealNumberCiphertext c = yashe->encrypt(yashe->encode(0.98, 64));
	RealNumberCiphertextVector vecSum = cVec - c;
	EXPECT_VECTOR_DOUBLE_EQ(vecMinusPoint98, yashe->decode(yashe->decrypt(vecSum)));

	c = yashe->encrypt(yashe->encode(23.182, 64));
	vecSum = cVec - c;
	EXPECT_VECTOR_DOUBLE_EQ(vecMinus23Point182, yashe->decode(yashe->decrypt(vecSum)));

	RealNumberCiphertext c6 = yashe->encrypt(yashe->encode(6, 64));
	RealNumberCiphertext c3point14 = yashe->encrypt(yashe->encode(3.14, 64));
	vecSum = cVec - c6;
	vecSum = vecSum - c3point14;
	EXPECT_VECTOR_DOUBLE_EQ(vecMinus6Minus3Point14, yashe->decode(yashe->decrypt(vecSum)));

	// subtracting zero
	c = yashe->encrypt(yashe->encode(0, 64));
	vecSum = cVec - c;
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(vecSum)));
}

TEST (YASHE_subtraction_vector, vector_and_vector) {
	vector<double> u = {99.9123, 182.1023, 0.00, -1239.1, -432};
	vector<double> v = {-9.321, 82.1023, -0.928, -9.1, -432};
	vector<double> w = {1, -922.3, 184, -1, 9.217};

	vector<double> vMinusW = {-10.321, 1004.4023, -184.928, -8.1, -441.2170};
	vector<double> wMinusV = {10.321, -1004.4023, 184.928, 8.1, 441.2170};
	vector<double> uMinusV = {109.2333, 100, 0.928, -1230, 0};
	vector<double> uMinusW = {98.9123, 1104.4023, -184, -1238.1, -441.217};

	RealNumberCiphertextVector cipherU = yashe->encrypt(yashe->encode(u));
	RealNumberCiphertextVector cipherV = yashe->encrypt(yashe->encode(v));
	RealNumberCiphertextVector cipherW = yashe->encrypt(yashe->encode(w));

	RealNumberCiphertextVector vecSub = cipherV - cipherW;
	EXPECT_VECTOR_DOUBLE_EQ(vMinusW, yashe->decode(yashe->decrypt(vecSub)), "v - w");;

	vecSub = cipherW - cipherV;
	EXPECT_VECTOR_DOUBLE_EQ(wMinusV, yashe->decode(yashe->decrypt(vecSub)), "w - v");

	vecSub = cipherU - cipherV;
	EXPECT_VECTOR_DOUBLE_EQ(uMinusV, yashe->decode(yashe->decrypt(vecSub)), "u - v");

	vecSub = cipherU - cipherW;
	EXPECT_VECTOR_DOUBLE_EQ(uMinusW, yashe->decode(yashe->decrypt(vecSub)), "u - w");
}

TEST (YASHE_multiplication_vector, vector_and_scalar) {
	vector<double> vec = {8.2, 1.273, 0.0, -918.921, -1, 10248283.919882};
	vector<double> vecZero = {0, 0, 0.0, 0, .0, 0.0};
	vector<double> vecTimesMinusOne = {-8.2, -1.273, 0.0, 918.921, 1, -10248283.919882};
	vector<double> vecTimesOnePointThree = {10.66, 1.6549, 0.0, -1194.5973, -1.3, 13322769.0958465989679098129272461};

	RealNumberCiphertextVector cipherVec = yashe->encrypt(yashe->encode(vec));
	RealNumberCiphertext c = yashe->encrypt(yashe->encode(0.0, 64));
	RealNumberCiphertextVector vecProd = c * cipherVec;
	EXPECT_VECTOR_DOUBLE_EQ(vecZero, yashe->decode(yashe->decrypt(vecProd)), "0.0 * vec");

	c = yashe->encrypt(yashe->encode(-1.0, 64));
	vecProd = c * cipherVec;
	EXPECT_VECTOR_DOUBLE_EQ(vecTimesMinusOne, yashe->decode(yashe->decrypt(vecProd)), "-1.0 * vec");
	vecProd = cipherVec * c;
	EXPECT_VECTOR_DOUBLE_EQ(vecTimesMinusOne, yashe->decode(yashe->decrypt(vecProd)), "vec * -1.0");

	c = yashe->encrypt(yashe->encode(1.3, 64));
	vecProd = c * cipherVec;
	EXPECT_VECTOR_DOUBLE_EQ(vecTimesOnePointThree, yashe->decode(yashe->decrypt(vecProd)), "1.3 * vec");

	c = yashe->encrypt(yashe->encode(1.0, 64));
	vecProd = c * cipherVec;
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(vecProd)), "1.0 * vec");

	vec = {918271.23234};
	vecTimesMinusOne = {-918271.23234};
	c = yashe->encrypt(yashe->encode(-1.0, 64));
	cipherVec = yashe->encrypt(yashe->encode(vec));
	vecProd = c * cipherVec;
	EXPECT_VECTOR_DOUBLE_EQ(vecTimesMinusOne, yashe->decode(yashe->decrypt(vecProd)), "-1.0 * vec");
}

TEST (YASHE_multiplication_vector, inner_product) {
	vector<double> u = {0.3219, 1023.281, 0.00, -1.1239, -432};
	vector<double> v = {-918.3, 183.27182, -0.928, -.0001, -432};
	vector<double> w = {-4.3123, -3.2, 92.3481, 0.0, 1};
	vector<double> zero = {0.0, 0, 0, 0, 0};

	RealNumberCiphertextVector cipherU = yashe->encrypt(yashe->encode(u));
	RealNumberCiphertextVector cipherV = yashe->encrypt(yashe->encode(v));
	RealNumberCiphertextVector cipherW = yashe->encrypt(yashe->encode(w));
	RealNumberCiphertextVector cipherZero = yashe->encrypt(yashe->encode(zero));

	double uTimesV = 373866.9705838100053;
	double uTimesW = -3707.8873293700003;
	double vTimesW = 2855.8162291999988;
	double vTimesV = 1063488.31119012227281928;

	RealNumberCiphertext innerProd = cipherU * cipherV;
	EXPECT_FLOAT_EQ(uTimesV, yashe->decode(yashe->decrypt(innerProd)));
	
	innerProd = cipherU * cipherW;
	EXPECT_FLOAT_EQ(uTimesW, yashe->decode(yashe->decrypt(innerProd)));

	innerProd = cipherV * cipherW;
	EXPECT_FLOAT_EQ(vTimesW, yashe->decode(yashe->decrypt(innerProd)));
	innerProd = cipherW * cipherV;
	EXPECT_FLOAT_EQ(vTimesW, yashe->decode(yashe->decrypt(innerProd)));

	innerProd = cipherZero * cipherU;
	EXPECT_FLOAT_EQ(0.0, yashe->decode(yashe->decrypt(innerProd)));

	innerProd = cipherV * cipherV;
	EXPECT_FLOAT_EQ(vTimesV, yashe->decode(yashe->decrypt(innerProd)));

}
	
TEST (YASHE_division_vector, vector_and_scalar) {
	vector<double> vec = {0.0, -9.102, 83.2143, 23.1828, -19.2934};
	vector<double> vecOver2 = {0.0, -4.551, 41.60715, 11.59140, -9.6467};
	vector<double> vecOver1Point43 = {0.0, -6.3650349650, 58.19182,  16.21175,-13.49189};
	vector<double> vecOverPoint71 = {0.0, -12.81972, 117.20324, 32.65183,-27.17380};
	vector<double> vecOverPoint00001 = {0, -910200, 8321429.99999999906867742538, 2318280.0, -1929339.99999999976716935635};
	vector<double> vecOver918347182993 = {0, -0.00000000000991128428, 0.00000000009061311620, 0.00000000002524404760, -0.00000000002100883016};
	
	RealNumberCiphertextVector cipherVec = yashe->encrypt(yashe->encode(vec));
	RealNumberCiphertextVector vecDiv = cipherVec / 2.0;
	EXPECT_VECTOR_DOUBLE_EQ(vecOver2, yashe->decode(yashe->decrypt(vecDiv)), "vec / 2.0");

	vecDiv = cipherVec / 1.43;
	EXPECT_VECTOR_DOUBLE_EQ(vecOver1Point43, yashe->decode(yashe->decrypt(vecDiv)), "vec / 1.43");

	vecDiv = cipherVec / 0.71;
	EXPECT_VECTOR_DOUBLE_EQ(vecOverPoint71, yashe->decode(yashe->decrypt(vecDiv)), "vec / 0.71");

	vecDiv = cipherVec / 0.00001;
	EXPECT_VECTOR_DOUBLE_EQ(vecOverPoint00001, yashe->decode(yashe->decrypt(vecDiv)), "vec / 0.00001");

	vecDiv = cipherVec / 918347182993;
	EXPECT_VECTOR_DOUBLE_EQ(vecOver918347182993, yashe->decode(yashe->decrypt(vecDiv)), "vec / 918347182993");

	vec = {91.29, -19238, 0.12, 0.1232};
	cipherVec = yashe->encrypt(yashe->encode(vec));
	vecDiv = cipherVec / 1.0;
	EXPECT_VECTOR_DOUBLE_EQ(vec, yashe->decode(yashe->decrypt(vecDiv)), "vec / 1");
}

GTEST_API_ int main(int argc, char **argv) {
	string fileName("../keys/L4_t2to32_yashe.keys");
	std::ifstream inFile(fileName);
	if (inFile.good()){
		cout << "Loading yashe." << endl;
		yashe = new Yashe(fileName);
	}else{
		cout << "Generating keys." << endl;
		yashe = new Yashe(params_L4_t2to32);
		yashe->serialize(fileName);
	}
	cout << "Starting tests." << endl;
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
