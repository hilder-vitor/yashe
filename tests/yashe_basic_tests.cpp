#include "yashe_basic_tests.h"
#include <random>
#include<vector>
#include <stdio.h>

#define EPSILON 0.00001

using namespace std;

Yashe* yashe;

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

TEST (YASHE_encoding_decoding, big_double_value) {
	double a = 1.7e+308 + .0959798;
	RealNumberPlaintext plain = yashe->encode(a, 64);
	double value = yashe->decode(plain);
	EXPECT_FLOAT_EQ(a, value);
}

TEST (YASHE_encoding_decoding, fixed_values) {
	vector<double> vec = {0.0, 1.2, -1, -1.232, -19221.28314, 827182.17261999997, -1.3864014009099695e-05, 1.0/(-72129.182742000005)};
	for (unsigned int i = 0; i < vec.size(); i++){
		RealNumberPlaintext plain = yashe->encode(vec[i], 64);
		double value = yashe->decode(plain);
		EXPECT_FLOAT_EQ(vec[i], value);
	}
}


TEST (YASHE_encoding_decoding, random_values) {
	unsigned int N = 25;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 0; i < N; i++){
		RealNumberPlaintext plain = yashe->encode(vec[i], 32);
		double value = yashe->decode(plain);
		EXPECT_FLOAT_EQ(vec[i], value);
	}
}

TEST (YASHE_encrypt_decrypt_with_encoding, fixed_values) {
	vector<double> vecDouble = {0.0, 0.5, 1, -2.5, 123, 123.123, 819231, 819231.2938192};
	vector<unsigned int> vecPrecisions = {16, 32, 64, 70, 50};
	for (unsigned int i = 0; i < vecDouble.size(); i++){
		for (unsigned int j = 0; j < vecPrecisions.size(); j++){
			RealNumberPlaintext plain = yashe->encode(vecDouble[i], vecPrecisions[j]);
			RealNumberCiphertext c = yashe->encrypt(plain);
			double value = yashe->decode(yashe->decrypt(c));
			EXPECT_FLOAT_EQ(vecDouble[i], value) << "(precision = " << vecPrecisions[j] << ")" << endl;
		}
	}
}

TEST (YASHE_encrypt_decrypt_with_encoding, random_values) {
	unsigned int N = 5;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 0; i < N; i++){
		RealNumberPlaintext plain = yashe->encode(vec[i], 32);
		RealNumberCiphertext c = yashe->encrypt(plain);
		double value = yashe->decode(yashe->decrypt(c));
		EXPECT_FLOAT_EQ(vec[i], value);
	}
}

TEST (YASHE_recrypt, fixed_values) {
	vector<double> vecDouble = {0.0, -0.5, 1, -2.5, 123, 123.123, 819231, 819231.2938192, 1.7e+308, 1.6e+308 + .1931235979};
	unsigned int N = vecDouble.size();
	for (unsigned int i = 0; i < N; i++){
		RealNumberCiphertext c = yashe->encrypt(yashe->encode(vecDouble[i]));
		yashe->recrypt(c);
		double value = yashe->decode(yashe->decrypt(c));
		EXPECT_FLOAT_EQ(vecDouble[i], value);
	}
}

void testSingleAddition(double a, unsigned int precA, double b, unsigned int precB){
	RealNumberPlaintext m1 = yashe->encode(a, precA);
	RealNumberPlaintext m2 = yashe->encode(b, precB);
	RealNumberCiphertext c1 = yashe->encrypt(m1);
	RealNumberCiphertext c2 = yashe->encrypt(m2);
	RealNumberCiphertext cSum = c1 + c2;
	double calculated = yashe->decode(yashe->decrypt(cSum));
	double expected = a + b;
	EXPECT_FLOAT_EQ(expected, calculated);
}

TEST (YASHE_additions, fixed_values) {
	vector<double> vecDouble = {0.0, -0.5, 0.5, 1, 1283, -9283.182742, 5, 827182.17262};
	vector<unsigned int> vecPrecisions = {16, 32, 64, 70};
	for (unsigned int i = 1; i < vecDouble.size(); i++){
		for (unsigned int j = 0; j < vecPrecisions.size(); j++){
			testSingleAddition(vecDouble[i-1], vecPrecisions[j], vecDouble[i], vecPrecisions[j]);
		}
	}
}


TEST (YASHE_additions, random_values) {
	unsigned int N = 5;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 1; i < N; i++){
		testSingleAddition(vec[i-1], 32, vec[i], 32);
	}
}

TEST (YASHE_additions, random_values_different_exp_shifts) {
	unsigned int N = 3;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 1; i < N; i++){
		testSingleAddition(vec[i-1], 128, vec[i], 64);
		testSingleAddition(vec[i-1], 128, vec[i], 150);
	}
}

TEST (YASHE_additions, plaintext_ciphertext) {
	unsigned int N = 3;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++){
			RealNumberCiphertext c = yashe->encrypt(yashe->encode(vec[i])) + vec[j];
			double calculated = yashe->decode(yashe->decrypt(c));
			EXPECT_FLOAT_EQ(vec[i] + vec[j], calculated);
		}
	}
}

void testSingleSubtraction(double a, unsigned int precA, double b, unsigned int precB){
	RealNumberPlaintext m1 = yashe->encode(a, precA);
	RealNumberPlaintext m2 = yashe->encode(b, precB);
	RealNumberCiphertext c1 = yashe->encrypt(m1);
	RealNumberCiphertext c2 = yashe->encrypt(m2);
	RealNumberCiphertext cSub = c1 - c2;
	double calculated = yashe->decode(yashe->decrypt(cSub));
	double expected = a - b;
	EXPECT_FLOAT_EQ(expected, calculated);
}

TEST (YASHE_subtractions, fixed_values) {
	vector<double> vecDouble = {0.0, -0.5, 0.5, 1, 11283, -72129.182742, 5, 827182.17262};
	vector<unsigned int> vecPrecisions = {16, 32, 64, 70};
	for (unsigned int i = 1; i < vecDouble.size(); i++){
		for (unsigned int j = 0; j < vecPrecisions.size(); j++){
			testSingleSubtraction(vecDouble[i-1], vecPrecisions[j], vecDouble[i], vecPrecisions[j]);
		}
	}
}

TEST (YASHE_subtractions, fixed_values_different_precisions) {
	testSingleSubtraction(0.0, 32, 1.5, 16);
	testSingleSubtraction(18273, 64, 1928817.29381, 63);
	testSingleSubtraction(18.273, 128, 91.1, 32);
	testSingleSubtraction(0.0, 32, 0.0, 0);
}



TEST (YASHE_subtractions, random_values) {
	unsigned int N = 10;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 1; i < N; i++){
		testSingleSubtraction(vec[i-1], 32, vec[i], 32);
	}
}

TEST (YASHE_subtractions, plaintext_ciphertext) {
	unsigned int N = 3;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++){
			RealNumberCiphertext c = yashe->encrypt(yashe->encode(vec[i])) - vec[j];
			double calculated = yashe->decode(yashe->decrypt(c));
			EXPECT_FLOAT_EQ(vec[i] - vec[j], calculated);
		}
	}
}


TEST (YASHE_products, random_values) {
	unsigned int N = 10;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 1; i < N; i++){
		RealNumberPlaintext m1 = yashe->encode(vec[i-1], 32);
		RealNumberPlaintext m2 = yashe->encode(vec[i], 32);
		RealNumberCiphertext c1 = yashe->encrypt(m1);
		RealNumberCiphertext c2 = yashe->encrypt(m2);
		RealNumberCiphertext cProd = c1 * c2;
		double calculated = yashe->decode(yashe->decrypt(cProd));
		double expected = vec[i-1] * vec[i];
		EXPECT_FLOAT_EQ(expected, calculated);
	}
}

TEST (YASHE_products, plaintext_ciphertext) {
	unsigned int N = 3;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++){
			RealNumberCiphertext c = yashe->encrypt(yashe->encode(vec[i])) * vec[j];
			double calculated = yashe->decode(yashe->decrypt(c));
			EXPECT_FLOAT_EQ(vec[i] * vec[j], calculated);
		}
	}
}

void testSingleDivision(double a, unsigned int precA, double N){
	RealNumberPlaintext m1 = yashe->encode(a, precA);
	RealNumberCiphertext c1 = yashe->encrypt(m1);
	RealNumberCiphertext cDiv = c1 / N;
	double calculated = yashe->decode(yashe->decrypt(cDiv));
	double expected = a / N;
	EXPECT_FLOAT_EQ(expected, calculated) << "a / N = " << a << " / " << N << endl;
}

TEST (YASHE_divisions, fixed_values) {
	testSingleDivision(0.0, 64, -0.5);
	vector<double> vecDouble = {0.0, -0.5, 0.5, 1, 2, 4, 8, 16, 32, -72129.182742, 5, 827182.17262};
	vector<unsigned int> vecPrecisions = {16, 32, 64, 70};
	for (unsigned int i = 1; i < vecDouble.size(); i++){
		for (unsigned int j = 0; j < vecPrecisions.size(); j++){
			testSingleDivision(vecDouble[i-1], vecPrecisions[j], vecDouble[i]);
		}
	}
}


TEST (YASHE_divisions, random_values) {
	unsigned int N = 10;
	vector<double> vec = createRandomFloatVector(N);
	for (unsigned int i = 1; i < N; i++){
		testSingleDivision(vec[i-1] + 0.1, 32, vec[i] + 0.1);
	}
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
