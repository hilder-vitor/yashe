#include "yashe_matrix_tests.h"
#include <random>
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

void EXPECT_MATRIX_DOUBLE_EQ(DoubleMatrix A, DoubleMatrix B, string message = ""){
	EXPECT_EQ(A.size(), B.size()) << "Number of lines is different. " << message;
	EXPECT_EQ(A[0].size(), B[0].size()) << "Number of columns is different. " << message;
	unsigned int N = A.size();
	unsigned int M = A[0].size();
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < M; j++)
			EXPECT_FLOAT_EQ(A[i][j], B[i][j]) << message;
	}
}

TEST (YASHE_encoding_decoding_matrix, random_and_fixed_values) {
	unsigned int N = 3;
	DoubleMatrix A(N);
	A[0] = {0.0, 1.2, -1, -1.232, -19221.28314, 827182.17261999997, -1.3864014009099695e-05, 1.0/(-72129.182742000005)};
	A[1] = {-92.13, 1, -1, -32, 1.2314, -812.1999, 1.3864014009099695e-05, -1.0/(-72129.182742000005)};
	A[2] = {0.0, 0, 0, 0, 0, 0.0, 0, .0};
	RealNumberPlaintextMatrix pA = yashe->encode(A);
	EXPECT_MATRIX_DOUBLE_EQ(A, yashe->decode(pA));

	A = createRandomFloatMatrix(10, 4);
	pA = yashe->encode(A);
	EXPECT_MATRIX_DOUBLE_EQ(A, yashe->decode(pA));
}

TEST (YASHE_encoding_encrypting_decrypting_decoding_matrix, random_and_fixed_values) {
	unsigned int N = 3;
	DoubleMatrix A(N);
	A[0] = {0.0, 1.2, -1, -1.232, -19221.28314, 827182.17261999997, -1.3864014009099695e-05, 1.0/(-72129.182742000005)};
	A[1] = {-92.13, 1, -1, -32, 1.2314, -812.1999, 1.3864014009099695e-05, -1.0/(-72129.182742000005)};
	A[2] = {0.0, 0, 0, 0, 0, 0.0, 0, .0};
	RealNumberPlaintextMatrix pA = yashe->encode(A);
	RealNumberCiphertextMatrix cA = yashe->encrypt(pA);
	EXPECT_MATRIX_DOUBLE_EQ(A, yashe->decode(yashe->decrypt(cA)));

	A = createRandomFloatMatrix(10, 4);
	pA = yashe->encode(A);
	cA = yashe->encrypt(pA);
	EXPECT_MATRIX_DOUBLE_EQ(A, yashe->decode(yashe->decrypt(cA)));
}

TEST (YASHE_addition_matrix, matrix_and_scalar) {
	unsigned int N = 2;
	DoubleMatrix A(N);
	A[0] = {0.0, -1.22, 92.381, 82.177, -231.3, 9.1};
	A[1] = {-91.23, 8.1, 0, 1, -2, 8};

	DoubleMatrix Aplus9Point7(N);
	Aplus9Point7[0] = {9.7, 8.48,   102.081, 91.877, -221.6, 18.8};
	Aplus9Point7[1] = {-81.53, 17.8, 9.7,    10.7,   7.7,   17.7};

	DoubleMatrix Aplus28133Point1897182(N);
	Aplus28133Point1897182[0] = {28133.1897182, 28131.9697182, 28225.5707182, 28215.3667182, 27901.8897182, 28142.2897182};
	Aplus28133Point1897182[1] = {28041.9597182, 28141.2897182, 28133.1897182, 28134.1897182, 28131.1897182, 28141.1897182};
	RealNumberPlaintextMatrix pA = yashe->encode(A);
	RealNumberCiphertextMatrix cA = yashe->encrypt(pA);
	RealNumberCiphertext c = yashe->encrypt(yashe->encode(9.7));
	RealNumberCiphertextMatrix matSum = cA + c;
	EXPECT_MATRIX_DOUBLE_EQ(Aplus9Point7, yashe->decode(yashe->decrypt(matSum)), "A + 9.7");

	c = yashe->encrypt(yashe->encode(28133.1897182));
	matSum = cA + c;
	EXPECT_MATRIX_DOUBLE_EQ(Aplus28133Point1897182, yashe->decode(yashe->decrypt(matSum)));
}


TEST (YASHE_addition_matrix, matrix_and_matrix) {
	unsigned int N = 4;
	DoubleMatrix A(N);
	A[0] = {-211.5853,  -626.0846,  -299.7301};
	A[1] = {-702.9763, -7.5688, 403.3649};
	A[2] = {-205.1212, 339.4169, -614.8392};
	A[3] = {895.3354, 809.7368, -28.6340};

	DoubleMatrix B(N);
	B[0] = {211.5853, 0.0846, 299};
	B[1] = {0, 1, 3649};
	B[2] = {3819.1, 183, 49.219};
	B[3] = {809.7368, -24.6340, 895.3354};

	DoubleMatrix AplusB(N);
	AplusB[0] = {0.0000, -626, -.7301};
	AplusB[1] = {-702.9763, -6.5688, 4052.3649};
	AplusB[2] = {3613.978799999999865, 522.416899999999941,  -565.620199999999954};
	AplusB[3] = {1705.0722, 785.1028, 866.7014};

	DoubleMatrix AplusA(N);
	AplusA[0] = {-423.170599999999979, -1252.1692, -599.460199999999986};
	AplusA[1] = {-1405.9526, -15.1376, 806.72979999999995};
	AplusA[2] = {-410.242399999999975,  678.833799999999997, -1229.6784};
	AplusA[3] = {1790.6708, 1619.4736, -57.2680};
		
	RealNumberCiphertextMatrix cA = yashe->encrypt(yashe->encode(A));
	RealNumberCiphertextMatrix cB = yashe->encrypt(yashe->encode(B));

	RealNumberCiphertextMatrix matSum = cA + cB;
	EXPECT_MATRIX_DOUBLE_EQ(AplusB, yashe->decode(yashe->decrypt(matSum)), "A + B");

	matSum = cA + cA;
	EXPECT_MATRIX_DOUBLE_EQ(AplusA, yashe->decode(yashe->decrypt(matSum)), "A + A");
}

TEST (YASHE_subtraction_matrix, matrix_and_scalar) {
	unsigned int N = 3;
	DoubleMatrix A(N);
	A[0] = {0.0, -22.1, 82.177, -231.3, 9.1};
	A[1] = {-19.32, 8.1, 0, 1, 8};
	A[2] = {-19.32, 0.0, 0, 1, 81923.123};

	DoubleMatrix Aminus9Point7(N);
	Aminus9Point7[0] = {-9.7, -31.8, 72.477, -241, -0.6};
	Aminus9Point7[1] = {-29.02, -1.6, -9.7, -8.7, -1.7};
	Aminus9Point7[2] = {-29.02, -9.7, -9.7, -8.7, 81913.423};

	DoubleMatrix Aminus847193Point1938591(N);
	Aminus847193Point1938591[0] = {-847193.1938591, -847215.2938591, -847111.0168591, -847424.4938591, -847184.0938591};
	Aminus847193Point1938591[1]= {-847212.5138590999,-847185.0938591, -847193.1938591, -847192.1938591, -847185.1938591};
	Aminus847193Point1938591[2]= {-847212.5138590999, -847193.1938591, -847193.1938591, -847192.1938591,-765270.0708591};
	RealNumberPlaintextMatrix pA = yashe->encode(A);
	RealNumberCiphertextMatrix cA = yashe->encrypt(pA);
	RealNumberCiphertext c = yashe->encrypt(yashe->encode(9.7));
	RealNumberCiphertextMatrix matSub = cA - c;
	EXPECT_MATRIX_DOUBLE_EQ(Aminus9Point7, yashe->decode(yashe->decrypt(matSub)), "A - 9.7");

	c = yashe->encrypt(yashe->encode(847193.1938591));
	matSub = cA - c;
	EXPECT_MATRIX_DOUBLE_EQ(Aminus847193Point1938591, yashe->decode(yashe->decrypt(matSub)));
}

TEST (YASHE_multiplication_matrix, matrix_and_scalar) {
	unsigned int N = 3;
	DoubleMatrix A(N);
	A[0]  = {0.0, 1.292, -19284.12, 9.1023, -8193.29481, 1};
	A[1]  = {0.1, 292, -19284.12, 8471, 1, 92491.1938};
	A[2]  = {2891.3, -19284.12, 47183,  17.48, 11, 3019.9421};

	RealNumberCiphertextMatrix cA = yashe->encrypt(yashe->encode(A));

	DoubleMatrix Atimes29point12(N);
	Atimes29point12[0] = {0.0, 37.62304, -561553.5744, 265.058975999999973, -238588.7448672, 29.12};
	Atimes29point12[1] = {2.912, 8503.04, -561553.5744, 246675.52, 29.12, 2693343.563455999828875};
	Atimes29point12[2] = {84194.656, -561553.5744, 1373968.959999999962747, 509.0176, 320.319999999999993, 87940.713952000005520};

	DoubleMatrix AtimesPoint293(N);
	AtimesPoint293[0] = {0, 0.378556, -5650.247159999998985, 2.6669739, -2400.635379329999523,0.293};
	AtimesPoint293[1] = {0.0293, 85.555999999999997, -5650.247159999998985, 2482.002999999999702, 0.293, 27099.919783399996959};
	AtimesPoint293[2] = {847.150899999999979, -5650.247159999998985, 13824.618999999998778, 5.12164,  3.223, 884.843035299999997};


	RealNumberCiphertext c = yashe->encrypt(yashe->encode(29.12));
	RealNumberCiphertextMatrix cProd = cA * c;
	EXPECT_MATRIX_DOUBLE_EQ(Atimes29point12, yashe->decode(yashe->decrypt(cProd)), "A * 29.12");

	c = yashe->encrypt(yashe->encode(.293));
	cProd = cA * c;
	EXPECT_MATRIX_DOUBLE_EQ(AtimesPoint293, yashe->decode(yashe->decrypt(cProd)), "A * 0.293");
	cProd = c * cA;
	EXPECT_MATRIX_DOUBLE_EQ(AtimesPoint293, yashe->decode(yashe->decrypt(cProd)), "A * 0.293");
}

TEST (YASHE_multiplication_matrix, matrix_and_matrix) {
	unsigned int N = 4;
	unsigned int P = 3;
	DoubleMatrix A(N);
	A[0] = {-211.5853, -626.0846,  -299.7301};
	A[1] = {-702.9763, -7.5688, 403.3649};
	A[2] = {-205.1212, 339.4169, -614.8392};
	A[3] = {895.3354, 809.7368, -28.6340};

	DoubleMatrix B(P);
	B[0] = {211.5853, -0.0846, 299, 0.21};
	B[1] = {0, 1, 3649, 0.0};
	B[2] = {3819.1, 183, 81.23, 49.219};

	DoubleMatrix AtimesB(N);
	AtimesB[0] = {-1189467.5640860898420214653015137, -55458.7927836200033198110759258, -2372193.7861229996196925640106201, -14796.8487049000013939803466201};
	AtimesB[1] = {1391751.4382616097573190927505493, 73867.6796949800045695155858994, -205043.1340730000229086726903915, 19705.5919900999979290645569563};
	AtimesB[2] = {-2391533.0193583602085709571838379, -112158.8034464800002751871943474, 1127257.6410840000025928020477295, -30304.8460368000014568679034710};
	AtimesB[3] = {80083.6998096199968131259083748, -4506.0305748399996446096338332, 3220108.9279800001531839370727539, -1221.3164120000001275911927223};

	RealNumberCiphertextMatrix cA = yashe->encrypt(yashe->encode(A));
	RealNumberCiphertextMatrix cB = yashe->encrypt(yashe->encode(B));

	RealNumberCiphertextMatrix cProd = cA * cB;
	EXPECT_MATRIX_DOUBLE_EQ(AtimesB, yashe->decode(yashe->decrypt(cProd)), "A * B");
}

TEST (YASHE_multiplication_matrix, matrix_and_vector) {
	unsigned int N = 4;
	unsigned int M = 3;
	DoubleMatrix A(N);
	A[0] = {-211.5853, -626.0846,  -299.7301};
	A[1] = {-702.9763, -7.5688, 403.3649};
	A[2] = {-205.1212, 339.4169, -614.8392};
	A[3] = {895.3354, 809.7368, -28.6340};

	DoubleVector v = {-92.1, 0.0192, -0.00682};

	DoubleVector AtimesV = {19477.0294649619972915388643742, 64741.2209604219970060512423515, 18902.3725278239980980288237333, -82444.6481095599883701652288437};

	RealNumberCiphertextMatrix cA = yashe->encrypt(yashe->encode(A));
	RealNumberCiphertextVector cv = yashe->encrypt(yashe->encode(v));
	RealNumberCiphertextVector cProd = cA * cv;

	DoubleVector calculatedAtimesV = yashe->decode(yashe->decrypt(cProd));

	for (unsigned int i = 0; i < M; i++){
		EXPECT_FLOAT_EQ(AtimesV[i], calculatedAtimesV[i]) << "A*v [" << i << "]";
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
