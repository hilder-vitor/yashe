#include "yashe_symmetric_matrix_tests.h"
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

SymmetricMatrix<double> createRandomFloatSymmetricMatrix(unsigned int N){

	vector<vector<double> > A = createRandomFloatMatrix(N, N);
	SymmetricMatrix<double> S(N, 0.0);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			S.set(A[i][j], i, j);
		}
	}
	return S;
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

void EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(const SymmetricMatrix<double>& A, const SymmetricMatrix<double>& B, string message = ""){
	EXPECT_EQ(A.size(), B.size()) << "Dimentions are different. " << message;
	unsigned int N = A.size();
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++)
			EXPECT_FLOAT_EQ(A.get(i,j), B.get(i,j)) << message;
	}
}

TEST (YASHE_encoding_decoding, symmetric_matrices_N6){
	unsigned int N = 6;
	SymmetricMatrix<double> DS = createRandomFloatSymmetricMatrix(N);
	SymmetricMatrix<RealNumberPlaintext> PS = yashe->encode(DS);
	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(DS, yashe->decode(PS));
}

TEST (YASHE_encoding_decoding, symmetric_matrices_N169){
	unsigned int N = 169;
	SymmetricMatrix<double> DS = createRandomFloatSymmetricMatrix(N);
	DS.set(0.0, 0, 0);
	SymmetricMatrix<RealNumberPlaintext> PS = yashe->encode(DS);
	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(DS, yashe->decode(PS));
}

TEST (YASHE_encoding_encrypting_decrypting_decoding, symmetric_matrices_N7){
	unsigned int N = 7;
	SymmetricMatrix<double> DS = createRandomFloatSymmetricMatrix(N);
	SymmetricMatrix<RealNumberCiphertext> CS = yashe->encrypt(yashe->encode(DS));
	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(DS, yashe->decode(yashe->decrypt(CS)));
}

TEST (YASHE_outer_product, double_values_u_and_v){
	DoubleVector u = {0.9, 918.3, -19.3134, .264};
	DoubleVector v = {12, -3.918, 3.294, 12384.64};

	DoubleMatrix outer(4);
	outer[0] = {10.8, -3.5262, 2.9646, 11146.176};
	outer[1] = {11019.600, -3597.8994, 3024.8802, 11372814.912};
	outer[2] = {-231.76080, 75.66990, -63.61834, -239189.50618};
	outer[3] = {3.16800, -1.034352, 0.869616, 3269.54496};

	DoubleMatrix calculated = outerProduct(u, v);
	EXPECT_MATRIX_DOUBLE_EQ(outer, calculated, "outer product of u and v.");
}

TEST (YASHE_outer_product, double_values_u_and_u){
	DoubleVector u = {0.9, 918.3, -19.3134, .264};

	SymmetricMatrix<double> sym = outerProduct(u);
	SymmetricMatrix<double> expected(4, 0.0);
	expected.set(0.8100000, 0, 0);
	expected.set(826.47, 1, 0); expected.set(843274.8900000, 1, 1);
	expected.set(-17.3820600, 2, 0); expected.set(-17735.4952200, 2, 1); expected.set(373.0074196, 2, 2);
	expected.set(0.2376000, 3, 0); expected.set(242.4312000, 3, 1); expected.set(-5.0987376, 3, 2); expected.set(0.0696960, 3, 3);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, sym, "outer product of u and u." );
}

TEST (YASHE_outer_product, encrypted_values_u_and_v){
	DoubleVector u = {0.9, 918.3, -19.3134, .264};
	DoubleVector v = {12, -3.918, 3.294, 12384.64};

	DoubleMatrix outer(4);
	outer[0] = {10.8, -3.5262, 2.9646, 11146.176};
	outer[1] = {11019.600, -3597.8994, 3024.8802, 11372814.912};
	outer[2] = {-231.76080, 75.66990, -63.61834, -239189.50618};
	outer[3] = {3.16800, -1.034352, 0.869616, 3269.54496};

	RealNumberCiphertextVector encU = yashe->encrypt(yashe->encode(u));
	RealNumberCiphertextVector encV = yashe->encrypt(yashe->encode(v));

	RealNumberCiphertextMatrix calculated = outerProduct(encU, encV);
	EXPECT_MATRIX_DOUBLE_EQ(outer, yashe->decode(yashe->decrypt(calculated)), "outer product of u and v.");
}

TEST (YASHE_outer_product, encrypted_values_u_and_u){
	DoubleVector u = {0.9, 918.3, -19.3134, .264};

	SymmetricMatrix<double> expected(4, 0.0);
	expected.set(0.8100000, 0, 0);
	expected.set(826.47, 1, 0); expected.set(843274.8900000, 1, 1);
	expected.set(-17.3820600, 2, 0); expected.set(-17735.4952200, 2, 1); expected.set(373.0074196, 2, 2);
	expected.set(0.2376000, 3, 0); expected.set(242.4312000, 3, 1); expected.set(-5.0987376, 3, 2); expected.set(0.0696960, 3, 3);

	RealNumberCiphertextVector encU = yashe->encrypt(yashe->encode(u));
	SymmetricMatrix<RealNumberCiphertext> sym = outerProduct(encU);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, yashe->decode(yashe->decrypt(sym)), "outer product of u and u." );
}


TEST(YASHE_multiply_transpose, double_values_N5_P3){
	unsigned int N = 5;
	unsigned int P = 3;
	DoubleMatrix A(N);
	A[0] = {0.9, 9182.2, -91.4};
	A[1] = {9.183, -0.25, 0.183};
	A[2] = {1831.9, 2.2, 53.2};
	A[3] = {-0.8, 0.0, -4};
	A[4] = {23, 1, 28471.1};
	
	SymmetricMatrix<double> expected(P, 0.0);
	expected.set(3356472.3874890003, 0, 0);
	expected.set(12314.8642500, 1, 0);expected.set(84312802.7425000221, 1, 1);
	expected.set(752215.0004890, 2, 0);expected.set(-810664.9857500001, 2, 1);expected.set(810614735.4434889555, 2, 2);

	SymmetricMatrix<double> calculated = multiplyTransposeMatrixByMatrix(A);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, calculated, "A' * A" );
}

TEST(YASHE_multiply_transpose, double_values_N3_P6){
	unsigned int N = 3;
	unsigned int P = 6;
	DoubleMatrix A(N);
	A[0] = {0.0, 8.21, 0.9, 77.32, 9182.2, -91.4};
	A[1] = {9.183, 81.29, -0.25, 947.281, 0.0123, 0.183};
	A[2] = {0.192, 1831.9, 75.19, 2.2, 53.2, 184.8};
	
	SymmetricMatrix<double> E(P, 0.0);
	E.set(84.3643530000, 0, 0);
	E.set(1098.21087, 1, 0); E.set(3362533.0782000003, 1, 1);
   	E.set(12.140730, 2, 0); E.set(137727.62750, 2, 1); E.set(5654.4086, 2, 2);
	E.set(8699.3038230, 3, 0); E.set(81669.449690, 3, 1); E.set(-1.814250, 3, 2); E.set(903324.5153609999, 3, 3);
	E.set(10.3273509, 4, 0); E.set(172843.9418670, 4, 1); E.set(12264.084925, 4, 2); E.set(710096.3955563001, 4, 3); E.set(84315627.0801513046, 4, 4);
	E.set(37.162089, 5, 0); E.set(337799.6020700001, 5,1); E.set(13812.80625, 5,2); E.set(-6487.135577, 5, 3); E.set(-829421.7177491001, 5, 4); E.set(42505.033489, 5, 5);

	SymmetricMatrix<double> calculated = multiplyTransposeMatrixByMatrix(A);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(E, calculated, "A' * A" );
}

TEST(YASHE_symmetric_multiplication, double_values_N4){
	unsigned int N = 4;
	SymmetricMatrix<double> A(N, 0.0);
	A.set(92.182, 0, 0);
	A.set(-381.94, 1, 0); A.set(-82.4, 1, 1);
	A.set(49, 2, 0); A.set(-.02, 2, 1); A.set(0.18, 2, 2);
	A.set(419, 3, 0); A.set(92.12, 3, 1); A.set(83.273, 3, 2); A.set(64.1, 3, 3);

	SymmetricMatrix<double> B(N, 0.0);
	B.set(289, 0, 0);
	B.set(-31.4, 1, 0); B.set(0.192, 1, 1);
	B.set(973.3, 2, 0); B.set(963.41, 2, 1); B.set(-0.873, 2, 2);
	B.set(419, 3, 0); B.set(0.0, 3, 1); B.set(8.84, 3, 2); B.set(64.1, 3, 3);

	SymmetricMatrix<double> AtimesB(N, 0.0);
	AtimesB.set(261886.214, 0, 0);
	AtimesB.set(-69214.486, 1, 0); AtimesB.set(11957.827, 1, 1);
	AtimesB.set(49228.2090, 2, 0); AtimesB.set(-1365.19004, 2, 1); AtimesB.set(48408.40798, 2, 2);
	AtimesB.set(226105.9429, 3, 0); AtimesB.set(67087.12797, 3, 1); AtimesB.set(497055.975871, 3, 2); AtimesB.set(180405.94332, 3, 3);

	SymmetricMatrix<double> calculated = A * B;

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(AtimesB, calculated, "A * B" );
}

TEST(YASHE_symmetric_multiplication, double_values_N4_square){
	unsigned int N = 4;
	SymmetricMatrix<double> A(N, 0.0);
	A.set(92.182, 0, 0);
	A.set(-381.94, 1, 0); A.set(-82.4, 1, 1);
	A.set(49, 2, 0); A.set(-.02, 2, 1); A.set(0.18, 2, 2);
	A.set(419, 3, 0); A.set(92.12, 3, 1); A.set(83.273, 3, 2); A.set(64.1, 3, 3);

	SymmetricMatrix<double> AtimesA(N, 0.0);
	AtimesA.set(332337.684724, 0, 0);
	AtimesA.set(34861.1629200000, 1, 0); AtimesA.set(161154.0184, 1, 1);
	AtimesA.set(39424.7638000000, 2, 0); AtimesA.set(-11042.3068400000, 2, 1); AtimesA.set(9335.4253290000, 2, 2);
	AtimesA.set(34378.22220, 3, 0); AtimesA.set(-161720.32146, 3, 1); AtimesA.set(25881.94604, 3, 2); AtimesA.set(195090.296929, 3, 3);

	SymmetricMatrix<double> calculated = A * A;

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(AtimesA, calculated, "A * A" );
}


TEST(YASHE_symmetric_addition, double_values_scalar_N3){
	unsigned int N = 3;
	SymmetricMatrix<double> A(N, 0.0);
	A.set(332337.684724, 0, 0);
	A.set(34861.162920, 1, 0); A.set(161154.0184, 1, 1);
	A.set(39424.7638000, 2, 0); A.set(-11042.3068400, 2, 1); A.set(9335.4253290000, 2, 2);

	double scalar = 736843.13931;

	SymmetricMatrix<double> E(N, 0.0);
	E.set(1069180.824034000, 0, 0);
	E.set(771704.30223, 1, 0); E.set(897997.1577099999, 1, 1);
	E.set(776267.90311, 2, 0); E.set(725800.83247, 2, 1); E.set(746178.564639, 2, 2);

	SymmetricMatrix<double> sum = A + scalar;

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(E, sum, "A + lambda" );
}

TEST(YASHE_symmetric_subtraction, double_values_scalar_N3){
	unsigned int N = 3;
	SymmetricMatrix<double> A(N, 0.0);
	A.set(37.6724, 0, 0);
	A.set(61.120, 1, 0); A.set(161154.0184, 1, 1);
	A.set(-4.7007, 2, 0); A.set(0, 2, 1); A.set(-95.453290000, 2, 2);

	double scalar = 73.13981831;

	SymmetricMatrix<double> E(N, 0.0);
	E.set(-35.467418309999992, 0, 0);
	E.set( -12.019818309999998, 1, 0); E.set(161080.878581689990824, 1, 1);
	E.set(-77.840518309999993, 2, 0); E.set(-73.139818309999995, 2, 1); E.set(-168.593108309999991, 2, 2);

	SymmetricMatrix<double> sub = A - scalar;

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(E, sub, "A - lambda" );
}

TEST(YASHE_symmetric_subtraction, double_values_matrix_N3){
	unsigned int N = 3;
	SymmetricMatrix<double> A(N, 0.0);
	A.set(37.6724, 0, 0);
	A.set(61.120, 1, 0); A.set(161154.0184, 1, 1);
	A.set(-4.7007, 2, 0); A.set(0, 2, 1); A.set(-95.453290000, 2, 2);

	SymmetricMatrix<double> B(N, 0.0);
	B.set(1145.7400, 0, 0);
	B.set(490.0100, 1, 0); B.set(4806.2900, 1, 1);
	B.set(-2197.6100, 2, 0); B.set(-2900.4100, 2, 1); B.set(8.8320, 2, 2);

	SymmetricMatrix<double> AminusB(N, 0.0);
	AminusB.set(-1108.06760, 0, 0);
	AminusB.set(-428.889999999999986, 1, 0); AminusB.set(156347.728399999992689, 1, 1);
	AminusB.set(2192.90930, 2, 0); AminusB.set(2900.409999999999854, 2, 1); AminusB.set(-104.28529, 2, 2);

	SymmetricMatrix<double> sub = A - B;

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(AminusB, sub, "A - B" );

}

TEST(YASHE_symmetric_covariance, double_values_matrix_N6_P3){
	unsigned int N = 6;
	unsigned int P = 3;
	DoubleMatrix X(N);
	X[0] = {29.1679, 0.0, -26.3484};
	X[1] = {16.5752, -68.6748, 8.5019};
	X[2] = {51.4491, -36.0149, 32.9470};
	X[3] = {-34.4962, 4.3080, -24.8359};
	X[4] = {0.928, 6.83580, 9.7333};
	X[5] = {60.9385, -47.7381, 18.8163};
		
	SymmetricMatrix<double> expected(P, 0.0);
	expected.set(6194.33269195, 0, 0);
	expected.set(-3226.79683060000, 1, 0); expected.set(5163.72740934, 1, 1);
	expected.set(2704.93346175000,2,0); expected.set(-2283.85849926200, 2, 1); expected.set(2861.0016390976, 2, 2);

	SymmetricMatrix<double> C = calculateCovariance(X);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, C);
}

TEST(YASHE_symmetric_covariance, encrypted_values_matrix_N6_P3){
	unsigned int N = 6;
	unsigned int P = 3;
	DoubleMatrix X(N);
	X[0] = {29.1679, 0.0, -26.3484};
	X[1] = {16.5752, -68.6748, 8.5019};
	X[2] = {51.4491, -36.0149, 32.9470};
	X[3] = {-34.4962, 4.3080, -24.8359};
	X[4] = {0.928, 6.83580, 9.7333};
	X[5] = {60.9385, -47.7381, 18.8163};
		
	SymmetricMatrix<double> expected(P, 0.0);
	expected.set(6194.33269195, 0, 0);
	expected.set(-3226.79683060000, 1, 0); expected.set(5163.72740934, 1, 1);
	expected.set(2704.93346175000,2,0); expected.set(-2283.85849926200, 2, 1); expected.set(2861.0016390976, 2, 2);

	RealNumberCiphertextMatrix cX = yashe->encrypt(yashe->encode(X));

	SymmetricMatrix<RealNumberCiphertext> C = calculateCovariance(cX);
	SymmetricMatrix<RealNumberPlaintext> S = yashe->decrypt(C);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, yashe->decode(yashe->decrypt(C)));
}


TEST(YASHE_symmetric_covariance, double_values_matrix_N8_P4){
	unsigned int N = 8;
	unsigned int P = 4;
	DoubleMatrix X(N);
	X[0] = {-58.5193, 37.3509, 84.6241, 15.5061};
	X[1] = {53.7900, 2.6724, 27.5779, 61.9787};
	X[2] = {-65.0143, 38.5376, -70.3401, 43.1870};
	X[3] = {35.8606, 47.0229,-35.2771, -76.6108};
	X[4] = {-5.9548, -22.6844, 72.8573, -43.5289};
	X[5] = {-52.5733, 11.0885,-13.6153,-38.4474};
	X[6] = {5.3933, -28.9780, 1.4551, 67.8734};
	X[7] = {-31.0055, 58.2462, 36.3037, 7.1579};

	SymmetricMatrix<double> expected(P, 0.0);
	expected.set(13914.90165964796, 0, 0);
	expected.set(-3201.03019600980,1, 0); expected.set(7455.47638349898, 1, 1);
	expected.set(500.51231543694, 2, 0); expected.set(-2683.44876950469, 2, 1); expected.set(19613.74615190041,2, 2);
	expected.set(-167.61952464898, 3, 0); expected.set(-2833.22457586510, 3, 1); expected.set(-73.82082653347,3,2); expected.set(19678.44629570449, 3,3);

	SymmetricMatrix<double> C = calculateCovariance(X);

	EXPECT_SYMMETRIC_MATRIX_DOUBLE_EQ(expected, C);
}

GTEST_API_ int main(int argc, char **argv) {
	string baseFileName("../keys/1455824090_yashe.keys");
	cout << "Loading yashe." << endl;
	yashe = new Yashe(baseFileName);
	cout << "Starting tests." << endl;

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
