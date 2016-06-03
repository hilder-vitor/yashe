#include "yashe_coefficient_crt_tests.h"
#include <random>
#include<vector>
#include <stdio.h>

#define EPSILON 0.00001

using namespace std;

Yashe* yashe;

RealNumberPlaintextVector plain_to_real_number_plain(PlaintextVector vec, int precision){
	RealNumberPlaintextVector r_vec;
	unsigned int N = vec.size();
	for (unsigned int i = 0; i < N; i++){
		RealNumberPlaintext r(vec[i], precision);
		r_vec.push_back(r);
	}
	return r_vec;
}

PlaintextVector real_number_plain_to_plain(RealNumberPlaintextVector r_vec){
	PlaintextVector vec;
	unsigned int N = r_vec.size();
	for (unsigned int i = 0; i < N; i++){
		Plaintext p(r_vec[i]);
		vec.push_back(p);
	}
	return vec;
}

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

TEST (YASHE_pack_unpacking, addition_packed_plaintext) {
	vector<fmpzxx> t = {fmpzxx("2"), fmpzxx("5"), fmpzxx("39")};
	Plaintext p1(t[1], yashe->get_phi());
	p1.set(2, 1); // p1 = x^2

	Plaintext p2(t[2], yashe->get_phi());
	p2.set(2, 1);
	p2.set(3, 1); // p2 = x^3 + x^2

	Plaintext p3(t[3], yashe->get_phi());
	p3.set(8, 8);
	p3.set(7, 7);
	p3.set(10, 3); // p2 = 3*x^10 + 8*x^8 + 7*x^7

	CoefficientwiseCRT crt(t, yashe->get_phi()); 

	Plaintext p = crt.pack({p1, p2, p3});
	Plaintext q = crt.pack({p2, p3, p1});

	Plaintext r = p + q;
	vector<Plaintext> additions = crt.unpack(r);

	EXPECT_EQ(additions[0].get(0), ((p1+p2).get(0))%t[0]);
	EXPECT_EQ(additions[0].get(1), ((p1+p2).get(1))%t[0]);
	EXPECT_EQ(additions[0].get(2), ((p1+p2).get(2))%t[0]);

	EXPECT_EQ(additions[1].get(0), ((p2+p3).get(0))%t[1]);
	EXPECT_EQ(additions[1].get(1), ((p2+p3).get(1))%t[1]);
	EXPECT_EQ(additions[1].get(2), ((p2+p3).get(2))%t[1]);
	EXPECT_EQ(additions[1].get(7), ((p2+p3).get(7))%t[1]);
	EXPECT_EQ(additions[1].get(8), ((p2+p3).get(8))%t[1]);
	EXPECT_EQ(additions[1].get(10), ((p2+p3).get(10))%t[1]);

	EXPECT_EQ(additions[2].get(0), ((p3+p1).get(0))%t[2]);
	EXPECT_EQ(additions[2].get(1), ((p3+p1).get(1))%t[2]);
	EXPECT_EQ(additions[2].get(2), ((p3+p1).get(2))%t[2]);
	EXPECT_EQ(additions[2].get(7), ((p3+p1).get(7))%t[2]);
	EXPECT_EQ(additions[2].get(8), ((p3+p1).get(8))%t[2]);
	EXPECT_EQ(additions[2].get(10), ((p3+p1).get(10))%t[2]);
	EXPECT_EQ(10, additions[2].degree());
}

TEST (YASHE_pack_unpacking, addition_packed_realnumberplaintext) {
	vector<fmpzxx> t = {fmpzxx("9"), fmpzxx("10"), fmpzxx("31")};
	CoefficientwiseCRT crt(t, yashe->get_phi()); 

	RealNumberPlaintext p1 = yashe->encode(0.9281);
	RealNumberPlaintext p2 = yashe->encode(28.911);
	RealNumberPlaintext p3 = yashe->encode(10.1809);
	RealNumberPlaintext p = crt.pack({p1, p2, p3});

	RealNumberPlaintext q1 = yashe->encode(1.887);
	RealNumberPlaintext q2 = yashe->encode(8.911);
	RealNumberPlaintext q3 = yashe->encode(10.0);
	RealNumberPlaintext q = crt.pack({q1, q2, q3});

	RealNumberPlaintext sum = p + q;
	
	PlaintextVector u = crt.unpack(sum);
	RealNumberPlaintextVector vec;
	for (unsigned int i = 0; i < u.size(); i++)
		vec.push_back(RealNumberPlaintext(u[i], 64));

	EXPECT_EQ((p1 + q1).double_value(), vec[0].double_value());
	EXPECT_EQ((p2 + q2).double_value(), vec[1].double_value());
	EXPECT_EQ((p3 + q3).double_value(), vec[2].double_value());
}

TEST (YASHE_unpacking_packing, addition_packed_plaintext) {
	vector<fmpzxx> coprimes = {fmpzxx("2"), fmpzxx("3"), fmpzxx("5"), fmpzxx("7"), fmpzxx("11")};
	CoefficientwiseCRT crt(coprimes, yashe->get_phi()); 

	fmpzxx t(coprimes[0]*coprimes[1]*coprimes[2]*coprimes[3]*coprimes[4]);
	fmpzxx original_yashes_t = yashe->t;
	yashe->t = t;

	RealNumberPlaintext packed1 = yashe->encode(-0.9281);
	RealNumberPlaintext packed2 = yashe->encode(9.3512);
	RealNumberPlaintext packed3 = yashe->encode(21.827);

	RealNumberPlaintextVector vec1 = plain_to_real_number_plain(crt.unpack(packed1), 64);
	RealNumberPlaintextVector vec2 = plain_to_real_number_plain(crt.unpack(packed2), 64);
	RealNumberPlaintextVector vec3 = plain_to_real_number_plain(crt.unpack(packed3), 64);

	RealNumberPlaintextVector vec_add = vec1 + vec2;

	RealNumberPlaintext packed_add(crt.pack(real_number_plain_to_plain(vec_add)), 64);
	EXPECT_EQ((packed1 + packed2).double_value(), packed_add.double_value());

	vec_add = vec1 + vec2 + vec3;
	packed_add = RealNumberPlaintext(crt.pack(real_number_plain_to_plain(vec_add)), 64);
	EXPECT_EQ((packed1 + packed2 + packed3).double_value(), packed_add.double_value());
	yashe->t = original_yashes_t;
}


TEST (YASHE_unpacking_encrypting_packing, addition_packed_plaintext) {
	vector<fmpzxx> coprimes = {fmpzxx("2"), fmpzxx("3"), fmpzxx("5"), fmpzxx("7"), fmpzxx("11")};
	CoefficientwiseCRT crt(coprimes, yashe->get_phi()); 

	fmpzxx t(coprimes[0]*coprimes[1]*coprimes[2]*coprimes[3]*coprimes[4]);

	RealNumberPlaintext packed1 = yashe->encode(1.9281);
	RealNumberPlaintext packed2 = yashe->encode(9.3512);

	RealNumberCiphertextVector vec1 = yashe->encrypt(plain_to_real_number_plain(crt.unpack(packed1), 64));
	RealNumberCiphertextVector vec2 = yashe->encrypt(plain_to_real_number_plain(crt.unpack(packed2), 64));

	RealNumberCiphertextVector vec_add = vec1 + vec2;

	RealNumberPlaintext packed_add(crt.pack(real_number_plain_to_plain(yashe->decrypt(vec_add))), 64);
	EXPECT_EQ((packed1 + packed2).double_value(), packed_add.double_value());

	vec_add = vec1 + vec1 + vec2;
	packed_add = RealNumberPlaintext(crt.pack(real_number_plain_to_plain(yashe->decrypt(vec_add))), 64);
	EXPECT_EQ((packed1 + packed1 + packed2).double_value(), packed_add.double_value());
}

TEST (YASHE_pack_unpacking, multiplication_packed_realnumberplaintext) {
	vector<fmpzxx> t = {fmpzxx("27"), fmpzxx("32"), fmpzxx("31")};
	CoefficientwiseCRT crt(t, yashe->get_phi()); 

	RealNumberPlaintext p1 = yashe->encode(0.9281);
	RealNumberPlaintext p2 = yashe->encode(28.911);
	RealNumberPlaintext p3 = yashe->encode(10.1809);
	RealNumberPlaintext p = crt.pack({p1, p2, p3});

	RealNumberPlaintext q1 = yashe->encode(0.9281);
	RealNumberPlaintext q2 = yashe->encode(28.911);
	RealNumberPlaintext q3 = yashe->encode(10.1809);
	RealNumberPlaintext q = crt.pack({q1, q2, q3});

	RealNumberPlaintext sum = p * q;
	
	PlaintextVector u = crt.unpack(sum);
	RealNumberPlaintextVector vec;
	for (unsigned int i = 0; i < u.size(); i++)
		vec.push_back(RealNumberPlaintext(u[i], 128));

	EXPECT_FLOAT_EQ((p1 * q1).double_value(), vec[0].double_value());
	EXPECT_FLOAT_EQ((p2 * q2).double_value(), vec[1].double_value());
	EXPECT_FLOAT_EQ((p3 * q3).double_value(), vec[2].double_value());
}

TEST (YASHE_pack_unpacking, multiplication_packed_plaintext) {
	vector<fmpzxx> t = {fmpzxx("2"), fmpzxx("5"), fmpzxx("39")};
	Plaintext p1(t[1], yashe->get_phi());
	p1.set(2, 1); // p1 = x^2

	Plaintext p2(t[2], yashe->get_phi());
	p2.set(2, 1);
	p2.set(3, 1); // p2 = x^3 + x^2

	Plaintext p3(t[3], yashe->get_phi());
	p3.set(8, 8);
	p3.set(7, 7);
	p3.set(10, 3); // p2 = 3*x^10 + 8*x^8 + 7*x^7

	CoefficientwiseCRT crt(t, yashe->get_phi()); 

	Plaintext p = crt.pack({p1, p2, p3});
	Plaintext q = crt.pack({p2, p3, p1});

	Plaintext r = p * q;
	vector<Plaintext> multiplications = crt.unpack(r);

	EXPECT_EQ(multiplications[0].get(0), ((p1*p2).get(0))%t[0]);
	EXPECT_EQ(multiplications[0].get(1), ((p1*p2).get(1))%t[0]);
	EXPECT_EQ(multiplications[0].get(2), ((p1*p2).get(2))%t[0]);

	EXPECT_EQ(multiplications[1].get(0), ((p2*p3).get(0))%t[1]);
	EXPECT_EQ(multiplications[1].get(1), ((p2*p3).get(1))%t[1]);
	EXPECT_EQ(multiplications[1].get(2), ((p2*p3).get(2))%t[1]);
	EXPECT_EQ(multiplications[1].get(7), ((p2*p3).get(7))%t[1]);
	EXPECT_EQ(multiplications[1].get(8), ((p2*p3).get(8))%t[1]);
	EXPECT_EQ(multiplications[1].get(10), ((p2*p3).get(10))%t[1]);

	EXPECT_EQ(multiplications[2].get(0), ((p3*p1).get(0))%t[2]);
	EXPECT_EQ(multiplications[2].get(1), ((p3*p1).get(1))%t[2]);
	EXPECT_EQ(multiplications[2].get(2), ((p3*p1).get(2))%t[2]);
	EXPECT_EQ(multiplications[2].get(7), ((p3*p1).get(7))%t[2]);
	EXPECT_EQ(multiplications[2].get(8), ((p3*p1).get(8))%t[2]);
	EXPECT_EQ(multiplications[2].get(10), ((p3*p1).get(10))%t[2]);
}

TEST (YASHE_unpacking_packing, multiplication_packed_plaintext) {
	vector<fmpzxx> coprimes = {fmpzxx("2"), fmpzxx("3"), fmpzxx("5"), fmpzxx("7"), fmpzxx("11")};
	CoefficientwiseCRT crt(coprimes, yashe->get_phi()); 

	fmpzxx t(coprimes[0]*coprimes[1]*coprimes[2]*coprimes[3]*coprimes[4]);

	RealNumberPlaintext packed1 = yashe->encode(0.9281);
	RealNumberPlaintext packed2 = yashe->encode(9.3512);

	RealNumberPlaintextVector vec1 = plain_to_real_number_plain(crt.unpack(packed1), 64);
	RealNumberPlaintextVector vec2 = plain_to_real_number_plain(crt.unpack(packed2), 64);

	RealNumberPlaintextVector vec_prod = product_component_by_component (vec1 , vec2);

	RealNumberPlaintext packed_prod(crt.pack(real_number_plain_to_plain(vec_prod)), 128);
	EXPECT_EQ((packed1 * packed2).double_value(), packed_prod.double_value());

	vec_prod = product_component_by_component(product_component_by_component(vec1, vec1), vec2);
	packed_prod = RealNumberPlaintext(crt.pack(real_number_plain_to_plain(vec_prod)), 3*64);
	EXPECT_EQ((packed1 * packed1 * packed2).double_value(), packed_prod.double_value());
}


TEST (YASHE_unpacking_encrypting_packing, multiplication_packed_plaintext) {
	vector<fmpzxx> coprimes = {fmpzxx("2"), fmpzxx("3"), fmpzxx("5"), fmpzxx("7"), fmpzxx("11")};
	CoefficientwiseCRT crt(coprimes, yashe->get_phi()); 

	fmpzxx t(coprimes[0]*coprimes[1]*coprimes[2]*coprimes[3]*coprimes[4]);

	RealNumberPlaintext packed1 = yashe->encode(1.9281);
	RealNumberPlaintext packed2 = yashe->encode(9.3512);

	RealNumberCiphertextVector vec1 = yashe->encrypt(plain_to_real_number_plain(crt.unpack(packed1), 64));
	RealNumberCiphertextVector vec2 = yashe->encrypt(plain_to_real_number_plain(crt.unpack(packed2), 64));

	RealNumberCiphertextVector vec_prod = product_component_by_component(vec1, vec2);

	RealNumberPlaintext packed_prod(crt.pack(real_number_plain_to_plain(yashe->decrypt(vec_prod))), 128);
	EXPECT_EQ((packed1 * packed2).double_value(), packed_prod.double_value());

	vec_prod = product_component_by_component(product_component_by_component(vec1, vec1), vec2);
	packed_prod = RealNumberPlaintext(crt.pack(real_number_plain_to_plain(yashe->decrypt(vec_prod))), 3*64);
	EXPECT_EQ((packed1 * packed1 * packed2).double_value(), packed_prod.double_value());
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

