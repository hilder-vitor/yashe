#ifndef __YASHE__MATRIX_UTILS
#define __YASHE__MATRIX_UTILS

#include <iostream>
#include <vector>

using namespace std;


template <typename ELEMENT>
void operator+=(vector<vector<ELEMENT> >& A, const vector<vector<ELEMENT> >& B){
	if (A.size() != B.size() || A[0].size() != B[0].size())
		throw std::invalid_argument("It is impossible to add matrix with different sizes.");
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = A[i][j] + B[i][j];
		}
	}
}

template <typename ELEMENT>
void operator+=(vector<vector<ELEMENT> >& A, const ELEMENT& c){
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = A[i][j] + c;
		}
	}
}

template <typename ELEMENT>
void operator-=(const vector<vector<ELEMENT> >& A, const vector<vector<ELEMENT> >& B){
	if (A.size() != B.size() || A[0].size() != B[0].size())
		throw std::invalid_argument("It is impossible to subtract matrix with different sizes.");
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = A[i][j] - B[i][j];
		}
	}
}

template <typename ELEMENT>
void operator-=(vector<vector<ELEMENT> >& A, const ELEMENT& c){
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = A[i][j] - c;
		}
	}
}

template <typename ELEMENT>
void operator*=(vector<vector<ELEMENT> >& A, const ELEMENT& c){
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = A[i][j] * c;
		}
	}
}

template <typename ELEMENT>
void operator*=(ELEMENT& c, const vector<vector<ELEMENT> >& A){
	unsigned int n = A.size();
	unsigned int m = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < m; j++){
			A[i][j] = c * A[i][j];
		}
	}
}

template <typename ELEMENT>
vector<ELEMENT> operator*(const vector<vector<ELEMENT> >& A, const vector<ELEMENT>& v){
	if (A[0].size() != v.size())
		throw std::invalid_argument("A and v don't have conformant dimensions to be multiplied.");

	unsigned int N = A.size();
	unsigned int M = A[0].size();

	vector<ELEMENT> resp;
	resp.reserve(M);

	#pragma omp parallel for
	for (unsigned int i = 0; i < N; i++){
		ELEMENT innerProd = A[i][0] * v[0];
		for (unsigned int j = 1; j < M; j++){
			innerProd = innerProd + A[i][j] * v[j];
		}
		resp.push_back(innerProd);
	}
	return resp;
}


template <typename ELEMENT, typename T2>
void operator/=(vector<vector<ELEMENT> >& A, const T2& m){
	unsigned int N = A.size();
	unsigned int M = A[0].size();
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < M; j++){
			A[i][j] = A[i][j] / m;
		}
	}
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator+(const vector<vector<ELEMENT> >& A, const vector<vector<ELEMENT> >& B){
	vector<vector<ELEMENT> > mat(A);
	mat += B;
	return mat;
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator+(const vector<vector<ELEMENT> >& A, const ELEMENT& c){
	vector<vector<ELEMENT> > B(A);
	B += c;
	return B;
}
template <typename ELEMENT>
vector<vector<ELEMENT> > operator+(const ELEMENT& c, const vector<vector<ELEMENT> >& A){
	return A + c;
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator-(const vector<vector<ELEMENT> >& A, const vector<vector<ELEMENT> >& B){
	vector<ELEMENT> C(A);
	C -= B;
	return C;
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator-(const vector<vector<ELEMENT> >& A, const ELEMENT& c){
	vector<vector<ELEMENT> > B(A);
	B -= c;
	return B;
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator*(const vector<vector<ELEMENT> >& A, const vector<vector<ELEMENT> >& B){
	if (A[0].size() != B.size())
		throw std::invalid_argument("A and B don't have conformant dimensions to be multiplied.");

	unsigned int N = A.size();
	unsigned int P = A[0].size();
	unsigned int M = B[0].size();

	vector<vector<ELEMENT> > C(N);
	
	#pragma omp parallel for
	for (unsigned int i = 0; i < N; i++){
		C[i].reserve(M);
		for (unsigned int j = 0; j < M; j++){
			ELEMENT innerProduct = A[i][0] * B[0][j];
			for (unsigned int k = 1; k < P; k++){
				innerProduct = innerProduct + A[i][k] * B[k][j];
			}
			C[i] = innerProduct;
		}
	}
	return C;
}

template <typename ELEMENT>
vector<vector<ELEMENT> > operator*(const vector<vector<ELEMENT> >& A, const ELEMENT& c){
	vector<vector<ELEMENT> > B(A);
	B *= c;
	return B;
}


template <typename ELEMENT>
vector<vector<ELEMENT> > operator*(const ELEMENT& c, const vector<vector<ELEMENT> >& A){
	return A * c;
}

template <typename ELEMENT, typename T2>
vector<vector<ELEMENT> > operator/(const vector<vector<ELEMENT> >& A, const T2& m){
	vector<vector<ELEMENT> > B(A);
	B /= m;
	return B;
}

template <typename ELEMENT>
std::ostream& operator<<(std::ostream& os, const vector<vector<ELEMENT> >& A){
	unsigned int N = A.size();
	unsigned int M = A[0].size();
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < M; j++){
			os << A[i][j] << " ";
		}
		os << endl;
	}
	return os;
}

template <class ELEMENT> class SymmetricMatrix{
	private:
		unsigned int N;
		vector<vector<ELEMENT> > A;

	public:

		SymmetricMatrix(unsigned int __N, const ELEMENT& zero) : N(__N), A(__N) {
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i].push_back(zero);
				}
			}
		}

		SymmetricMatrix(const vector<vector<ELEMENT> >& matrix) : N(matrix.size()), A(matrix.size()) {
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i].push_back(matrix[i][j]);
				}
			}
		}



		ELEMENT get(unsigned int i, unsigned int j) const{
			if (j > i)
				return A[j][i];
			return A[i][j];
		}

		void set(const ELEMENT& x, unsigned int i, unsigned int j){
			if (j > i){
				A[j][i] = ELEMENT(x);
			}else{
				A[i][j] = ELEMENT(x);
			}
		}

		unsigned int size() const{
			return N;
		}

		SymmetricMatrix<ELEMENT> operator*(const SymmetricMatrix<ELEMENT>& S){
			if (S.size() != N){
				throw std::invalid_argument("It is impossible to multiply symmetric matrices with different dimensions");
			}
			SymmetricMatrix<ELEMENT> prod(S.size(), S.get(0, 0));
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					ELEMENT innerProd = get(0, i) * S.get(0, j);
					for (unsigned int k = 1; k < N; k++){
						innerProd = innerProd + get(k, i) * S.get(k, j);
					}
					prod.set(innerProd, i, j);
				}
			}
			return prod; 
		}


		template <typename T>
		vector<ELEMENT> operator*(const vector<T>& v){
			if (v.size() != N){
				throw std::invalid_argument("Dimensions of symmetric matrices and vectorwith different dimensions");
			}
			vector<ELEMENT> prod(N, get(0, 0)); // create a N-dimensional vector with dummy values
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				ELEMENT innerProd = get(i, 0) * v[0];
				for (unsigned int j = 1; j < N; j++){
					innerProd = innerProd + get(i, j) * v[j];
				}
				prod[i] = innerProd;
			}
			return prod; 
		}

		void operator*=(const ELEMENT& c){
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] * c;
				}
			}
		}
		SymmetricMatrix<ELEMENT> operator*(const ELEMENT& c) const{
			SymmetricMatrix<ELEMENT> prod(*this);
			prod *= c;
			return prod;
		}

		void operator+=(const SymmetricMatrix<ELEMENT>& S){
			if (S.size() != N){
				throw std::invalid_argument("It is impossible to sum symmetric matrices with different dimensions.");
			}
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] + S.get(i, j);
				}
			}
		}
	
		void operator+=(const ELEMENT& c){
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] + c;
				}
			}
		}
		
		SymmetricMatrix<ELEMENT> operator+(const SymmetricMatrix<ELEMENT>& S) const{
			SymmetricMatrix<ELEMENT> sum(*this);
			sum += S;
			return sum;
		}

		SymmetricMatrix<ELEMENT> operator+(const ELEMENT& c){
			SymmetricMatrix<ELEMENT> sum(*this);
			sum += c;
			return sum;
		}

		void operator-=(const SymmetricMatrix<ELEMENT>& S){
			if (S.size() != N){
				throw std::invalid_argument("It is impossible to subtract matrices with different dimensions.");
			}
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] - S.get(i, j);
				}
			}
		}

		void operator-=(const ELEMENT& c){
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] - c;
				}
			}
		}

		SymmetricMatrix<ELEMENT> operator-(const ELEMENT& c) const{
			SymmetricMatrix<ELEMENT> sub(*this);
			sub -= c;
			return sub;
		}

		SymmetricMatrix<ELEMENT> operator-(const SymmetricMatrix<ELEMENT>& S) const{
			SymmetricMatrix<ELEMENT> sub(*this);
			sub -= S;
			return sub;
		}

		template <typename T2>
		void operator/=(const T2& scalar){
			#pragma omp parallel for
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j <= i; j++){
					A[i][j] = A[i][j] / scalar;
				}
			}
		}
		
		template <typename T2>
		SymmetricMatrix<ELEMENT> operator/(const T2& scalar){
			SymmetricMatrix<ELEMENT> div(*this);
			div /= scalar;
			return div;
		}

		vector<vector<ELEMENT> > to_usual_matrix() const{
			vector<vector<ELEMENT> > mat(N);
			for (unsigned int i = 0; i < N; i++){
				for (unsigned int j = 0; j < N; j++){
					mat[i].push_back(get(i, j));
				}
			}
			cout << "returning from to_usual_matrix" << endl;
			return mat;
		}
};

template <typename ELEMENT>
std::ostream& operator<<(std::ostream& os, const SymmetricMatrix<ELEMENT>& S){
	unsigned int N = S.size();
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++){
			os << S.get(i, j) << " ";
		}
		os << endl;
	}
	return os;
}

/* This function returns the inner product of the i-th column of A by the j-th column of B */
template <typename ELEMENT>
ELEMENT inner_product_column_by_column(const vector<vector< ELEMENT> >& A, unsigned int i, const vector<vector<ELEMENT> >& B, unsigned int j){
	unsigned int N = A.size();
	ELEMENT inner_product = A[0][i] * B[0][j];

	for (unsigned int k = 1; k < N; k++){
		inner_product = inner_product + A[k][i] * B[k][j];
	}
	return inner_product;
}


/**
 *	Receives a NxP matrix A and returns a PxP matrix equal to transpose(A) * A.
 *
*/
template <typename ELEMENT>
SymmetricMatrix<ELEMENT> multiply_transpose_matrix_by_matrix(const vector<vector< ELEMENT> >& A){
	unsigned int P = A[0].size();

	SymmetricMatrix<ELEMENT> C(P, A[0][0]);
	#pragma omp parallel for
	for (unsigned int i = 0; i < P; i++){
		for (unsigned int j = 0; j <= i; j++){
			C.set(inner_product_column_by_column(A, i, A, j), i, j);
		}
	}
	return C;
}


template <typename ELEMENT>
vector<vector<ELEMENT> > outer_product(const vector<ELEMENT>& u, const vector<ELEMENT>& v){
	if (u.size() != v.size()){
		throw std::invalid_argument("Outer product should be done with vectors of same dimension.");
	}
	unsigned int N = u.size();
	vector<vector<ELEMENT> > O(N);
	#pragma omp parallel for
	for (unsigned int i = 0; i < N; i++){
		O[i].reserve(N);
		for (unsigned int j = 0; j < N; j++){
			O[i].push_back(u[i] * v[j]);
		}
	}
	return O;
}
template <typename ELEMENT>
SymmetricMatrix<ELEMENT> outer_product(const vector<ELEMENT>& u){
	unsigned int N = u.size();
	SymmetricMatrix<ELEMENT> O(N, u[0]);
	#pragma omp parallel for
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j <= i; j++){
			O.set(u[i] * u[j], i, j);
		}
	}
	return O;
}

template <typename ELEMENT>
ELEMENT mean_of_column(const vector<vector<ELEMENT> >& A, unsigned int j){
	unsigned int N = A.size();
	ELEMENT mean = A[0][j];
	for (unsigned int i = 1; i < N; i++){
		mean = mean + A[i][j];
	}
	return mean / (N - 1);
}

template <typename ELEMENT>
void subtract_mean_column_wise(vector<vector<ELEMENT> >& A){
	unsigned int N = A.size();
	unsigned int P = A[0].size();

	#pragma omp parallel for
	for (unsigned int j = 0; j < P; j++){
		ELEMENT mean = mean_of_column(A, j);
		for (unsigned int i = 0; i < N; i++){
			A[i][j] =  mean - A[i][j];
		}
	}
}

template <typename ELEMENT>
SymmetricMatrix<ELEMENT> calculate_covariance(vector<vector<ELEMENT> >& data){
	subtract_mean_column_wise(data);
	return multiply_transpose_matrix_by_matrix(data);
}

template <typename ELEMENT>
vector<vector<ELEMENT> > get_lines(unsigned int line, unsigned int qnt, const vector<vector<ELEMENT> >& A){
	unsigned int N = A.size();
	vector<vector<ELEMENT> > B;
	for (unsigned int i = line; i < line + qnt && i < N; i++){
		B.push_back(A[i]);
	}
	return B;
}


#endif
