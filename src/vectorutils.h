#ifndef __YASHE__VECTOR_UTILS
#define __YASHE__VECTOR_UTILS

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

template <typename ELEMENT>
void operator+=(vector<ELEMENT>& u, const vector<ELEMENT>& v){
	if (v.size() != u.size())
		throw std::invalid_argument("It is impossible to add vectors of different sizes.");
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] + v[i];
	}
}

template <typename ELEMENT>
void operator+=(vector<ELEMENT>& u, const ELEMENT& c){
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] + c;
	}
}

template <typename ELEMENT>
void operator-=(vector<ELEMENT>& u, const vector<ELEMENT>& v){
	if (v.size() != u.size())
		throw std::invalid_argument("It is impossible to subtract vectors of different sizes.");
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] - v[i];
	}
}
template <typename ELEMENT>
void operator-=(vector<ELEMENT>& u, const ELEMENT& c){
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] - c;
	}
}

template <typename ELEMENT>
void operator*=(vector<ELEMENT>& u, const ELEMENT& c){
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] * c;
	}
}

template <typename T1, typename T2>
void operator/=(vector<T1>& u, const T2& m){
	unsigned int n = u.size();
	for (unsigned int i = 0; i < n; i++){
		u[i] = u[i] / m;
	}
}

template <typename ELEMENT>
vector<ELEMENT> operator+(const vector<ELEMENT>& u, const vector<ELEMENT>& v){
	vector<ELEMENT> vec(u);
	vec += v;
	return vec;
}

template <typename ELEMENT>
vector<ELEMENT> operator+(const vector<ELEMENT>& u, const ELEMENT& c){
	vector<ELEMENT> vec(u);
	vec += c;
	return vec;
}
template <typename ELEMENT>
vector<ELEMENT> operator+(const ELEMENT& c, const vector<ELEMENT>& v){
	return v + c; 
}

template <typename ELEMENT>
vector<ELEMENT> operator-(const vector<ELEMENT>& u, const vector<ELEMENT>& v){
	vector<ELEMENT> vec(u);
	vec -= v;
	return vec;
}

template <typename ELEMENT>
vector<ELEMENT> operator-(const vector<ELEMENT>& u, const ELEMENT& c){
	vector<ELEMENT> vec(u);
	vec -= c;
	return vec;
}

template <typename ELEMENT>
ELEMENT operator*(const vector<ELEMENT>& u, const vector<ELEMENT>& v){
	if (u.size() != v.size())
		throw std::invalid_argument("It is impossible to multiply vectors of different sizes.");
	unsigned int n = u.size();
	ELEMENT innerProduct(u[0] * v[0]);
	for (unsigned int i = 1; i < n; i++){
		innerProduct = innerProduct + u[i] * v[i];
	}

	return innerProduct;
}

template <typename ELEMENT>
vector<ELEMENT> operator*(const vector<ELEMENT>& u, const ELEMENT& c){
	vector<ELEMENT> vec(u);
	vec *= c;
	return vec;
}


template <typename ELEMENT>
vector<ELEMENT> operator*(const ELEMENT& c, const vector<ELEMENT>& v){
	return v * c;
}

template <typename ELEMENT>
vector<ELEMENT> product_component_by_component(const vector<ELEMENT>& u, const vector<ELEMENT>& v){
	if (u.size() != v.size())
		throw std::invalid_argument("It is impossible to multiply vectors of different sizes.");
	unsigned int n = u.size();
	vector<ELEMENT> resp;
	for (unsigned int i = 0; i < n; i++){
		resp.push_back(u[i] * v[i]);
	}
	return resp;
}

template <typename T1, typename T2>
vector<T1> operator/(const vector<T1>& u, const T2& m){
	vector<T1> vec(u);
	vec /= m;
	return vec;
}

template <typename ELEMENT>
std::ostream& operator<<(std::ostream& os, const vector<ELEMENT>& u){
	unsigned int lastPosition = u.size() - 1;
	for (unsigned int i = 0; i < lastPosition; i++){
		os << u[i] << ", ";
	}
	os << u[lastPosition];
	return os;
}


#endif
