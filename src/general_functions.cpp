#include "general_functions.h"

using namespace std;

void printMessageWithTime(string message){
	time_t moment = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << message << " | " << ctime(&moment)  << endl;
}

double random_double(){
	static std::default_random_engine e{((unsigned int)time(NULL))};
	static std::uniform_real_distribution<float> d{0, 30};	
	return d(e) + d(e);
}

inline bool equal(double a, double b){
	return a - EPSILON <= b && b <= a + EPSILON;
}

bool equal(const vector<double>& a, const vector<double>& b){
	unsigned int N = a.size();
	if (N != b.size())
		return false;
	for (unsigned int i = 0; i < N; i++){
		if (!equal(a[i], b[i]))
			return false;
	}
	return true;
}

bool equal(const vector<vector<double> >& a, const vector<vector<double> >& b){
	unsigned int N = a.size();
	if (N != b.size())
		return false;
	for (unsigned int i = 0; i < N; i++){
		if (!equal(a[i], b[i]))
			return false;
	}
	return true;
}
