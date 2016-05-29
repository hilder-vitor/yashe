#ifndef YASHE__GENERAL_FUNCTIONS
#define YASHE__GENERAL_FUNCTIONS


#include <ctime>
#include <chrono>
#include <random>
#include <vector>
#include <string>
#include <iostream>

#define EPSILON 0.0001

void printMessageWithTime(std::string message);

double random_double();

inline bool equal(double a, double b);

bool equal(const std::vector<double>& a, const std::vector<double>& b);

bool equal(const std::vector<std::vector<double> >& a, const std::vector<std::vector<double> >& b);

#endif
