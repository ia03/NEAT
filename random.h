// random.h
// I. Ahmed
//
// This declares the random functions, which use the C++ STL to
// generate pseudo-random numbers.

#pragma once

// Generates a random real number between the low and high values.
double random_real_num(const double low, double high);

// Generates a random real number between 0 and 1.
double random_real_num();

// Generates a random integer between the low and high values.
int random_int(const int low, const int high);

// Returns a random index. The index can be any number from 0 to size - 1.
template<typename T>
int random_idx(const T &container)
{
    return random_int(0, container.size() - 1);
}

// Uses a random number generator to return true or false based on the
// inputted probability
bool true_or_false(const double probability);

// Generates a normally distributed real number using the mean and variance.
double random_distribution(const double mean, const double variance);