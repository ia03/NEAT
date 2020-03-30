// random.cpp
// I. Ahmed
//
// This defines the random functions, which use the C++ STL to
// generate pseudo-random numbers.

#include <random>
#include "random.h"

// Random number generator.
typedef std::mt19937 Rand_Num_Gen;

std::random_device random_device;

Rand_Num_Gen rand_num_gen(random_device());

// Generates a random real number between the low and high values.
double random_real_num(const double low, double high)
{
    std::uniform_real_distribution<> generator(low, high);
    return generator(rand_num_gen);
}

// Generates a random real number between 0 and 1.
double random_real_num()
{
    return random_real_num(0.0, 1.0);
}

// Generates a random integer between the low and high values.
int random_int(const int low, const int high)
{
    std::uniform_int_distribution<> generator(low, high);
    return generator(rand_num_gen);
}

// Uses a random number generator to determine whether the answer is "yes" or
// "no."
bool true_or_false(const double probability)
{
    if (random_real_num() < probability)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Generates a normally distributed real number using the mean and variance.
double random_distribution(const double mean, const double variance)
{
    std::normal_distribution<> normal_distribution;
    return variance * normal_distribution(rand_num_gen) + mean;
}