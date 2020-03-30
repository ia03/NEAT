// Neuron_gene.cpp
// I. Ahmed
//
// This defines the Neuron_gene class, which represents a
// neuron gene in the NEAT algorithm.

#include <cmath>
#include "Neuron_gene.h"

// Calculates the sigmoid of x.
double sigmoid(double x)
{
    return x / (1 + std::abs(x));
}

// Calculates the activation.
void Neuron_gene::calc_activation()
{
    // Calculate the sigmoid of the input.
    activation = sigmoid(input);
}

// Set the input and activation to 0.
void Neuron_gene::flush()
{
    input = 0.0;
    activation = 0.0;
}