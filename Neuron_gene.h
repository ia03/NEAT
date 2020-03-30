// Neuron_gene.h
// I. Ahmed
//
// This declares the Neuron_gene class, which represents a
// neuron gene in the NEAT algorithm.

#pragma once

// The type of neuron.
enum class Neuron_type
{
    INPUT,
    BIAS,
    OUTPUT,
    HIDDEN
};

// Neuron gene class.
class Neuron_gene
{
public:
    // Neuron ID.
    int id;

    // Neuron type.
    Neuron_type type;

    // Input sum.
    double input;

    // Activation.
    double activation;

    // Depth. The depth for input and bias neurons is 0, whereas the depth for
    // output neurons is 1. The depth for hidden nodes is in between.
    double depth;

    // Constructor
    Neuron_gene() = default;

    // Calculate the activation.
    void calc_activation();

    // Set the input and activation to 0.
    void flush();
};