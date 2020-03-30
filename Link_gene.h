// Link_gene.h
// I. Ahmed
//
// This declares the Link_gene struct, which represents a link gene
// in the NEAT algorithm.

#pragma once

enum class Link_type
{
    FORWARD,
    BIAS
};

struct Link_gene
{
    // Constructor.
    Link_gene() = default;
    
    // Constructor.
    Link_gene(int from_neuron_id, int to_neuron_id, double weight,
        int innov_id, bool enabled) : from_neuron_id(from_neuron_id),
        to_neuron_id(to_neuron_id), weight(weight), innov_id(innov_id),
        enabled(enabled) {};

    // From neuron ID.
    int from_neuron_id;
    
    // To neuron ID.
    int to_neuron_id;

    // Weight.
    double weight;

    // Innovation ID.
    int innov_id;

    // Used to determine whether or not the link is enabled.
    bool enabled;
};