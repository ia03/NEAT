// Neuron_innovation.h
// I. Ahmed
//
// This declares the Neuron_innovation struct, which represents a
// neuron innovation in the NEAT algorithm.

#pragma once

struct Neuron_innovation
{
    // ID.
    int id;

    // The innovation ID of the link that was split to create this neuron.
    int link_split_innov_id;

    // The innovation ID of the output link when this neuron was created.
    int out_link_innov_id;

    // The innovation ID of the input link when this neuron was created.
    int in_link_innov_id;
};