// Link_innovation.h
// I. Ahmed
//
// This declares the Link_innovation struct, which represents a link
// innovation in the NEAT algorithm.

#pragma once

struct Link_innovation
{
    // Innovation ID.
    int innov_id;

    // The ID of the neuron this link innovation goes from.
    int from_neuron_id;

    // The ID of the neuron this link innovation goes to.
    int to_neuron_id;
};