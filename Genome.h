// Genome.h
// I. Ahmed
//
// This declares the Genome class, which represents a genome (neural
// network) in the NEAT algorithm.

#pragma once
#include <vector>
#include "Neuron_gene.h"
#include "Link_gene.h"

class Genome
{
public:
    // ID.
    int id;

    // Input length.
    int input_length = 0;

    // Output length.
    int output_length = 0;

    // Fitness.
    double fitness = 0;

    // Neuron genes.
    std::vector<Neuron_gene> neuron_genes;

    // Link genes.
    std::vector<Link_gene> link_genes;

    // Constructor.
    Genome() = default;

    // Finds a neuron.
    std::vector<Neuron_gene>::const_iterator find_neuron_const(const int id)
        const;
    std::vector<Neuron_gene>::iterator find_neuron(const int id);
    
    // Finds a link.
    std::vector<Link_gene>::const_iterator find_link_const(const int innov_id)
        const;
    std::vector<Link_gene>::iterator find_link(const int innov_id);
    std::vector<Link_gene>::const_iterator find_link_const(
        const Link_gene &link_gene) const;

    // Determines if a neuron is in this genome.
    bool contains_neuron(const int id) const;
    bool contains_neuron(Neuron_gene &neuron_gene) const;

    // Determines if a link is in this genome.
    bool contains_link(const int innov_id) const;
    bool contains_link(const Link_gene &link_gene) const;

    // Determines if a link is in this genome and is enabled.
    bool contains_enabled_link(const int innov_id) const;

    // Sets the neurons' inputs and activations to 0.
    void flush();

    // Sorts link genes by the to neuron's depth.
    void sort_link_genes();

    // Calculates the output of the genome given an input.
    void activate(const std::vector<double>& input);

    // Returns the output of the genome.
    std::vector<double> output();

    // Randomizes link weights between a low value and a high value.
    void randomize_weights(const double low, const double high);

    // Compares the fitnesses of 2 genomes (a higher fitness is better).
    // If the fitnesses are the same, compare the sizes (fewer links is better).
    bool operator<(const Genome &right_val) const;
    bool operator>(const Genome &right_val) const;

    // Adds a link gene.
    void add_link_gene(const Link_gene &link_gene);

    // Adds a neuron gene.
    void add_neuron_gene(const Neuron_gene &neuron_gene);

    // Inherits a link gene from another genome.
    void inherit_link_gene(const Genome &genome, const Link_gene &link_gene);

    // Counts the number of links present in this genome that are not present
    // in another genome. This means that excess links are included.
    int links_not_present_in(const Genome &genome) const;
};

// Finds a genome.
std::vector<Genome>::const_iterator find_genome(
    const std::vector<Genome> &genomes, const Genome &genome);