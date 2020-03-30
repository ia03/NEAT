// NEAT.h
// I. Ahmed
//
// This declares the NEAT class, which manages an implementation of
// the NEAT algorithm.

#pragma once
#include <random>
#include <vector>
#include "Population.h"
#include "Species.h"
#include "Genome.h"
#include "Neuron_innovation.h"
#include "Link_innovation.h"

double sigmoid(double x);

class NEAT
{
public:
    // The ID of the next neuron to be created.
    int next_neuron_id;

    // The innovation ID of the next link to be created.
    int next_link_innov_id;

    // The genome used to initialize the population.
    Genome seed_genome;

    // Link innovations.
    std::vector<Link_innovation> link_innovations;

    // Neuron innovations.
    std::vector<Neuron_innovation> neuron_innovations;

    // The population.
    Population population;

    // The fitness calculation function.
    double (*calc_fitness)(Genome &genome);

    // Parameters.

    // The fitness required to consider the network successful.
    double goal_fitness = -0.1;

    // Population size.
    int population_size = 300;

    // Number of generations.
    int max_generation = 500;

    // Coefficient for the disjoint in the equation that determines whether
    // two genomes belong to the same species.
    double coeff_disjoint = 2.0;

    // Coefficient for the weight difference in the equation that determines
    // whether two genomes belong to the same species.
    double coeff_weight_diff = 0.4;

    // Threshold for determining whether 2 genomes belong to the same species.
    double same_species_thresh = 1.0;

    // Threshold for a species' stale age.
    int stale_age_thresh = 15;

    // Crossover rate.
    double crossover_rate = 0.75;

    // Percentage of each species to be removed each generation.
    double cull_species_percent = 0.5;

    // Probability to mutate a genome's weight.
    double mutate_weight_prob = 0.2;

    // Probability to mutate a genome's weight in a biased way (using Gaussian
    // perturb noise).
    double bias_weight_prob = 0.9;

    // The Gaussian noise variance to use when mutating genome weights in 
    // biased way.
    double mutate_weight_size = 0.1;

    // Probability to add a forward link.
    double mutate_add_forward_link_prob = 2.0;

    // Probability to add a bias link.
    double mutate_add_bias_link_prob = 0.4;

    // Probability to add a neuron.
    double mutate_add_neuron_prob = 0.5;

    // Probability to disable an enabled link.
    double mutate_enabled_prob = 0.2;

    // Probability to enable a disabled link.
    double mutate_disabled_prob = 0.2;

    // Species number threshold.
    int num_species_thresh = 10;

    // Constructor.
    NEAT() = default;

    // Finds a link innovation based on the from and to neuron IDs.
    std::vector<Link_innovation>::iterator find_link_innovation(
        const int from_neuron_id, const int to_neuron_id);

    // Finds a neuron innovation based on the innovation ID of the link to be split.
    std::vector<Neuron_innovation>::iterator find_neuron_innovation(
        const int link_split_innov_id);

    // Adds a link innovation based on the to and from neuron IDs.
    void add_link_innovation(const int from_neuron_id, const int to_neuron_id,
        Link_innovation &link_innovation);

    // Adds a neuron innovation based on the innovation ID of the link to be
    // split.
    void add_neuron_innovation(const int link_split_innov_id,
        Neuron_innovation &neuron_innovation);

    // Mutates a genome by adding a link.
    void mutate_add_link(Genome &genome, const Link_type link_type,
        const double probability);
    
    // Mutates a genome by adding a neuron.
    void mutate_add_neuron(Genome &genome, const double probability);

    // Toggles a link in a genome.
    void mutate_toggle(Genome &genome, const bool enabled,
        const double probability);

    // Mutates the weights of a genome.
    void mutate_weights(Genome &genome, const double probability,
        const double bias_probability, const double mutate_size);

    // Mutates a genome.
    void mutate(Genome &genome);

    // Performs a crossover on 2 genomes to obtain a child genome.
    Genome crossover(Genome &first_genome, Genome &second_genome);

    // Measures the disjoint between 2 genomes (includes exceed). The result is
    // the ratio of disjoint links to the number of links in the genome that has
    // the most links.
    double disjoint(const Genome &first_genome, const Genome &second_genome)
        const;

    // Measures the ratio of the weight difference between 2 genomes to the number
    // of times a link is present in both genomes.
    double weight_diff(const Genome &first_genome,
        const Genome &second_genome) const;

    // Determine whether two genomes should be considered of the same species.
    bool same_species(const Genome &first_genome, const Genome &second_genome)
        const;

    // Adds the genome to an existing species or creates a new species.
    void add_genome(Genome &genome);

    // Removes stale species.
    void remove_stale_species();

    // Returns a list of all the genomes in the population.
    std::vector<Genome> aggregate_genomes() const;

    // Calculates species' average ranks in a population. A higher result is
    // better.
    std::vector<double> calc_species_avg_rank() const;

    // Removes weak species. Weak species are those with an average rank lower
    // than the average average rank.
    void remove_weak_species();

    // Removes empty species.
    void remove_empty_species();

    // Removes the weakest x% of a species. This is determined based on the
    // cull_species_percent parameter.
    void cull_species();

    // Only keeps the best genome in each species.
    void cull_species_to_one();

    // Breeds a child for a species. This can be either sexual or asexual
    // reproduction. Returns true if succeeded.
    bool breed_child(Species &species, Genome &child);

    // Initializes the population by cloning the seed genome.
    void init_population();

    // Reproduce the next generation of the population.
    void reproduce();

    // Evaluate the population.
    void evaluate();

    // Evolve the population over multiple generations.
    void evolve();
};

