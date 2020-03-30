// Population.h
// I. Ahmed
//
// This declares the Population class, which represents a population
// of species in the NEAT algorithm.

#pragma once
#include <vector>
#include <cfloat>
#include "Species.h"

class Population
{
public:
    // Species.
    std::vector<Species> species_list;

    // Best fitness.
    double best_fitness = -DBL_MAX;

    // Best genome.
    Genome best_genome;

    // Next genome ID.
    int next_genome_id = 0;

    // Next species ID.
    int next_species_id = 0;

    // Default constructor.
    Population() = default;

    // Creates 1 species and clones the seed genome.
    Population(const Genome &seed_genome, int population_size);

    // Set the best genome and save its fitness.
    void set_best_fitness_and_genome();

    // Adds a species.
    void add_species(Species &species);

    // Removes a species.
    void remove_species(const int index);

    // Returns a list of all the genomes in the population.
    std::vector<Genome> aggregate_genomes();

    // Reassigns genomes' IDs from 0 to the population size.
    void reassign_genome_ids();
};

