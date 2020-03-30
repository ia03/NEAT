// Species.h
// I. Ahmed
//
// This declares the Species class, which represents a species
// in the NEAT algorithm.

#pragma once
#include <vector>
#include "Genome.h"

class Species
{
public:
    // Genomes.
    std::vector<Genome> genomes;

    // ID.
    int id = 0;
    
    // Age.
    int age = 0;

    // Best fitness.
    double best_fitness = 0;

    // Best genome.
    Genome best_genome;

    // Stale age (how many generations in which the best fitness hasn't
    // improved).
    int stale_age = 0;

    // Next genome ID.
    int next_genome_id = 0;

    // Default constructor.
    Species() = default;

    // Clones the seed genome and randomizes its weights.
    Species(const Genome &seed_genome, int species_size);

    // Sets the best fitness to the maximum of all genomes' fitnesses and
    // updates the best genome correspondingly.
    void set_best_fitness_and_genome();

    // Sorts genomes by fitness (higher fitness is better) in descending
    // order.
    void sort_genomes();

    // Adds a genome.
    void add_genome(Genome &genome);
};

