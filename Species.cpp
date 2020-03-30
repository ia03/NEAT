// Species.cpp
// I. Ahmed
//
// This defines the Species class, which represents a species
// in the NEAT algorithm.

#include "Species.h"
#include "Genome.h"
#include <algorithm>
#include <cfloat>
#include <functional>

// Clones the seed genome and randomizes its weights.
Species::Species(const Genome &seed_genome, int species_size)
{
    next_genome_id = species_size;

    // Clone the genome, set member variables, randomize weights, and add it.
    for (int i = 0; i < species_size; i++)
    {
        Genome genome = seed_genome;
        genome.id = i;
        genome.randomize_weights(-1.0, 1.0);
        genomes.push_back(genome);
    }
}

// Sets the best fitness to the maximum of all genomes' fitnesses.
void Species::set_best_fitness_and_genome()
{
    double best_fitness = -DBL_MAX;

    // Go through each genome and update the best fitness accordingly.
    for (auto &genome : genomes)
    {
        if (genome.fitness > best_fitness)
        {
            best_fitness = genome.fitness;
            best_genome = genome;
        }
    }
}

// Sorts genomes by fitness (higher fitness is better) in descending
// order.
void Species::sort_genomes()
{
    std::sort(genomes.begin(), genomes.end(), std::greater<Genome>());
}

// Add genome.
void Species::add_genome(Genome &genome)
{
    genome.id = next_genome_id++;
    genomes.push_back(genome);
}