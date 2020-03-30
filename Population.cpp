// Population.cpp
// I. Ahmed
//
// This defines the Population class, which represents a population
// of species in the NEAT algorithm.

#include "Population.h"
#include "Species.h"

// Creates 1 species and clones the seed genome.
Population::Population(const Genome &seed_genome, int population_size)
{
    Species species = Species(seed_genome, population_size);
    species_list.push_back(species);
    next_species_id = 1;
    next_genome_id = population_size;
}

// Set the best genome and save its fitness.
void Population::set_best_fitness_and_genome()
{
    double current_best_fitness = -DBL_MAX;

    // Iterate through all the genomes in each species and update the best
    // fitness and genome if necessary.
    for (auto &species : species_list)
    {
        for (auto &genome : species.genomes)
        {
            if (genome.fitness > current_best_fitness)
            {
                current_best_fitness = genome.fitness;
                best_genome = genome;
            }
        }
    }
    best_fitness = current_best_fitness;
}

// Add a species.
void Population::add_species(Species &species)
{
    species.id = next_species_id++;
    species_list.push_back(species);
}

// Returns a list of all the genomes in the population.
std::vector<Genome> Population::aggregate_genomes()
{
    std::vector<Genome> result;

    // Loop through the species and add each genome in the species.
    for (auto &species : species_list)
    {
        for (auto &genome : species.genomes)
        {
            result.push_back(genome);
        }
    }
    return result;
}

// Removes a species.
void Population::remove_species(const int index)
{
    species_list.erase(species_list.begin() + index);
}

// Reassigns genomes' IDs from 0 to the population size.
void Population::reassign_genome_ids()
{
    int id = 0;

    // Go through each species.
    for (auto &species : species_list)
    {
        // Go through each genome.
        for (auto &genome : species.genomes)
        {
            // Set the genome's ID.
            genome.id = id;
            id++;
        }
    }

    next_genome_id = id + 1;
}