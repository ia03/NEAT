// NEAT.cpp
// I. Ahmed
//
// This defines the NEAT class, which manages an implementation of
// the NEAT algorithm.

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include "NEAT.h"
#include "Genome.h"
#include "random.h"

// Finds a link innovation based on the from and to neuron IDs.
std::vector<Link_innovation>::iterator NEAT::find_link_innovation(
    const int from_neuron_id, const int to_neuron_id)
{
    return std::find_if(link_innovations.begin(), link_innovations.end(),
        [from_neuron_id, to_neuron_id](
            const Link_innovation &link_innovation)
    {
        return link_innovation.from_neuron_id == from_neuron_id &&
            link_innovation.to_neuron_id == to_neuron_id;
    });
}

// Finds a neuron innovation based on the innovation ID of the link to be
// split.
std::vector<Neuron_innovation>::iterator NEAT::find_neuron_innovation(
    const int link_split_innov_id)
{
    return std::find_if(neuron_innovations.begin(), neuron_innovations.end(),
        [link_split_innov_id](const Neuron_innovation & neuron_innovation)
    {
        return neuron_innovation.link_split_innov_id == link_split_innov_id;
    });
}

// Adds a link innovation based on the to and from neuron IDs.
void NEAT::add_link_innovation(const int from_neuron_id,
    const int to_neuron_id, Link_innovation &link_innovation)
{
    link_innovation.from_neuron_id = from_neuron_id;
    link_innovation.to_neuron_id = to_neuron_id;
    link_innovation.innov_id = next_link_innov_id++;
    link_innovations.push_back(link_innovation);
}

// Adds a neuron innovation based on the innovation ID of the link to be
// split.
void NEAT::add_neuron_innovation(const int link_split_innov_id,
    Neuron_innovation &neuron_innovation)
{
    neuron_innovation.link_split_innov_id = link_split_innov_id;
    neuron_innovation.id = next_neuron_id++;
    neuron_innovation.in_link_innov_id = next_link_innov_id++;
    neuron_innovation.out_link_innov_id = next_link_innov_id++;
    neuron_innovations.push_back(neuron_innovation);
}

// Mutates a genome by adding a link.
void NEAT::mutate_add_link(Genome &genome, const Link_type link_type,
    const double probability)
{
    Link_gene link_gene;

    // Use a random number generator to determine whether or not this
    // mutation should occur.
    if (!true_or_false(probability))
    {
        return;
    }

    int from_neuron_idx;
    Neuron_gene *from_neuron;

    int to_neuron_idx;
    Neuron_gene *to_neuron;


    switch (link_type)
    {
    case Link_type::FORWARD:
        // Randomly pick from and to neurons.
        from_neuron_idx = random_idx(genome.neuron_genes);
        from_neuron = &genome.neuron_genes[from_neuron_idx];

        to_neuron_idx = random_idx(genome.neuron_genes);
        to_neuron = &genome.neuron_genes[to_neuron_idx];

        // Check for a same-depth connection.
        if (from_neuron->depth == to_neuron->depth)
        {
            return;
        }
        
        // Swap if the depths are reversed.
        if (from_neuron->depth > to_neuron->depth)
        {
            std::swap(from_neuron_idx, to_neuron_idx);
            std::swap(from_neuron, to_neuron);
        }
        break;
    case Link_type::BIAS:
        // Set the bias neuron as the from neuron.
        from_neuron_idx = genome.input_length - 1;
        from_neuron = &genome.neuron_genes[from_neuron_idx];

        // Set the to neuron (which can't be an input neuron).
        to_neuron_idx = random_int(genome.input_length, genome.neuron_genes.size() - 1);
        to_neuron = &genome.neuron_genes[to_neuron_idx];
        break;
    default:
        return;
    }

    auto found_link_innovation = find_link_innovation(
        from_neuron->id, to_neuron->id);

    // Check if the link innovation already exists.
    if (found_link_innovation != link_innovations.end())
    {
        // Create a link gene using the existing link innovation.
        link_gene.from_neuron_id = from_neuron->id;
        link_gene.to_neuron_id = to_neuron->id;
        link_gene.innov_id = found_link_innovation->innov_id;
        link_gene.weight = random_real_num();
        link_gene.enabled = true;
        genome.add_link_gene(link_gene);
        return;
    }
    // Create a new link innovation.
    Link_innovation new_link_innovation;
    add_link_innovation(from_neuron->id, to_neuron->id, new_link_innovation);

    // Create a link gene using the new link innovation.
    link_gene.from_neuron_id = from_neuron->id;
    link_gene.to_neuron_id = to_neuron->id;
    link_gene.innov_id = new_link_innovation.innov_id;
    link_gene.weight = random_real_num();
    link_gene.enabled = true;
    genome.add_link_gene(link_gene);
}

// Mutates a genome by adding a neuron.
void NEAT::mutate_add_neuron(Genome &genome, const double probability)
{
    Neuron_gene neuron_gene;
    Link_gene input_link_gene;
    Link_gene output_link_gene;

    // Use a random number generator to determine whether or not this
    // mutation should occur.
    if (!true_or_false(probability))
    {
        return;
    }

    // Make sure there are links.
    if (genome.link_genes.size() == 0)
    {
        return;
    }

    // Randomly select a link to split.
    int split_link_idx = random_idx(genome.link_genes);
    Link_gene *split_link = &genome.link_genes[split_link_idx];

    auto from_neuron = genome.find_neuron(split_link->from_neuron_id);
    auto to_neuron = genome.find_neuron(split_link->to_neuron_id);

    // Make sure the link is enabled.
    if (!split_link->enabled)
    {
        return;
    }

    // Disable the link.
    split_link->enabled = false;

    // Check if the neuron innovation already exists.
    auto found_neuron_innovation = find_neuron_innovation(split_link->innov_id);
    
    // If the innovation already exists.
    if (found_neuron_innovation != neuron_innovations.end())
    {
        // The neuron already exists.
        if (genome.contains_neuron(found_neuron_innovation->id))
        {
            // Enable the input and output links
            auto input_link = genome.find_link(
                found_neuron_innovation->in_link_innov_id);
            auto output_link = genome.find_link(
                found_neuron_innovation->out_link_innov_id);

            input_link->enabled = true;
            output_link->enabled = true;
        }
        else  // The neuron doesn't already exist.
        {
            // Create the neuron.
            neuron_gene.id = found_neuron_innovation->id;
            neuron_gene.depth = (from_neuron->depth + to_neuron->depth) / 2;
            neuron_gene.type = Neuron_type::HIDDEN;
            neuron_gene.input = 0.0;
            neuron_gene.activation = 0.0;

            // Create the input link.
            input_link_gene.from_neuron_id = from_neuron->id;
            input_link_gene.to_neuron_id = found_neuron_innovation->id;
            input_link_gene.innov_id = found_neuron_innovation->in_link_innov_id;
            input_link_gene.weight = 1.0;
            input_link_gene.enabled = true;

            // Create the output link.
            output_link_gene.from_neuron_id = found_neuron_innovation->id;
            output_link_gene.to_neuron_id = to_neuron->id;
            output_link_gene.innov_id = found_neuron_innovation->out_link_innov_id;
            output_link_gene.weight = split_link->weight;
            output_link_gene.enabled = true;

            // Add the neuron and links.
            genome.add_neuron_gene(neuron_gene);
            genome.add_link_gene(input_link_gene);
            genome.add_link_gene(output_link_gene);
        }
        return;
    }
    
    // Create neuron innovation, input link innovation, output link
    // innovation.
    Neuron_innovation new_neuron_innovation;
    add_neuron_innovation(split_link->innov_id, new_neuron_innovation);

    Link_innovation input_link_innovation;
    input_link_innovation.from_neuron_id = split_link->from_neuron_id;
    input_link_innovation.to_neuron_id = new_neuron_innovation.id;
    input_link_innovation.innov_id = new_neuron_innovation.in_link_innov_id;
    link_innovations.push_back(input_link_innovation);

    Link_innovation output_link_innovation;
    output_link_innovation.from_neuron_id = new_neuron_innovation.id;
    output_link_innovation.to_neuron_id = split_link->to_neuron_id;
    output_link_innovation.innov_id = new_neuron_innovation.out_link_innov_id;
    link_innovations.push_back(output_link_innovation);

    // Create neuron gene, input link gene, and output link gene.
    neuron_gene.id = new_neuron_innovation.id;
    neuron_gene.depth = (from_neuron->depth + to_neuron->depth) / 2;
    neuron_gene.type = Neuron_type::HIDDEN;
    neuron_gene.input = 0.0;
    neuron_gene.activation = 0.0;

    input_link_gene.from_neuron_id = split_link->from_neuron_id;
    input_link_gene.to_neuron_id = new_neuron_innovation.id;
    input_link_gene.innov_id = new_neuron_innovation.in_link_innov_id;
    input_link_gene.weight = 1.0;
    input_link_gene.enabled = true;

    output_link_gene.from_neuron_id = new_neuron_innovation.id;
    output_link_gene.to_neuron_id = split_link->to_neuron_id;
    output_link_gene.innov_id = new_neuron_innovation.out_link_innov_id;
    output_link_gene.weight = split_link->weight;
    output_link_gene.enabled = true;

    genome.add_neuron_gene(neuron_gene);
    genome.add_link_gene(input_link_gene);
    genome.add_link_gene(output_link_gene);
    
}

// Toggles a link.
void NEAT::mutate_toggle(Genome &genome, const bool enabled,
    const double probability)
{
    std::vector<std::vector<Link_gene>::iterator> matching_link_genes;

    // Use a random number generator to determine whether or not this
    // mutation should occur.
    if (!true_or_false(probability))
    {
        return;
    }

    // Loop through the links to find the matching ones.
    for (auto link_gene = genome.link_genes.begin();
        link_gene != genome.link_genes.end(); link_gene++)
    {
        if (link_gene->enabled == enabled)
        {
            matching_link_genes.push_back(link_gene);
        }
    }

    // Make sure at least 1 link was found.
    if (matching_link_genes.size() == 0)
    {
        return;
    }

    // Randomly pick a link to toggle.
    int picked_link_num = random_idx(matching_link_genes);
    matching_link_genes[picked_link_num]->enabled = !enabled;
}

// Mutates the weights of a genome.
void NEAT::mutate_weights(Genome &genome, const double probability,
    const double bias_probability, const double mutate_size)
{
    // Use a random number generator to determine whether or not this
    // mutation should occur.
    if (!true_or_false(probability))
    {
        return;
    }

    // Go through each link and mutate the weight.
    for (auto &link_gene : genome.link_genes)
    {
        // Biased mutation (takes current weight into account).
        if (true_or_false(bias_probability))
        {
            double delta_weight = random_distribution(0, mutate_size);
            link_gene.weight += delta_weight;
        }
        else  // Unbiased mutation.
        {
            link_gene.weight = random_distribution(0, mutate_size);
        }
    }
}

// Mutates a genome.
void NEAT::mutate(Genome &genome)
{
    // Mutate its weights.
    mutate_weights(genome, mutate_weight_prob, bias_weight_prob,
        mutate_weight_size);
    
    // Use while loops so the probability can be over 1.0.

    // Add forward links.
    double current_add_forward_link_prob = mutate_add_forward_link_prob;
    while (current_add_forward_link_prob > 0.0)
    {
        mutate_add_link(genome, Link_type::FORWARD, current_add_forward_link_prob);
        current_add_forward_link_prob--;
    }

    // Add bias links.
    double current_add_bias_link_prob = mutate_add_bias_link_prob;
    while (current_add_bias_link_prob > 0.0)
    {
        mutate_add_link(genome, Link_type::BIAS, current_add_bias_link_prob);
        current_add_bias_link_prob--;
    }

    // Add neurons.
    double current_add_neuron_prob = mutate_add_neuron_prob;
    while (current_add_neuron_prob > 0.0)
    {
        mutate_add_neuron(genome, current_add_neuron_prob);
        current_add_neuron_prob--;
    }

    // Toggle on links.
    double current_mutate_enabled_prob = mutate_enabled_prob;
    while (current_mutate_enabled_prob > 0.0)
    {
        mutate_toggle(genome, true, current_mutate_enabled_prob);
        current_mutate_enabled_prob--;
    }

    // Toggle off links.
    double current_mutate_disabled_prob = mutate_disabled_prob;
    while (current_mutate_disabled_prob > 0.0)
    {
        mutate_toggle(genome, false, current_mutate_disabled_prob);
        current_mutate_disabled_prob--;
    }
}

// Performs a crossover on 2 genomes to obtain a child genome.
Genome NEAT::crossover(Genome &first_genome, Genome &second_genome)
{
    // Make sure the first genome is better.
    if (first_genome < second_genome)
    {
        return crossover(second_genome, first_genome);
    }

    Genome child_genome;

    // Set the input and output lengths.
    child_genome.input_length = first_genome.input_length;
    child_genome.output_length = first_genome.output_length;

    // Add the input and output neurons.
    for (int i = 0; i < (first_genome.input_length + first_genome.output_length);
        i++)
    {
        child_genome.add_neuron_gene(first_genome.neuron_genes[i]);
    }

    // Go through the other neurons and links in the first genome.
    for (auto &link_gene : first_genome.link_genes)
    {
        // An excess or disjoint link from the better genome should be added
        // to the child.
        if (!second_genome.contains_link(link_gene))
        {
            child_genome.inherit_link_gene(first_genome, link_gene);
        }
        else  // Link shared by both parents.
        {
            // Randomly pick whose parent's link should be inherited.
            if (true_or_false(0.5))
            {
                child_genome.inherit_link_gene(first_genome, link_gene);
            }
            else
            {
                child_genome.inherit_link_gene(second_genome, link_gene);
            }
        }
    }

    return child_genome;
}

// Measures the disjoint between 2 genomes (includes exceed). The result is
// the ratio of disjoint links to the number of links in the genome that has
// the most links.
double NEAT::disjoint(const Genome &first_genome, const Genome &second_genome)
    const
{
    // Count the number of disjoint links.
    int num_of_disjoint_links = first_genome.links_not_present_in(
        second_genome) + second_genome.links_not_present_in(first_genome);

    // Get the larger genome size.
    int larger_genome_size = std::max(first_genome.link_genes.size(),
        second_genome.link_genes.size());

    // Divide them.
    return static_cast<double>(num_of_disjoint_links) /
        static_cast<double>(larger_genome_size);
}

// Measures the ratio of the weight difference between 2 genomes to the number
// of times a link is present in both genomes.
double NEAT::weight_diff(const Genome &first_genome,
    const Genome &second_genome) const
{
    double weight_difference = 0;
    int non_disjoint_links = 0;

    for (const auto &link_gene : first_genome.link_genes)
    {
        auto link_in_second_genome = second_genome.find_link_const(link_gene);
        if (second_genome.contains_link(link_gene))
        {
            // Consider the weight 0 if the link is disabled.
            weight_difference += std::abs((link_gene.weight * link_gene.enabled) -
                link_in_second_genome->weight * link_in_second_genome->enabled);

            non_disjoint_links++;
        }
    }

    return weight_difference / static_cast<double>(non_disjoint_links);
}

// Determine whether two genomes should be considered of the same species.
bool NEAT::same_species(const Genome &first_genome, const Genome &second_genome)
    const
{
    double disjoint_val = disjoint(first_genome, second_genome);
    double weight_diff_val = weight_diff(first_genome, second_genome);

    // Use the parameter coefficients to get a value which represents how
    // different the 2 genomes are.
    double combined_diff = (coeff_disjoint * disjoint_val) + (
        coeff_weight_diff * weight_diff_val);
    
    // Compare the value to the threshold.
    return combined_diff < same_species_thresh;
}

// Adds the genome to an existing species or creates a new species.
void NEAT::add_genome(Genome &genome)
{
    // Check if the genome can be added to any existing species.
    for (auto &species : population.species_list)
    {
        // Make sure the species is not empty.
        if (species.genomes.size() != 0)
        {
            // If the species' first genome is similar enough, add the genome
            // to that species.
            if (same_species(species.genomes[0], genome))
            {
                species.add_genome(genome);
                return;
            }
        }
    }

    // No matching species found, so create a new one.
    Species new_species;
    new_species.add_genome(genome);
    population.add_species(new_species);
}

// Removes stale species.
void NEAT::remove_stale_species()
{
    std::vector<Species> &species_list = population.species_list;

    // Go through all the species and find ones that are stale.
    int num_species = species_list.size();

    // Go through each species and determine if it is stale.
    for (int species = 0; species < num_species; species++)
    {
        if (species_list[species].stale_age > stale_age_thresh)
        {
            population.remove_species(species);
            num_species--;
        }
        else
        {
            species++;
        }
    }
}

// Returns a list of all the genomes in the population.
std::vector<Genome> NEAT::aggregate_genomes() const
{
    std::vector<Genome> result;

    // Iterate through every species and add each of its genomes to the
    // result.
    for (auto const &species : population.species_list)
    {
        for (auto const &genome : species.genomes)
        {
            result.push_back(genome);
        }
    }

    return result;
}

// Calculates species' average ranks in a population. A higher result is
// better.
std::vector<double> NEAT::calc_species_avg_rank() const
{
    std::vector<double> result;
    std::vector<Genome> genomes = aggregate_genomes();
    std::sort(genomes.begin(), genomes.end(), std::greater<Genome>());

    // Loop through each species
    for (auto &species : population.species_list)
    {
        int total_rank = 0;

        // Add up the ranks.
        for (auto &genome : species.genomes)
        {
            total_rank += genomes.cend() - find_genome(genomes, genome);
        }

        // Divide the total rank by the total number of genomes.
        result.push_back(static_cast<double>(total_rank) / static_cast<double>(
            species.genomes.size()));
    }

    return result;
}

// Removes weak species. Weak species are those with an average rank lower
// than the average average rank.
void NEAT::remove_weak_species()
{
    std::vector<double> species_avg_rank = calc_species_avg_rank();

    double total_avg_rank = std::accumulate(species_avg_rank.begin(),
        species_avg_rank.end(), 0);

    int current_species = 0;
    int num_species = population.species_list.size();

    // Go through each species and determine if it is weak.
    for (int species = 0; species < num_species; species++)
    {
        // If the species' rank is lower than the average, erase it.
        bool weak = (std::floor(species_avg_rank[species] *
            static_cast<double>(num_species) / total_avg_rank) < 1);

        // Erasing the species messes with the iteration, so an
        // increment is only done when the species is not being erased.
        if (weak)
        {
            population.remove_species(species);
            num_species--;
        }
        else
        {
            species++;
        }
    }
}

// Removes empty species.
void NEAT::remove_empty_species()
{
    // Go through each species and determine if it is empty.
    for (int species = 0; species < population.species_list.size(); species++)
    {
        // Erasing the species messes with the iteration, so an
        // increment is only done when the species is not being erased.
        if (population.species_list[species].genomes.size() == 0)
        {
            population.remove_species(species);
        }
    }
}

// Removes the weakest x% of a species. This is determined based on the
// cull_species_percent parameter.
void NEAT::cull_species()
{
    // Go through each species.
    for (auto &species : population.species_list)
    {
        // Sort the species.
        species.sort_genomes();

        // Figure out how many need to be removed.
        int num_to_remove = std::floor(
            species.genomes.size() * cull_species_percent);

        // After the genomes are sorted, remove the appropriate number of
        // genomes from the species.
        for (; num_to_remove > 0; num_to_remove--)
        {
            species.genomes.pop_back();
        }
    }
}

// Only keeps the best genome in each species.
void NEAT::cull_species_to_one()
{
    // Go through each species.
    for (auto &species : population.species_list)
    {
        // Make sure the species is not empty.
        if (species.genomes.size() > 0)
        {
            // Set best genome.
            species.set_best_fitness_and_genome();

            // Keep only the best genome.
            Genome best_genome = species.best_genome;
            species.genomes.clear();
            species.genomes.push_back(best_genome);
        }
    }

    // Remove empty species.
    remove_empty_species();
}

// Breeds a child for a species. This can be either sexual or asexual
// reproduction. Returns true if succeeded.
bool NEAT::breed_child(Species &species, Genome &child)
{
    // Make sure the species is not empty.
    if (species.genomes.size() == 0)
    {
        return false;
    }

    // Choose between a crossover (sexual reproduction) and a clone (asexual
    // reproduction).
    if (true_or_false(crossover_rate))  // Crossover.
    {
        int first_parent_idx = random_idx(species.genomes);
        int second_parent_idx = random_idx(species.genomes);

        // Don't breed a genome with itself.
        if (first_parent_idx != second_parent_idx)
        {
            child = crossover(species.genomes[first_parent_idx],
                species.genomes[second_parent_idx]);
        }
        else
        {
            return false;
        }
    }
    else  // Clone.
    {
        int parent_idx = random_idx(species.genomes);
        child = species.genomes[parent_idx];
    }

    // Mutate the child.
    mutate(child);

    return true;
}

// Initializes the population by cloning the seed genome.
void NEAT::init_population()
{
    population = Population(seed_genome, population_size);
}

// Reproduce the next generation of the population.
void NEAT::reproduce()
{
    // Keep the previous best genome.
    std::vector<Genome> children;
    children.push_back(population.best_genome);

    // Remove weak genomes in each species.
    cull_species();

    // Remove stale and weak species.
    if (population.species_list.size() > num_species_thresh)
    {
        remove_stale_species();
        remove_weak_species();
    }

    // Breed children in each species.
    std::vector<double> average_ranks = calc_species_avg_rank();
    double total_avg_rank = std::accumulate(average_ranks.begin(),
        average_ranks.end(), 0);

    // Go through each species.
    for (int i = 0; i < population.species_list.size(); i++)
    {
        // Calculate the number of children to breed based on the species'
        // average rank.
        int num_to_breed = std::floor(average_ranks[i] * static_cast<double>(
            population_size) / total_avg_rank) - 1;
        int num_successfully_bred = 0;

        // Breed until the number to breed has been reached.
        while (num_successfully_bred < num_to_breed)
        {
            Genome child;
            bool successful_birth = breed_child(population.species_list[i],
                child);
            if (successful_birth)
            {
                children.push_back(child);
                num_successfully_bred++;
            }
        }
    }

    // Keep the best member of each species.
    cull_species_to_one();

    // Randomly choose species and breed children until the population size
    // has been reached.
    int num_genomes = children.size() + population.species_list.size();
    while (num_genomes < population_size)
    {
        // Randomly select a species.
        int species_idx = random_idx(population.species_list);
        Genome child;

        // Breed a new child.
        bool successful_birth = breed_child(
            population.species_list[species_idx], child);

        if (successful_birth)
        {
            children.push_back(child);
            num_genomes++;
        }
    }

    // Speciate the new children into species.
    for (auto &child : children)
    {
        add_genome(child);
    }

    // Reassign genome IDs.
    population.reassign_genome_ids();
}

// Evaluates the population.
void NEAT::evaluate()
{
    // Loop through each genome in each species.
    for (auto &species : population.species_list)
    {
        for (auto &genome : species.genomes)
        {
            // Calculate the fitness.
            genome.flush();
            double fitness = calc_fitness(genome);
            genome.fitness = fitness;
        }
        // Set the best fitness and genome of the species.
        double old_best_fitness = species.best_fitness;
        species.set_best_fitness_and_genome();
        double new_best_fitness = species.best_fitness;

        // Check if the species is stale.
        if (new_best_fitness > old_best_fitness)
        {
            species.stale_age = 0;
        }
        else
        {
            species.stale_age++;
        }
    }

    // Set the best fitness and genome of the population.
    population.set_best_fitness_and_genome();
}

// Evolve the population over multiple generations.
void NEAT::evolve()
{
    // Generate initial species at random.
    int generation = 0;
    init_population();

    // Speciate the genome into species.
    std::vector<Genome> genomes = aggregate_genomes();
    population.species_list.clear();
    for (auto &genome : genomes)
    {
        add_genome(genome);
    }

    // Loop through successive generations.
    while (generation < max_generation)
    {
        evaluate();
        const double best_fitness = population.best_fitness;
        const Genome &best_genome = population.best_genome;
        std::cout << "Best fitness of generation " << generation << " is "
            << best_fitness << " for genome " << best_genome.id << ". There "
            "are " << population.species_list.size() << " species.\n";
        if (best_fitness >= goal_fitness)
        {
            std::cout << "Task succeeded.\n";
            return;
        }
        reproduce();
        generation++;
    }
}