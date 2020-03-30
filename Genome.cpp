// Genome.cpp
// I. Ahmed
//
// This defines the Genome class, which represents a genome (neural
// network) in the NEAT algorithm.

#include <algorithm>
#include <vector>
#include "Genome.h"
#include "random.h"

// Finds a neuron.
std::vector<Neuron_gene>::const_iterator Genome::find_neuron_const(
    const int id) const
{
    return std::find_if(neuron_genes.begin(), neuron_genes.end(),
        [id](const Neuron_gene &neuron_gene)
    {
        return neuron_gene.id == id;
    });
}

// Finds a neuron.
std::vector<Neuron_gene>::iterator Genome::find_neuron(
    const int id)
{
    return std::find_if(neuron_genes.begin(), neuron_genes.end(),
        [id](const Neuron_gene &neuron_gene)
    {
        return neuron_gene.id == id;
    });
}

// Finds a link.
std::vector<Link_gene>::const_iterator Genome::find_link_const(
    const int innov_id) const
{
    return std::find_if(link_genes.begin(), link_genes.end(),
        [innov_id](const Link_gene &link_gene)
    {
        return link_gene.innov_id == innov_id;
    });
}

// Finds a link.
std::vector<Link_gene>::iterator Genome::find_link(const int innov_id)
{
    return std::find_if(link_genes.begin(), link_genes.end(),
        [innov_id](const Link_gene &link_gene)
    {
        return link_gene.innov_id == innov_id;
    });
}

// Finds a link.
std::vector<Link_gene>::const_iterator Genome::find_link_const(
    const Link_gene &link_gene) const
{
    return find_link_const(link_gene.innov_id);
}

// Determines if a neuron is in this genome.
bool Genome::contains_neuron(const int id) const
{
    return find_neuron_const(id) != neuron_genes.end();
}
bool Genome::contains_neuron(Neuron_gene &neuron_gene) const
{
    return contains_neuron(neuron_gene.id);
}

// Determines if a link is in this genome.
bool Genome::contains_link(const int innov_id) const
{
    return find_link_const(innov_id) != link_genes.end();
}

// Determines if a link is in this genome.
bool Genome::contains_link(const Link_gene &link_gene) const
{
    return contains_link(link_gene.innov_id);
}

// Determines if a link is in this genome and is enabled.
bool Genome::contains_enabled_link(const int innov_id) const
{
    // Find the link.
    auto found_link = find_link_const(innov_id);

    // Make sure the link has been found and is enabled.
    return found_link != link_genes.end() && found_link->enabled;
}

// Sets the neurons' inputs and activations to 0.
void Genome::flush()
{
    // Go through each neuron gene and flush.
    for (auto &neuron_gene : neuron_genes)
    {
        neuron_gene.flush();
    }
}

// Sorts link genes by the to neuron's depth.
void Genome::sort_link_genes()
{
    struct Depth_and_link_gene
    {
        Depth_and_link_gene() = default;

        double depth;
        Link_gene link_gene;

        // Overload the < operator for std::sort.
        bool operator<(const Depth_and_link_gene& right_val) const
        {
            return depth < right_val.depth;
        }
    };

    std::vector<Depth_and_link_gene> depths_and_links;

    // Create a list of the links and their corresponding to neuron depths.
    for (const auto &link_gene : link_genes)
    {
        Depth_and_link_gene depth_and_link;
        depth_and_link.link_gene = link_gene;
        depth_and_link.depth = find_neuron_const(
            link_gene.to_neuron_id)->depth;
        depths_and_links.push_back(depth_and_link);
    }

    // Sort the list.
    std::sort(depths_and_links.begin(), depths_and_links.end());

    link_genes.clear();

    // Extract the link genes from the list to set the link gene list.
    for (const auto &depth_and_link : depths_and_links)
    {
        link_genes.push_back(depth_and_link.link_gene);
    }
}

// Calculate the output of the genome given an input.
void Genome::activate(const std::vector<double> &input)
{
    // Sort the link genes in order of the to neuron's depth.
    sort_link_genes();

    // Set neurons' inputs to 0.
    for (auto &neuron_gene : neuron_genes)
    {
        neuron_gene.input = 0.0;
    }

    // Set the input neurons.
    for (int i = 0; i < input_length; i++)
    {
        neuron_genes[i].activation = input[i];
    }

    // Activate hidden and output neurons.
    for (auto link = link_genes.begin(); link != link_genes.end(); link++)
    {
        auto to_neuron = find_neuron(link->to_neuron_id);
        auto from_neuron = find_neuron(link->from_neuron_id);

        // If the link is enabled, update the to neuron's input.
        if (link->enabled)
        {
            to_neuron->input += (from_neuron->activation * link->weight);
        }

        // If the next link does not exist or has a different to neuron,
        // calculate the activation of the to neuron of this link and set it.
        auto next_link = link + 1;

        if (next_link == link_genes.end() || 
            next_link->to_neuron_id != link->to_neuron_id)
        {
            to_neuron->calc_activation();
        }

    }
}

// Returns the output of the genome.
std::vector<double> Genome::output()
{
    std::vector<double> result;

    // Loop through the output genomes, adding the activation of each into
    // a vector.
    for (int i = input_length; i < (input_length + output_length); i++)
    {
        result.push_back(neuron_genes[i].activation);
    }

    return result;
}

// Randomizes link weights between a low value and a high value.
void Genome::randomize_weights(const double low, const double high)
{
    // Assign a random weight between the 2 values for each link gene.
    for (auto &link_gene : link_genes)
    {
        link_gene.weight = random_real_num(low, high);
    }
}

// Compares the fitnesses of 2 genomes (a higher fitness is better).
// If the fitnesses are the same, compare the sizes (fewer links is better).
bool Genome::operator<(const Genome &right_val) const
{
    // Higher fitness is better.
    if (fitness != right_val.fitness)
    {
        return fitness < right_val.fitness;
    }
    // Fewer links is better.
    else if (link_genes.size() != right_val.link_genes.size())
    {
        return link_genes.size() > right_val.link_genes.size();
    }
    else  // Pick randomly.
    {
        return true_or_false(0.5);
    }
}

// Compares the fitnesses of 2 genomes (a higher fitness is better).
    // If the fitnesses are the same, compare the sizes (fewer links is better).
bool Genome::operator>(const Genome &right_val) const
{
    return right_val < *this;
}

// Adds a link gene.
void Genome::add_link_gene(const Link_gene &link_gene)
{
    link_genes.push_back(link_gene);
}

// Adds a neuron gene.
void Genome::add_neuron_gene(const Neuron_gene &neuron_gene)
{
    neuron_genes.push_back(neuron_gene);
}

// Inherits a link gene from another genome.
void Genome::inherit_link_gene(const Genome &genome,
    const Link_gene &link_gene)
{
    add_link_gene(link_gene);

    // Make sure the from and to neurons are present in the child.
    if (!contains_neuron(link_gene.from_neuron_id))
    {
        // Add the from neuron gene to the child.
        auto &from_neuron_gene = *(genome.find_neuron_const(
            link_gene.from_neuron_id));

        add_neuron_gene(from_neuron_gene);
    }

    if (!contains_neuron(link_gene.to_neuron_id))
    {
        // Add the to neuron gene to the child.
        auto &to_neuron_gene = *(genome.find_neuron_const(
            link_gene.to_neuron_id));

        add_neuron_gene(to_neuron_gene);
    }
}

// Counts the number of links present in this genome that are not present
// in another genome. This means that excess links are included.
int Genome::links_not_present_in(const Genome &genome) const
{
    int result = 0;

    // Go through each link and check if it is present in the other genome.
    for (const auto &link_gene : link_genes)
    {
        if (!genome.contains_link(link_gene))
        {
            result++;
        }
    }

    return result;
}

// Finds a genome.
std::vector<Genome>::const_iterator find_genome(
    const std::vector<Genome> &genomes, const Genome &genome)
{
    return std::find_if(genomes.begin(), genomes.end(),
        [genome](const Genome &current_genome)
    {
        return current_genome.id == genome.id;
    });
}