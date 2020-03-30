// main.cpp
// I. Ahmed
//
// This program evolves neural networks to perform an XOR
// operation using the NEAT algorithm.

#include <iostream>
#include "NEAT.h"
#include "Genome.h"
#include "random.h"

// Calculates the fitness of a genome using the XOR function.
double calc_fitness(Genome &genome)
{
    const int num_inputs = 4;
    double fitness = 0;
    std::vector<std::vector<double>> inputs;
    std::vector<double> outputs;

    // Set up the expected inputs.
    // Note that the last argument is the bias.
    inputs.push_back({0.0, 0.0, 1.0});
    inputs.push_back({0.0, 1.0, 1.0});
    inputs.push_back({1.0, 0.0, 1.0});
    inputs.push_back({1.0, 1.0, 1.0});

    // Set up the expected outputs.
    outputs.push_back(0.0);
    outputs.push_back(1.0);
    outputs.push_back(1.0);
    outputs.push_back(0.0);

    // For every combination of input/expected output,
    // input the input into the genome. Square the
    // difference between the actual output and the
    // expected output. Subtract the result from the
    // fitness.
    for (int i = 0; i < num_inputs; i++)
    {
        const double expected_output = outputs[i];
        genome.activate(inputs[i]);
        const double output = genome.output()[0];
        const double output_difference = std::abs(
            output - expected_output);
        fitness -= std::pow(output_difference, 2);
    }

    return fitness;
}

// Sets up the network to perform the XOR function.
NEAT setup_network()
{
    // Initialize the network.
    NEAT network;
    Genome &seed_genome = network.seed_genome;
    network.calc_fitness = calc_fitness;

    // Initialize the seed genome.
    seed_genome.id = 0;
    seed_genome.input_length = 2;
    seed_genome.output_length = 1;

    // Create the neurons.
    Neuron_gene first_input_neuron;
    first_input_neuron.id = 0;
    first_input_neuron.type = Neuron_type::INPUT;
    first_input_neuron.depth = 0.0;
    Neuron_gene second_input_neuron;
    second_input_neuron.id = 1;
    second_input_neuron.type = Neuron_type::INPUT;
    second_input_neuron.depth = 0.0;
    Neuron_gene bias_neuron;
    bias_neuron.id = 2;
    bias_neuron.type = Neuron_type::BIAS;
    bias_neuron.depth = 0.0;
    Neuron_gene output_neuron;
    output_neuron.id = 3;
    output_neuron.type = Neuron_type::OUTPUT;
    output_neuron.depth = 1.0;
    Neuron_gene hidden_neuron;
    hidden_neuron.id = 4;
    hidden_neuron.type = Neuron_type::HIDDEN;
    hidden_neuron.depth = 0.5;

    // Add the neurons.
    seed_genome.neuron_genes.push_back(first_input_neuron);
    seed_genome.neuron_genes.push_back(second_input_neuron);
    seed_genome.neuron_genes.push_back(bias_neuron);
    seed_genome.neuron_genes.push_back(output_neuron);
    seed_genome.neuron_genes.push_back(hidden_neuron);

    // Create the links.
    // 1. Links from input neurons to hidden neuron.
    // 2. Links from input neurons and hidden neuron
    // to output neuron.
    Link_gene link_1(0, 3, 0.0, 0, true);
    Link_gene link_2(1, 3, 0.0, 0, true);
    Link_gene link_3(2, 3, 0.0, 1, true);
    Link_gene link_4(0, 4, 0.0, 2, true);
    Link_gene link_5(1, 4, 0.0, 2, true);
    Link_gene link_6(2, 4, 0.0, 3, true);
    Link_gene link_7(4, 3, 0.0, 4, true);

    // Add the links.
    seed_genome.link_genes.push_back(link_1);
    seed_genome.link_genes.push_back(link_2);
    seed_genome.link_genes.push_back(link_3);
    seed_genome.link_genes.push_back(link_4);
    seed_genome.link_genes.push_back(link_5);
    seed_genome.link_genes.push_back(link_6);
    seed_genome.link_genes.push_back(link_7);

    // Set the next link and neuron innovation IDs.
    network.next_link_innov_id = seed_genome.link_genes.size();
    network.next_neuron_id = seed_genome.neuron_genes.size();
    return network;
}

// Sets up the network and evolves it.
int main()
{
    std::cout << "This program evolves neural networks to perform "
        << "the XOR operation.\n";
    NEAT network = setup_network();
    network.evolve();
}