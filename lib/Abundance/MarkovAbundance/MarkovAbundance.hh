/*
 * File:   MarkovAbundance.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef MARKOV_ABUNDANCE_HH
#define MARKOV_ABUNDANCE_HH

#include "Abundance/Abundance.hh"

#include "Logging/Logging.hh"

class MarkovAbundance : public Abundance {
public:
	/**
	 * Compute abundances for the paths provided.
	 */
	virtual boost::unordered_map<std::string, double> computeAbundances() = 0;
protected:
	/**
	 * Constructor.
	 * @param graph the graph to use to construct paths
	 * @param paths the paths to compute abundances for
	 */
	MarkovAbundance(boost::shared_ptr<HeftyGraph> graph, boost::unordered_set<std::string> paths);

	/**
	 * Get the sum of the edges that transition from the 'start' state to a state that
	 * emits symbols. The total edge weight is computed by determining total number of instances
	 * of the first kmer in the vertices which have no incoming edges.
	 * @return the total sum of 'edge weights' that have vertices that have no incoming edges.
	 */
	std::size_t getBeginStateTransitionSum();

private:
	std::size_t beginStateTransitionSum;

	/**
	 * Compute the total sum of transitions from the start state to a state in the model.
	 * @return the toal sum of edges.
	 */
	std::size_t computeBeginStateTransitionSum();
};

#endif // MARKOV_ABUNDANCE_HH
