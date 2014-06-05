/* 
 * File:   WeightedEdge.hh
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef WEIGHTED_EDGE_HH
#define WEIGHTED_EDGE_HH

#include <boost/cstdint.hpp>

class WeightedEdge {
public:
	/**
	 *  Default constructor. Initial edge weight is 1.
	 */
	WeightedEdge();
	/**
	 * Constructor specifying edge weight.
	 * @param weight the weight of this edge
	 */
	WeightedEdge(std::size_t weight);

	/**
	 * The current weight of this edge.
	 * @return the weight of the edge
	 */
	std::size_t getWeight();
	/**
	 * Set the current weight of this edge.
	 * @param weight the weight to set
	 */
	void setWeight(std::size_t weight);
	/**
	 * Overloaded ++ operator adds 1 to edge weight.
	 */
	WeightedEdge operator++(int);

	/**
	 * Increase the weight of this edge by the specified amount.
	 * @param amount the amount to increase the weight by.
	 */
	void increaseWeight(std::size_t amount);

	/**
	 * Decrease the weight of this edge by the specified weight.
	 * @param amount the amount to decrease the weight by.
	 */
	void decreaseWeight(std::size_t amount);

	/**
	 * Copy the old edge weight back.
	 */
	void resetWeight();

	/**
	 * Lock an edge weight. When reset is called later, the current edge
	 * weight when the edge was locked is set as the current weight.
	 */
	void lockWeight();

	/**
	 * If the weight of this edge is below 1, then the edge should be 
	 * treated as though it were removed from the graph.
	 * @return true if the weight is less than 1, false otherwise.
	 */
	bool removed();
private:
	/** the current edge weight */
	std::size_t weight;
	/** the original weight when constructed */
	std::size_t lockedWeight;
};

#endif // WEIGHTED_EDGE_HH
