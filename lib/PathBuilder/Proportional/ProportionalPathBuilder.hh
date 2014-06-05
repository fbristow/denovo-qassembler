/* 
 * File:   ProportionalPathBuilder.hh
 * Author: fbristow
 *
 * Created on February 8, 2012
 */
#ifndef PROPORTIONAL_PATH_BUILDER_HH
#define PROPORTIONAL_PATH_BUILDER_HH

#include "PathBuilder/PathBuilder.hh"
#include "Logging/Logging.hh"

class ProportionalPathBuilder: public PathBuilder {
public:
	/**
	 * Constructor.
	 * @param graph the graph to use to construct paths
	 * @param epsilon the epsilon to use when constructing paths
	 */
	ProportionalPathBuilder(boost::shared_ptr<SkinnyGraph> graph, double epsilon);

	/**
	 * Construct paths from the supplied graph using a proportional approach.
	 * @return a set of proportional paths from the graph.
	 */
	boost::unordered_set<PathBuilder::Path> buildPaths();

	/**
	 * Get the epsilon used to build paths
	 * @return epsilon used to build paths
	 */
	double getEpsilon();

	/**
	 * Set the epsilon used to build paths
	 * @param epsilon the epsilon used to build paths
	 */
	void setEpsilon(double epsilon);
private:
	/** epsilon used to construct paths */
	double epsilon;
	
	/**
	 * Sum the weights of the edges incoming to a certain vertex
	 * @param v the vertex in question
	 */
	double sumIncomingEdges(SkinnyGraph::Vertex v);

	/**
	 * Sum the weights of the edges leaving a certain vertex
	 * @param v the vertex in question
	 */
	double sumOutgoingEdges(SkinnyGraph::Vertex v);

	/**
	 * Find the edge with the largest weight exiting a node
	 * @param v the vertex in question
	 * @return a tuple in which first is the edge descriptor and second is its weight
	 */
	boost::tuple<SkinnyGraph::Edge, double> findMaxEdge(SkinnyGraph::Vertex v);

	/**
	 * Find the edge exiting a ndoe that has a proportion closest to that specified
	 * @param v the vertex in question
	 * @param p the proportion to look for
	 * @return a tuple in which first is an edge descriptor and second is a boolean
	 * 	   indicating that (true) an edge with a proportion within epsilon WAS
	 * 	   found, or (false) no such edge was found.
	 */
	boost::tuple<SkinnyGraph::Edge, bool> findClosestEdge(SkinnyGraph::Vertex v, double p);
};

#endif // PROPORTIONAL_PATH_BUILDER_HH
