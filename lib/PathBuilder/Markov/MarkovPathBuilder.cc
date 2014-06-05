/*
 * MarkovPathBuilder.cc
 *
 *  Created on: 2012-06-29
 *      Author: fbristow
 */

#ifndef MARKOV_PATH_BUILDER_CC
#define MARKOV_PATH_BUILDER_CC

#include <algorithm>
#include <vector>
#include <boost/foreach.hpp>
#include "PathBuilder/Markov/MarkovPathBuilder.hh"

DECLARE_LOG(logger, "qassembler.MarkovPathBuilder");

MarkovPathBuilder::MarkovPathBuilder(boost::shared_ptr<SkinnyGraph> graph) : PathBuilder(graph),
		uniform_distribution(0, 1), generator(42u), rng(generator, uniform_distribution) {}

boost::unordered_set<PathBuilder::Path>
MarkovPathBuilder::buildPaths() {
	boost::unordered_set<PathBuilder::Path> paths;
	std::vector<SkinnyGraph::Vertex> startingPoints = getStartingPoints();
	boost::shared_ptr<SkinnyGraph::Graph> g = this->graph->graph();

	TRACE(logger, "Generating paths.");
	while (!startingPoints.empty()) {
		std::vector<SkinnyGraph::Edge> edgesFollowed;
		std::vector<SkinnyGraph::Vertex> verticesFollowed;
		SkinnyGraph::Vertex v = startingPoints.front();
		SkinnyGraph::Edge e;
		std::size_t smallestEdge = (std::size_t) -1;

		verticesFollowed.push_back(v);

		if (getOutgoingEdges(v).size() == 0) {
			startingPoints.erase(startingPoints.begin());
		} else {
			while (getOutgoingEdges(v).size() > 0) {
				std::string vertexName = this->graph->node(v)->getName();
				// if the vertex has only one outgoing edge, then follow it:
				TRACE(logger, "Vertex [" << vertexName << "].");
				if (getOutgoingEdges(v).size() == 1) {
					TRACE(logger, "Vertex [" << vertexName << "] has one outgoing edge, following that edge.");
					SkinnyGraph::OutgoingEdges edges;
					boost::tie(edges, boost::tuples::ignore) = boost::out_edges(v, *g);
					v = boost::target(*edges, *g);
					e = *edges;
				} else {
					TRACE(logger, "Vertex [" << vertexName << "] has multiple outgoing edges, picking edge.");
					// the vertex has multiple outgoing edges. going to select an edge by markov process:
					// first get the sum of weights of outgoing edges:
					double sum = sumOutgoingEdges(v);
					// create a list of all outgoing edges and their respective weights
					std::vector<EdgeWeightPair> outgoingEdges;
					TRACE(logger, "Constructing list of outgoing edge weights.");
					BOOST_FOREACH (SkinnyGraph::Edge e, getOutgoingEdges(v)) {
						outgoingEdges.push_back(std::make_pair(e, this->graph->edge(e)->getWeight() / sum));
					}
					TRACE(logger, "Sorting edge weights.");
					std::sort(outgoingEdges.begin(), outgoingEdges.end(), compareEdgeWeightPairs);
					// generate a random number:
					TRACE(logger, "Generating random number.");
					double markov = rng();
					TRACE(logger, "Random number is [" << markov << "]");
					// select an edge using the random number
					SkinnyGraph::Edge selectedEdge;
					bool selected = false;
					double total = 0.;
					TRACE(logger, "Selecting edge.");
					BOOST_FOREACH (EdgeWeightPair ewp, outgoingEdges) {
						TRACE(logger, "Testing [" << total << "] > [" << markov << "]");
						total += ewp.second;
						if (total > markov) {
							selectedEdge = ewp.first;
							selected = true;
							break;
						}
					}
					assert(selected);
					TRACE(logger, "Selected edge: [" << std::hex << selectedEdge << std::dec << "]");
					v = boost::target(selectedEdge, *g);
					TRACE(logger, "Target: [" << this->graph->node(v)->getName() << "]");
					e = selectedEdge;
				}

				verticesFollowed.push_back(v);
				edgesFollowed.push_back(e);

				TRACE(logger, "Determining if edge [" << std::hex << e << std::dec << "] is smaller than smallest edge.");
				std::size_t followedWeight = this->graph->edge(e)->getWeight();
				if (followedWeight > 1 && followedWeight < smallestEdge) {
					TRACE(logger, "New smallest edge weight: [" << followedWeight << "]");
					smallestEdge = followedWeight;
				} else {
					TRACE(logger, "Smallest edge weight remains the same [" << smallestEdge << "]");
				}
			}

			bool allPositive = true;
			TRACE(logger, "Reducing edge weights by [" << smallestEdge << "]");
			BOOST_FOREACH (SkinnyGraph::Edge e, edgesFollowed) {
				boost::shared_ptr<WeightedEdge> edge = this->graph->edge(e);
				std::size_t weight = edge->getWeight() - smallestEdge;
				edge->setWeight(weight);
				allPositive &= !edge->removed();
			}

			if (!allPositive) {
				TRACE(logger, "Removing path.");
				startingPoints.erase(startingPoints.begin());
			}
		}

		PathBuilder::Path path = verticesToSequenceNodes(verticesFollowed);
		paths.insert(path);


	}

	return paths;
}

double
MarkovPathBuilder::sumOutgoingEdges(SkinnyGraph::Vertex v) {
	double sum = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, getOutgoingEdges(v)) {
		sum += this->graph->edge(e)->getWeight();
	}

	return sum;
}

bool
compareEdgeWeightPairs (EdgeWeightPair a, EdgeWeightPair b) {
	return a.second < b.second;
}

#endif // MARKOV_PATH_BUILDER_CC
