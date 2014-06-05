/* 
 * File:   ProportionalPathBuilder.cc
 * Author: fbristow
 *
 * Created on February 8, 2012
 */
#ifndef PROPORTIONAL_PATH_BUILDER_CC
#define PROPORTIONAL_PATH_BUILDER_CC

#include <boost/foreach.hpp>

#include "PathBuilder/Proportional/ProportionalPathBuilder.hh"

DECLARE_LOG(logger, "qassembler.ProportionalPathBuilder");

ProportionalPathBuilder::ProportionalPathBuilder(boost::shared_ptr<SkinnyGraph> graph, double epsilon) : PathBuilder(graph) {
	this->epsilon = epsilon;
}

boost::unordered_set<PathBuilder::Path>
ProportionalPathBuilder::buildPaths() {
	boost::unordered_set<PathBuilder::Path> paths;
	std::vector<SkinnyGraph::Vertex> startingPoints = getStartingPoints();
	boost::shared_ptr<SkinnyGraph::Graph> g = this->graph->graph();


	while (!startingPoints.empty()) {
		// we haven't yet selected a proportion, start with a sentinel value of -1
		double p = -1.;
		std::vector<SkinnyGraph::Edge> followed;
		SkinnyGraph::Vertex v = startingPoints.front();
		SkinnyGraph::Edge lastEdge;
		std::vector<SkinnyGraph::Vertex> verticesFollowed;
		PathBuilder::Path path;
		std::string constructed = this->graph->node(v)->fullSequence();
		std::size_t smallestEdge = (std::size_t) -1; // cast -1 to std::size_t (which is unsigned), so we get the maximal value

		TRACE(logger, "Initial sequence from node [" << this->graph->node(v)->getName() << "] is: [" << constructed << "]");

		// keep track of all of the vertices we've already covered, including the
		// starting point
		verticesFollowed.push_back(v);

		// for each starting point, we're going to continually follow a path until we
		// arrive at a node that has no more outgoing edges. When a node has no outgoing
		// edges, then we've reached the end of the possible path that we're searching.
		while (getOutgoingEdges(v).size() > 0) {
			std::string vertexName = this->graph->node(v)->getName();
			TRACE(logger, "This node has outgoing edges, going to pick which edge to follow.");
			double sum = 0;
			// if v has more than one incoming edge and we haven't yet selected a proportion, then we should
			// select a proportion based upon the weight of the edges incoming to this node and the edge we took
			// to get here.
			if (getIncomingEdges(v).size() > 1 && p < 0) {
				TRACE(logger, "Vertex " << vertexName << " has more than one incoming edge [" << getIncomingEdges(v).size() << "], but we haven't selected a proportion. Defining proportion.");
				sum = sumIncomingEdges(v);
				boost::shared_ptr<WeightedEdge> e = this->graph->edge(lastEdge);
				p = e->getWeight() / sum;
				TRACE(logger, "Selected proportion: [" << p << "].");
			}

			// if v only has one outgoing edge, then we only have one possible edge to take, so take
			// that edge.
			if (getOutgoingEdges(v).size() == 1) {
				TRACE(logger, "Vertex " << vertexName << " only has one outgoing edge, following that edge.");
				SkinnyGraph::OutgoingEdges e;
				boost::tie(e, boost::tuples::ignore) = boost::out_edges(v, *g);
				v = boost::target(*e, *g);
				lastEdge = *e;
			} else {
				TRACE(logger, "Vertex " << vertexName << " has more than one outgoing edge, going to decide which to use.");
				// v doesn't have only one outgoing edge. We're going to start by finding the edge with
				// the largest outgoing weight to see if that's pretty close to the proportion we previously
				// selected.
				SkinnyGraph::Edge maxEdge; std::size_t maxWeight;
				boost::tie(maxEdge, maxWeight) = findMaxEdge(v);
				sum = sumOutgoingEdges(v);

				if (p < 0) {
					TRACE(logger, "We haven't selected a proportion yet, so we're just going to pick the edge with the largest weight.");
					// we haven't selected a proportion yet, so we can define it using the max edge
					// weight and the sum of the outgoing edges.
					p = maxWeight / sum;
					v = boost::target(maxEdge, *g);
					lastEdge = maxEdge;

					TRACE(logger, "Selected proportion [" << p << "], following vertex " << this->graph->node(v)->getName() << ".");
				} else {
					TRACE(logger, "We've already selected a proportion, going to try finding a similar edge.");
					// we've already selected a proportion. Try to find an edge exiting this node that
					// has a proportion similar to what we selected already.
					SkinnyGraph::Edge closest; bool found;
					boost::tie(closest, found) = findClosestEdge(v, p);

					if (found) {
						// we found a proportion similar to the one we'd defined already,
						// we don't need to modify the proportion we've selected, just follow
						// that edge to the next node.
						v = boost::target(closest, *g);
						lastEdge = closest;

						TRACE(logger, "Found a similar proportional edge, following to " << this->graph->node(v)->getName() << ".");
					} else {
						TRACE(logger, "We couldn't find a similar edge, picking the edge with the largest outgoing weight.");
						// we couldn't find a proportion similar to the one we defined already, so
						// now we need to decide whether or not the largest outgoing edge has a proportion
						// that's larger than the proportion we selected previously.
						double maxEdgeP = maxWeight / sum;
						// we're definitely going to follow the edge:
						v = boost::target(maxEdge, *g);
						lastEdge = maxEdge;
						// now decide if we need to redefine the proportion:
						if (maxEdgeP < p) {
							TRACE(logger, "The edge with the largest weight has a proportion that's less than what we selected before: [" << maxEdgeP << " < " << p << "], using that as a guiding proportion.");
							// if the edge with the largest weight has a proportion smaller than what
							// we selected previously, we should redefine the proportion to be that proportion.
							p = maxEdgeP;
						} else {
							TRACE(logger, "The edge with the largest weight has a proportion that's larger than what we selected before: [" << maxEdgeP << " > " << p << "], not selecting new guide proportion.");
						}
					}
				}
			}

			if (std::find(verticesFollowed.begin(), verticesFollowed.end(), v) != verticesFollowed.end()) {
				TRACE(logger, "Encountered cycle, breaking from path.");
				break;
			}

			// keep track of the edges we followed to get to the vertex we're visiting.
			followed.push_back(lastEdge);

			std::size_t followedWeight = this->graph->edge(lastEdge)->getWeight();
			if (followedWeight > 1 && followedWeight < smallestEdge) {
				smallestEdge = this->graph->edge(lastEdge)->getWeight();
			}
			// keep track of the vertices that we've visited on this path
			verticesFollowed.push_back(v);

			constructed += this->graph->node(v)->sequence();

			TRACE(logger, "After adding node [" << this->graph->node(v)->getName() << "], sequence is: [" << constructed << "]");

//			if (this->graph->edge(lastEdge)->removed()) {
//				TRACE(logger, "Last edge is removed from the graph, hopping out of loop.");
//				break; // we've completely consumed the last edge, stop following it.
//			}
		}
		TRACE(logger, "Finished creating path in graph [" << this->graph->getId() << "]. Final sequence is [" << constructed << "] which is [" << constructed.size() << "] characters long. Followed [" << followed.size() << "] edges, smallest edge was: [" << smallestEdge << "]");

		BOOST_FOREACH (SkinnyGraph::Edge e, followed) {
			this->graph->edge(e)->decreaseWeight(smallestEdge);
		}

		// keep track of the paths that we've followed
		path = verticesToSequenceNodes(verticesFollowed);
		paths.insert(path);

		if (p < 0) {
			// if p is less than 0, then we never encountered a node with more than one edge exiting. If we didn't
			// make any decisions, then we can cleanly remove this node.
			startingPoints.erase(startingPoints.begin());
		}
	}

	return paths;
}

double
ProportionalPathBuilder::getEpsilon() {
	return this->epsilon;
}

void
ProportionalPathBuilder::setEpsilon(double epsilon) {
	this->epsilon = epsilon;
}

double
ProportionalPathBuilder::sumIncomingEdges(SkinnyGraph::Vertex v) {
	double sum = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, getIncomingEdges(v)) {
		sum += this->graph->edge(e)->getWeight();
	}

	return sum;
}

double
ProportionalPathBuilder::sumOutgoingEdges(SkinnyGraph::Vertex v) {
	double sum = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, getOutgoingEdges(v)) {
		sum += this->graph->edge(e)->getWeight();
	}

	return sum;
}

boost::tuple<SkinnyGraph::Edge, double>
ProportionalPathBuilder::findMaxEdge(SkinnyGraph::Vertex v) {
	SkinnyGraph::Edge maxEdge;
	double maxWeight = 0;

	BOOST_FOREACH (SkinnyGraph::Edge e, getOutgoingEdges(v)) {
		double weight = this->graph->edge(e)->getWeight();
		if (weight > maxWeight) {
			maxWeight = weight;
			maxEdge = e;
		}
	}

	return boost::make_tuple(maxEdge, maxWeight);
}

boost::tuple<SkinnyGraph::Edge, bool>
ProportionalPathBuilder::findClosestEdge(SkinnyGraph::Vertex v, double p) {
	SkinnyGraph::Edge closest;
	bool found = false;

	double outgoingSum = sumOutgoingEdges(v);

	BOOST_FOREACH (SkinnyGraph::Edge e, getOutgoingEdges(v)) {
		boost::shared_ptr<WeightedEdge> edge = this->graph->edge(e);
		double edgeP = edge->getWeight() / outgoingSum;

		if (edgeP > p - this->epsilon && edgeP < p + this->epsilon) {
			closest = e;
			found = true;
			break;
		}
	}

	return boost::make_tuple(closest, found);
}

#endif // PROPORTIONAL_PATH_BUILDER_CC
