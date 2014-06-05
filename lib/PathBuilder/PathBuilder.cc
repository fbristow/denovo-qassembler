/* 
 * File:   PathBuilder.cc
 * Author: fbristow
 *
 * Created on February 8, 2012
 */
#ifndef PATH_BUILDER_CC
#define PATH_BUILDER_CC

#include <boost/foreach.hpp>
#include <algorithm>

#include "PathBuilder.hh"

DECLARE_LOG(logger, "qassembler.PathBuilder");

PathBuilder::PathBuilder(boost::shared_ptr<SkinnyGraph> graph) {
	this->graph = graph;
}

PathBuilder::~PathBuilder() {
	this->graph->resetEdgeWeights();
}

boost::shared_ptr<SkinnyGraph>
PathBuilder::getGraph() {
	return this->graph;
}

void
PathBuilder::setGraph(boost::shared_ptr<SkinnyGraph> graph) {
	this->graph = graph;
}

std::vector<SkinnyGraph::Vertex>
PathBuilder::getStartingPoints() {
	std::vector<SkinnyGraph::Vertex> startingPoints;

	BOOST_FOREACH (SkinnyGraph::Vertex v, this->graph->getVertexIterators()) {
		if (boost::in_degree(v, *this->graph->graph()) == 0) {
			startingPoints.push_back(v);
		}
	}

	std::sort(startingPoints.begin(), startingPoints.end(), compareVertexIdentifiers);
	return startingPoints;
}

PathBuilder::Path
PathBuilder::verticesToSequenceNodes(std::vector<SkinnyGraph::Vertex> vertices) {
	PathBuilder::Path p;

	BOOST_FOREACH (SkinnyGraph::Vertex v, vertices) {
		p.push_back(this->graph->node(v));
	}

	return p;
}

std::vector<SkinnyGraph::Edge>
PathBuilder::getOutgoingEdges(SkinnyGraph::Vertex vertex) {
	std::vector<SkinnyGraph::Edge> outgoingEdges;

	BOOST_FOREACH (SkinnyGraph::Edge e, boost::out_edges(vertex, *this->graph->graph())) {
		if (!this->graph->edge(e)->removed()) {
			outgoingEdges.push_back(e);
		}
	}

	return outgoingEdges;
}

std::vector<SkinnyGraph::Edge>
PathBuilder::getIncomingEdges(SkinnyGraph::Vertex vertex) {
	std::vector<SkinnyGraph::Edge> incomingEdges;

	BOOST_FOREACH (SkinnyGraph::Edge e, boost::in_edges(vertex, *this->graph->graph())) {
		TRACE (logger, "Current edge weight: [" << this->graph->edge(e)->getWeight() << "]");
		if (!this->graph->edge(e)->removed()) {
			TRACE (logger, "adding edge to collection");
			incomingEdges.push_back(e);
		}
	}

	return incomingEdges;
}

bool
compareVertexIdentifiers(SkinnyGraph::Vertex v1, SkinnyGraph::Vertex v2) {
	return v1 < v2;
}

#endif // PATH_BUILDER_CC
