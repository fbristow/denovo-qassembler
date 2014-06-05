/* 
 * File:   GraphWriter.cc
 * Author: fbristow
 *
 * Created on January 17, 2012
 */
#ifndef GRAPH_WRITER_CC
#define GRAPH_WRITER_CC

#include "GraphWriter.hh"
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>

GraphWriter::GraphWriter(boost::shared_ptr<SkinnyGraph> g, std::string filename, std::string directory) {
	this->g = g;
	this->filename = filename;
	this->directory = directory;
}

void
GraphWriter::write() {
	boost::filesystem::create_directory(directory);
	this->filename = directory + "/" + this->filename;
	std::ofstream file(this->filename.c_str());
	// write the graph header:
	file << "digraph G {" << std::endl;
	file << "\trankdir=LR;" << std::endl;

	// write out all of the nodes:
	BOOST_FOREACH(SkinnyGraph::Vertex v, this->g->getVertexIterators()) {
		boost::shared_ptr<SequenceNode> n = this->g->node(v);
		std::size_t id = n->getId(); 
		std::string name = n->getName();
		std::size_t coverageSum = 0;
		BOOST_FOREACH(boost::shared_ptr<Kmer> k, n->getKmers()) {
			coverageSum += k->getCount();
		}
		std::size_t kmers = n->kmerCount();
		file << id << " [label=\"" << name << ": kmers(" << kmers << "), avg coverage("<< coverageSum / (double)kmers << ")\"];" << std::endl; 
/*		std::string sequence;
	       if (boost::in_degree(v, *this->g->graph()) == 0) {
			sequence = n->fullSequence();
	       } else {
		       sequence = n->sequence();
	       }

		file << id << " [label=\"" << sequence << "\"];" << std::endl;*/
	}

	// write out all of the edges:
	BOOST_FOREACH(SkinnyGraph::Edge e, this->g->edges()) {
		boost::shared_ptr<WeightedEdge> w = this->g->edge(e);
		std::size_t weight = w->getWeight();
		SkinnyGraph::Vertex source, target;
		std::size_t sourceId, targetId;
		source = boost::source(e, *this->g->graph());
		target = boost::target(e, *this->g->graph());
		sourceId = this->g->node(source)->getId();
		targetId = this->g->node(target)->getId();

		file << sourceId << "->" << targetId << " [label=\"" << weight << "\"];" << std::endl;
	}

	// write the graph footer:
	file << "}" << std::endl;
}

#endif // GRAPH_WRITER_CC
