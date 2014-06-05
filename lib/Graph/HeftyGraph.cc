/* 
 * File:   HeftyGraph.cc
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef HEFTY_GRAPH_CC
#define HEFTY_GRAPH_CC

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>


#include <list>
#include "Graph/HeftyGraph.hh"
#include "Util/Util.hh"
#include "Exception/ReadSizeException.hh"

DECLARE_LOG(logger, "qassembler.HeftyGraph");

HeftyGraph::HeftyGraph(uint16_t kmerLength) {
	this->nextGraphId = 0;
	this->kmerLength = kmerLength;
	this->trackReads = HeftyGraph::DONT_TRACK_READS;
}

HeftyGraph::HeftyGraph(uint16_t kmerLength, HeftyGraph::TrackReads trackReads) {
	this->nextGraphId = 0;
	this->kmerLength = kmerLength;
	this->trackReads = trackReads;
}

HeftyGraph::HeftyGraph(uint16_t kmerLength, HeftyGraph::TrackReads trackReads, boost::shared_ptr<PreHash> guide, std::size_t minEdgeWeight) {
	this->nextGraphId = 0;
	this->kmerLength = kmerLength;
	this->trackReads = trackReads;
	this->guide = guide;
	this->minEdgeWeight = minEdgeWeight;
}

HeftyGraph::~HeftyGraph() {

}

HeftyGraph::ReadLookup
HeftyGraph::getForwardReads() {
	return this->read2graphForward;
}

HeftyGraph::ReadLookup
HeftyGraph::getReverseReads() {
	return this->read2graphReverse;
}

void
HeftyGraph::addReadToGraph(boost::shared_ptr<Sequence> read) {
	std::string sequence = read->getSequence();
	std::string name = read->getName();
	std::size_t id = read->getID();


	if (read->getLength() < kmerLength) {
		throw ReadSizeException("Read is too short.");
	}

	if (guide) {
		TRACE(logger, "Adding read with guide.");
		DEBUG(logger, "Adding read " << name << " (forward) to graph.");
		addReadToGraphWithGuide(sequence, id, name, Kmer::FORWARD);
		sequence = read->getReverseComplement();
		DEBUG(logger, "Adding read " << name << " (reverse) to graph.");
		addReadToGraphWithGuide(sequence, id, name, Kmer::REVERSE);
	} else {
		// if we have no guide, blindly add the read to the graph.
		TRACE(logger, "Adding read without guide.");
		DEBUG(logger, "Adding read " << name << " (forward) to graph.");
		addReadToGraph(sequence, id, name, Kmer::FORWARD);
		sequence = read->getReverseComplement();
		DEBUG(logger, "Adding read " << name << " (reverse complement) read to graph.");
		addReadToGraph(sequence, id, name, Kmer::REVERSE);
		DEBUG(logger, "Finished adding read " << name << " to graph.");
	}
}

void
HeftyGraph::addReadToGraphWithGuide(std::string sequence, std::size_t source, std::string sourceName, Kmer::Strand direction) {
	std::string upperSequence = boost::to_upper_copy(sequence);
	if (upperSequence.size() == kmerLength) {
		if (guide->kmerCount(upperSequence) > minEdgeWeight) {
			addSingleKmerToGraph(upperSequence, source, sourceName, direction);
		} else {
			TRACE (logger, "Not adding single k-mer [" << upperSequence << "]");
		}
	} else {
		// identify k-mer pairs that both have hash counts greater than minEdge weight.
		for (std::size_t i = 0; i < upperSequence.size() - kmerLength; i++) {
			std::string firstKmer = upperSequence.substr(i, kmerLength);
			std::string secondKmer = upperSequence.substr(i + 1, kmerLength);

			if (guide->kmerCount(firstKmer) > minEdgeWeight && guide->kmerCount(secondKmer) > minEdgeWeight) {
				addKmerPairToGraph(std::make_pair(firstKmer, secondKmer), source, sourceName, direction);
				TRACE(logger, "Adding [" << firstKmer << "] and [" << secondKmer << "] to graph as pair with guide.");
			} else if (guide->kmerCount(firstKmer) > minEdgeWeight) {
				addSingleKmerToGraph(firstKmer, source, sourceName, direction);
				TRACE(logger, "Adding [" << firstKmer << "] to graph as single with guide.");
			} else if (guide->kmerCount(secondKmer) > minEdgeWeight) {
				addSingleKmerToGraph(secondKmer, source, sourceName, direction);
				TRACE(logger, "Adding [" << secondKmer << "] to graph as single with guide.");
			} else {
				TRACE (logger, "Not adding kmer pair [" << firstKmer << "] -> [" << secondKmer << "] to graph");
			}
		}
	}
}

void
HeftyGraph::addReadToGraph(std::string sequence, std::size_t source, std::string sourceName, 
			   Kmer::Strand direction) {

	std::string upperSequence = boost::to_upper_copy(sequence);
	TRACE(logger, "Adding [" << sequence << "] to graph (kmerLength = " << kmerLength << ").");
	if (upperSequence.size() == kmerLength) {
		TRACE(logger, "Adding single kmer to graph.");
		addSingleKmerToGraph(upperSequence, source, sourceName, direction);	
	} else {
		for (std::size_t i = 0; i < upperSequence.size() - kmerLength; i++) {
			std::string firstKmer = upperSequence.substr(i, kmerLength);
			std::string secondKmer = upperSequence.substr(i + 1, kmerLength);

			TRACE(logger, "Adding [" << firstKmer << "] and [" << secondKmer << "] to graph.");

			addKmerPairToGraph(std::make_pair(firstKmer, secondKmer), source, sourceName, direction);
		}
	}
}

void
HeftyGraph::addSingleKmerToGraph(std::string kmer, std::size_t source, std::string sourceName, Kmer::Strand direction) {
	std::size_t hash = qassembler::hash(kmer);
	boost::shared_ptr<SequenceNode> node;
	boost::shared_ptr<SkinnyGraph> graph;
	SkinnyGraph::Vertex vertex;

	if (!hashExists(hash)) {
		boost::tie(graph, vertex) = createGraphWithVertex(hash, kmer, sourceName, source, 0, direction);
	} else {
		boost::tie(graph, vertex) = getGraphAndVertexForHash(hash);

		node = graph->node(vertex);
		node->addKmerSourceAt(node->findKmer(hash), source, 0, direction);
	}
	addReference(hash, source, direction, graph);
}

boost::shared_ptr<SkinnyGraph>
HeftyGraph::findOrCreateGraph(std::size_t hash, std::string kmer, std::string sourceName, std::size_t source, std::size_t position,
			      Kmer::Strand direction) {
	boost::shared_ptr<SkinnyGraph> graph;
	boost::shared_ptr<SequenceNode> node;
	SkinnyGraph::Vertex vertex;

	if (!hashExists(hash)) {
		boost::tie(graph, boost::tuples::ignore) = 
			createGraphWithVertex(hash, kmer, sourceName, source, position, direction);
	} else {
		boost::tie(graph, vertex) = getGraphAndVertexForHash(hash);

		node = graph->node(vertex);
		node->addKmerSourceAt(node->findKmer(hash), source, position, direction);
	}

	return graph;
}

void
HeftyGraph::addKmerPairToGraph(std::pair<std::string, std::string> kmers, std::size_t source, std::string sourceName,
			       Kmer::Strand direction) {
	std::size_t hash1, hash2;
	std::size_t kmer1Pos, kmer2Pos;
	SkinnyGraph::Vertex kmer1Vertex, kmer2Vertex;
	boost::shared_ptr<SequenceNode> node;
	boost::shared_ptr<SkinnyGraph> kmer1Graph, kmer2Graph, graph;

	char significant;

	significant = kmers.second[kmers.second.size() - 1];

	hash1 = qassembler::hash(kmers.first);
	hash2 = qassembler::hash(kmers.second);

	// check to see if either of the two k-mers have been created. if any hasn't been previously
	// created, then add a new vertex for that k-mer.
	kmer1Graph = findOrCreateGraph(hash1, kmers.first, sourceName, source, 0, direction);
	kmer2Graph = findOrCreateGraph(hash2, kmers.second, sourceName, source, 0, direction);

	// at this point, the two k-mers must exist. identify whether or not the two k-mers belong to
	// the same graph. If they do not, then merge the two graphs that they belong to.
	if (kmer1Graph != kmer2Graph) {
		graph = mergeGraphs(kmer1Graph, kmer2Graph);
	} else {
		graph = kmer1Graph; // they're both the same, it's irrelevant which one we use
	}

	// 1) identify the two vertices:
	kmer1Vertex = graph->getVertexForHash(hash1);
	kmer2Vertex = graph->getVertexForHash(hash2);

	// 2) identify the location of the hashes within those vertices:
	kmer1Pos = graph->node(kmer1Vertex)->findKmer(hash1);
	kmer2Pos = graph->node(kmer2Vertex)->findKmer(hash2);

	// 2.5) add a reference to the transition that took place between this k-mer pair:
	TRACE(logger, "adding a transition reference from kmer1 [" << kmers.first[kmers.first.size() - 1] << "] to kmer2 [" << significant << "]");
	graph->node(kmer1Vertex)->getKmer(kmer1Pos)->addTransition(significant);

	// if the two kmers occupy the same vertex, then we don't need to add an edge between them and we don't need
	// to add references because those should have been added already.
	if (kmer1Vertex == kmer2Vertex) {
		return;
	}

	// 3) add the edges:
	graph->addEdgeBetweenNodes(kmer1Vertex, kmer2Vertex, kmer1Pos + 1, kmer2Pos);

	// 4) update references to graphs:
	addReference(hash1, source, direction, graph);
	addReference(hash2, source, direction, graph);
}

void
HeftyGraph::addReference(std::size_t hash, std::size_t source, Kmer::Strand direction, boost::shared_ptr<SkinnyGraph> g) {
	if (this->trackReads == HeftyGraph::TRACK_READS) {
		if (direction == Kmer::FORWARD) {
			this->read2graphForward.put(source, g);
		} else {
			this->read2graphReverse.put(source, g);
		}
	}
	this->biGraphs.put(hash, g);
}

boost::tuple<boost::shared_ptr<SkinnyGraph>, SkinnyGraph::Vertex>
HeftyGraph::getGraphAndVertexForHash(std::size_t hash) {
	boost::shared_ptr<SkinnyGraph> graph = getGraphForHash(hash);
	TRACE(logger, "Found [" << std::hex << hash << std::dec << "] in graph [" << std::hex << graph << std::dec << "]");
	SkinnyGraph::Vertex vertex = graph->getVertexForHash(hash);
	return boost::make_tuple(graph, vertex);
}

boost::tuple<boost::shared_ptr<SkinnyGraph>, SkinnyGraph::Vertex>
HeftyGraph::createGraphWithVertex(std::size_t hash, std::string sequence, std::string sourceName, std::size_t sourceId, std::size_t position, Kmer::Strand direction) {
	boost::shared_ptr<SkinnyGraph> g = boost::make_shared<SkinnyGraph>(getNextGraphId());
	SkinnyGraph::Vertex v = g->createFirstSequenceNode(hash, sequence, sourceName, sourceId, position, direction);
	
	// return the tuple
	return boost::make_tuple(g, v);
}

boost::shared_ptr<SkinnyGraph>
HeftyGraph::getGraphForHash(std::size_t hash) {
	return this->biGraphs.get(hash);
}

bool
HeftyGraph::hashExists(std::size_t hash) {
	return this->biGraphs.count(hash) > 0;
}

std::set<boost::shared_ptr<SkinnyGraph> >
HeftyGraph::getGraphs() {
	std::set<boost::shared_ptr<SkinnyGraph> > graphSet;

	BOOST_FOREACH (boost::shared_ptr<SkinnyGraph> g, this->biGraphs.getValues()) {
		graphSet.insert(g);
	}

	return graphSet;
}

int
HeftyGraph::numGraphs() {
	return getGraphs().size();
}

boost::shared_ptr<SkinnyGraph>
HeftyGraph::mergeGraphs(boost::shared_ptr<SkinnyGraph> g1, boost::shared_ptr<SkinnyGraph> g2) {
	boost::shared_ptr<SkinnyGraph> from, to;

	// 1) Figure out which graph has fewer nodes.
	if (g1->numVertices() > g2->numVertices()) {
		to = g1;
		from = g2;
	} else {
		to = g2;
		from = g1;
	}

	TRACE(logger, "Merging [" << std::hex << from << "] into [" << to << std::dec << "]");

	to->merge(from);

	// copy all of the references for hashes to the old graph to the new one
	updateReferences(from, to);
	return to;
}

void
HeftyGraph::updateReferences(boost::shared_ptr<SkinnyGraph> from, boost::shared_ptr<SkinnyGraph> to) {
	TRACE(logger, "Updating references of hashes to graphs.");
	// update hashes:
	boost::unordered_set<std::size_t> keys = this->biGraphs.get(from);
	biGraphs.rehash(keys.size());
	BOOST_FOREACH(std::size_t k, keys) {
		this->biGraphs.put(k, to);
	}
	this->biGraphs.clear(from);

	TRACE(logger, "Updating references of forward read identifiers to graphs.");
	// update forward reads:
	BOOST_FOREACH (std::size_t read, this->read2graphForward.get(from)) {
		this->read2graphForward.put(read, to);
	}
	this->read2graphForward.clear(from);

	TRACE(logger, "Updating references of reverse read identifiers to graphs.");
	// update reverse reads:
	BOOST_FOREACH (std::size_t read, this->read2graphReverse.get(from)) {
		this->read2graphReverse.put(read, to);
	}
	this->read2graphReverse.clear(from);
}

std::size_t
HeftyGraph::getNextGraphId() {
	return this->nextGraphId++;
}

void
HeftyGraph::removeEdgesBelowThreshold(std::size_t threshold) {
	BOOST_FOREACH(boost::shared_ptr<SkinnyGraph> g, this->getGraphs()) {
		g->removeSmallEdges(threshold);
	}
}

void
HeftyGraph::removeGraphsShorterThan(std::size_t threshold) {
	BOOST_FOREACH(boost::shared_ptr<SkinnyGraph> g, this->getGraphs()) {
		if (g->numVertices() == 1) {
			SkinnyGraph::Vertices iterator;
			boost::tie(iterator, boost::tuples::ignore) = g->getVertexIterators();
			boost::shared_ptr<SequenceNode> n = g->node(*iterator);
			if (n->kmerCount() + kmerLength < threshold) {
				this->biGraphs.clear(g);
			}
		}
	}
}

uint16_t
HeftyGraph::getKmerLength() {
	return this->kmerLength;
}

void
HeftyGraph::lockEdgeWeights() {
	BOOST_FOREACH(boost::shared_ptr<SkinnyGraph> g, this->getGraphs()) {
		g->lockEdgeWeights();
	}
}

void
HeftyGraph::resetEdgeWeights() {
	BOOST_FOREACH(boost::shared_ptr<SkinnyGraph> g, this->getGraphs()) {
		g->resetEdgeWeights();
	}
}
#endif // HEFTY_GRAPH
