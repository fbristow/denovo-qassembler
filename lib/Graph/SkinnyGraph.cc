/* 
 * File:   SkinnyGraph.cc
 * Author: fbristow
 *
 * Created on January 5, 2012
 */
#ifndef SKINNY_GRAPH_CC
#define SKINNY_GRAPH_CC

#define DEFAULT_HASH_TABLE_SIZE 10000

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include "SkinnyGraph.hh"
#include "Kmer/FirstKmer.hh"
#include "Exception/InvalidGraphStateException.hh"

DECLARE_LOG(logger, "qassembler.SkinnyGraph");

SkinnyGraph::SkinnyGraph(std::size_t identifier) : g(boost::make_shared<Graph>()){
	this->id = identifier;
	this->nextVertexId = 0;
	this->hash2vertex.rehash(DEFAULT_HASH_TABLE_SIZE);
}

SkinnyGraph::~SkinnyGraph() {
	this->hash2vertex.clear();
}

std::size_t
SkinnyGraph::getId() {
	return this->id;
}

SkinnyGraph::Vertex
SkinnyGraph::getVertexForHash(std::size_t hash) {
	return this->hash2vertex[hash];
}

void
SkinnyGraph::setVertexForHash(std::size_t hash, Vertex v) {
	this->hash2vertex[hash] = v;
}

SkinnyGraph::Vertex
SkinnyGraph::createSequenceNode(std::size_t hash, char nucleotide, std::string sourceName, std::size_t sourceId, std::size_t position, Kmer::Strand direction) {
	SkinnyGraph::Vertex v;
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>(getNextVertexId(), sourceName);

	n->addKmer(hash, nucleotide, sourceId, position, direction);
	v = boost::add_vertex(n, *this->g);
	hash2vertex[hash] = v;

	return v;
}

SkinnyGraph::Vertex
SkinnyGraph::createSequenceNode(boost::shared_ptr<SequenceNode> node) {
	SkinnyGraph::Vertex v;
	std::vector<boost::shared_ptr<Kmer> > kmers;

	// add a new vertex with a reference to this sequence node
	v = boost::add_vertex(node, *this->g);
	// we can't let the sequence node use the old identifier that it has, it probably collides with
	// one of the vertices in the graph that we have already.
	node->setId(getNextVertexId());

	// copy all of the hashes contained in this sequence node and update
	// all of the references
	kmers = node->getKmers();

	TRACE(logger, "Copying kmer [" << kmers.size() << "] references.");
	for (std::size_t i = 0; i < kmers.size(); i++) {
		TRACE(logger, "Copying hash [" << std::hex << kmers[i]->getHash() << std::dec << "]");
		hash2vertex[kmers[i]->getHash()] = v;
	}
	TRACE(logger, "Finished copying kmer references.");

	return v;
}

SkinnyGraph::Vertex
SkinnyGraph::createFirstSequenceNode(std::size_t hash, std::string sequence, std::string sourceName, std::size_t sourceId, std::size_t position, Kmer::Strand direction) {
	SkinnyGraph::Vertex v;
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>(getNextVertexId(), sourceName);
	boost::shared_ptr<FirstKmer> mer = boost::make_shared<FirstKmer>(hash, sequence, sourceId, position, direction);

	n->addKmer(mer);
	v = boost::add_vertex(n, *this->g);
	hash2vertex[hash] = v;

	return v;
}


boost::tuple<SkinnyGraph::Vertex, SkinnyGraph::Vertex>
SkinnyGraph::split(SkinnyGraph::Vertex v, std::size_t position) {
	TRACE(logger, "checking whether or not we actually need to split a vertex.");
	if (position == 0 || position == node(v)->kmerCount()) {
		return boost::make_tuple(v, v);
	}

	TRACE(logger, "going to split a vertex.");
	boost::shared_ptr<WeightedEdge> backingEdge = boost::make_shared<WeightedEdge>();
	// assign or create new nodes for the two parts of the node to be split
	TRACE(logger, "getting the front half.");
	SkinnyGraph::Vertex frontHalf = getFrontDestination(v);
	TRACE(logger, "getting the back half.");
	SkinnyGraph::Vertex backHalf = getBackDestination(v);
	TRACE(logger, "copying the k-mer count as the new edge weight.");
	TRACE(logger, "\tk-mer count: [" << node(v)->kmerCount() << "], requested position: [" << std::dec << position << "].");
	boost::shared_ptr<Kmer> k = node(v)->getKmer(position);
	TRACE(logger, "copying from k-mer [" << std::hex << k << std::dec << "].");
	// don't want a count here, rather we want to know how many times the kmer at this position transitioned to it's neighbour
	// so find out what the most significant base of the neighbour is, then find out how many times the transition was made
	boost::shared_ptr<Kmer> prev = node(v)->getKmer(position - 1);
	std::size_t count = prev->getTransitionCount(k->getBase());
	TRACE(logger, "transition count from prev with base [" << k->getBase() << "] is [" << count << "]");
	// copy the number of k-mers at the place we're splitting as the new edge weight
	backingEdge->setWeight(count);
	// add an edge between the two nodes that we're working with
	TRACE(logger, "adding an edge between the two halves. [" << std::hex << frontHalf << std::dec << "] and [" << std::hex << backHalf << std::dec << "] with weight [" << count << "]");
	boost::add_edge(frontHalf, backHalf, backingEdge, *g);

	TRACE(logger, "copying the k-mers to each respective half.");
	// move the k-mers from the node we're looking at to each respective node
	for (std::size_t i = 0; i < node(v)->kmerCount(); i++) {
		boost::shared_ptr<Kmer> mer = node(v)->getKmer(i);
		SkinnyGraph::Vertex append;
		if (i < position) {
			append = frontHalf;
			node(append)->addKmer(mer);
		} else {
			append = backHalf;
			node(append)->addKmerAt(mer, i - position);
		}
		TRACE(logger, "Copying: [" << mer->getBase() << "] from [" << std::hex << v << std::dec << "] to [" << std::hex << append << std::dec << "]");
		setVertexForHash(mer->getHash(), append);
	}

	TRACE(logger, "removing the vertex from the graph.");
	// remove the original vertex from the graph
	boost::clear_vertex(v, *this->g);
	boost::remove_vertex(v, *this->g);

	TRACE(logger, "Front half has sequence [" << node(frontHalf)->fullSequence() << "], back half has sequence [" << node(backHalf)->sequence() << "]");

	return boost::make_tuple(frontHalf, backHalf);
}

SkinnyGraph::Vertex
SkinnyGraph::getFrontDestination(SkinnyGraph::Vertex vertex) {
	SkinnyGraph::Vertex frontHalf = NULL;
	// we can reduce the necessity of creating and adding new nodes to the graph
	// by checking to see how many incoming neighbours and how many outgoing
	// neighbours the current node has. If the node has one incoming neighbour
	// AND this node is the only neighbour of the incoming vertex, then we're
	// just going to move the leading half of the vertex we're splitting to the
	// end of the single incoming neighbour.
	if (boost::in_degree(vertex, *g) == 1) {
		// we have only one incoming neighbour, make sure that we're the only
		// neighbour of that node
		SkinnyGraph::IncomingEdges in;
		boost::tie(in, boost::tuples::ignore) = boost::in_edges(vertex, *g);
		SkinnyGraph::Vertex source = boost::source(*in, *g);
		if (boost::out_degree(source, *g) == 1) {
			TRACE(logger, "Front destination neighbour only has one outgoing edge, going to merge nodes.");
			// we're the only neighbour of our neighbour, we're going to
			// dump all of our stuff into their house
			frontHalf = source;
		}
	}

	// in the event that we can't keep all of our stuff at our neighbours house
	// (the above constraints aren't true), then we should add a new node to the
	// graph and copy all our incoming edges to the new node:
	if (frontHalf == NULL) {
		TRACE(logger, "Constraints not satisfied for merging nodes, creating new node and copying edges.");
		// 1) add the new node:
		boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>(getNextVertexId(), node(vertex)->getName());
		frontHalf = boost::add_vertex(n, *g);
		// 2) copy all of the incoming edges to the new node:
		SkinnyGraph::IncomingEdges vIn, vInEnd;
		boost::tie(vIn, vInEnd) = boost::in_edges(vertex, *g);
		for (; vIn != vInEnd; vIn++) {
			boost::shared_ptr<WeightedEdge> backingEdge = edge(*vIn);
			SkinnyGraph::Vertex source = boost::source(*vIn, *g);
			boost::add_edge(source, frontHalf, backingEdge, *g);
		}
	}

	return frontHalf;
}

SkinnyGraph::Vertex
SkinnyGraph::getBackDestination(SkinnyGraph::Vertex vertex) {
	SkinnyGraph::Vertex backHalf = NULL;
	// exactly the same idea as getFrontDestination, except we're applying the concept
	// to the outgoing edges of the vertex:
	if (boost::out_degree(vertex, *g) == 1) {
		// we only have one outgoing neighbour. Check to see that we're the only
		// incoming neighbour of that node:
		SkinnyGraph::OutgoingEdges out;
		boost::tie(out, boost::tuples::ignore) = boost::out_edges(vertex, *g);
		SkinnyGraph::Vertex target = boost::target(*out, *g);

		if (boost::in_degree(target, *g) == 1) {
			TRACE(logger, "Front destination neighbour only has one outgoing edge, going to merge nodes.");
			backHalf = target;
		}
	}

	// Again, if we can't stuff the second half into its neighbour, create a new node
	// and copy all of the outgoing edges from vertex into the new node:
	if (backHalf == NULL) {
		TRACE(logger, "Constraints not satisfied for merging nodes, creating new node and copying edges.");
		// 1) add the new vertex:
		boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>(getNextVertexId(), node(vertex)->getName());
		backHalf = boost::add_vertex(n, *g);
		// 2) copy all of the outgoing edges to the new vertex:
		SkinnyGraph::OutgoingEdges vOut, vOutEnd;
		boost::tie(vOut, vOutEnd) = boost::out_edges(vertex, *g);
		for (; vOut != vOutEnd; vOut++) {
			boost::shared_ptr<WeightedEdge> backingEdge = edge(*vOut);
			SkinnyGraph::Vertex target = boost::target(*vOut, *g);
			boost::add_edge(backHalf, target, backingEdge, *g);
		}
	}

	return backHalf;
}

boost::shared_ptr<WeightedEdge>
SkinnyGraph::edge(SkinnyGraph::Edge e) {
	return (*this->g)[e];
}

boost::shared_ptr<SequenceNode>
SkinnyGraph::node(SkinnyGraph::Vertex v) {
	return (*this->g)[v];
}

boost::shared_ptr<SkinnyGraph::Graph>
SkinnyGraph::graph() {
	return this->g;
}

std::size_t
SkinnyGraph::getNextVertexId() {
	nextVertexId++;
	return nextVertexId;
}

std::size_t
SkinnyGraph::numVertices() {
	return boost::num_vertices(*this->g);
}

boost::unordered_map<std::size_t, SkinnyGraph::Vertex>
SkinnyGraph::getVertices() {
	return this->hash2vertex;
}

SkinnyGraph::Vertex
SkinnyGraph::addEdge(SkinnyGraph::Vertex source, SkinnyGraph::Vertex dest) {
	// otherwise, we have to add an edge between the two nodes, so add it:
	boost::shared_ptr<WeightedEdge> backingEdge = boost::make_shared<WeightedEdge>();

	return addEdge(source, dest, backingEdge);
}

SkinnyGraph::Vertex
SkinnyGraph::addEdge(SkinnyGraph::Vertex source, SkinnyGraph::Vertex dest, boost::shared_ptr<WeightedEdge> backingEdge) {
	SkinnyGraph::Edge e;
	bool added;

	// add an edge between the two nodes
	boost::tie(e, added) = boost::add_edge(source, dest, backingEdge, *this->g);

	// boost will say that the edge wasn't added if the edge already exists
	if (!added) {
		// so if the edge wasn't added, increment the weight of the edge
		edge(e)->increaseWeight(1);
	}

	return dest;
}

SkinnyGraph::Vertex
SkinnyGraph::addEdgeOrMerge(SkinnyGraph::Vertex source, SkinnyGraph::Vertex dest) {
	SkinnyGraph::Vertex next;
	// find out if these two vertices are candidates for merging. Two nodes can be merged
	// if source has one outgoing edge and dest has one incoming edge.
	bool neighbours = (boost::out_degree(source, *this->g) == 1 && boost::in_degree(dest, *this->g) == 1);
	bool noNeighbours = (boost::out_degree(source, *this->g) == 0 && boost::in_degree(dest, *this->g) == 0);
	if (neighbours) {
		// make sure that the outgoing/incoming edge that these nodes have is shared as an edge:
		boost::tie(boost::tuples::ignore, neighbours) = boost::edge(source, dest, *this->g);
	}
	if (neighbours || noNeighbours) { 
		TRACE(logger, "Conditions are valid to merge source and dest, going to merge.");
		// merge the k-mers from source into dest (the merge operation will insert
		// kmers at the beginning of the node).
		node(dest)->merge(node(source));
		// redirect all hashes that pointed at source to point at dest:
		BOOST_FOREACH(boost::shared_ptr<Kmer> k, node(source)->getKmers()) {
			setVertexForHash(k->getHash(), dest);
		}
		// copy the incoming edges from source:
		BOOST_FOREACH(SkinnyGraph::Edge e, boost::in_edges(source, *this->g)) {
			SkinnyGraph::Vertex edgeSource = boost::source(e, *this->g);
			boost::add_edge(edgeSource, dest, edge(e), *this->g);
		}
		// now clear that node and remove the node from the graph:
		boost::clear_vertex(source, *this->g);
		boost::remove_vertex(source, *this->g);

		next = dest;
	} else {
		TRACE(logger, "Cannot merge nodes, going to add an edge.");
		// the two nodes are not candidates for merging. Just add an edge between them
		next = addEdge(source, dest);
	}

	return next;
}

void
SkinnyGraph::addEdgeBetweenNodes(SkinnyGraph::Vertex source, SkinnyGraph::Vertex dest, std::size_t sourcePos, std::size_t destPos) {
	SkinnyGraph::Vertex frontSource, backDest;

	boost::tie(frontSource, boost::tuples::ignore) = split(source, sourcePos);
	boost::tie(boost::tuples::ignore, backDest) = split(dest, destPos);

	addEdgeOrMerge(frontSource, backDest);
}

std::pair<SkinnyGraph::Edges, SkinnyGraph::Edges>
SkinnyGraph::edges() {
	return boost::edges(*this->g);
}

std::size_t
SkinnyGraph::numEdges() {
	return boost::num_edges(*this->g);
}

void
SkinnyGraph::merge(boost::shared_ptr<SkinnyGraph> from) {
	if (from->numEdges() > 0) {
		boost::unordered_map<SkinnyGraph::Vertex, SkinnyGraph::Vertex> oldToNew;
		// if the graph that we're copying from has more than one edge, then we should just copy
		// source and destinations of edges from that graph
		TRACE(logger, "The graph to merge from has edges, copying all nodes and edges.");
		BOOST_FOREACH(SkinnyGraph::Edge e, from->edges()) {
			SkinnyGraph::Vertex source, target;
			boost::shared_ptr<SequenceNode> nodeSource, nodeTarget;
			boost::shared_ptr<WeightedEdge> backingEdge = from->edge(e);

			source = boost::source(e, *from->graph());
			target = boost::target(e, *from->graph());
			
			nodeSource = from->node(source);
			nodeTarget = from->node(target);

			// if we haven't already added this node to the graph, then add it
			if (oldToNew.count(source) == 0) {
				oldToNew[source] = createSequenceNode(nodeSource);
			}

			if (oldToNew.count(target) == 0) {
				oldToNew[target] = createSequenceNode(nodeTarget);
			}

			// finally, add an edge between those two nodes and set the weight:
			addEdge(oldToNew[source], oldToNew[target], backingEdge);
		}
	} else {
		// otherwise, just copy the one vertex that's here
		TRACE(logger, "The graph to merge from has no edges, going to copy one vertex.");
		if (from->numVertices() == 1) {
			SkinnyGraph::Vertex v = from->getVertices().begin()->second;
			boost::shared_ptr<SequenceNode> vNode = from->node(v);
			TRACE(logger, "source node has [" << from->node(v)->getKmers().size() << "] kmers.");
			TRACE(logger, "source node has id [" << vNode->getId() << "]");
			createSequenceNode(vNode);
		} else {
			TRACE(logger, "The graph to merge from has no edges, but more than one vertex. Invalid state, bailing.");
			throw InvalidGraphStateException("The graph to merge from has no edges, but more than one vertex.");
		}
	}
}

std::pair<SkinnyGraph::Vertices, SkinnyGraph::Vertices>
SkinnyGraph::getVertexIterators() {
	return boost::vertices(*this->graph());
}

std::size_t
SkinnyGraph::removeSmallEdges(std::size_t weightThreshold) {
	std::size_t edgesBeforeRemoval = boost::num_edges(*this->graph());
	SkinnyGraph::EdgeRemovalFilter filter(*this->graph(), weightThreshold);
	TRACE(logger, "Removing edges with weight below threshold.");
	boost::remove_edge_if(filter, *this->graph());
	//TRACE(logger, "Removing disconnected vertices.");
	// now remove any disconnected vertices:
	//BOOST_FOREACH (SkinnyGraph::Vertex v, getVertexIterators()) {
	//	if (boost::in_degree(v, *this->graph()) == 0 && boost::out_degree(v, *this->graph()) == 0) {
	//		boost::remove_vertex(v, *this->graph());
	//	}
	//}

	return edgesBeforeRemoval - boost::num_edges(*this->graph());
}

void
SkinnyGraph::lockEdgeWeights() {
	BOOST_FOREACH(SkinnyGraph::Edge e, this->edges()) {
		edge(e)->lockWeight();
	}
}

void
SkinnyGraph::resetEdgeWeights() {
	BOOST_FOREACH(SkinnyGraph::Edge e, this->edges()) {
		edge(e)->resetWeight();
	}
}

#endif // SKINNY_GRAPH_CC
