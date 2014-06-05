/* 
 * File:   SkinnyGraph.hh
 * Author: fbristow
 *
 * Created on January 5, 2012
 */
#ifndef SKINNY_GRAPH_HH 
#define SKINNY_GRAPH_HH

#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>

#include "Graph/Edge/WeightedEdge.hh"
#include "Graph/Node/SequenceNode.hh"

#include "Logging/Logging.hh"

class SkinnyGraph {
public:
	/** constructor, specifying an identifier */
	SkinnyGraph(std::size_t);
	/** destructor */
	~SkinnyGraph();

	/** get the graph identifier */
	std::size_t getId();

	// define a graph as an adjacency list that uses vectors to store edges and vertices, is directed, and uses
	// the definitions for vertex and edge as described above
	typedef boost::adjacency_list<boost::setS, // use a Set to store the edge list
				      boost::setS, // use a Set to store the vertex list
				      boost::bidirectionalS, // I have this set as bidirectional, but I think this could be just directedS instead.
				      boost::shared_ptr<SequenceNode>, // vertex properties are represented by a SequenceNode.
				      boost::shared_ptr<WeightedEdge> // edge properties are represented by a WeightedEdge.
				      > Graph;
	// get a type from the graph defined above for vertex
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	// get a type from the graph defined above for edge
	typedef boost::graph_traits<Graph>::edge_descriptor Edge;
	// typedef an iterator for incoming edges to a node
	typedef boost::graph_traits<Graph>::in_edge_iterator IncomingEdges;
	// typedef an iterator for outgoing edges from a node
	typedef boost::graph_traits<Graph>::out_edge_iterator OutgoingEdges;
	// typedef an iterator for vertices in the graph
	typedef boost::graph_traits<Graph>::vertex_iterator Vertices;
	// typedef an iterator for edges in the graph
	typedef boost::graph_traits<Graph>::edge_iterator Edges;

	/** 
	 * Add a new vertex to the graph with a sequence node.
	 * @param hash the hash of the first kmer added to this node.
	 * @param base the most significant nucleotide of the first kmer added to this node.
	 * @param sourceName the name of the read for the first kmer added to this node.
	 * @param position the position in the read where the first kmer was generated.
	 * @param direction the strand from the read where the first kmer was generated.
	 * @return the Vertex descriptor for the new sequence node.
	 */
	Vertex createSequenceNode(std::size_t hash, char base, std::string sourceName,
		       std::size_t sourceId, std::size_t position, Kmer::Strand direction);
	/** 
	 * Add a new vertex to the graph with the existing sequence node.
	 * @param node the node to copy kmers from
	 * @return the Vertex descriptor for the new node in the graph.
	 */
	Vertex createSequenceNode(boost::shared_ptr<SequenceNode> node);
	/** 
	 * Add a new vertex that holds the first k-mer.
	 * @param hash the hash of the kmer added to this node.
	 * @param sequence the complete sequence of the kmer added to this node.
	 * @param sourceName the name of the read for the kmer added to this node.
	 * @param sourceId the identifier of the read for the kmer added to this node.
	 * @param position the place in the read where the kmer was found.
	 * @param strand the direction of the read when this kmer was generated.
	 */
	Vertex createFirstSequenceNode(std::size_t hash, std::string sequence, std::string sourceName, 
			std::size_t sourceId, std::size_t position, Kmer::Strand strand);
	/** find out which vertex belongs to the specified hash */
	Vertex getVertexForHash(std::size_t);
	/** set which vertex belongs to the specified hash */
	void setVertexForHash(std::size_t, Vertex);
	/** split the specified vertex at a particular position */
	boost::tuple<Vertex, Vertex> split(Vertex, std::size_t);
	/** get a weighted edge for the specified edge */
	boost::shared_ptr<WeightedEdge> edge(Edge);
	/** get a sequence node for the specified vertex */
	boost::shared_ptr<SequenceNode> node(Vertex);
	/** get the backing graph for this skinny graph */
	boost::shared_ptr<Graph> graph();
	/** get all the edges from the graph */
	std::pair<Edges, Edges> edges();
	/** how many vertices are in this graph? */
	std::size_t numVertices();
	/** how many edges are in this graph? */
	std::size_t numEdges();
	/** get the vertices in this graph */
	boost::unordered_map<std::size_t, Vertex> getVertices();
	/** get a unique list of vertices in this graph */
	std::pair<Vertices, Vertices> getVertexIterators();
	/** add an edge between two vertices in this graph */
	Vertex addEdge(Vertex, Vertex);
	/** add an edge between two vertices in this graph (or merge if the two vertices share a single edge) */
	Vertex addEdgeOrMerge(Vertex, Vertex);
	/** add an edge between two vertices in this graph (using an existing weighted edge) */
	Vertex addEdge(Vertex, Vertex, boost::shared_ptr<WeightedEdge>);
	/** add an edge between two existing nodes, may result in splitting or merging of existing nodes */
	void addEdgeBetweenNodes(Vertex source, Vertex dest, std::size_t sourcePos, std::size_t destPos);
	/** merge the nodes an edges from the passed graph into this one */
	void merge(boost::shared_ptr<SkinnyGraph>);
	/** remove edges below a certain threshold */
	std::size_t removeSmallEdges(std::size_t);

	/**
	 * Lock a snapshot of all edge weights currently in the graph.
	 */
	void lockEdgeWeights();

	/**
	 * Reset all edge weights to the locked snapshot.
	 */
	void resetEdgeWeights();
private:
	/** default constructor (shouldn't be called directly) */
	SkinnyGraph();
	/** a unique identifier for this graph object */
	std::size_t id;
	/** a reference to a boost graph */
	boost::shared_ptr<Graph> g;
	/** references from hash values to vertices */
	boost::unordered_map<std::size_t, Vertex> hash2vertex;
	/** the next vertex id to use when creating a new vertex */
	std::size_t nextVertexId;
	/** get the next vertex identifier and increment */
	std::size_t getNextVertexId();

	/** get the front half of a node when splitting */
	Vertex getFrontDestination(Vertex);
	/** get the back half of a node for splitting */
	Vertex getBackDestination(Vertex);

	/** a filter for removing edges from the graph */
	struct EdgeRemovalFilter {
		EdgeRemovalFilter(Graph &g_, std::size_t weight_) : g(g_), weight(weight_) {}
		template <class Edge>
		bool operator() (const Edge e) {
			return g[e]->getWeight() <= weight;
		}
		
		Graph& g;
		std::size_t weight;
	};
};

#endif // SKINNY_GRAPH_HH
