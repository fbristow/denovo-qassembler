/* 
 * File:   PathBuilder.hh
 * Author: fbristow
 *
 * Created on February 8, 2012
 */
#ifndef PATH_BUILDER_HH
#define PATH_BUILDER_HH

#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>

#include "Graph/SkinnyGraph.hh"
#include "Graph/Node/SequenceNode.hh"
#include "Logging/Logging.hh"

class PathBuilder {
public:
	/** to be used as a way to describe a path through the graph */
	typedef std::vector<boost::shared_ptr<SequenceNode> > Path;

	/** get the graph that this path builder uses to construct paths */
	boost::shared_ptr<SkinnyGraph> getGraph();
	/** set the graph that this path builder uses to construct paths */
	void setGraph(boost::shared_ptr<SkinnyGraph> graph);

	/** 
	 * construct some set of paths from the graph.
	 * this is a pure-virtual function, sub-classes must
	 * implement this method, the notation for this in C++ is
	 * to append '= 0' to the end of a method definition
	 */
	virtual boost::unordered_set<Path> buildPaths() = 0;
protected:
	/** Default constructor. */
	PathBuilder(boost::shared_ptr<SkinnyGraph> graph);
	/** Default destructor. */
	~PathBuilder();
	/**
	 * Get the set of starting points to use for constructing paths
	 * (i.e., the vertices that have no incoming edges).
	 * @return the set of vertices with no incoming edges.
	 */
	std::vector<SkinnyGraph::Vertex> getStartingPoints();

	/**
	 * Convert an ordered list of vertices into an ordered list of
	 * sequence nodes that an implementation can use to construct
	 * proper sequences.
	 * @param vertices the set of vertices to convert
	 * @return an ordered list of sequence nodes represented by the vertices
	 */
	Path verticesToSequenceNodes(std::vector<SkinnyGraph::Vertex> vertices);

	/**
	 * Get all edges outgoing from a vertex which have not yet been completely
	 * consumed by previous path construction.
	 * @param vertex the vertex to retrieve outgoing edges for.
	 * @return the set of outgoing edges for that vertex.
	 */
	std::vector<SkinnyGraph::Edge> getOutgoingEdges(SkinnyGraph::Vertex vertex);

	/**
	 * Get all edges incoming to a vertex which have not yet been completely
	 * consumed by previous path construction.
	 * @param vertex the vertex to retrieve incoming edges for.
	 * @return the set of incoming edges for that vertex.
	 */
	std::vector<SkinnyGraph::Edge> getIncomingEdges(SkinnyGraph::Vertex vertex);

	/** the graph that we'll search for paths in */
	boost::shared_ptr<SkinnyGraph> graph;
};

bool compareVertexIdentifiers(SkinnyGraph::Vertex v1, SkinnyGraph::Vertex v2);

#endif // PATH_BUILDER_HH
