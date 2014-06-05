/**
 * File:   GraphWriter.hh
 * Author: fbristow
 *
 * Created on January 17, 2012
 */
#ifndef GRAPH_WRITER_HH
#define GRAPH_WRITER_HH

#include <boost/shared_ptr.hpp>
#include "Graph/SkinnyGraph.hh"

class GraphWriter {
public:
	/** constructor with graph and filename specified */
	GraphWriter(boost::shared_ptr<SkinnyGraph> graph, std::string filename, std::string directory);
	/** write the graph out in dot format to the specified filename */
	void write();
private:
	/** a reference to the graph we're going to write */
	boost::shared_ptr<SkinnyGraph> g;
	/** a filename to write to */
	std::string filename;
	/** a directory name to store the files in */
	std::string directory;
};

#endif // GRAPH_WRITER_HH
