/* 
 * File:   HeftyGraph.hh
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef HEFTY_GRAPH_HH
#define HEFTY_GRAPH_HH

#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/cstdint.hpp>

#include "Graph/SkinnyGraph.hh"
#include "Lookup/GraphLookup.hh"
#include "PreHash/PreHash.hh"
#include "Sequence/Sequence.hh"

#include "Logging/Logging.hh"

class HeftyGraph {
public:

	enum TrackReads {
		TRACK_READS,
		DONT_TRACK_READS
	};

	/**
	 * Constructor where kmerLength is specified. Read tracking is disabled by default and
	 * graph is constructed without the assistance of a guide.
	 * @param kmerLength the value to use for k
	 */
	HeftyGraph(uint16_t kmerLength);
	/**
	 * Constructor specifying length and explicitly specifying whether or not to track reads.
	 * Graph is constructed without the assistance of a guide.
	 * @param kmerLength the value to use for k
	 * @param trackReads whether or not the placement of reads should be tracked.
	 */
	HeftyGraph(uint16_t kmerLength, TrackReads trackReads);
	/**
	 * Constructor specifying length, explicitly specifying whether or not to track reads and
	 * supplying a hash count for use in deciding which portions of reads to add to the graph.
	 * @param kmerLength the value to use for k
	 * @param trackReads whether or not the placement of reads should be tracked.
	 * @param guide the hash counter to use as a guide for graph construction
	 * @param minEdgeWeight the minimum edge weight allowed to be constructed for this graph.
	 */
	HeftyGraph(uint16_t kmerLength, TrackReads trackReads, boost::shared_ptr<PreHash> guide, std::size_t minEdgeWeight);
	/**
	 * Destructor
	 */
	~HeftyGraph();
	
	// typedef of a set of paths
	typedef std::set<std::vector<SkinnyGraph::Vertex> > Paths;
	// used for determining which graph can be found in which graph
	typedef GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > ReadLookup;
	// used for determining which hash is found in which graph
	typedef GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > HashLookup;

	/**
	 * Add the supplied read to the graph. 
	 * @param read the read to add to the graph.
	 */
	void addReadToGraph(boost::shared_ptr<Sequence> read);
	/** 
	 * Count the number of sub-graphs in the graph.
	 * @return the number of sub-graphs present in the graph.
	 */
	int numGraphs();
	/** 
	 * Get a collection of the graphs in this graph.
	 * @return references to all graphs found in this graph.
	 */
	std::set<boost::shared_ptr<SkinnyGraph> > getGraphs();
	/**
	 * Get the graphs in this graph by their read identifiers in the forward direction.
	 * @return a key-value store identifying which graph a read was placed in.
	 */
	ReadLookup getForwardReads();
	/**
	 * Get the graphs in this graph by their read identifiers in the reverse direction.
	 * @return a key-value store identifying which graph a read was placed in.
	 */
	ReadLookup getReverseReads();

	/**
	 * Remove edges from all graphs that are below a threshold value.
	 * @param threshold the edge weight threshold.
	 */
	void removeEdgesBelowThreshold(std::size_t threshold);
	/**
	 * Remove graphs with single nodes below a certain length.
	 * @param thresholdLength the minimum number of k-mers required to be present in a node.
	 */
	void removeGraphsShorterThan(std::size_t thresholdLength);
	/**
	 * Get the length of kmer used to construct this graph.
	 * @return kmerLength
	 */
	uint16_t getKmerLength();
	/**
	 * Get both the graph and vertex for a specified hash.
	 * @param hash the hash value to lookup.
	 * @return the graph and vertex for that hash.
	 */
	boost::tuple<boost::shared_ptr<SkinnyGraph>, SkinnyGraph::Vertex> getGraphAndVertexForHash(std::size_t hash);

	/**
	 * Optionally lock a snapshot of all edge weights. Should be called after all reads
	 * have been added to the graph so that subsequent algorithms can destructively modify 
	 * edge weights, then reset them to the locked snapshot.
	 */
	void lockEdgeWeights();

	/**
	 * Reset all edge weights to the locked snapshot.
	 */
	void resetEdgeWeights();
private:
	/** internal graph identifier */
	std::size_t nextGraphId;
	/** bidirectional lookup for graphs */
	HashLookup biGraphs;
	/** a mapping of reads to graphs */
	ReadLookup read2graphForward;
	/** a mapping of reads to graphs */
	ReadLookup read2graphReverse;
	/** the length of k-mers for this graph */
	uint16_t kmerLength;
	/** the minimum edge weight allowed to be created. */
	std::size_t minEdgeWeight;
	/** the guide */
	boost::shared_ptr<PreHash> guide;
	/** should we bother keeping track of where reads are being put? */
	TrackReads trackReads;
	/**
	 * Get the next unique identifier for graphs and increment.
	 * @return the next unique identifier for graphs.
	 */
	std::size_t getNextGraphId();

	/**
	 * Add a read to the graph by manually specifying all components instead of supplying
	 * an AMOS read.
	 * @param sequence the sequence of the read to add
	 * @param source the AMOS identifier where this read came from
	 * @param sourceName the external AMOS identifier where this read came from
	 * @param direction the orientation of the read when adding this sequence
	 */
	void addReadToGraph(std::string sequence, std::size_t source, std::string sourceName, Kmer::Strand direction);

	/**
	 * Add a sequence from a read to the graph using the guide.
	 * @param sequence the sequence of the read to add
	 * @param source the AMOS identifier where this read came from
	 * @param sourceName the external AMOS identifier where this read came from
	 * @param direction the orientation of the read when adding this sequence
	 */
	void addReadToGraphWithGuide(std::string sequence, std::size_t source, std::string sourceName, Kmer::Strand direction);

	/**
	 * Add a pair of overlapping k-mers to the graph.
	 * @param kmers the pair of overlapping k-mers in the pair to add
	 * @param source the AMOS identifier where these kmers came from
	 * @param sourceName the external identifier where these kmers came from
	 * @param direction the orientation of the read where these kmers came from
	 */
	void addKmerPairToGraph(std::pair<std::string, std::string> kmers, std::size_t source, std::string sourceName, 
				Kmer::Strand direction);

	/**
	 * Add a single kmer to the graph (should only be used when read length is equal to kmer length). No edges
	 * are added in this method.
	 * @param kmer the kmer to add to the graph
	 * @param source the AMOS identifier where this kmer came from
	 * @param sourceName the external identifier where this kmer came from
	 * @param direction the orientation of the read where this kmer came from.
	 */
	void addSingleKmerToGraph(std::string kmer, std::size_t source, std::string sourceName, Kmer::Strand direction);

	/** 
	 * Check to see if a hash exists in the current set of graphs.
	 * @param hash the hash to check for.
	 * @return whether or not the hash exists in any graph.
	 */
	bool hashExists(std::size_t hash);
	
	/**
	 * Get the graph that contains this hash value.
	 * @param hash the hash to get a graph for.
	 * @return the graph where this hash resides.
	 */
	boost::shared_ptr<SkinnyGraph> getGraphForHash(std::size_t hash);

	/**
	 * Find an existing graph, or construct one and store the supplied k-mer information.
	 * @param hash the hash value for this kmer
	 * @param kmer the actual kmer sequence
	 * @param sourceName the external identifier for this kmer
	 * @param source the internal AMOS identifier for this kmer
	 * @param location where was this kmer found in the original sequence?
	 * @param direction the orientation of the read that this kmer was constructed from
	 * @return the graph that has the kmer
	 */
	boost::shared_ptr<SkinnyGraph> findOrCreateGraph(std::size_t hash, std::string kmer, std::string sourceName, std::size_t source,
							 std::size_t location, Kmer::Strand direction);
	
	/**
	 * Set references for the hash and read in the specified graph.
	 * @param hash the hash to update references for
	 * @param source the read where this hash came from
	 * @param direction the strand the read was oriented in when generating the hash
	 * @param graph the graph where this hash is stored
	 */
	void addReference(std::size_t hash, std::size_t source, Kmer::Strand direction, boost::shared_ptr<SkinnyGraph> graph);

	/**
	 * Create a new graph containing a vertex with the specified params.
	 * @param hash the hash to put into the sequence node.
	 * @param sequence the sequence to use for the first kmer in the sequence node.
	 * @param sourceName the external AMOS identifier for the source of this kmer
	 * @param sourceId the internal AMOS identifier for the source of this kmer.
	 * @param position where was this sequence found in the source?
	 * @param direction the orientation of the read when this kmer was generated.
	 * @return the graph and vertex that was created to store this kmer.
	 */
	boost::tuple<boost::shared_ptr<SkinnyGraph>, SkinnyGraph::Vertex> 
		createGraphWithVertex(std::size_t hash, std::string sequence, std::string sourceName, std::size_t sourceId,
			              std::size_t position, Kmer::Strand direction);

	/**
	 * Merge two graphs together. Merges the smaller graph (the graph with fewer nodes) into the larger graph
	 * and returns the resulting merged graph.
	 * @param g1 the first graph to merge
	 * @param g2 the second graph to merge
	 * @return a graph containing all the nodes and edges in g1 and g2.
	 */
	boost::shared_ptr<SkinnyGraph> mergeGraphs(boost::shared_ptr<SkinnyGraph> g1, boost::shared_ptr<SkinnyGraph> g2);

	/**
	 * Update all references from the old graph to the new graph for hashes.
	 * @param source the graph where the hashes were previously located
	 * @param the graph where the hashes are currently located
	 */
	void updateReferences(boost::shared_ptr<SkinnyGraph> source, boost::shared_ptr<SkinnyGraph> dest);
};

#endif // HEFTY_GRAPH_HH

