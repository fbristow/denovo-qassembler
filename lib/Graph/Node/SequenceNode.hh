/* 
 * File:   SequenceNode.hh
 * Author: fbristow
 *
 * Created on December 14, 2011
 */
#ifndef SEQUENCE_NODE_HH
#define SEQUENCE_NODE_HH

#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <vector>

#include "Kmer/Kmer.hh"

class SequenceNode {
public:
	/**
	 * Default constructor.
	 */
	SequenceNode();
	/** 
	 * Copy constructor.
	 * @param copy the copy to use as a template.
	 */
	SequenceNode(SequenceNode *copy);
	/**
	 * constructor with values for all members.
	 * @param id the unique identifier to use for this sequence node.
	 * @param name the human-readable name for this sequence node.
	 */
	SequenceNode(std::size_t id, std::string name);
	/**
	 * Destructor.
	 */
	~SequenceNode();
		
	/**
	 * Merge the kmers from the supplied node into this one.
	 * @param source the node from which we can get kmers.
	 */
	void merge(boost::shared_ptr<SequenceNode> source);
	/**
	 * Get the sequence represented by this node. Does not include the complete sequence
	 * for the first kmer.
	 * @return the sequence represented by this node (without the complete first kmer).
	 */
	std::string sequence();
	/**
	 * Get the full sequence (including the first k-mer) represented by this node.
	 * @return the complete sequence represented by this node (including the complete first kmer).
	 */
	std::string fullSequence();
	/**
	 * Get the name for this node.
	 * @return the name for this node.
	 */
	std::string getName();
	/**
	 * Set the name for this node.
	 * @param name the new name for this node.
	 */
	void setName(std::string name);
	/**
	 * Get the identifier for this node.
	 * @return the identifier for this node.
	 */
	std::size_t getId();
	/**
	 * Set the identifier for this node.
	 * @param id the new identifier for this node.
	 */
	void setId(std::size_t id);
	
	/**
	 * Get the kmers for this node.
	 * @return the list of kmers for this node.
	 */
	std::vector<boost::shared_ptr<Kmer> > getKmers();
	/**
	 * Get a single k-mer at a position.
	 * @param position the requested position.
	 * @return the kmer at position.
	 */
	boost::shared_ptr<Kmer> getKmer(std::size_t position);
	/** 
	 * Set the kmers for this node.
	 * @param kmers the new list of kmers to use for this node.
	 */
	void setKmers(std::vector<boost::shared_ptr<Kmer> > kmers);
	/**
	 * Add a new kmer to this node by specifying the values.
	 * @param hash the hashed value for this kmer.
	 * @param base the most-significant base for this kmer.
	 * @param source the identifier of the read where this kmer came from.
	 * @param offset the position in the read where this kmer came from.
	 * @param strand the direction we were iterating over windows in the read when we generated this kmer.
	 */
	void addKmer(std::size_t hash, char base, std::size_t source, std::size_t offset, Kmer::Strand strand);
	/**
	 * Add a new kmer to the end of the lsit of kmers in this node by passing a kmer object.
	 * @param mer the kmer to add to this sequence node.
	 */
	void addKmer(boost::shared_ptr<Kmer> mer);
	/**
	 * Add a new kmer to this node at some specific location.
	 * @param mer the kmer to add to this sequence node.
	 * @param position the position that this kmer should be added to the node.
	 */
	void addKmerAt(boost::shared_ptr<Kmer> mer, std::size_t position);
	/**
	 * Add a new source for an existing kmer.
	 * @param position the kmer position that we're adding a new source to.
	 * @param source the read identifier that this kmer came from.
	 * @param offset the location in the read where this kmer came from.
	 * @param strand the direction of the read when we generated this kmer.
	 */
	void addKmerSourceAt(std::size_t position, std::size_t source, std::size_t offset, Kmer::Strand strand);
	/**
	 * Find out where a hash is in the list of kmers.
	 * @param hash the hash to search for.
	 * @return the location of the hash (or -1 if we didn't find it).
	 */
	int32_t findKmer(std::size_t hash);
	/**
	 * How many kmers are in this node?
	 * @return the number of kmers in this node.
	 */
	std::size_t kmerCount();
private:
	/** the human-readable name for this sequence node */
	std::string name;
	/** the ordered collection of kmers found in this node */
	std::vector<boost::shared_ptr<Kmer> > kmers;
	/** a cache of positions for fast lookup of kmer by hash */
	boost::unordered_map<std::size_t, std::size_t> kmerLocation;
	/** identifier for this sequence node */
	std::size_t id;

	/** 
	 * Update the cache of locations of all the k-mers.
	 */
	void updateKmerLocations();
};

#endif // SEQUENCE_NODE_HH
