/* 
 * File:   Kmer.hh
 * Author: fbristow
 *
 * Created on December 15, 2011
 */
#ifndef KMER_HH
#define KMER_HH

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

class Kmer {
public:
	enum Strand {
		FORWARD,
		REVERSE
	};
	
	/** typedef for kmer sources, tuple contains:
	 * [0] -> read identifier,
	 * [1] -> read offset,
	 * [2] -> strand.
	 */
	//typedef boost::tuple<std::size_t, std::size_t, Strand> Source;
	typedef std::pair<std::size_t /* offset */, Strand /* direction */> Source;
	/**
	 * Default constructor
	 */
	Kmer();
	/**
	 * Copy constructor
	 */
	Kmer(Kmer *);
	/**
	 * Specify all parameters of a Kmer.
	 * @param hash the hash generated by this kmer sequence
	 * @param base the most significant base (i.e., the last base) for this kmer
	 * @param source the read that this kmer came from
	 * @param position the place where this kmer came from in the source read
	 * @param strand the direction the read was in when the kmer was generated
	 */
	Kmer(std::size_t hash, char base, std::size_t source, std::size_t position, Strand strand);
	
	/**
	 * Get the hash generated by this kmer seqence.
	 * @return the hash generated by this kmer sequence
	 */
	std::size_t getHash();
	/**
	 * Set the hash for this kmer sequence.
	 * @param hash the hash to use for this kmer sequence.
	 */
	void setHash(std::size_t hash);
	/**
	 * Get the complete sequence for this kmer. Sub-classes may store more than a single character, particularly
	 * in the case where the kmer was the first generated for a specific read. This method always returns a string,
	 * even if the kmer only stores a single character.
	 * @return the complete sequence represented by this kmer
	 */
	virtual std::string getSequence();
	/**
	 * Get the most significant base for this kmer, even if the kmer stores more than one character.
	 * @return the most significant base for this kmer.
	 */
	virtual char getBase();
	/**
	 * Set the most significant base for this kmer.
	 * @param base the most significant base for this kmer.
	 */
	void setBase(char base);
	
	/**
	 * Add a new source to this kmer. We expect to see identical Kmers in many reads (that's how de Bruijn graph assembly works!),
	 * so we may want to keep track of where all the kmers are coming from.
	 * @param source the read where the kmer came from
	 * @param position the position where the kmer came from in the read
	 * @param strand the direction we were scanning the read when the kmer was generated
	 */
	void addSource(std::size_t source, std::size_t position, Strand strand);
	/**
	 * Get the list of all sources where this kmer can be found.
	 * @return the list of all sources where this kmer can be found.
	 */
	boost::unordered_map<std::size_t /* readIdentifier */, Source /* position */> getSources();
	/**
	 * Set the list of all sources where this kmer can be found.
	 * @param sources the new list of sources where this kmer can be found.
	 */
	void setSources(boost::unordered_map<std::size_t /* readIdentifier */, Source /* position */> sources);
	/**
	 * Find out how many times this kmer has been observed in the data set.
	 * @return the number of instances of this kmer in the data set.
	 */
	std::size_t getCount();
	/**
	 * Get the first read where this kmer was observed.
	 * @return the first read where this kmer was observed.
	 */
	Source getSource();

	/**
	 * Add a new transition that was observed in the data set (i.e., if you add a kmer pair where the next kmer has a
	 * most significant base of 'G', then you've observed a transition from this kmer to 'G').
	 * @param base the base where you transitioned to.
	 */
	void addTransition(char base);

	/**
	 * Determine how many transitions were made for a certain nucleotide.
	 * @param base the base to count transitions for
	 * @return the number of times the current kmer transitioned to that nucleotide
	 */
	std::size_t getTransitionCount(char base);
private:
	/** the most significant base in this kmer */
	char base;
	/** the hash generated by this kmer */
	std::size_t hash;
	/** the collection of places where this kmer was observed in the data set */
	boost::unordered_map<std::size_t /* readIdentifier */, Source /* position */> sources;
	/** the total number of times that this kmer transitioned to a different base */
	boost::unordered_map<char /* transition */, std::size_t /* transitionCount */> transitions;
};

#endif // KMER_HH
