/* 
 * File:   PreHash.hh
 * Author: fbristow
 *
 * Created on January 30, 2012
 */
#ifndef PRE_HASH_HH
#define PRE_HASH_HH

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <boost/shared_ptr.hpp>

#include <vector>

#include "Util/Util.hh"
#include "Kmer/Kmer.hh"
#include "Sequence/Sequence.hh"

#include "Logging/Logging.hh"

class PreHash {
public:
	/**
	 * Default constructor.
	 */
	PreHash();
	/**
	 * specify a k-mer size and initial size of hash to use for the pre-hash.
	 * @param k the k size for this pre-hash
	 * @param initialHashSize how many buckets should be used for the sets (performance related).
	 */
	PreHash(std::size_t k, std::size_t initialHashSize);

	/**
	 * add a read to this pre-hash
	 * @param read the read to use when generating k-mers
	 */
	void addRead(boost::shared_ptr<Sequence> read);

	/**
	 * get the hashes that belong to a specific read in a specific orientation
	 * @param readId the read to get hashes for
	 * @param direction the orientation of read to get hashes for.
	 */
	std::vector<std::size_t> getHashes(std::size_t readId, Kmer::Strand direction);

	/**
	 * Get all reads that contain a specific hash.
	 * @param hash the hash to find reads for
	 * @return the set of reads that contain this hash.
	 */
	boost::unordered_set<std::size_t> getReads(std::size_t hash);

	/** 
	 * How many times has this hash been generated in the set of reads?
	 * @param hash the hash to count
	 * @return the number of times the hash was seen in the data set
	 */
	std::size_t hashCount(std::size_t hash);

	/**
	 * How many times was the kmer seen in the set of reads? Convenience method
	 * so that we don't have to hash the kmer twice unnecessarily.
	 * @param the kmer to count
	 * @return the number of times the kmer was seen in the read set.
	 */
	std::size_t kmerCount(std::string kmer);

	/**
	 * Get all hashes generated in the data set.
	 * @return the set of all hashes in the data set
	 */
	boost::unordered_set<std::size_t> getAllHashes();
private:
	/**
	 * Add a read from a string.
	 * @param sequence the sequence to add.
	 * @param readId the read identifier.
	 * @param strand the strand used when the read is being added.
	 */
	void addRead(std::string sequence, std::size_t readId, Kmer::Strand strand);
	/**
	 * Add a kmer to this pre-hash.
	 * @param hash the hash to add for the kmer.
	 * @param readId the read identifier where this kmer was generated from.
	 * @param strand the strand that the kmer came from.
	 * @param position the position in the read where the kmer came from.
	 */
	void addKmer(std::size_t hash, std::size_t readId, Kmer::Strand strand, std::size_t position);
	/** the set of all records encountered */
	boost::unordered_map<std::size_t, boost::unordered_set<std::size_t> > hashes2reads;
	/** a record of hashes found in reads */
	boost::unordered_map<std::size_t,
		boost::unordered_map<Kmer::Strand, std::vector<std::size_t> > > reads2hashes;
	/** length of k-mers */
	std::size_t kmerLength;
};

#endif // PRE_HASH_HH
