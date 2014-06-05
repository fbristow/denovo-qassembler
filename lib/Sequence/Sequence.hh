/* 
 * File:   Sequence.hh
 * Author: fbristow
 *
 * Created on June 13, 2012 
 */
#ifndef SEQUENCE_HH
#define SEQUENCE_HH

#include <string>

class Sequence {
public:
	/**
	 * Default constructor.
	 */
	Sequence ();

	/**
	 * Copy constructor.
	 * @param seq the sequence to copy from.
	 */
	Sequence (Sequence *seq);

	/**
	 * Constructor specifying all parts.
	 * @param sequence the sequence.
	 * @param name the name of the sequence record.
	 * @param comment the comment for the sequence record.
	 * @param qual the quality string associated with the sequence record.
	 */
	Sequence (std::string sequence, std::string name, std::string comment, std::string qual);

	/**
	 * Get the sequence from this read.
	 * @return the sequence in this read.
	 */
	std::string getSequence();

	/**
	 * Setter for sequence.
	 * @param sequence the new sequence.
	 */
	void setSequence(std::string sequence);

	/**
	 * Get the reverse complement of the sequence in this read.
	 * @return the reverse complemented sequence.
	 */
	std::string getReverseComplement();

	/**
	 * Get the name of this sequence record.
	 * @return the human-readable name for this sequence record.
	 */
	std::string getName();

	/**
	 * Setter for name.
	 * @param the new name.
	 */
	void setName(std::string sequence);

	/**
	 * Get the comment of this sequence record.
	 * @return the human-readable comment for this sequence record.
	 */
	std::string getComment();

	/**
	 * Setter for comment.
	 * @param comment the new comment.
	 */
	void setComment(std::string comment);

	/**
	 * Get the quality string for this sequence record.
	 * @return the quality string of this sequence record.
	 */
	std::string getQual();

	/**
	 * Setter for qual.
	 * @param qual the new quality string.
	 */
	void setQual(std::string qual);

	/**
	 * Get the number of base-pairs in this sequence.
	 * @return the length of the sequence.
	 */
	std::size_t getLength();

	/**
	 * Get a unique, numeric identifier for this read for quick lookups.
	 * @return a unique, numeric identifier for this read.
	 */
	std::size_t getID();

	/**
	 * Override the generated identifier.
	 * @param id the new identifier for this read.
	 */
	void setID(std::size_t id);
private:
	std::string sequence;
	std::string reverse;
	std::string name;
	std::string comment;
	std::string qual;
	std::size_t id;

	void revcom();

};

#endif // SEQUENCE_HH
