/* 
 * File:   SequenceNodeTest.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef SEQUENCE_TEST_CC
#define SEQUENCE_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Sequence/Sequence.hh"
#include "Util/Util.hh"
#include "Exception/InvalidInputException.hh"

BOOST_AUTO_TEST_SUITE (SEQUENCE)

BOOST_AUTO_TEST_CASE (constructors_test) {
	Sequence seq ("ACGT", "read1", "comment", "++++");
	BOOST_REQUIRE_EQUAL("ACGT", seq.getSequence());
	BOOST_REQUIRE_EQUAL("read1", seq.getName());
	BOOST_REQUIRE_EQUAL("++++", seq.getQual());
	BOOST_REQUIRE_EQUAL("comment", seq.getComment());
	BOOST_REQUIRE_EQUAL("ACGT", seq.getReverseComplement());
	BOOST_REQUIRE_EQUAL(4, seq.getLength());

	Sequence seq2(&seq);
	BOOST_REQUIRE_EQUAL("ACGT", seq2.getSequence());
	BOOST_REQUIRE_EQUAL("read1", seq2.getName());
	BOOST_REQUIRE_EQUAL("++++", seq2.getQual());
	BOOST_REQUIRE_EQUAL("comment", seq2.getComment());
	BOOST_REQUIRE_EQUAL("ACGT", seq2.getReverseComplement());

}

BOOST_AUTO_TEST_CASE (reverse_complement_test) {
	Sequence seq ("AAAATCTG", "read1", "comment", "++++++++");
	BOOST_REQUIRE_EQUAL(seq.getReverseComplement(), "CAGATTTT");
}

BOOST_AUTO_TEST_CASE (numeric_identifier_generator) {
	Sequence seq ("AAAATCTG", "read1", "comment", "++++++++");
	BOOST_REQUIRE_EQUAL(seq.getID(), qassembler::hash("read1"));
}

// ambiguous nucleotides as per: http://www.bioinformatics.org/sms/iupac.html
BOOST_AUTO_TEST_CASE (ambiguous_reverse_complement) {
	Sequence seq ("ACGTRYSWKMBDHVN", "read1", "comment", "++++++++++++++");
	BOOST_REQUIRE_EQUAL(seq.getReverseComplement(), "NBDHVKMWSRYACGT");
}

BOOST_AUTO_TEST_CASE (non_dna_input) {
	try {
		Sequence seq("NOTDNAXXX", "read1", "comment", "++++");
		BOOST_REQUIRE_EQUAL(0, 1); // AUTOFAIL
	} catch (InvalidInputException &e) {

	}
	
}

BOOST_AUTO_TEST_SUITE_END ()
#endif // SEQUENCE_TEST_CC
