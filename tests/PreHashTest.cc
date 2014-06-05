/* 
 * File:   PreHashTest.cc
 * Author: fbristow
 *
 * Created on February 1, 2012
 */
#ifndef PRE_HASH_TEST_CC
#define PRE_HASH_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "PreHash/PreHash.hh"
#include "Sequence/Sequence.hh"

struct PreHashFixture {
	PreHashFixture() {
		read1 = boost::make_shared<Sequence>("AACT", "read1", "", "++++");
		read2 = boost::make_shared<Sequence>("ACTC", "read2", "", "++++");
		read3 = boost::make_shared<Sequence>("actc", "read3", "", "++++");

		AAC = qassembler::hash("AAC");
		ACT = qassembler::hash("ACT");
		AGT = qassembler::hash("AGT");
		GTT = qassembler::hash("GTT");
		GAG = qassembler::hash("GAG");
		CTC = qassembler::hash("CTC");
	}

	boost::shared_ptr<Sequence> read1;
	boost::shared_ptr<Sequence> read2;
	boost::shared_ptr<Sequence> read3;
	std::size_t AAC, ACT, AGT, GTT, GAG, CTC;
};

BOOST_FIXTURE_TEST_SUITE(pre_hash, PreHashFixture)

BOOST_AUTO_TEST_CASE (add_one_sequence) {
	BOOST_TEST_CHECKPOINT("Add single sequence to pre-hash.");
	PreHash p(3, 3);
	p.addRead(read1);

	BOOST_REQUIRE_EQUAL(p.hashCount(0xdeadbeef), 0);
	BOOST_REQUIRE_EQUAL(p.hashCount(AAC), 1); // AAC
	BOOST_REQUIRE_EQUAL(p.hashCount(ACT), 1); // ACT
	BOOST_REQUIRE_EQUAL(p.hashCount(AGT), 1); // AGT (reverse)
	BOOST_REQUIRE_EQUAL(p.hashCount(GTT), 1); // GTT (reverse)

	BOOST_REQUIRE_EQUAL(p.getHashes(read1->getID(), Kmer::FORWARD).size(), 2);
	BOOST_REQUIRE_EQUAL(p.getHashes(read1->getID(), Kmer::REVERSE).size(), 2);
}

BOOST_AUTO_TEST_CASE (add_overlapping_sequences) {
	PreHash p(3, 3);
	p.addRead(read1);
	p.addRead(read2);

	BOOST_REQUIRE_EQUAL(p.hashCount(0xdeadbeef), 0);
	BOOST_REQUIRE_EQUAL(p.hashCount(AAC), 1); // AAC
	BOOST_REQUIRE_EQUAL(p.hashCount(ACT), 2); // ACT
	BOOST_REQUIRE_EQUAL(p.hashCount(AGT), 2); // AGT (reverse)
	BOOST_REQUIRE_EQUAL(p.hashCount(GTT), 1); // GTT (reverse)

	BOOST_REQUIRE_EQUAL(p.getReads(ACT).size(), 2);
}

BOOST_AUTO_TEST_CASE (add_upper_lower_case_sequences) {
	PreHash p(3, 3);
	p.addRead(read2);
	p.addRead(read3);

	BOOST_REQUIRE_EQUAL(p.hashCount(0xdeadbeef), 0);
	BOOST_REQUIRE_EQUAL(p.hashCount(ACT), 2);
	BOOST_REQUIRE_EQUAL(p.hashCount(CTC), 2);
	BOOST_REQUIRE_EQUAL(p.hashCount(GAG), 2);
	BOOST_REQUIRE_EQUAL(p.hashCount(AGT), 2);

	BOOST_REQUIRE_EQUAL(p.getHashes(read2->getID(), Kmer::FORWARD).size(), 2);
	BOOST_REQUIRE_EQUAL(p.getHashes(read3->getID(), Kmer::FORWARD).size(), 2);
}

BOOST_AUTO_TEST_CASE (kmer_count) {
	PreHash p(3, 3);
	p.addRead(read2);
	p.addRead(read3);

	BOOST_REQUIRE_EQUAL(p.kmerCount("0xdeadbeef"), 0);
	BOOST_REQUIRE_EQUAL(p.kmerCount("ACT"), 2);
	BOOST_REQUIRE_EQUAL(p.kmerCount("CTC"), 2);
	BOOST_REQUIRE_EQUAL(p.kmerCount("GAG"), 2);
	BOOST_REQUIRE_EQUAL(p.kmerCount("AGT"), 2);
}

BOOST_AUTO_TEST_SUITE_END()
#endif // PRE_HASH_TEST_CC
