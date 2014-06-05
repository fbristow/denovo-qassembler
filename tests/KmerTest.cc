/* 
 * File:   KmerTest.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef KMER_TEST_CC
#define KMER_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "Kmer/Kmer.hh"

BOOST_AUTO_TEST_SUITE (kmer)

BOOST_AUTO_TEST_CASE (constructors_test) {
	std::size_t hash = 0x42;
	char base = 'c';
	std::size_t id = 0x9001;
	std::size_t pos = 0x1337;
	Kmer::Strand strand = Kmer::FORWARD;
	BOOST_TEST_CHECKPOINT("Creating new k-mer");
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>(hash, base, id, pos, strand);
	
	BOOST_REQUIRE_EQUAL(k->getHash(), hash);
	BOOST_REQUIRE_EQUAL(k->getBase(), base);
	
	BOOST_TEST_CHECKPOINT("getting the sources from a k-mer");
	// check that the first kmer was inserted correctly:
	boost::unordered_map<std::size_t, Kmer::Source> sources = k->getSources();
	BOOST_REQUIRE_EQUAL(sources.size(), 1);

	std::size_t sourceId = sources.begin()->first;
	std::size_t sourcePos;
	Kmer::Strand sourceStrand;
	
	boost::tie(sourcePos, sourceStrand) = sources[sourceId];
	
	BOOST_REQUIRE_EQUAL(sourceId, id);
	BOOST_REQUIRE_EQUAL(sourcePos, pos);
	BOOST_REQUIRE_EQUAL(sourceStrand, strand);

	BOOST_TEST_CHECKPOINT("calling copy constructor");
	Kmer k2(*k);
	BOOST_REQUIRE_EQUAL(k2.getHash(), hash);
	BOOST_REQUIRE_EQUAL(k2.getBase(), base);
	
	sources = k2.getSources();
	BOOST_REQUIRE_EQUAL(sources.size(), 1);
	sourceId = sources.begin()->first;
	
	boost::tie(sourcePos, sourceStrand) = sources[sourceId];
	
	BOOST_REQUIRE_EQUAL(sourceId, id);
	BOOST_REQUIRE_EQUAL(sourcePos, pos);
	BOOST_REQUIRE_EQUAL(sourceStrand, strand);
}

BOOST_AUTO_TEST_CASE (getters_setters_test) {
	std::size_t hash = 0x42;
	char base = 'c';
	std::size_t id = 0x9001;
	std::size_t pos = 0x1337;
	Kmer::Strand strand = Kmer::FORWARD;
	
	boost::unordered_map<std::size_t, Kmer::Source> sources;
	sources[id] = std::make_pair(pos, strand);
	
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>();
	
	k->setHash(hash);
	BOOST_REQUIRE_EQUAL(k->getHash(), hash);
	k->setBase(base);
	BOOST_REQUIRE_EQUAL(k->getBase(), base);
	k->setSources(sources);
	BOOST_REQUIRE_EQUAL(k->getCount(), 1);
	
	std::size_t sourceId = k->getSources().begin()->first;
	std::size_t sourcePos;
	Kmer::Strand sourceStrand;
	
	boost::tie(sourcePos, sourceStrand) = k->getSource();
	
	BOOST_REQUIRE_EQUAL(sourceId, id);
	BOOST_REQUIRE_EQUAL(sourcePos, pos);
	BOOST_REQUIRE_EQUAL(sourceStrand, strand);
}

BOOST_AUTO_TEST_CASE (add_source_test) {
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>();
	
	std::size_t id = 0x9001;
	std::size_t pos = 0x1337;
	Kmer::Strand strand = Kmer::FORWARD;
	
	BOOST_REQUIRE_EQUAL(k->getSources().size(), 0);
	
	k->addSource(id, pos, strand);

	BOOST_REQUIRE_EQUAL(k->getSources().size(), 1);
	
	std::size_t sourceId = k->getSources().begin()->first;
	std::size_t sourcePos;
	Kmer::Strand sourceStrand;
	
	boost::tie(sourcePos, sourceStrand) = k->getSources()[sourceId];
	
	BOOST_REQUIRE_EQUAL(sourceId, id);
	BOOST_REQUIRE_EQUAL(sourcePos, pos);
	BOOST_REQUIRE_EQUAL(sourceStrand, strand);
}

BOOST_AUTO_TEST_CASE (get_source_test) {
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>();
	
	std::size_t id = 0x9001;
	std::size_t pos = 0x1337;
	Kmer::Strand strand = Kmer::FORWARD;
	
	std::size_t id2 = 0x42;
	std::size_t pos2 = 0x42;
	Kmer::Strand strand2 = Kmer::REVERSE;
	
	BOOST_REQUIRE_EQUAL(k->getSources().size(), 0);
	
	k->addSource(id, pos, strand);
	k->addSource(id2, pos2, strand2);

	BOOST_REQUIRE_EQUAL(k->getSources().size(), 2);
	
	std::size_t sourceId = k->getSources().begin()->first;
	std::size_t sourcePos;
	Kmer::Strand sourceStrand;
	
	boost::tie(sourcePos, sourceStrand) = k->getSources()[sourceId];
	
	BOOST_REQUIRE_EQUAL(sourceId, id);
	BOOST_REQUIRE_EQUAL(sourcePos, pos);
	BOOST_REQUIRE_EQUAL(sourceStrand, strand);
}

BOOST_AUTO_TEST_CASE (get_count_test) {
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>();
	
	std::size_t id = 0x9001;
	std::size_t pos = 0x1337;
	Kmer::Strand strand = Kmer::FORWARD;
	
	std::size_t id2 = 0x42;
	std::size_t pos2 = 0x42;
	Kmer::Strand strand2 = Kmer::REVERSE;
	
	BOOST_REQUIRE_EQUAL(k->getCount(), 0);
	
	k->addSource(id, pos, strand);
	BOOST_REQUIRE_EQUAL(k->getCount(), 1);
	
	k->addSource(id2, pos2, strand2);
	BOOST_REQUIRE_EQUAL(k->getCount(), 2);
}

BOOST_AUTO_TEST_CASE (transition_count_test) {
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>();
	BOOST_REQUIRE_EQUAL(k->getTransitionCount('C'), 0);
	k->addTransition('C');
	BOOST_REQUIRE_EQUAL(k->getTransitionCount('C'), 1);

	BOOST_REQUIRE_EQUAL(k->getTransitionCount('G'), 0);
	for (int i = 0; i < 300; i++) {
		k->addTransition('G');
	}
	BOOST_REQUIRE_EQUAL(k->getTransitionCount('G'), 300);
	BOOST_REQUIRE_EQUAL(k->getTransitionCount('C'), 1);
}

BOOST_AUTO_TEST_SUITE_END()
#endif // KMER_TEST_CC
