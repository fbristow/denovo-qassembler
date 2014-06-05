/* 
 * File:   SequenceNodeTest.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef SEQUENCE_NODE_TEST_CC
#define SEQUENCE_NODE_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Graph/Node/SequenceNode.hh"

struct SequenceNodeFixture {
	SequenceNodeFixture() { 
		std::size_t hash = 0x42;
		char base = 'c';
		std::size_t id = 0x9001;
		std::size_t pos = 0x1337;
		Kmer::Strand strand = Kmer::FORWARD;
		boost::shared_ptr<Kmer> a = boost::make_shared<Kmer>(hash, base, id, pos, strand);
		
		kmers1.push_back(a);
		
		hash++; base = 'g'; id++; pos++; strand = Kmer::REVERSE;
		boost::shared_ptr<Kmer> b = boost::make_shared<Kmer>(hash, base, id, pos, strand);
		kmers1.push_back(b);
		
		hash++; base = 'a'; id++; pos++; strand = Kmer::FORWARD;
		boost::shared_ptr<Kmer> c = boost::make_shared<Kmer>(hash, base, id, pos, strand);
		kmers2.push_back(c);
		
		hash++; base = 't'; id++; pos++; strand = Kmer::REVERSE;
		boost::shared_ptr<Kmer> d = boost::make_shared<Kmer>(hash, base, id, pos, strand);
		kmers2.push_back(d);
	}
	
	~SequenceNodeFixture() { 
		kmers1.clear();
		kmers2.clear();
	}
	
	std::vector<boost::shared_ptr<Kmer> > kmers1;
	std::vector<boost::shared_ptr<Kmer> > kmers2;
};

BOOST_FIXTURE_TEST_SUITE (sequence_node, SequenceNodeFixture)

BOOST_AUTO_TEST_CASE (constructors_test) {
	std::size_t id = 0x30;
	std::string name = "nodeName";
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>(id, name);
	BOOST_REQUIRE_EQUAL(n->getId(), id);
	BOOST_REQUIRE_EQUAL(n->getName(), name);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 0);
	BOOST_REQUIRE_EQUAL(n->sequence(), "");
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 0);
	
	boost::shared_ptr<SequenceNode> n2 = boost::make_shared<SequenceNode>(n.get());
	BOOST_REQUIRE_EQUAL(n2->getId(), id);
	BOOST_REQUIRE_EQUAL(n2->getName(), name);
	BOOST_REQUIRE_EQUAL(n2->getKmers().size(), 0);
	BOOST_REQUIRE_EQUAL(n2->sequence(), "");
	BOOST_REQUIRE_EQUAL(n2->getKmers().size(), 0);
}

BOOST_AUTO_TEST_CASE (getter_setter_test) {
	std::string name = "nodeName";
	std::size_t id = 0x30;
	
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	n->setName(name);
	BOOST_REQUIRE_EQUAL(n->getName(), name);
	
	n->setId(id);
	BOOST_REQUIRE_EQUAL(n->getId(), id);
	
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 0);
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
}

BOOST_AUTO_TEST_CASE (get_sequence) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->sequence(), "cg");
}

BOOST_AUTO_TEST_CASE (merge) {
	boost::shared_ptr<SequenceNode> n1 = boost::make_shared<SequenceNode>();
	boost::shared_ptr<SequenceNode> n2 = boost::make_shared<SequenceNode>();
	n1->setKmers(kmers1);
	n2->setKmers(kmers2);
	
	n1->merge(n2);
	
	BOOST_REQUIRE_EQUAL(n1->getKmers().size(), 4);
	BOOST_REQUIRE_EQUAL(n1->sequence(), "atcg");
	
	std::vector<boost::shared_ptr<Kmer> > kmers = n1->getKmers();
	BOOST_REQUIRE_EQUAL(kmers[0]->getHash(), 0x44);
	BOOST_REQUIRE_EQUAL(kmers[1]->getHash(), 0x45);
	BOOST_REQUIRE_EQUAL(kmers[2]->getHash(), 0x42);
	BOOST_REQUIRE_EQUAL(kmers[3]->getHash(), 0x43);
}

BOOST_AUTO_TEST_CASE (add_kmer_params) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	n->setKmers(kmers1);
	
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	
	n->addKmer(0x42, 'c', 0x42, 0x42, Kmer::REVERSE);
	
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 3);
	BOOST_REQUIRE_EQUAL(n->getKmers()[2]->getHash(), 0x42);
	BOOST_REQUIRE_EQUAL(n->getKmers()[2]->getBase(), 'c');
}

BOOST_AUTO_TEST_CASE (add_kmer) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>(0x42, 'c', 0x42, 0x42, Kmer::REVERSE);
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	n->addKmer(k);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 3);
	BOOST_REQUIRE_EQUAL(n->getKmers()[2]->getHash(), k->getHash());
	BOOST_REQUIRE_EQUAL(n->getKmers()[2]->getBase(), k->getBase());
}

BOOST_AUTO_TEST_CASE (add_kmer_at) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	boost::shared_ptr<Kmer> k = boost::make_shared<Kmer>(0x9001, 'u', 0x9001, 0x9001, Kmer::REVERSE);
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	n->addKmerAt(k, 1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 3);
	BOOST_REQUIRE_EQUAL(n->getKmers()[1]->getHash(), 0x9001);
	BOOST_REQUIRE_EQUAL(n->getKmers()[1]->getBase(), 'u');
}

BOOST_AUTO_TEST_CASE (add_kmer_source_at) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	n->addKmerSourceAt(1, 0x9001, 0x9001, Kmer::REVERSE);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	BOOST_REQUIRE_EQUAL(n->getKmers()[1]->getCount(), 2);
	std::size_t id = 0x9001;
	std::size_t pos;
	Kmer::Strand strand;
	boost::tie(pos, strand) = n->getKmers()[1]->getSources()[id];
	BOOST_REQUIRE_EQUAL(id, 0x9001);
	BOOST_REQUIRE_EQUAL(pos, 0x9001);
	BOOST_REQUIRE_EQUAL(strand, Kmer::REVERSE);
}

BOOST_AUTO_TEST_CASE (find_kmer) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->getKmers().size(), 2);
	BOOST_REQUIRE_EQUAL(n->findKmer(0x42), 0);
	BOOST_REQUIRE_EQUAL(n->findKmer(0x9001), -1);
	n->addKmer(0x9001, 'u', 0x9001, 0x9001, Kmer::REVERSE);
	BOOST_REQUIRE_EQUAL(n->findKmer(0x9001), 2);
}

BOOST_AUTO_TEST_CASE (kmer_count) {
	boost::shared_ptr<SequenceNode> n = boost::make_shared<SequenceNode>();
	BOOST_REQUIRE_EQUAL(n->kmerCount(), 0);
	n->setKmers(kmers1);
	BOOST_REQUIRE_EQUAL(n->kmerCount(), 2);
	n->addKmer(0x9001, 'u', 0x9001, 0x9001, Kmer::REVERSE);
	BOOST_REQUIRE_EQUAL(n->kmerCount(), 3);
}

BOOST_AUTO_TEST_SUITE_END ()
#endif // SEQUENCE_NODE_TEST_CC
