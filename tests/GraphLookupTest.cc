/* 
 * File:   GraphLookupTest.cc
 * Author: fbristow
 *
 * Created on January 25, 2012
 */
#ifndef GRAPH_LOOKUP_TEST_CC
#define GRAPH_LOOKUP_TEST_CC

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include "Lookup/GraphLookup.hh"
#include "Graph/SkinnyGraph.hh"

BOOST_AUTO_TEST_SUITE (GRAPH_LOOKUP)

BOOST_AUTO_TEST_CASE (put_test) {
	GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > lookup;
	boost::shared_ptr<SkinnyGraph> g = boost::make_shared<SkinnyGraph>(1);
	lookup.put(0x42, g); 
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x42), g);
	BOOST_FOREACH(std::size_t gt, lookup.get(g)) {
		BOOST_REQUIRE_EQUAL(gt, 0x42);
	}
}

BOOST_AUTO_TEST_CASE (replace_test) {
	GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > lookup;
	boost::shared_ptr<SkinnyGraph> g = boost::make_shared<SkinnyGraph>(1);
	boost::shared_ptr<SkinnyGraph> g2 = boost::make_shared<SkinnyGraph>(2);
	lookup.put(0x42, g);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x42), g);
	lookup.put(0x42, g2);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g2), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x42), g2);
}

BOOST_AUTO_TEST_CASE (clear_hash_test) {
	GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > lookup;
	boost::shared_ptr<SkinnyGraph> g = boost::make_shared<SkinnyGraph>(1);
	boost::shared_ptr<SkinnyGraph> g2 = boost::make_shared<SkinnyGraph>(2);
	lookup.put(0x42, g);
	lookup.put(0x43, g);
	lookup.put(0x44, g2);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(0x43), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(0x44), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 2);
	BOOST_REQUIRE_EQUAL(lookup.count(g2), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x42), g);
	BOOST_REQUIRE_EQUAL(lookup.get(0x43), g);
	BOOST_REQUIRE_EQUAL(lookup.get(0x44), g2);

	lookup.clear(0x42);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 0);
	BOOST_REQUIRE_EQUAL(lookup.count(0x43), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(0x44), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g2), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x43), g);
	BOOST_REQUIRE_EQUAL(lookup.get(0x44), g2);
}

BOOST_AUTO_TEST_CASE (clear_graph_test) {
	GraphLookup<std::size_t, boost::shared_ptr<SkinnyGraph> > lookup;
	boost::shared_ptr<SkinnyGraph> g = boost::make_shared<SkinnyGraph>(1);
	boost::shared_ptr<SkinnyGraph> g2 = boost::make_shared<SkinnyGraph>(2);
	lookup.put(0x42, g);
	lookup.put(0x43, g);
	lookup.put(0x44, g2);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(0x43), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(0x44), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 2);
	BOOST_REQUIRE_EQUAL(lookup.count(g2), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x42), g);
	BOOST_REQUIRE_EQUAL(lookup.get(0x43), g);
	BOOST_REQUIRE_EQUAL(lookup.get(0x44), g2);

	lookup.clear(g);
	BOOST_REQUIRE_EQUAL(lookup.count(0x42), 0);
	BOOST_REQUIRE_EQUAL(lookup.count(0x43), 0);
	BOOST_REQUIRE_EQUAL(lookup.count(0x44), 1);
	BOOST_REQUIRE_EQUAL(lookup.count(g), 0);
	BOOST_REQUIRE_EQUAL(lookup.count(g2), 1);
	BOOST_REQUIRE_EQUAL(lookup.get(0x44), g2);
}

BOOST_AUTO_TEST_SUITE_END()
#endif // GRAPH_LOOKUP_TEST_CC
