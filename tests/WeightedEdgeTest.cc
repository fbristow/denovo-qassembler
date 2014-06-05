/* 
 * File:   WeightedEdgeTest.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef WEIGHTED_EDGE_TEST_CC
#define WEIGHTED_EDGE_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "Graph/Edge/WeightedEdge.hh"

BOOST_AUTO_TEST_SUITE (weighted_edge)

BOOST_AUTO_TEST_CASE (constructors_test) {
	WeightedEdge e;
	BOOST_REQUIRE_EQUAL(e.getWeight(), 1.);
	
	WeightedEdge e2(42.);
	BOOST_REQUIRE_EQUAL(e2.getWeight(), 42.);
}

BOOST_AUTO_TEST_CASE (getter_setter_test) {
	WeightedEdge e;

	e.setWeight(42.);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 42.);
}

BOOST_AUTO_TEST_CASE (increment_test) {
	WeightedEdge e;
	
	e++;
	BOOST_REQUIRE_EQUAL(e.getWeight(), 2.);
}

BOOST_AUTO_TEST_CASE (reset_test) {
	WeightedEdge e;
	// default edge weight is 1
	e.setWeight(30);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 30.);
	e.resetWeight();
	BOOST_REQUIRE_EQUAL(e.getWeight(), 1.);
}

BOOST_AUTO_TEST_CASE (lock_test) {
	WeightedEdge e;
	e.setWeight(30);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 30.);
	e.lockWeight();
	BOOST_REQUIRE_EQUAL(e.getWeight(), 30.);
	e.setWeight(42);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 42.);
	e.resetWeight(); 
	BOOST_REQUIRE_EQUAL(e.getWeight(), 30.);
}

BOOST_AUTO_TEST_CASE (decrease_by_large_amount) {
	// When decreasing an edge weight by an amount larger than the edge, just set the edge weight to 0:
	WeightedEdge e;
	e.setWeight(20);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 20);
	e.decreaseWeight(40);
	BOOST_REQUIRE_EQUAL(e.getWeight(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
#endif // WEIGHTED_EDGE_TEST_CC
