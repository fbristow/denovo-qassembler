/*
 * File: SkinnyGraphTest.cc
 * Author: fbristow
 *
 * Created on January 4, 2012
 */
#ifndef SKINNY_GRAPH_TEST_CC
#define SKINNY_GRAPH_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include "Graph/SkinnyGraph.hh"
#include "Exception/InvalidGraphStateException.hh"

struct SkinnyGraphFixture {
	SkinnyGraphFixture() : g(3), g2(boost::make_shared<SkinnyGraph>(4)), g3(boost::make_shared<SkinnyGraph>(5)) {
		char nuc[5] = {'a', 'c', 'g', 't', 'u'};
		std::size_t id = 0x9001;
		std::size_t hash = 0x0042;
		BOOST_TEST_CHECKPOINT("Going to create a new sequence node");
		v = g.createSequenceNode(hash, nuc[0], "read0001", id, 0, Kmer::FORWARD);

		for (int i = 1; i < 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g.node(v)->addKmer(++hash, nuc[i % 5], ++id, i, Kmer::FORWARD);
		}

		v2 = g.createSequenceNode(++hash, nuc[4], "read0002", ++id, 0, Kmer::FORWARD);
		for (int i = 8; i >= 0; i--) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g.node(v2)->addKmer(++hash, nuc[i % 5], ++id, i, Kmer::FORWARD);
		}

		v3 = g2->createSequenceNode(++hash, 'a', "read0003", ++id, 0, Kmer::FORWARD);
		for (int i = 1; i <= 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g2->node(v3)->addKmer(++hash, 'a', ++id, i, Kmer::FORWARD);
		}
		v4 = g2->createSequenceNode(++hash, 'g', "read0004", ++id, 0, Kmer::FORWARD);
		for (int i = 1; i <= 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g2->node(v4)->addKmer(++hash, 'g', ++id, i, Kmer::FORWARD);
		}
		v5 = g3->createSequenceNode(++hash, 't', "read0005", ++id, 0, Kmer::FORWARD);
		for (int i = 1; i <= 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g3->node(v5)->addKmer(++hash, 't', ++id, i, Kmer::FORWARD);
		}

		v6 = g3->createSequenceNode(++hash, 'a', "read0006", ++id, 0, Kmer::FORWARD);
		for (int i = 1; i <= 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g3->node(v6)->addKmer(++hash, 'a', ++id, i, Kmer::FORWARD);
		}
		v7 = g3->createSequenceNode(++hash, 'g', "read0007", ++id, 0, Kmer::FORWARD);
		for (int i = 1; i <= 10; i++) {
			BOOST_TEST_CHECKPOINT("Adding new k-mer to node: " << i);
			g3->node(v7)->addKmer(++hash, 'g', ++id, i, Kmer::FORWARD);
		}
	}

	~SkinnyGraphFixture() {
	}

	SkinnyGraph g;
	boost::shared_ptr<SkinnyGraph> g2;
	boost::shared_ptr<SkinnyGraph> g3;
	SkinnyGraph::Vertex v, v2, v3, v4, v5, v6, v7;
};

BOOST_FIXTURE_TEST_SUITE (skinny_graph, SkinnyGraphFixture)

BOOST_AUTO_TEST_CASE (constructor_test) {
	SkinnyGraph g(3);

	BOOST_REQUIRE_EQUAL(g.numVertices(), 0);
}

BOOST_AUTO_TEST_CASE (create_sequence_node) {
	SkinnyGraph g(3);

	BOOST_REQUIRE_EQUAL(g.numVertices(), 0);

	g.createSequenceNode(0x9001, 'c', "read0001", 0x9001, 0x9001, Kmer::FORWARD);

	BOOST_REQUIRE_EQUAL(g.numVertices(), 1);
}

BOOST_AUTO_TEST_CASE (get_vertex_for_hash) {
	SkinnyGraph g(2);
	SkinnyGraph::Vertex v;

	BOOST_REQUIRE_EQUAL(g.numVertices(), 0);

	v = g.createSequenceNode(0x9001, 'c', "read0001", 0x9001, 0x9001, Kmer::FORWARD);

	BOOST_REQUIRE_EQUAL(g.numVertices(), 1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x9001), v);
}

BOOST_AUTO_TEST_CASE (set_vertex_for_hash) {
	SkinnyGraph g(2);
	SkinnyGraph::Vertex v, v2;

	BOOST_REQUIRE_EQUAL(g.numVertices(), 0);
	v = g.createSequenceNode(0x9001, 'c', "read0001", 0x9001, 0x9001, Kmer::FORWARD);
	BOOST_REQUIRE_EQUAL(g.numVertices(), 1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x9001), v);
	v2 = g.createSequenceNode(0x42, 'c', "read0001", 0x9001, 0x9001, Kmer::REVERSE);
	BOOST_REQUIRE_EQUAL(g.numVertices(), 2);
	g.setVertexForHash(0x9001, v2);
	BOOST_REQUIRE_EQUAL(g.numVertices(), 2);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x9001), v2);
}

BOOST_AUTO_TEST_CASE (split) {
	SkinnyGraph::Vertex v1, v2;
	boost::tie(v1, v2) = g.split(v, 5);

	BOOST_REQUIRE_EQUAL(g.node(v1)->sequence(), "acgtu");
	BOOST_REQUIRE_EQUAL(g.node(v2)->sequence(), "acgtu");

	
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0042), v1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0043), v1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0044), v1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0045), v1);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0046), v1);

	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0047), v2);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0048), v2);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x0049), v2);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x004a), v2);
	BOOST_REQUIRE_EQUAL(g.getVertexForHash(0x004b), v2);

	SkinnyGraph::Edge e;
	boost::tie(e, boost::tuples::ignore) = boost::edge(v1, v2, *g.graph());
	// BOOST_REQUIRE_EQUAL(g.edge(e)->getWeight(), 1);
}

BOOST_AUTO_TEST_CASE (add_edge_between_nodes) {
	SkinnyGraph::Vertex next;
	next = g.addEdge(v, v2);
	BOOST_REQUIRE_EQUAL(g.numVertices(), 2);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g.graph()), 1);
	// the next node that we retrieve should be the same as the target of the edge
	BOOST_REQUIRE_EQUAL(next, v2);
}

BOOST_AUTO_TEST_CASE (merge_graphs_invalid_state) {
	BOOST_TEST_CHECKPOINT("going to merge the graphs.");
	try {
		g.merge(g2);
		BOOST_TEST_CHECKPOINT("automatically fail, this test merge shouldn't pass");
		BOOST_REQUIRE_EQUAL(0, 1);
	} catch (InvalidGraphStateException &e) {
		BOOST_TEST_CHECKPOINT("passed, we threw the required exception.");
		BOOST_REQUIRE_EQUAL(1, 1);
	}
}

BOOST_AUTO_TEST_CASE (merge_graphs) {
	BOOST_TEST_CHECKPOINT("adding an edge between nodes in g2 to make valid.");
	g2->addEdge(v3, v4);
	BOOST_TEST_CHECKPOINT("going to merge the graphs.");
	g.merge(g2);
	// now the graph should have 4 nodes with 1 edge
	BOOST_TEST_CHECKPOINT("validating the number of vertices");
	BOOST_REQUIRE_EQUAL(g.numVertices(), 4);
	// the graph should have 1 edge
	BOOST_TEST_CHECKPOINT("validating the number of edges");
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g.graph()), 1);
	// v and v2 shouldn't have been invalidated, so I should still be able to use
	// those handles to validate:
	BOOST_TEST_CHECKPOINT("validating that the two vertices originally in this graph were not modified");
	BOOST_REQUIRE_EQUAL(g.node(v)->sequence(), "acgtuacgtu");
	BOOST_REQUIRE_EQUAL(g.node(v2)->sequence(), "utgcautgca");
	// now I need to get a handle on v3 and v4. v3 should be the node that has one outgoing
	// edge and v4 should be the node that has one incoming edge:
	BOOST_TEST_CHECKPOINT("getting a handle on the two vertices merged into this graph");
	std::pair<std::size_t, SkinnyGraph::Vertex> p;
	BOOST_FOREACH(p, g.getVertices()) {
		if (boost::out_degree(p.second, *g.graph()) == 1) {
			v3 = p.second;
		} else if (boost::in_degree(p.second, *g.graph()) == 1) {
			v4 = p.second;
		}
	}
	BOOST_TEST_CHECKPOINT("validating that the sequence in these nodes wasn't changed");
	BOOST_REQUIRE_EQUAL(g.node(v3)->sequence(), "aaaaaaaaaaa");
	BOOST_REQUIRE_EQUAL(g.node(v4)->sequence(), "ggggggggggg");
	BOOST_TEST_CHECKPOINT("validating that the vertices still have an edge between them");
	// is there still an edge between the new v3 and v4?
	bool edgePresent;
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v3, v4, *g.graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	// check that the hashes have been updated correctly:
	for (std::size_t h = 0x0056; h <= 0x0060; h++) {
		BOOST_TEST_CHECKPOINT("Checking hash [" << std::hex << h << std::dec << "]");
		BOOST_REQUIRE_EQUAL(g.getVertices()[h], v3);
	}
	for (std::size_t h = 0x0061; h <= 0x006b; h++) {
		BOOST_TEST_CHECKPOINT("Checking hash [" << std::hex << h << std::dec << "]");
		BOOST_REQUIRE_EQUAL(g.getVertices()[h], v4);
	}
	BOOST_REQUIRE(v3 != v4);
}

BOOST_AUTO_TEST_CASE (merge_graphs_with_cycle) {
	BOOST_TEST_CHECKPOINT("adding an edge between v3 and v4, then v4 and v3 to make a cycle");
	g2->addEdge(v3, v4);
	g2->addEdge(v4, v3);
	BOOST_TEST_CHECKPOINT("going to merge the graphs");
	g.merge(g2);
	// the graph should still have 4 vertices and should have 2 edges:
	BOOST_TEST_CHECKPOINT("validating the number of vertices");
	BOOST_REQUIRE_EQUAL(g.numVertices(), 4);
	BOOST_TEST_CHECKPOINT("validating the number of edges");
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g.graph()), 2);
	// make sure that the two original vertices were left unchanged:
	BOOST_TEST_CHECKPOINT("checking original vertices");
	BOOST_REQUIRE_EQUAL(g.node(v)->sequence(), "acgtuacgtu");
	BOOST_REQUIRE_EQUAL(g.node(v2)->sequence(), "utgcautgca");
	// now we need to get a handle on v3 and v4. Since they're in a cycle we're going to have to
	// use their sequences to get the handle.
	std::pair<std::size_t, SkinnyGraph::Vertex> p;
	BOOST_FOREACH(p, g.getVertices()) {
		if (g.node(p.second)->sequence() == "aaaaaaaaaaa") {
			v3 = p.second;
		} else if (g.node(p.second)->sequence() == "ggggggggggg") {
			v4 = p.second;
		}
	}
	// even though we already sort of checked this in the above loop, check it again
	BOOST_REQUIRE_EQUAL(g.node(v3)->sequence(), "aaaaaaaaaaa");
	BOOST_REQUIRE_EQUAL(g.node(v4)->sequence(), "ggggggggggg");
	BOOST_TEST_CHECKPOINT("checking that the edges were carried over correctly");
	bool edgePresent;
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v3, v4, *g.graph());
	BOOST_REQUIRE(edgePresent);
	// check the other half of the cycle
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v4, v3, *g.graph());
	BOOST_REQUIRE(edgePresent);
	for (std::size_t h = 0x0056; h <= 0x0060; h++) {
		BOOST_REQUIRE_EQUAL(g.getVertices()[h], v3);
	}
	for (std::size_t h = 0x0061; h <= 0x006b; h++) {
		BOOST_REQUIRE_EQUAL(g.getVertices()[h], v4);
	}
	BOOST_REQUIRE(v3 != v4);
}

BOOST_AUTO_TEST_CASE (merge_graphs_one_node) {
	BOOST_TEST_CHECKPOINT("removing a vertex from g2");
	boost::remove_vertex(v4, *g2->graph());
	BOOST_REQUIRE_EQUAL(g2->numVertices(), 1);

	BOOST_TEST_CHECKPOINT("going to merge the graphs");
	g.merge(g2);
	BOOST_TEST_CHECKPOINT("validating the number of vertices");
	BOOST_REQUIRE_EQUAL(g.numVertices(), 3);
	// this graph should have no edges still:
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g.graph()), 0);
	// we should still have a handle on v and v2
	BOOST_TEST_CHECKPOINT("validating that the two vertices originally in this graph were not modified");
	BOOST_REQUIRE_EQUAL(g.node(v)->sequence(), "acgtuacgtu");
	BOOST_REQUIRE_EQUAL(g.node(v2)->sequence(), "utgcautgca");
	// we should be able to get v3 as being the only node that's not v or v2:
	std::pair<std::size_t, SkinnyGraph::Vertex> p;
	BOOST_FOREACH(p, g.getVertices()) {
		if (p.second != v && p.second != v2) {
			v3 = p.second;
		}
	}
	BOOST_REQUIRE_EQUAL(g.node(v3)->sequence(), "aaaaaaaaaaa");
	BOOST_TEST_CHECKPOINT("checking that all hash values were updated correctly");
	for (std::size_t h = 0x0056; h <= 0x0060; h++) {
		BOOST_REQUIRE_EQUAL(g.getVertices()[h], v3);
	}
}

BOOST_AUTO_TEST_CASE (merge_no_edges) {
	BOOST_TEST_CHECKPOINT("going to merge two nodes that have no edges");
	std::string mergedSequence = g3->node(v5)->sequence() + g3->node(v6)->sequence();
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 0);
	g3->addEdgeOrMerge(v5, v6);
	// should still have no edges
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 0);
	// should have two nodes:
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 2);
	BOOST_REQUIRE_EQUAL(g3->node(v6)->sequence(), mergedSequence);
}

BOOST_AUTO_TEST_CASE (merge_neighbours) {
	std::string mergedSequence = g3->node(v5)->sequence() + g3->node(v6)->sequence();
	boost::add_edge(v5, v6, *g3->graph());
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 1);

	BOOST_TEST_CHECKPOINT("going to merge two nodes that share a single edge");
	g3->addEdgeOrMerge(v5, v6);
	// should now have no edges and 2 nodes, v5 and v6 should have been merged
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 2);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 0);

	BOOST_REQUIRE_EQUAL(g3->node(v6)->sequence(), mergedSequence);
}

BOOST_AUTO_TEST_CASE (no_merge_add_edge_possible_neighbours) {
	// add edges such that v5 has one outgoing edge and v6 has one incoming edge, but
	// they do not share this edge.
	boost::add_edge(v5, v7, *g3->graph());
	boost::add_edge(v7, v6, *g3->graph());

	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 2);

	BOOST_TEST_CHECKPOINT("Adding an edge between nodes, no merging should take place.");
	g3->addEdgeOrMerge(v5, v6);
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 3);
}

BOOST_AUTO_TEST_CASE (no_merge_add_edge) {
	// add an edge to the graph such that one of the vertices already has an outgoing edge
	// and the other vertex does not.
	boost::add_edge(v5, v7, *g3->graph());
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 1);

	BOOST_TEST_CHECKPOINT("Adding an edge between nodes, no merging should take place.");
	g3->addEdgeOrMerge(v5, v6);
	// no merging, so 3 nodes and 2 edges should be present
	BOOST_REQUIRE_EQUAL(g3->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(boost::num_edges(*g3->graph()), 2);
}

BOOST_AUTO_TEST_SUITE_END()
#endif // SKINNY_GRAPH_TEST_CC
