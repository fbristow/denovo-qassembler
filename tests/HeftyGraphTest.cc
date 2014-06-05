/* 
 * File:   HeftyGraphTest.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef HEFTY_GRAPH_TEST_CC
#define HEFTY_GRAPH_TEST_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include "Graph/HeftyGraph.hh"

struct HeftyGraphFixture {
	HeftyGraphFixture() : read1(boost::make_shared<Sequence>()),
			      read2(boost::make_shared<Sequence>()),
			      read3(boost::make_shared<Sequence>()),
			      read4(boost::make_shared<Sequence>())
	{
		read1->setSequence("ACCT");
		read1->setQual("++++");
		read1->setID(0x9001);
		read1->setName("read1");

		read2->setSequence("CCTA");
		read2->setQual("++++");
		read2->setID(0x0042);
		read2->setName("read2");

		read3->setSequence("CCTT");
		read3->setQual("++++");
		read3->setID(0x0043);
		read3->setName("read3");

		read4->setSequence("TTGG");
		read4->setQual("++++");
		read4->setID(0x0044);
		read4->setName("read4");
	}

	boost::shared_ptr<Sequence> read1;
	boost::shared_ptr<Sequence> read2;
	boost::shared_ptr<Sequence> read3;
	boost::shared_ptr<Sequence> read4;
};

BOOST_FIXTURE_TEST_SUITE (hefty_graph, HeftyGraphFixture)

BOOST_AUTO_TEST_CASE (constructor_test) {
	HeftyGraph hg(3);
	
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 0);
}

BOOST_AUTO_TEST_CASE (add_read_to_graph_num_graphs) {
	HeftyGraph hg(3);
	
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 0);
	BOOST_TEST_CHECKPOINT("Adding read1 to graph");
	hg.addReadToGraph(read1);
	// after adding the read to the graph, we should have 2 graphs, one for the
	// forward direction and one for the reverse direction.
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);
}

BOOST_AUTO_TEST_CASE (add_distinct_reads_to_graph) {
	HeftyGraph hg (3);

	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 0);
	BOOST_TEST_CHECKPOINT("Adding read3 to graph");
	hg.addReadToGraph(read3);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);
	BOOST_TEST_CHECKPOINT("Adding read4 to graph");
	hg.addReadToGraph(read4);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 4);
}

BOOST_AUTO_TEST_CASE (add_overlapping_reads_in_same_graph) {
	HeftyGraph hg(5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();
       	boost::shared_ptr<Sequence> readC = boost::make_shared<Sequence>();
	readA->setSequence("ACTGGTAAATGTATG");
	readA->setQual("+++++++++++++++");
	readA->setID(0x0042);
	readA->setName("readA");

	readB->setSequence("ACTGGTAATG");
	readB->setQual("++++++++++");
	readB->setID(0x0043);
	readB->setName("readB");

	readC->setSequence("TAATGCGTAAA");
	readC->setQual("+++++++++++");
	readC->setID(0x0044);
	readC->setName("readC");

	BOOST_TEST_CHECKPOINT("loading readA into graph");
	hg.addReadToGraph(readA);
	BOOST_TEST_CHECKPOINT("loading readB into graph");
	hg.addReadToGraph(readB);

	/**
+----------+  1   +---------+
| ACTGGTAA | ---> | ATGTATG |
+----------+      +---------+
  |
  | 1
  v
+----------+
|    TG    |
+----------+
	 */
	// this should leave us with 1 graph, 3 nodes in the graph (the start of a bulge):
	// I don't care about reverse complement right now
	BOOST_TEST_CHECKPOINT("getting the constructed graph");
	HeftyGraph::ReadLookup graphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = graphs.get(0x0042);
	BOOST_REQUIRE_EQUAL(forward->numVertices(), 3);
	// this graph should only have 2 edges between those nodes
	BOOST_REQUIRE_EQUAL(boost::num_edges(*forward->graph()), 2);
	BOOST_FOREACH(SkinnyGraph::Edge e, boost::edges(*forward->graph())) {
		BOOST_REQUIRE_EQUAL(forward->edge(e)->getWeight(), 1);
	}

	BOOST_TEST_CHECKPOINT("loading readC into graph");
	hg.addReadToGraph(readC);

	// after adding readC, the beginning of the bulge should have been completed, so
	/**
             1
  +-----------------------------+
  |                             v
+----------+  1   +----+  1   +---------+
| ACTGGTAA | ---> | TG | ---> | ATGTATG |
+----------+      +----+      +---------+
	 */	
	// we should not have added any more nodes to the graph:
	BOOST_REQUIRE_EQUAL(forward->numVertices(), 3);

	// now can validate the structure of the graph
	// the graph should now have 3 edges joining the vertices.
	BOOST_REQUIRE_EQUAL(boost::num_edges(*forward->graph()), 3);
}

BOOST_AUTO_TEST_CASE (add_read_to_graph_overlapping_reads) {
	HeftyGraph hg(3, HeftyGraph::TRACK_READS);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 0);
	BOOST_TEST_CHECKPOINT("Adding read1 to graph");
	hg.addReadToGraph(read1);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);
	BOOST_TEST_CHECKPOINT("Adding read2 to graph");
	hg.addReadToGraph(read2);
	// we should still have two graphs, the pair of reads are overlapping
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);
	// now let's make sure that the graphs only have one node
	BOOST_TEST_CHECKPOINT("Getting graphs");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	HeftyGraph::ReadLookup reverseGraphs = hg.getReverseReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x9001);
	boost::shared_ptr<SkinnyGraph> reverse = reverseGraphs.get(0x9001);

	BOOST_TEST_CHECKPOINT("checking that both the forward and reverse graphs only have one vertex");
	// both graphs should have only one vertex
	BOOST_REQUIRE_EQUAL(forward->numVertices(), 1);
	BOOST_TEST_CHECKPOINT("checking that both the forward and reverse graphs only have one vertex");
	BOOST_REQUIRE_EQUAL(reverse->numVertices(), 1);

	// get the vertex and sequence node from the forward graph
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v = vertices.begin()->second;
	boost::shared_ptr<SequenceNode> n = forward->node(v);
	BOOST_TEST_CHECKPOINT("Checking that the forward sequence is correct");
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "ACCTA");

	// get the vertex and sequence node from the reverse graph
	SkinnyGraph::Vertex v2 = reverse->getVertices().begin()->second;
	n = reverse->node(v2);
	// reverse complement of the above sequence:
	BOOST_TEST_CHECKPOINT("Checking that the reverse sequence is correct");
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "TAGGT");
}

BOOST_AUTO_TEST_CASE (merge_graphs_and_nodes) {
	// this test case is going to cover the situation where we first add two reads that
	// belong in separate graphs, then add a third read that overlaps with both of those
	// reads such that the graphs will be merged, then the nodes should also be merged.
	
	HeftyGraph hg (5, HeftyGraph::TRACK_READS);
	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readC = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAAAA");
	readA->setQual("+++++");
	readA->setID(0x42);

	readB->setName("readB");
	readB->setSequence("CCCCC");
	readB->setQual("+++++");
	readB->setID(0x43);

	readC->setName("readC");
	readC->setSequence("AAAAACCCCC");
	readC->setQual("++++++++++");
	readC->setID(0x9001);

	BOOST_TEST_CHECKPOINT("Adding first read to graph");
	hg.addReadToGraph(readA);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);
	BOOST_TEST_CHECKPOINT("Adding second read to graph");
	hg.addReadToGraph(readB);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 4);
	BOOST_TEST_CHECKPOINT("Adding joining read to graph");
	hg.addReadToGraph(readC);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2); // the joining read should have merged both
						// the forward and reverse graphs.
						
	BOOST_TEST_CHECKPOINT("Getting forward graph.");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);
	
	BOOST_REQUIRE_EQUAL(forward->numVertices(), 1);

	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v = vertices[qassembler::hash("AAAAA")];
	boost::shared_ptr<SequenceNode> n = forward->node(v);
	BOOST_TEST_CHECKPOINT("Checking that the forward sequence is correct");
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAAAACCCCC");
}

BOOST_AUTO_TEST_CASE (add_edge_at_middle_of_dest) {
	// this test covers the situation where we add a read that requires splitting the target
	// node and adding an edge.
	HeftyGraph hg (5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAAACCCTT");
	readA->setQual("+++++++++");
	readA->setID(0x42);

	readB->setSequence("GACCCTT");
	readB->setQual("+++++++");
	readB->setID(0x43);

	BOOST_TEST_CHECKPOINT("Adding first read to graph");
	hg.addReadToGraph(readA);
	BOOST_TEST_CHECKPOINT("Adding second read to graph");
	hg.addReadToGraph(readB);
	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2); // graphs should have merged.

	BOOST_TEST_CHECKPOINT("Getting forward graph");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);

	BOOST_REQUIRE_EQUAL(forward->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(forward->numEdges(), 2);

	BOOST_TEST_CHECKPOINT("Checking validity of node sequences.");
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v1 = vertices[qassembler::hash("AAAAC")];
	boost::shared_ptr<SequenceNode> n = forward->node(v1);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAAACCC");

	SkinnyGraph::Vertex v2 = vertices[qassembler::hash("GACCC")];
	n = forward->node(v2);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "GACCC");

	SkinnyGraph::Vertex v3 = vertices[qassembler::hash("ACCCT")];
	n = forward->node(v3);
	BOOST_REQUIRE_EQUAL(n->sequence(), "TT");

	BOOST_TEST_CHECKPOINT("Validating edges");
	bool edgePresent;
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v2, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
}

BOOST_AUTO_TEST_CASE (add_kmer_at_middle) {
	// this test validates what happens when a new hash being added to the graph
	// is a variant that deviates from an existing node in the middle of the node.
	// this situation can occur during the very initial construction phases of the graph
	// when we encounter a homopolymer region.
	HeftyGraph hg(5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAACCCCGT");
	readA->setQual("+++++++++");
	readA->setID(0x42);

	readB->setName("readB");
	readB->setSequence("AAACCCGA");
	readB->setQual("++++++++");
	readB->setID(0x43);

	BOOST_TEST_CHECKPOINT("Adding first read to graph");
	hg.addReadToGraph(readA);
	BOOST_TEST_CHECKPOINT("Adding second read to graph");
	hg.addReadToGraph(readB);

	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2); // should never have added any other graphs

	BOOST_TEST_CHECKPOINT("Getting forward graph");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);

	BOOST_REQUIRE_EQUAL(forward->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(forward->numEdges(), 2);

	BOOST_TEST_CHECKPOINT("Checking validity of node sequences");
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v1 = vertices[qassembler::hash("AAACC")];
	boost::shared_ptr<SequenceNode> n = forward->node(v1);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAACCC");

	SkinnyGraph::Vertex v2 = vertices[qassembler::hash("CCCGT")];
	n = forward->node(v2);
	BOOST_REQUIRE_EQUAL(n->sequence(), "CGT");

	SkinnyGraph::Vertex v3 = vertices[qassembler::hash("CCCGA")];
	n = forward->node(v3);
	BOOST_REQUIRE_EQUAL(n->sequence(), "GA");

	BOOST_TEST_CHECKPOINT("Checking edges");
	bool edgePresent;
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v2, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
}

BOOST_AUTO_TEST_CASE (add_kmer_at_middle_with_cycle) {
	HeftyGraph hg(5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAACCCCGT");
	readA->setQual("+++++++++");
	readA->setID(0x42);

	readB->setName("readB");
	readB->setSequence("AAACCCGT");
	readB->setQual("++++++++");
	readB->setID(0x43);

	BOOST_TEST_CHECKPOINT("Adding first read to graph");
	hg.addReadToGraph(readA);
	BOOST_TEST_CHECKPOINT("Adding second read to graph");
	hg.addReadToGraph(readB);

	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);

	BOOST_TEST_CHECKPOINT("Getting forward graph");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);

	BOOST_REQUIRE_EQUAL(forward->numVertices(), 4);
	BOOST_REQUIRE_EQUAL(forward->numEdges(), 4);

	BOOST_TEST_CHECKPOINT("Checking validity of node sequences");
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v1 = vertices[qassembler::hash("AAACC")];
	boost::shared_ptr<SequenceNode> n = forward->node(v1);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAACCC");

	SkinnyGraph::Vertex v2 = vertices[qassembler::hash("ACCCC")];
	n = forward->node(v2);
	BOOST_REQUIRE_EQUAL(n->sequence(), "CG");

	SkinnyGraph::Vertex v3 = vertices[qassembler::hash("CCCGT")];
	n = forward->node(v3);
	BOOST_REQUIRE_EQUAL(n->sequence(), "T");

	SkinnyGraph::Vertex v4 = vertices[qassembler::hash("ACCCG")];
	n = forward->node(v4);
	BOOST_REQUIRE_EQUAL(n->sequence(), "G");

	BOOST_TEST_CHECKPOINT("Checking edges");
	bool edgePresent;

	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v2, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v4, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v2, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v4, v3, *forward->graph());
}

BOOST_AUTO_TEST_CASE (add_kmer_at_middle_of_both) {
	HeftyGraph hg(5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readC = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAAATTCCCC");
	readA->setQual("++++++++++");
	readA->setID(0x42);

	readB->setName("readB");
	readB->setSequence("CCAATTAAAA");
	readB->setQual("++++++++++");
	readB->setID(0x43);

	readC->setName("readC");
	readC->setSequence("AAAATTAAA");
	readC->setQual("+++++++++");
	readC->setID(0x44);

	BOOST_TEST_CHECKPOINT("Adding first read to graph.");
	hg.addReadToGraph(readA);

	BOOST_TEST_CHECKPOINT("Adding second read to graph.");
	hg.addReadToGraph(readB);

	BOOST_TEST_CHECKPOINT("Adding third read to graph.");
	hg.addReadToGraph(readC);

	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2); // I expect this to be 2 graphs after all is said and done

	BOOST_TEST_CHECKPOINT("Getting forward graph");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);

	BOOST_REQUIRE_EQUAL(forward->numVertices(), 4);
	BOOST_REQUIRE_EQUAL(forward->numEdges(), 3);
	
	BOOST_TEST_CHECKPOINT("Checking validity of node sequences");
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v1 = vertices[qassembler::hash("AAAAT")];
	boost::shared_ptr<SequenceNode> n = forward->node(v1);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAAATT");

	SkinnyGraph::Vertex v2 = vertices[qassembler::hash("CCAAT")];
	n = forward->node(v2);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "CCAATT");

	SkinnyGraph::Vertex v3 = vertices[qassembler::hash("AATTC")];
	n = forward->node(v3);
	BOOST_REQUIRE_EQUAL(n->sequence(), "CCCC");

	SkinnyGraph::Vertex v4 = vertices[qassembler::hash("AATTA")];
	n = forward->node(v4);
	BOOST_REQUIRE_EQUAL(n->sequence(), "AAAA");

	BOOST_TEST_CHECKPOINT("Checking edges");
	bool edgePresent;

	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v2, v4, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v4, *forward->graph());
}

BOOST_AUTO_TEST_CASE (add_kmer_at_middle_of_source) {
	HeftyGraph hg (5, HeftyGraph::TRACK_READS);

	boost::shared_ptr<Sequence> readA = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readB = boost::make_shared<Sequence>();
	boost::shared_ptr<Sequence> readC = boost::make_shared<Sequence>();

	readA->setName("readA");
	readA->setSequence("AAAAATTCCCC");
	readA->setQual("+++++++++++");
	readA->setID(0x42);

	readB->setName("readB");
	readB->setSequence("AAAAC");
	readB->setQual("+++++");
	readB->setID(0x43);
	
	readC->setName("readC");
	readC->setSequence("AAAAAC");
	readC->setQual("++++++");
	readC->setID(0x44);

	BOOST_TEST_CHECKPOINT("Adding first read to graph.");
	hg.addReadToGraph(readA);

	BOOST_TEST_CHECKPOINT("Adding second read to graph.");
	hg.addReadToGraph(readB);

	BOOST_TEST_CHECKPOINT("Adding third read to graph.");
	hg.addReadToGraph(readC);

	BOOST_REQUIRE_EQUAL(hg.numGraphs(), 2);

	BOOST_TEST_CHECKPOINT("Getting forward graph");
	HeftyGraph::ReadLookup forwardGraphs = hg.getForwardReads();
	boost::shared_ptr<SkinnyGraph> forward = forwardGraphs.get(0x42);

	BOOST_REQUIRE_EQUAL(forward->numVertices(), 3);
	BOOST_REQUIRE_EQUAL(forward->numEdges(), 2);

	BOOST_TEST_CHECKPOINT("Checking validity of node sequences");
	boost::unordered_map<std::size_t, SkinnyGraph::Vertex> vertices = forward->getVertices();
	SkinnyGraph::Vertex v1 = vertices[qassembler::hash("AAAAA")];
	boost::shared_ptr<SequenceNode> n = forward->node(v1);
	BOOST_REQUIRE_EQUAL(n->fullSequence(), "AAAAA");

	SkinnyGraph::Vertex v2 = vertices[qassembler::hash("AAAAC")];
	n = forward->node(v2);
	BOOST_REQUIRE_EQUAL(n->sequence(), "C");

	SkinnyGraph::Vertex v3 = vertices[qassembler::hash("AAATT")];
	n = forward->node(v3);
	BOOST_REQUIRE_EQUAL(n->sequence(), "TTCCCC");

	BOOST_TEST_CHECKPOINT("Checking edges.");
	bool edgePresent;

	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v3, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
	boost::tie(boost::tuples::ignore, edgePresent) = boost::edge(v1, v2, *forward->graph());
	BOOST_REQUIRE_EQUAL(edgePresent, true);
}


BOOST_AUTO_TEST_SUITE_END()

#endif // HEFTY_GRAPH_TEST
