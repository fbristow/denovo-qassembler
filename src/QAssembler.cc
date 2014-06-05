/**
 * File:   BuildGraph.cc
 * Author: fbristow
 *
 * Created on January 17, 2012
 */

#include <boost/cstdint.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <zlib.h>
#include <stdio.h>

#include "QAssemblerConfig.h"

#include "Graph/HeftyGraph.hh"
#include "IO/GraphWriter.hh"
#include "IO/FastaStream.hh"
#include "PreHash/PreHash.hh"
#include "PathBuilder/Proportional/ProportionalPathBuilder.hh"
#include "PathBuilder/Markov/MarkovPathBuilder.hh"
#include "Abundance/MarkovAbundance/ForwardAlgorithmAbundance.hh"
#include "Abundance/MarkovAbundance/MarkovChainAbundance.hh"
#include "Sequence/Sequence.hh"
#include "Exception/KmerLengthException.hh"
#include "QAssemblerParameterException.hh"

#include "Logging/Logging.hh"

/** input parameters */
std::string inputSequences;
/** graph construction parameters */
std::size_t kmerLength = 31;
bool preHash = false;
HeftyGraph::TrackReads trackReads = HeftyGraph::DONT_TRACK_READS;
/** graph modification parameters */
std::size_t aggressiveLength = 0;
std::size_t aggressiveEdgeWeight = 0;
/** sequence construction parameters */
double epsilon = 0.1;
std::string pathMethod = "proportional";
/** abundance estimation parameters */
std::string abundanceMethod = "";
/** output parameters */
bool printGraph = false;
bool printSequences = false;
std::size_t minimumLength = 0;
std::string configFile = "./log.config";
std::string sequenceDir = "./sequences";
std::string graphDir = "./graphs";
/*********************/

DECLARE_LOG(logger, "qassembler.QAssembler");

int parseArgs(int, char**);
boost::shared_ptr<HeftyGraph> buildGraph(boost::shared_ptr<PreHash>);

int main(int argc, char **argv) {
	if (parseArgs(argc, argv)) {
		exit(0);
	}
	CONFIGURE_LOG(configFile);

	boost::shared_ptr<PreHash> preHasher;
	preHasher = boost::make_shared<PreHash>(kmerLength, 20000);
	std::size_t totalReadsProcessed = 0;

	if (preHash) {
		INFO(logger, "Pre-hashing reads.");
		FastaStream fastaStream(inputSequences);
		boost::progress_display progress(fastaStream.seqCount());
		while (boost::shared_ptr<Sequence> read = fastaStream.nextSeq()) {
			preHasher->addRead(read);
			totalReadsProcessed++;
			++progress;
		}
		INFO(logger, "Pre-hashed [" << totalReadsProcessed << "] reads.");
	} else {
		preHasher.reset();
	}

	boost::shared_ptr<HeftyGraph> g = buildGraph(preHasher);

	if (aggressiveLength) {
		INFO(logger, "Removing edges from all graphs with single nodes with length less than [" << aggressiveLength << "]");
		g->removeGraphsShorterThan(aggressiveLength);
		INFO(logger, "[" << g->numGraphs() << "] graphs remain after removal.");
	}

	if (aggressiveEdgeWeight && !preHash) {
		INFO(logger, "Removing edges from all graphs below threshold weight [" << aggressiveEdgeWeight << "].");
		g->removeEdgesBelowThreshold(aggressiveEdgeWeight);
	}

	if (printGraph) {
		INFO(logger, "Writing graphs to files...");
		std::set<boost::shared_ptr<SkinnyGraph> > graphs = g->getGraphs();
		BOOST_FOREACH (boost::shared_ptr<SkinnyGraph> graph, graphs) {
			std::string filename = boost::lexical_cast<std::string>(graph->getId());
		        filename += ".dot";
			GraphWriter gw(graph, filename, graphDir);
			gw.write();
		}
	}

	if (printSequences) {
		boost::shared_ptr<Abundance> abundanceEstimator;

		INFO(logger, "Generating sequences into directory [" << sequenceDir << "]");
		std::set<boost::shared_ptr<SkinnyGraph> > graphs = g->getGraphs();
		boost::filesystem::create_directory(sequenceDir);
		boost::unordered_map<boost::shared_ptr<SkinnyGraph>, boost::unordered_set<std::string> > paths;
		//boost::unordered_set<std::string> paths;
		std::size_t sequenceCount = 0;

		boost::progress_display progress(graphs.size());
		BOOST_FOREACH (boost::shared_ptr<SkinnyGraph> graph, graphs) {
			DEBUG(logger, "Generating paths for graph [" << graph->getId() << "]");
			boost::shared_ptr<PathBuilder> pathBuilder;
			if (pathMethod == "proportional") {
				pathBuilder = boost::make_shared<ProportionalPathBuilder>(graph, epsilon);
			} else if (pathMethod == "markov") {
				pathBuilder = boost::make_shared<MarkovPathBuilder>(graph);
			} else if (pathMethod == "random") {
				throw QAssemblerParameterException ("random path builder is unimplemented.");
			}
			std::string filename = sequenceDir + "/" + boost::lexical_cast<std::string>(graph->getId()) + ".fna";
			std::ofstream sequenceFile(filename.c_str());
			boost::unordered_set<PathBuilder::Path> paths = pathBuilder->buildPaths();
			boost::unordered_map<std::string, double> abundances;
			boost::unordered_set<std::string> sequencePaths;
			BOOST_FOREACH (PathBuilder::Path p, paths) {
				std::string sequence = p[0]->fullSequence();
				for (std::size_t i = 1; i < p.size(); i++) {
					sequence += p[i]->sequence();
				}
				TRACE(logger, "Constructed sequence [" << sequence << "] for graph [" << graph->getId() << "].");
				// don't bother reporting sequences less than k
				if (sequence.size() > kmerLength && (!minimumLength || sequence.size() >= minimumLength)) {
					sequencePaths.insert(sequence);
				}
			}
			graph->resetEdgeWeights();
			if (abundanceMethod != "") {
				if (abundanceMethod == "markov-chain") {
					abundanceEstimator = boost::make_shared<MarkovChainAbundance>(g, sequencePaths);
				} else if (abundanceMethod == "forward-algorithm") {
					abundanceEstimator = boost::make_shared<ForwardAlgorithmAbundance>(g, sequencePaths);
				}

				abundances = abundanceEstimator->computeAbundances();
			}
			BOOST_FOREACH (std::string sequence, sequencePaths) {
				std::string abundance = "";
				if (abundanceEstimator) {
					abundance = " (" + abundanceMethod + ": " +
						boost::lexical_cast<std::string>(abundances[sequence]) + ")";
				}
				sequenceFile << ">" << ++sequenceCount << "(" << sequence.size() << "bp)" << abundance << std::endl;
				sequenceFile << sequence << std::endl << std::endl;
			}
			INFO(logger, "Generated [" << paths.size() << "] sequences for graph [" << graph->getId() << "].");
			++progress;
		}
	}

	return 0;
}

boost::shared_ptr<HeftyGraph> buildGraph(boost::shared_ptr<PreHash> preHasher) {
	std::size_t totalReadsProcessed = 1;
	INFO(logger, "Constructing graph...");
	boost::shared_ptr<HeftyGraph> g = boost::make_shared<HeftyGraph>(kmerLength, trackReads, preHasher, aggressiveEdgeWeight);
	try {
		FastaStream fastaStream(inputSequences);
		boost::progress_display progress(fastaStream.seqCount());
		while (boost::shared_ptr<Sequence> read = fastaStream.nextSeq()) {
			if (read->getLength() >= kmerLength) {
				g->addReadToGraph(read);
			}
			totalReadsProcessed++;
			++progress;
		}
	} catch (std::exception &e) {
		FATAL(logger, e.what());
	}

	g->lockEdgeWeights();

	INFO(logger, "Created [" << g->numGraphs() << "] graphs.");
	return g;
}

int parseArgs(int argc, char **argv) {
	namespace boost_po = boost::program_options;
	boost_po::options_description desc("Options");
	bool trackReadsBool = false;
	desc.add_options()
		("help,h", "list all options.")
		("input-sequences,i", boost_po::value<std::string>(&inputSequences),
		 	 "fasta/fastq file containing reads (required).")
		("kmer-size,k", boost_po::value<std::size_t>(&kmerLength)->default_value(31),
			 "set the k-mer size.")
		("pre-hash,p", boost_po::value<bool>(&preHash)->default_value(false)->zero_tokens(),
		 	 "pre-hash the reads to guide graph construction.")
		("aggressive-edge-removal,a", boost_po::value<std::size_t>(&aggressiveEdgeWeight)->default_value(0),
		 	 "remove edges from graphs where the edge weight is below a specified threshold.")
		("print-graphs,g", boost_po::value<bool>(&printGraph)->default_value(false)->zero_tokens(),
		 	 "print the graphs that were generated.")
		("graph-dir", boost_po::value<std::string>(&graphDir)->default_value("graphs"),
		 	 "directory to dump DOT formatted graphs.")
		("minimum-bases,m", boost_po::value<std::size_t>(&aggressiveLength)->default_value(0),
		 	 "remove graphs that have only one node with a length less than specified.")
		("sequences,s", boost_po::value<bool>(&printSequences)->default_value(false)->zero_tokens(),
		 	 "generate sequences and print to files.")
		("sequence-dir", boost_po::value<std::string>(&sequenceDir)->default_value("sequences"),
		 	 "directory to dump reconstructed sequences.")
		("path-method", boost_po::value<std::string>(&pathMethod)->default_value("proportional"),
			 "method used to generate paths through the graph (one of proportional, markov or random)")
		("epsilon,e", boost_po::value<double>(&epsilon)->default_value(0.01),
		 	 "allowable difference between paths during path generation.")
		("minimum-length,l", boost_po::value<std::size_t>(&minimumLength)->default_value(0),
		 	 "only report sequence longer than the specified length.")
		("abundance-method", boost_po::value<std::string>(&abundanceMethod)->default_value(""),
		 	 "the abundance estimation method (one of markov-chain or forward-algorithm)")
#ifdef USE_LOG4CXX
		("log-config", boost_po::value<std::string>(&configFile)->default_value("log.config"),
		 	 "location of config file for logging (log4cxx).")
#endif
		;

	boost_po::variables_map vm;
	try {
		boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
		boost_po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << std::endl;
		}

		if (trackReadsBool) {
			trackReads = HeftyGraph::TRACK_READS;
		} else {
			trackReads = HeftyGraph::DONT_TRACK_READS;
		}

		if (kmerLength % 2 == 0) {
			throw KmerLengthException("k-mer length must be odd.");
		}

		if (inputSequences == "") {
			throw QAssemblerParameterException("input-sequences is a required option.");
		}

		if (abundanceMethod != "" && abundanceMethod != "markov-chain" && abundanceMethod != "forward-algorithm") {
			throw QAssemblerParameterException("abundance method must be one of markov-chain or forward-algorithm");
		}

		if (pathMethod != "" && pathMethod != "proportional" && pathMethod != "markov" && pathMethod != "random") {
			throw QAssemblerParameterException("path method must be one of proportional, markov or random");
		}
	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << desc << std::endl;
		return 1;
	}

	return 0;
}
