/* 
 * File:   QAssemblerTestSuite.cc
 * Author: fbristow
 *
 * Created on December 21, 2011
 */
#ifndef QASSEMBLER_TEST_SUITE_CC
#define QASSEMBLER_TEST_SUITE_CC

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#ifdef USE_LOG4CXX
#include <log4cxx/propertyconfigurator.h>
#endif

struct QAssemblerTestConfig {
	QAssemblerTestConfig() {
#ifdef USE_LOG4CXX
		log4cxx::PropertyConfigurator::configure("log.config");
#endif
	}

};

BOOST_GLOBAL_FIXTURE(QAssemblerTestConfig)

#endif // QASSEMBLER_TEST_SUITE_CC
