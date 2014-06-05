/**
 * File:   Logging.hh
 * Author: fbristow
 *
 * Created on June 16, 2012
 */

#ifndef LOGGING_HH
#define LOGGING_HH

#include "QAssemblerConfig.h"

#ifdef USE_LOG4CXX

#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#define DECLARE_LOG(var, name) static log4cxx::LoggerPtr var(log4cxx::Logger::getLogger(name))
#define CONFIGURE_LOG(configFile) log4cxx::PropertyConfigurator::configure(configFile)

#define TRACE(logger, msg) LOG4CXX_TRACE(logger, msg)
#define DEBUG(logger, msg) LOG4CXX_DEBUG(logger, msg)
#define INFO(logger, msg)  LOG4CXX_INFO (logger, msg)
#define WARN(logger, msg)  LOG4CXX_WARN (logger, msg)
#define ERROR(logger, msg) LOG4CXX_ERROR(logger, msg)
#define FATAL(logger, msg) LOG4CXX_FATAL(logger, msg)

#else
// logging should be conditionally enabled at compile time.The following re-define the log4cxx logging
// macros to be no-ops or print to stderr statements. The do {} while(0) part enforces a semi-colon after
// the logging line (in case someone does something tricky like trying to put
// a block after the logging line, see: http://stackoverflow.com/questions/1306611/how-do-i-implement-no-op-macro-or-template-in-c#1306690)
#define DECLARE_LOG(var, name) ; 
#define CONFIGURE_LOG(configFile) do {} while(0)

#define FATAL(L, X) std::cerr << "FATAL: " << X << std::endl
#define ERROR(L, X) std::cerr << "ERROR: " << X << std::endl
#define WARN(L, X)  std::cerr << "WARN:  " << X << std::endl
#if LOGGING_SEVERITY >= 0
#define INFO(L, X)  std::cerr << "INFO:  " << X << std::endl
#else
#define INFO(L, X)  do {} while(0)
#endif

#if LOGGING_SEVERITY >= 1
#define DEBUG(L, X) std::cerr << "DEBUG: " << X << std::endl
#else
#define DEBUG(L, X) do {} while(0)
#endif

#if LOGGING_SEVERITY >= 2
#define TRACE(L, X) std::cerr << "TRACE: " << X << std::endl
#else
#define TRACE(L, X) do {} while(0)
#endif

#endif // USE_LOG4CXX


#endif // LOGGING_HH
