/* 
 * File:   InvalidGraphStateException.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef INVALID_GRAPH_STATE_EXCEPTION_HH
#define INVALID_GRAPH_STATE_EXCEPTION_HH

#include <stdexcept>

class InvalidGraphStateException : public std::runtime_error {
public:
	explicit InvalidGraphStateException(const std::string& m) : std::runtime_error(m) {}
	virtual ~InvalidGraphStateException() throw() {};
};

#endif // INVALID_GRAPH_STATE_EXCEPTION_HH

