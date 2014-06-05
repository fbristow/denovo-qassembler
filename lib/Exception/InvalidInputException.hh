/* 
 * File:   InvalidInputException.hh
 * Author: fbristow
 *
 * Created on June 28, 2012
 */
#ifndef INVALID_INPUT_EXCEPTION_HH
#define INVALID_INPUT_EXCEPTION_HH

#include <stdexcept>

class InvalidInputException : public std::runtime_error {
public:
	explicit InvalidInputException(const std::string& m) : std::runtime_error(m) {}
	virtual ~InvalidInputException() throw() {};
};

#endif // INVALID_INPUT_EXCEPTION_HH

