/* 
 * File:   ReadSizeException.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef READ_SIZE_EXCEPTION_HH
#define READ_SIZE_EXCEPTION_HH

#include <stdexcept>

class ReadSizeException : public std::runtime_error {
public:
	explicit ReadSizeException(const std::string& m) : std::runtime_error(m) {}
	virtual ~ReadSizeException() throw() {};
};

#endif // READ_SIZE_EXCEPTION_HH

