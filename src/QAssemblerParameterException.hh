/* 
 * File:   QAssemblerParameterException.hh
 * Author: fbristow
 *
 * Created on June 27, 2012
 */
#ifndef QASSEMBLER_PARAMETER_EXCEPTION_HH
#define QASSEMBLER_PARAMETER_EXCEPTION_HH

#include <stdexcept>

class QAssemblerParameterException : public std::runtime_error {
public:
	explicit QAssemblerParameterException(const std::string& m) : std::runtime_error(m) {}
	virtual ~QAssemblerParameterException() throw() {};
};

#endif // QASSEMBLER_PARAMETER_EXCEPTION_HH

