/* 
 * File:   Util.hh
 * Author: fbristow
 *
 * Created on Nov 2, 2011
 */
#ifndef UTIL_HH
#define UTIL_HH

#include <boost/functional/hash.hpp>
#include <boost/cstdint.hpp>

namespace qassembler {

	inline std::size_t hash(std::string value) {
		return boost::hash_value(value);
	}

} // namespace

#endif
