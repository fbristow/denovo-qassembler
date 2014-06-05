/* 
 * File:   GraphLookup.hh
 * Author: fbristow
 *
 * Created on January 25, 2012
 */
#ifndef GRAPH_LOOKUP_HH
#define GRAPH_LOOKUP_HH

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

template <class Key, class Value> class GraphLookup {
public:
	GraphLookup() {};

	typedef boost::unordered_map<Key, Value> Key2Value;
	typedef boost::unordered_map<Value, boost::unordered_set<Key> > Value2Key;

	typedef typename Key2Value::iterator KeyIterator;
	typedef typename Value2Key::iterator ValueIterator;
	typedef typename Key2Value::value_type KeyType;
	typedef typename Value2Key::value_type ValueType;

	void put(Key k, Value v){
		this->key2value[k] = v;
		this->value2key[v].insert(k);
	}
	
	Value get(Key k) {
		return this->key2value[k];
	}
	
	boost::unordered_set<Key> get(Value v) {
		return this->value2key[v];
	}
	
	std::size_t count(Key k){
		return this->key2value.count(k);
	}
	
	std::size_t count(Value v){
		return this->value2key[v].size();
	}
	
	void clear(Key k) {
		Value v = this->key2value[k];
		this->key2value.erase(k);
		this->value2key[v].erase(k);
	}
	
	void clear(Value v) {
		BOOST_FOREACH(Key k, this->value2key[v]) {
			if (this->key2value[k] == v) {
				this->key2value.erase(k);
			}
		}
		this->value2key.erase(v);
	}
	
	void rehash(std::size_t size) {
		this->key2value.rehash(size);
		this->value2key.rehash(size);
	}

	boost::unordered_set<Key> getKeys() {
		boost::unordered_set<Key> keys;

		BOOST_FOREACH (typename GraphLookup::KeyType k, this->key2value) {
			keys.insert(k.first);
		}

		return keys;
	}
	
	boost::unordered_set<Value> getValues() {
		boost::unordered_set<Value> values;

		BOOST_FOREACH (typename GraphLookup::ValueType v, this->value2key) {
			values.insert(v.first);
		}

		return values;
	}
private:
	Key2Value key2value;
	Value2Key value2key;
};

#endif // GRAPH_LOOKUP_HH
