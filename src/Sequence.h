#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include <vector>
#include "globals.h"

using namespace std;

class Sequence {
public:
	Sequence(string name, string seq);
	virtual ~Sequence();
	void translateToNum(int dataType, Options *options);
	string& getName() { return _name; };
	string& getSequence() { return _sequence; };
	string getColumns(vector<int> cols);
	unsigned int getNumerical(int pos);
	size_t getLength() { return _sequence.length(); };
	bool charIsUnambiguous(unsigned int n);
private:
	int _dataType;
	string _name;
	string _sequence;
	vector<unsigned int> _numericSeq;
	vector<int> _grouping;
	int _unambiguousThreshold;
};

#endif /* SEQUENCE_H_ */
