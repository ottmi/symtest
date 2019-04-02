#include "Sequence.h"
#include "helper.h"

Sequence::Sequence(string name, string seq)
{
	int n = name.find_first_of(" \t");
	if (n == -1)
	    _name = name;
	else
	    _name = name.substr(0, n);
	_longName = name;
	_sequence = seq;
}

Sequence::~Sequence()
{
	// TODO Auto-generated destructor stub
}

void Sequence::translateToNum(charMap_t& map, Options *options)
{
	_unambiguousThreshold = options->unambigousChars.size()-1;
	unsigned int numOfSites = _sequence.length() / options->groupLength;
	for (unsigned int i = 0; i < numOfSites; i++) {
		vector<int> cols;
		for (unsigned int j = 0; j < options->grouping.size(); j++)
			cols.push_back(options->groupLength * i + options->grouping[j] - 1);
		_numericSeq.push_back(mapCharToNum(getColumns(cols), map));
	}
	_groupSize = options->grouping.size();

}

string Sequence::getColumns(vector<int> cols)
{
	string s;

	for (unsigned int i = 0; i < cols.size(); i++)
		s += _sequence[cols[i]];

	return s;
}

unsigned int Sequence::getNumerical(unsigned int pos)
{
	if (pos < _numericSeq.size())
		return _numericSeq[pos];
	else
		return 0;
}

bool Sequence::charIsUnambiguous(unsigned int n)
{
	for (int i = 0; i < _groupSize; i++)
	{
		if ((int) (n & 255) > _unambiguousThreshold) return false;
		n = n >> 8;
	}
	return true;
}

