#include "Sequence.h"
#include "helper.h"

Sequence::Sequence(string name, string seq)
{
	_name = name;
	_sequence = seq;
}

Sequence::~Sequence()
{
	// TODO Auto-generated destructor stub
}

void Sequence::translateToNum(int dataType, Options *options)
{
	_dataType = dataType;
	if (dataType == _DNA_DATA)
		_unambiguousThreshold = _DNA_UNAMB_THRES;
	else if (dataType == _AA_DATA)
		_unambiguousThreshold = _AA_UNAMB_THRES;
	else
		_unambiguousThreshold = _ALPHANUM_UNAMB_THRES;

	unsigned int numOfSites = _sequence.length() / options->groupLength;
	for (unsigned int i = 0; i < numOfSites; i++)
	{
		vector<int> cols;
		for (unsigned int j = 0; j < options->grouping.size(); j++)
			cols.push_back(options->groupLength * i + options->grouping[j] - 1);
		_numericSeq.push_back(mapCharToNum(getColumns(cols), dataType));
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

unsigned int Sequence::getNumerical(int pos)
{
	if (pos < _numericSeq.size())
		return _numericSeq[pos];
	else
		return 0;
}
;

bool Sequence::charIsUnambiguous(unsigned int n)
{
	for (unsigned int i = 0; i < _groupSize; i++)
	{
		if ((n & 255) > _unambiguousThreshold) return false;
		n = n >> 8;
	}
	return true;
}

