#include "Sequence.h"
#include "globals.h"
#include "helper.h"

Sequence::Sequence(string name, string seq)
{
    _name = name;
    _sequence = seq;
    _grouping.push_back(1);
}

Sequence::~Sequence()
{
    // TODO Auto-generated destructor stub
}

void Sequence::translateToNum(int dataType)
{
    _dataType = dataType;
    if (dataType == _DNA_DATA)
	_unambiguousThreshold = _DNA_UNAMB_THRES;
    else if (dataType == _AA_DATA)
	_unambiguousThreshold = _AA_UNAMB_THRES;
    else
	_unambiguousThreshold = _ALPHANUM_UNAMB_THRES;

    for (unsigned int i = 0; i < _sequence.size(); i++)
	_numericSeq.push_back(mapCharToNum(_sequence.substr(i, 1), dataType));
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
};

bool Sequence::charIsUnambiguous(unsigned int n)
{
    for (unsigned int i = 0; i < _grouping.size(); i++)
    {
	if ((n & 255) > _unambiguousThreshold)
	    return false;
	n = n >> 8;
    }
    return true;
}

