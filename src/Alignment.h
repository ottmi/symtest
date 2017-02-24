#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include <vector>
#include <map>
#include <string>
#include "globals.h"
#include "Sequence.h"

//using namespace std;

class Site;
class Alignment
{
public:
	Alignment();
	Alignment(Options *options);
	virtual ~Alignment();

	void addSequence(Sequence sequence);

	void testSymmetry(string prefix, int windowSize, int windowStep);

	void printStatistics(vector<double> &pValues, char id);

	void writeResults(Options* options);
	void writeExtendedResult(string title, string baseName, string ext, unsigned int windowSize, unsigned int windowStep, vector< vector<double> >& matrix);
	void writeExtendedDistances(string title, string baseName, string ext, unsigned int windowSize, unsigned int windowStep, vector< vector<double> >& matrix);

	unsigned int getNumOfRows() { return _alignment.size(); };
	unsigned int getNumOfCols() { return _cols; };

private:
	int _dim;
	int _dataType;
	int _groupSize;
	int _format;
	unsigned int _cols;
	vector<Sequence> _alignment;
	map<unsigned int, unsigned int> _mapEnumeration;

	vector< vector<unsigned int> > _df;
	vector< vector<double> > _bowker;
	vector< vector<double> > _pBowker;
	vector< vector<double> > _stuart;
	vector< vector<double> > _pStuart;
	vector< vector<double> > _pAbabneh;
	vector< vector<double> > _aitchisonMarg;
	vector< vector<double> > _aitchisonFull;
	vector< vector<double> > _ds;
	vector< vector<double> > _dms;
};

#endif /* ALIGNMENT_H_ */
