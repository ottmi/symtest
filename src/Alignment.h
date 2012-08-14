#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include <vector>
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

	void writeSummary(string prefix, int windowSize, int windowStep);
	void writeExtendedResults(string prefix, int windowSize, int windowStep);

	unsigned int getNumOfRows() { return _alignment.size(); };
	unsigned int getNumOfCols() { return _cols; };

private:
	int _dim;
	int _dataType;
	int _format;
	unsigned int _cols;
	vector<Sequence> _alignment;

	vector< vector<unsigned int> > _df;
	vector< vector<double> > _bowker;
	vector< vector<double> > _ds;
	vector< vector<double> > _dms;
	vector< vector<double> > _stuart;
	vector< vector<double> > _aitchison;
};

#endif /* ALIGNMENT_H_ */
