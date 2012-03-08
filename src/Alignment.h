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
	void testSymmetry(string prefix, bool extended, int windowSize, int windowStep);

private:
	int _dataType;
	int _format;
	unsigned int _cols;
	vector<Sequence> _alignment;
};

#endif /* ALIGNMENT_H_ */
