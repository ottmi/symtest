#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <vector>
#include <set>
using namespace std;

#define PROGNAME "symtest"
#define VERSION "2.0.41"
#define PROGDATE "2015-12-17"

#define _DNA_DATA	0
#define	_AA_DATA	1
#define	_ALPHANUM_DATA	2

#define _DNA_MAP      "ACGTURYKMSWBDHVN?-"
#define _AA_MAP       "ACDEFGHIKLMNPQRSTVWYBJXZ?-"
#define _ALPHANUM_MAP "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?-"

#define _DNA_UNAMB_THRES 	4
#define _AA_UNAMB_THRES 	19
#define _ALPHANUM_UNAMB_THRES 	35

#define	_FASTA_FORMAT	0
#define	_PHYLIP_FORMAT	1

#define CSV_SEPARATOR ","

extern int verbose;

typedef struct opt_struct
{
    string inputAlignment;
    set<string> listOfSequences;
    int dataType;
    string prefix;
    vector<int> grouping;
    int groupLength;
    int columnFrom;
    int columnTo;
    bool writeExtendedTestResults;
    bool writeBowkerFile;
    bool writeStuartFile;
    bool writeAbabnehFile;
    bool writeAmsFile;
    bool writeAfsFile;
    bool writeEmsFile;
    bool writeEfsFile;
    int windowSize;
    int windowStep;
    int help;
} Options;

#endif /* GLOBALS_H_ */
