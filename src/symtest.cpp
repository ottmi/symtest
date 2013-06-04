#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstdio>
#include <cctype>
#include <cstring>
#include <unistd.h>
#include "globals.h"
#include "Alignment.h"
#include "helper.h"

using namespace std;

int verbose = 0;

int parseArguments(int argc, char** argv, Options *options)
{
	char c;

	if (argc < 2) return 0;

	if (argc == 2 && argv[1][0] == '-')
	{
		options->help = true;
		return 0;
	}

	options->inputAlignment = string(argv[--argc]);

	options->dataType = -1;
	options->help = 0;
	options->grouping.push_back(1);
	options->writeExtendedTestResults = false;
	options->writeBowkerFile = false;
	options->writeStuartFile = false;
	options->writeAbabnehFile = false;
	options->writeAitchisonFile = false;
	options->writeDelta_sFile = false;
	options->writeDelta_msFile = false;
	options->windowSize = -1;
	options->windowStep = -1;
	options->columnFrom = -1;
	options->columnTo = -1;

	int minGroup = 0;
	int maxGroup = 0;
	string listOfSeq;

	while ((c = getopt(argc, argv, "t:s:p:g:c:x::w:n:v::h")) != -1)
	{
		switch (c)
		{
			case 't':
			{
				char type = optarg[0];
				switch (type)
				{
					case 'a':
					case 'A':
						options->dataType = _AA_DATA;
						break;
					case 'd':
					case 'D':
						options->dataType = _DNA_DATA;
						break;
					case 'n':
					case 'N':
						options->dataType = _ALPHANUM_DATA;
						break;
					default:
						cerr << "Unknown data type: " << optarg << endl;
						return 2;
						break;
				}
				break;
			}
			case 's':
				listOfSeq = optarg;
				break;
			case 'p':
				options->prefix = optarg;
				break;
			case 'g':
			{
				options->grouping.clear();
				minGroup = 255;
				int i;
				stringstream ss(optarg);
				while (ss >> i)
				{
					if (i == 0)
					{
						cerr << "The site numbers used in groupings start with 1." << endl;
						return 3;
					}
					if (i > 0)
						options->grouping.push_back(i);
					else
						i = -i;
					if (i > maxGroup) maxGroup = i;
					if (i < minGroup) minGroup = i;

					if (ss.peek() == ',') ss.ignore();
				}
				break;
			}
			case 'c':
			{
				stringstream ss(optarg);
				ss >> options->columnFrom;
				if (options->columnFrom == 0)
				{
					cerr << "Alignment column enumeration starts with 1" << endl;
					return 4;
				}
				if (ss.peek() == '-')
				{
					ss.ignore();
					ss >> options->columnTo;
				}
				break;
			}
			case 'x':
				options->writeExtendedTestResults = true;
				if (optarg && strlen(optarg) != 0) {
					string s(optarg);
					transform(s.begin(), s.end(), s.begin(), ::tolower);
					cout << "Parsing " << s << endl;
					if (s.find("bowker") != string::npos)
						options->writeBowkerFile = true;
					if (s.find("stuart") != string::npos)
						options->writeStuartFile = true;
					if (s.find("ababneh") != string::npos)
						options->writeAbabnehFile = true;
					if (s.find("aitchison") != string::npos)
						options->writeAitchisonFile= true;
					if (s.find("delta_s") != string::npos)
						options->writeDelta_sFile = true;
					if (s.find("delta_ms") != string::npos)
						options->writeDelta_msFile = true;
				} else {
					options->writeBowkerFile = true;
					options->writeStuartFile = true;
					options->writeAbabnehFile = true;
					options->writeAitchisonFile= true;
					options->writeDelta_sFile = true;
					options->writeDelta_msFile = true;
				}
				break;
			case 'w':
			{
				int i;
				stringstream ss(optarg);
				while (ss >> i)
				{
					if (options->windowSize < 0) options->windowSize = i;
					options->windowStep = i;
					if (ss.peek() == ',') ss.ignore();
				}
				break;
			}
#ifdef _OPENMP
				case 'n':
				omp_set_num_threads(atoi(optarg));
				break;
#endif
			case 'v':
				if (optarg)
					verbose = atoi(optarg);
				else
					verbose = 1;
				break;
			case 'h':
				options->help = 1;
				break;
			default:
				if (c != '?') cerr << "Unknown parameter: " << c << endl;
				return 1;
		}
	}

	options->groupLength = maxGroup - minGroup + 1;
	if (verbose)
	{
		cout << "Grouping: length=" << options->groupLength << " grouping=";
		for (unsigned int i = 0; i < options->grouping.size(); i++)
		{
			if (i > 0) cout << ",";
			cout << options->grouping[i];
		}
		cout << endl;
	}

	if (options->columnFrom >= 0) cout << "Columns: " << options->columnFrom << " - " << options->columnTo << endl;

	if (options->prefix.length() == 0)
	{
		int m = options->inputAlignment.find_last_of('/') + 1;
		int n = options->inputAlignment.find_last_of('.');
		if (n > -1) n = n - m;
		options->prefix = options->inputAlignment.substr(m, n);
	}

	if (listOfSeq.length() > 0) {
		if (listOfSeq.find_first_of(',') == string::npos) {
			ifstream _fileReader;
			_fileReader.open(listOfSeq.c_str());
			if (!_fileReader.is_open()) {
				cerr << "Parameter -s does neither contain a list of sequences nor point to a file." << endl;
				cerr << "Please specify either a comma-separated list of sequence names or a text file" << endl;
				cerr << "containing one sequence name per line." << endl;
				return 5;
			}
			string seq;
			while (!_fileReader.eof()) {
				safeGetline(_fileReader, seq);
				if (seq.length()) {
					options->listOfSequences.insert(seq);
				}
			}
			_fileReader.close();
		} else {
			size_t from, to;
			from = 0;
			while ((to = listOfSeq.find_first_of(',', from)) != string::npos) {
				string seq = listOfSeq.substr(from, to-from);
				if (seq.length()) {
					options->listOfSequences.insert(seq);
				}
				from = to+1;
			}
			string seq = listOfSeq.substr(from);
			if (seq.length()) {
				options->listOfSequences.insert(seq);
			}
		}

		cout << "Sequences to be considered:" << endl;
		set<string>::iterator it;
		for (it = options->listOfSequences.begin(); it != options->listOfSequences.end(); it++) {
			cout << *it << endl;
		}
		cout << endl;
	}

	return 0;
}

void printSyntax()
{
/*
                       1         2         3         4         5         6         7         8
             012345678901234567890123456789012345678901234567890123456789012345678901234567890
*/
	cout << "Usage:" << endl;
	cout << "  symtest [options] <alignment>" << endl;
	cout << "  symtest -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
    cout << "  -t<a|d|n>      Data types: a = AA; d = DNA; n = Alphanumeric" << endl;
    cout << "                 [default: auto-detect]" << endl;
	cout << "  -s<FILE|LIST>  Only consider sequences in FILE or comma-separated LIST" << endl;
    cout << "  -p<STRING>     Prefix for output files [default: name of alignment w/o .ext]" << endl;
    cout << endl;
    cout << "  -c<n1-n2>      Limited window analysis: from column n1 to column n2," << endl;
    cout << "                 enumeration starts with 1" << endl;
    cout << "  -w<n3,n4>      Sliding window analysis: window size = n3; step size = n4" << endl;
    cout << "  -g<LIST>       Grouping of codon sites: comma-separated, pos/neg to use/skip" << endl;
    cout << "                 e.g., 1,2,-3 for duplets comprising 1st and 2nd position" << endl;
    cout << "                       1,2,3  for triplets comprising 1st, 2nd, and 3rd position" << endl;
    cout << "  -x[LIST]       Write extended output files [default: all available]" << endl;
    cout << "                 Optional: restrict to files in LIST (comma-separated)" << endl;
    cout << "                 Available: bowker, stuart, ababneh, aitchison, delta_s, delta_ms" << endl;
    cout << endl;
    cout << "  -v<n5>         Be increasingly verbose [n5 = 0|1|2]" << endl;
#ifdef _OPENMP
    cout << "  -n<n6>         Number of threads [default: " << omp_get_max_threads() << "]" << endl;
#endif
	cout << "  -h             This help page" << endl;
	cout << endl;
}

int main(int argc, char** argv)
{
	Options options;

	cout << PROGNAME << " " << VERSION << "|";
#ifdef _OPENMP
	cout << "OpenMP|";
#endif

	cout << PROGDATE << endl << endl;

	int ret = parseArguments(argc, argv, &options);
	if (ret) return ret;

	if (!options.inputAlignment.length() || options.help)
	{
		printSyntax();
		return 254;
	}

#ifdef _OPENMP
	cout << "Parallel execution with " << omp_get_max_threads() << " threads." << endl << endl;
#endif

	try
	{
		Alignment alignment = Alignment(&options);
		alignment.testSymmetry(options.prefix, options.windowSize, options.windowStep);
		alignment.writeResults(&options);
	} catch (string& s)
	{
		cerr << s << endl;
		return (255);
	}

	return 0;
}
