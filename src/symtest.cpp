#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstdio>
#include <unistd.h>
#include "globals.h"
#include "Alignment.h"
#include "helper.h"
#ifdef _MPFR
#include "mpfr.h"
#endif

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
	options->windowSize = -1;
	options->windowStep = -1;
	options->columnFrom = -1;
	options->columnTo = -1;

	int minGroup = 0;
	int maxGroup = 0;
	string listOfSeq;

	while ((c = getopt(argc, argv, "t:s:p:g:c:xw:n:v::h")) != -1)
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
				cerr << "Parameter -s did not contain a list of sequences but also does not point to a file." << endl;
				cerr << "Please specify either a comma-separated list of sequences or a text file with the sequence names." << endl;
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
	cout << "Usage:" << endl;
	cout << "  symtest [options] <alignment>" << endl;
	cout << "  symtest -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -t<a|d|n>      Data type a=AA, d=DNA, n=Alphanumeric [default: auto-detect]" << endl;
	cout << "  -s<LIST|FILE>  Only consider sequences in comma-separated LIST or FILE" << endl;
	cout << "  -p<STRING>     Prefix for output files [default: name of alignment w/o .ext]" << endl;
	cout << endl;
	cout << "  -c<from-to>    Only consider alignment columns from-to, enumeration starts with 1" << endl;
	cout << "  -g<LIST>       Grouping of sites, e.g. 1,2,-3 for duplets, 1,2,3 for codons" << endl;
	cout << endl;
	cout << "  -x             Write extended output files with test results" << endl;
	cout << "  -w<NUM>,[NUM]  Window size and step width for symmetry test" << endl;
	cout << endl;
#ifdef _OPENMP
	cout << "  -n<NUM>        Number of threads [default: " << omp_get_max_threads() << "]" << endl;
#endif
	cout << "  -v[NUM]        Be increasingly verbose" << endl;
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
		alignment.writeSummary(options.prefix, options.windowSize, options.windowStep);
		if (options.writeExtendedTestResults)
			alignment.writeExtendedResults(options.prefix, options.windowSize, options.windowStep);
	} catch (string& s)
	{
		cerr << s << endl;
		return (255);
	}

	return 0;
}
