#include <iostream>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <sstream>
#include <cstdio>
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

	if (argc < 2)
		return 0;

	if (argc == 2 && argv[1][0] == '-')
	{
		options->help = true;
		return 0;
	}

	options->inputAlignment = string(argv[--argc]);

	options->dataType = -1;
	options->help = 0;
	options->grouping.push_back(0);
	options->writeExtendedTestResults = false;
	options->windowSize = -1;
	options->windowStep = -1;

	int minGroup = 0;
	int maxGroup = 0;

	while ((c = getopt(argc, argv, "t:p:g:xw:n:v::h")) != -1)
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
					if (i >= 0)
						options->grouping.push_back(i);
					else
						i = -i;
					if (i > maxGroup)
						maxGroup = i;
					if (i < minGroup)
						minGroup = i;

					if (ss.peek() == ',')
						ss.ignore();
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
					if (options->windowSize < 0)
						options->windowSize = i;
					options->windowStep = i;
					if (ss.peek() == ',')
						ss.ignore();
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
				if (c != '?')
					cerr << "Unknown parameter: " << c << endl;
				return 1;
		}
	}

	options->groupOffset = minGroup;
	options->groupLength = maxGroup - minGroup + 1;

	if (options->prefix.length() == 0)
	{
		int m = options->inputAlignment.find_last_of('/') + 1;
		int n = options->inputAlignment.find_last_of('.');
		if (n > -1)
			n = n - m;
		options->prefix = options->inputAlignment.substr(m, n);
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
	cout << "  -p<STRING>     Prefix for output files [default: name of alignment w/o .ext]" << endl;
	cout << "  -g<LIST>       Grouping of sites, e.g. 0,1,-2 for duplets, 0,1,2 for codons" << endl;
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
	if (ret)
		return ret;

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
		alignment.testSymmetry(options.prefix, options.writeExtendedTestResults, options.windowSize, options.windowStep);
	} catch (string& s)
	{
		cerr << s << endl;
		return(255);
	}

	return 0;
}
