#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cfloat>
#include <sstream>
#include <map>
#include <utility>
#include "AlignmentReader.h"
#include "Alignment.h"
#include "helper.h"

Alignment::Alignment()
{
    _cols = 0;
}

Alignment::Alignment(Options *options)
{
    AlignmentReader alignmentReader(options->inputAlignment);
    _alignment = alignmentReader.getSequences(options->columnFrom, options->columnTo);

    _cols = alignmentReader.getCols();
    if (options->columnTo > (int) _cols)
    {
    	stringstream ss;
    	ss << options->columnTo << " has been specified as upper boundary for the alignment columns to be read, but the alignment contains only " << _cols << " columns.";
    	throw(ss.str());
    }
    if (options->columnFrom > 0 && options->columnTo > 0)
    	_cols = options->columnTo - options->columnFrom + 1;
    else if (options->columnFrom > 0)
    	_cols-= (options->columnFrom - 1);
    cout << _cols << endl;
    _cols/= options->groupLength;

    string dataTypeDesc[] = { "DNA", "AA", "alphanumeric" };
    if (options->dataType < 0)
    {
	map<char, unsigned long> baseOccurences;
	for (unsigned int i = 0; i < _alignment.size(); i++)
	{
	    string s = _alignment[i].getSequence();
	    for (unsigned int j = 0; j < s.length(); j++)
		baseOccurences[s[j]]++;
	}

	string maps[] = { _DNA_MAP, _AA_MAP, _ALPHANUM_MAP };
	unsigned long counts[3];
	for (unsigned int i = 0; i < 3; i++)
	{
	    counts[i] = 0;
	    string map = maps[i];
	    for (unsigned j = 0; j < map.length(); j++)
		counts[i] += baseOccurences[map[j]];
	}

	if (verbose)
	    cout << counts[0] << " DNA characters, " << counts[1] << " AA characters, " << counts[2] << " Alphanum characters." << endl;
	int dataTypeGuess = _ALPHANUM_DATA;
	if (counts[2] == counts[0])
	    dataTypeGuess = _DNA_DATA;
	else if (counts[2] == counts[1])
	    dataTypeGuess = _AA_DATA;
	_dataType = dataTypeGuess;

	cout << "It contains " << getNumOfRows() << " sequences " << "which appear to be " << dataTypeDesc[_dataType] << "." << endl;
    } else
    {
	_dataType = options->dataType;
	cout << "It contains " << getNumOfRows() << " sequences which have been defined to be " << dataTypeDesc[_dataType] << "." << endl;
    }

    for (unsigned int i = 0; i < _alignment.size(); i++)
    {
	_alignment[i].translateToNum(_dataType, options);
	if (verbose >= 2)
	{
	    cout << i << ": ";
	    for (unsigned int j = 0; j < _alignment[i].getLength() / options->groupLength; j++)
	      cout << mapNumToChar(_alignment[i].getNumerical(j), _dataType, options->grouping.size()) << " ";
	    cout << endl;
	}
    }
}

Alignment::~Alignment()
{
}

void Alignment::testSymmetry(string prefix, bool extended, int windowSize, int windowStep)
{
    string resultsFileName = prefix + ".symmetry.csv";
    string bowkerFileName = prefix + ".bowker.csv";
    string delta_sFileName = prefix + ".delta_s.csv";
    string delta_msFileName = prefix + ".delta_ms.csv";
    ofstream resultsFile, bowkerFile, delta_sFile, delta_msFile;
    cout << endl << "Performing tests of pairwise symmetry, writing results to: ";
    if (extended)
    {
	cout << endl;
	cout << "  comprehensive spreadsheet: " << resultsFileName << endl;
	cout << "  Bowker matrix:             " << bowkerFileName << endl;
	cout << "  delta_s distance matrix:   " << delta_sFileName << endl;
	cout << "  delta_ms distance matrix:  " << delta_msFileName << endl;

	bowkerFile.open(bowkerFileName.c_str(), ifstream::trunc);
	if (!bowkerFile.is_open())
	    throw("Error, cannot open file " + bowkerFileName);

	delta_sFile.open(delta_sFileName.c_str(), ifstream::trunc);
	if (!delta_sFile.is_open())
	    throw("Error, cannot open file " + delta_sFileName);

	delta_msFile.open(delta_msFileName.c_str(), ifstream::trunc);
	if (!delta_msFile.is_open())
	    throw("Error, cannot open file " + delta_msFileName);
    } else
    {
	cout << resultsFileName << endl;
    }

    resultsFile.open(resultsFileName.c_str(), ifstream::trunc);
    if (!resultsFile.is_open())
	throw("Error, cannot open file " + resultsFileName);

    unsigned int n = _alignment.size();
    if (windowSize <= 0)
	windowSize = _cols;
    if (windowStep <= 0)
	windowStep = windowSize;

    if (windowSize < (int) _cols)
	cout << "  WindowSize=" << windowSize << " StepWidth=" << windowStep << " _cols=" << _cols << endl;

    int dim;
    if (_dataType == _DNA_DATA)
	dim = 5;
    else if (_dataType == _AA_DATA)
	dim = 20;
    else
	dim = 36;

    resultsFile << "Seq1\tSeq2\tChi-square\tdf\tp-value\tDelta_s\tDelta_ms\tSites\tStart\tEnd" << endl;
    for (unsigned int windowStart = 0; windowStart + windowSize <= _cols; windowStart += windowStep)
    {
	vector<double> bowker_mat(n * n, 0);
	vector<double> ds_mat(n * n, 0);
	vector<double> dms_mat(n * n, 0);
	for (unsigned int k = 0; k < n; k++) // 1st sequence
	{
	    Sequence s1 = _alignment[k];
	    for (unsigned int l = k + 1; l < n; l++) // 2nd sequence
	    {
		Sequence s2 = _alignment[l];
		unsigned int sum = 0;
		map<pair<unsigned int, unsigned int>, unsigned int> dm;
		map<pair<unsigned int, unsigned int>, unsigned int>::iterator it;
		unsigned long id;
		for (unsigned int m = windowStart; m < windowStart + windowSize; m++)
		{
		    unsigned int c1 = s1.getNumerical(m);
		    unsigned int c2 = s2.getNumerical(m);
		    if (s1.charIsUnambiguous(c1) && s2.charIsUnambiguous(c2))
		    {
		        pair<unsigned int, unsigned int> p(c1, c2);
			if ((it = dm.find(p)) != dm.end())
			    it->second = it->second + 1;
			else
			    dm.insert(pair<pair<unsigned int, unsigned int>, unsigned int>(p, 1));
			sum++;
		    }
		}

		unsigned int df = 0;
		double bowker = .0;
		double delta_s = 0.0;
		for (int i = 0; i < dim; i++)
		{
		    for (int j = i + 1; j < dim; j++)
		    {
			double dm_ij = 0;
			pair<unsigned int, unsigned int> p1(i, j);
			if ((it = dm.find(p1)) != dm.end())
			    dm_ij = it->second;

			double dm_ji = 0;
			pair<unsigned int, unsigned int> p2(j, i);
			if ((it = dm.find(p2)) != dm.end())
			    dm_ji = it->second;

			if (dm_ij + dm_ji > 0)
			{
			    df++;
			    bowker += ((dm_ij - dm_ji) * (dm_ij - dm_ji)) / (dm_ij + dm_ji);
			}

			double x = dm_ji - dm_ij;
			delta_s += (x / sum) * (x / sum);
		    }
		}

		double Q = 1.0;
		if (df > 0)
		    Q = gammq((df / 2.0), (bowker / 2.0));
		bowker_mat[k * n + l] = Q;
		bowker_mat[l * n + k] = Q;

		delta_s = sqrt(delta_s);
		ds_mat[k * n + l] = delta_s;
		ds_mat[l * n + k] = delta_s;

		double delta_ms = 0.0;
		for (int i = 0; i < dim; i++)
		{
		    double row = 0.0;
		    double col = 0.0;
		    for (int j = 0; j < dim; j++)
		    {
			double dm_ij = 0;
			pair<unsigned int, unsigned int> p1(i, j);
			if ((it = dm.find(p1)) != dm.end())
			    dm_ij = it->second;

			double dm_ji = 0;
			pair<unsigned int, unsigned int> p2(j, i);
			if ((it = dm.find(p2)) != dm.end())
			    dm_ji = it->second;

			row += dm_ij;
			col += dm_ji;
		    }
		    delta_ms += ((row - col) / sum) * ((row - col) / sum);
		}
		delta_ms = sqrt(delta_ms) / sqrt(2.0);
		dms_mat[k * n + l] = delta_ms;
		dms_mat[l * n + k] = delta_ms;

		resultsFile << _alignment[k].getName() << "\t" << _alignment[l].getName() << "\t" << scientific << bowker << "\t" << df << "\t" << Q << "\t" << delta_s << "\t" << delta_ms << "\t"
			<< sum << "\t" << windowStart << "\t" << windowStart + windowSize - 1 << endl;
	    }

	}

	if (extended)
	{
	    bowkerFile.flags(ios::left);
	    delta_sFile.flags(ios::left);
	    delta_msFile.flags(ios::left);
	    bowkerFile << windowStart << "-" << windowStart + windowSize - 1;
	    delta_sFile << windowStart << "-" << windowStart + windowSize - 1;
	    delta_msFile << windowStart << "-" << windowStart + windowSize - 1;

	    for (unsigned int l = 0; l < n; l++)
	    {
		bowkerFile << "\t" << setw(12) << _alignment[l].getName();
		delta_sFile << "\t" << setw(12) << _alignment[l].getName();
		delta_msFile << "\t" << setw(12) << _alignment[l].getName();
	    }
	    bowkerFile << endl;
	    delta_sFile << endl;
	    delta_msFile << endl;

	    for (unsigned int k = 0; k < n; k++)
	    {
		bowkerFile << setw(12) << _alignment[k].getName();
		delta_sFile << setw(12) << _alignment[k].getName();
		delta_msFile << setw(12) << _alignment[k].getName();
		for (unsigned int l = 0; l < n; l++)
		{
		    bowkerFile << "\t" << scientific << bowker_mat[k * n + l];
		    delta_sFile << "\t" << scientific << ds_mat[k * n + l];
		    delta_msFile << "\t" << scientific << dms_mat[k * n + l];
		}
		bowkerFile << endl;
		delta_sFile << endl;
		delta_msFile << endl;
	    }
	}
    }
    resultsFile.close();
    if (extended)
    {
	bowkerFile.close();
	delta_sFile.close();
	delta_msFile.close();
    }
}
