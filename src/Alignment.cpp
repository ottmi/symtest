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
#include "Matrix.h"
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
	_cols /= options->groupLength;

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

		if (verbose) cout << counts[0] << " DNA characters, " << counts[1] << " AA characters, " << counts[2] << " Alphanum characters." << endl;
		int dataTypeGuess = _ALPHANUM_DATA;
		if (counts[2] == counts[0])
			dataTypeGuess = _DNA_DATA;
		else if (counts[2] == counts[1]) dataTypeGuess = _AA_DATA;
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
	string stuartFileName = prefix + ".stuart.csv";
	string internalFileName = prefix + ".internal.csv";
	ofstream resultsFile, bowkerFile, delta_sFile, delta_msFile, stuartFile, internalFile;
	cout << endl << "Performing tests of pairwise symmetry, writing results to: ";
	if (extended)
	{
		cout << endl;
		cout << "  comprehensive spreadsheet: " << resultsFileName << endl;
		cout << "  Bowker matrix:             " << bowkerFileName << endl;
		cout << "  delta_s distance matrix:   " << delta_sFileName << endl;
		cout << "  delta_ms distance matrix:  " << delta_msFileName << endl;
		cout << "  Stuart matrix:             " << stuartFileName << endl;
		cout << "  Internal Symmetry:         " << internalFileName << endl;

		bowkerFile.open(bowkerFileName.c_str(), ifstream::trunc);
		if (!bowkerFile.is_open()) throw("Error, cannot open file " + bowkerFileName);

		delta_sFile.open(delta_sFileName.c_str(), ifstream::trunc);
		if (!delta_sFile.is_open()) throw("Error, cannot open file " + delta_sFileName);

		delta_msFile.open(delta_msFileName.c_str(), ifstream::trunc);
		if (!delta_msFile.is_open()) throw("Error, cannot open file " + delta_msFileName);

		stuartFile.open(stuartFileName.c_str(), ifstream::trunc);
		if (!stuartFile.is_open()) throw("Error, cannot open file " + stuartFileName);

		internalFile.open(internalFileName.c_str(), ifstream::trunc);
		if (!internalFile.is_open()) throw("Error, cannot open file " + internalFileName);
} else
	{
		cout << resultsFileName << endl;
	}

	resultsFile.open(resultsFileName.c_str(), ifstream::trunc);
	if (!resultsFile.is_open()) throw("Error, cannot open file " + resultsFileName);

	unsigned int len = _alignment.size();
	if (windowSize <= 0) windowSize = _cols;
	if (windowStep <= 0) windowStep = windowSize;

	if (windowSize < (int) _cols) cout << "  WindowSize=" << windowSize << " StepWidth=" << windowStep << " _cols=" << _cols << endl;

	int dim;
	if (_dataType == _DNA_DATA)
		dim = 4;
	else if (_dataType == _AA_DATA)
		dim = 20;
	else
		dim = 36;

	cout.precision(6);
	resultsFile << "Seq1\tSeq2\tChi-square\tdf\tp-value\tDelta_s\tDelta_ms\tStuart\tInternal\tSites\tStart\tEnd" << endl;
	for (unsigned int windowStart = 0; windowStart + windowSize <= _cols; windowStart += windowStep)
	{
		double bowker_mat[len][len];
		double ds_mat[len][len];
		double dms_mat[len][len];
		double stuart_mat[len][len];
		double internal_mat[len][len];
		vector<unsigned int> count(10, 0);
		unsigned int counter = 0;
		double minQ = 1.0;
		for (unsigned int k = 0; k < len; k++) // 1st sequence
		{
			bowker_mat[k][k] = 0;
			ds_mat[k][k] = 0;
			dms_mat[k][k] = 0;
			stuart_mat[k][k] = 0;
			Sequence s1 = _alignment[k];
			for (unsigned int l = k + 1; l < len; l++) // 2nd sequence
			{
				Sequence s2 = _alignment[l];
				unsigned int sum = 0;
				map<pair<unsigned int, unsigned int> , unsigned int> nmap;
				map<pair<unsigned int, unsigned int> , unsigned int>::iterator it;
				unsigned long id;
				for (unsigned int m = windowStart; m < windowStart + windowSize; m++)
				{
					unsigned int c1 = s1.getNumerical(m);
					unsigned int c2 = s2.getNumerical(m);
					if (s1.charIsUnambiguous(c1) && s2.charIsUnambiguous(c2))
					{
						pair<unsigned int, unsigned int> p(c1, c2);
						if ((it = nmap.find(p)) != nmap.end())
							it->second = it->second + 1;
						else
							nmap.insert(pair<pair<unsigned int, unsigned int> , unsigned int>(p, 1));
						sum++;
					}
				}

				Matrix n(dim);
				for (int i = 0; i < dim; i++)
					for (int j = 0; j < dim; j++)
					{
						pair<unsigned int, unsigned int> p(i, j);
						if ((it = nmap.find(p)) != nmap.end())
							n(i,j) = it->second;
						else
							n(i,j) = 0;
					}

				unsigned int df = 0;
				double bowker = .0;
				double delta_s = 0.0;
				for (int i = 0; i < dim; i++)
				{
					for (int j = i + 1; j < dim; j++)
					{
						if (n(i,j) + n(j,i) > 0)
						{
							df++;
							bowker += (double) ((n(i,j) - n(j,i)) * (n(i,j) - n(j,i))) / (n(i,j) + n(j,i));
						}

						double x = n(j,i) - n(i,j);
						delta_s += (x / sum) * (x / sum);
					}
				}

				double Q = 1.0;
				if (df > 0) Q = gammq((df / 2.0), (bowker / 2.0));
				bowker_mat[k][l] = Q;
				bowker_mat[l][k] = Q;

				if(Q < minQ) minQ = Q;
				if(Q < 0.05) count[0]++;
				if(Q < 0.01) count[1]++;
				if(Q < 0.005) count[2]++;
				if(Q < 0.001) count[3]++;
				if(Q < 0.0005) count[4]++;
				if(Q < 0.0001) count[5]++;
				if(Q < 0.00005) count[6]++;
				if(Q < 0.00001) count[7]++;
				if(Q < 0.000005) count[8]++;
				if(Q < 0.000001) count[9]++;

				delta_s = sqrt(delta_s);
				ds_mat[k][l] = delta_s;
				ds_mat[l][k] = delta_s;

				double delta_ms = 0.0;
				for (int i = 0; i < dim; i++)
				{
					double rowSum = n.getRowSum(i);
					double colSum = n.getColSum(i);
					delta_ms += ((rowSum - colSum) / sum) * ((rowSum - colSum) / sum);
				}
				delta_ms = sqrt(delta_ms) / sqrt(2.0);
				dms_mat[k][l] = delta_ms;
				dms_mat[l][k] = delta_ms;

				Matrix V(dim-1);
				for (unsigned int i = 0; i < dim - 1; i++)
					for (unsigned int j = 0; j < dim - 1; j++)
						if (i == j)
							V(i, j) = n.getRowSum(i) + n.getColSum(i) - 2 * n(i, i);
						else
							V(i, j) = -(n(i, j) + n(j, i));
				V.inverse();

				double stuart = 0;
				for (unsigned int i = 0; i < dim - 1; i++)
				{
					double d_i = n.getRowSum(i) - n.getColSum(i);
					for (unsigned int j = 0; j < dim - 1; j++)
					{
						double d_j = n.getRowSum(j) - n.getColSum(j);
						stuart+= V(i, j) * d_i * d_j;
					}
				}
				stuart_mat[k][l] = gammq(1.5, (stuart / 2.0));
				stuart_mat[l][k] = stuart_mat[k][l];

				double internal = bowker - stuart;
				internal_mat[k][l] = gammq(1.5, (internal / 2.0));
				internal_mat[l][k] = internal_mat[k][l];


				resultsFile << _alignment[k].getName() << "\t" << _alignment[l].getName() << "\t" << scientific << bowker << "\t" << df << "\t" << Q << "\t" << delta_s
						<< "\t" << delta_ms << "\t" << stuart_mat[k][l] << "\t" << internal_mat[k][l] << "\t" << sum << "\t" << windowStart << "\t" << windowStart + windowSize - 1 << endl;

				counter++;
			}
		}

		if (extended)
		{
			bowkerFile.flags(ios::left);
			delta_sFile.flags(ios::left);
			delta_msFile.flags(ios::left);
			stuartFile.flags(ios::left);
			internalFile.flags(ios::left);
			bowkerFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
			delta_sFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
			delta_msFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
			stuartFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
			internalFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;

			for (unsigned int l = 0; l < len; l++)
			{
				bowkerFile << "\t" << setw(12) << _alignment[l].getName();
				delta_sFile << "\t" << setw(12) << _alignment[l].getName();
				delta_msFile << "\t" << setw(12) << _alignment[l].getName();
				stuartFile << "\t" << setw(12) << _alignment[l].getName();
				internalFile << "\t" << setw(12) << _alignment[l].getName();
			}
			bowkerFile << endl;
			delta_sFile << endl;
			delta_msFile << endl;
			stuartFile << endl;
			internalFile << endl;

			for (unsigned int k = 0; k < len; k++)
			{
				bowkerFile.flags(ios::left);
				delta_sFile.flags(ios::left);
				delta_msFile.flags(ios::left);
				stuartFile.flags(ios::left);
				internalFile.flags(ios::left);
				bowkerFile << setw(12) << _alignment[k].getName();
				delta_sFile << setw(12) << _alignment[k].getName();
				delta_msFile << setw(12) << _alignment[k].getName();
				stuartFile << setw(12) << _alignment[k].getName();
				internalFile << setw(12) << _alignment[k].getName();
				bowkerFile.flags(ios::right);
				delta_sFile.flags(ios::right);
				delta_msFile.flags(ios::right);
				stuartFile.flags(ios::right);
				internalFile.flags(ios::right);
				for (unsigned int l = 0; l < len; l++)
				{
					if (k == l)
					{
						bowkerFile   << "\t      -      ";
						delta_sFile  << "\t      -      ";
						delta_msFile << "\t      -      ";
						stuartFile   << "\t      -      ";
						internalFile   << "\t      -      ";
					}
					else
					{
						bowkerFile << "\t" << scientific << bowker_mat[k][l];
						delta_sFile << "\t" << scientific << ds_mat[k][l];
						delta_msFile << "\t" << scientific << dms_mat[k][l];
						stuartFile << "\t" << scientific << bowker_mat[k][l];
						internalFile << "\t" << scientific << internal_mat[k][l];
					}
				}
				bowkerFile << endl;
				delta_sFile << endl;
				delta_msFile << endl;
				stuartFile << endl;
				internalFile << endl;
			}
		}

		cout << endl << "Highlights from the analysis (window "<< windowStart << "-" << windowStart + windowSize - 1 << "):" << endl;
		cout.precision(2);
		if(minQ < 0.05)     cout << "P-values < 0.05                 " << setw(8) << count[0] << " (" << fixed << (double) count[0]*100/counter << "%)" << endl;
		if(minQ < 0.01)     cout << "P-values < 0.01                 " << setw(8) << count[1] << " (" << fixed << (double) count[1]*100/counter << "%)" << endl;
		if(minQ < 0.005)    cout << "P-values < 0.005                " << setw(8) << count[2] << " (" << fixed << (double) count[2]*100/counter << "%)" << endl;
		if(minQ < 0.001)    cout << "P-values < 0.001                " << setw(8) << count[3] << " (" << fixed << (double) count[3]*100/counter << "%)" << endl;
		if(minQ < 0.0005)   cout << "P-values < 0.0005               " << setw(8) << count[4] << " (" << fixed << (double) count[4]*100/counter << "%)" << endl;
		if(minQ < 0.0001)   cout << "P-values < 0.0001               " << setw(8) << count[5] << " (" << fixed << (double) count[5]*100/counter << "%)" << endl;
		if(minQ < 0.00005)  cout << "P-values < 0.00005              " << setw(8) << count[6] << " (" << fixed << (double) count[6]*100/counter << "%)" << endl;
		if(minQ < 0.00001)  cout << "P-values < 0.00001              " << setw(8) << count[7] << " (" << fixed << (double) count[7]*100/counter << "%)" << endl;
		if(minQ < 0.000005) cout << "P-values < 0.000005             " << setw(8) << count[8] << " (" << fixed << (double) count[8]*100/counter << "%)" << endl;
		if(minQ < 0.000001) cout << "P-values < 0.000001             " << setw(8) << count[9] << " (" << fixed << (double) count[9]*100/counter << "%)" << endl;
		cout << "Number of tests                   " << setw(15) << counter << endl;
		cout.precision(8);
		cout << "Smallest P-value                  " << setw(15) << scientific << minQ << endl;
	}
	resultsFile.close();
	if (extended)
	{
		bowkerFile.close();
		delta_sFile.close();
		delta_msFile.close();
	}
}
