#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>
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

Alignment::Alignment() {
	_cols = 0;
}

Alignment::Alignment(Options *options) {
	AlignmentReader alignmentReader(options->inputAlignment);
	_alignment = alignmentReader.getSequences(options->columnFrom, options->columnTo, options->listOfSequences);
	_cols = alignmentReader.getCols();
	_cols /= options->groupLength;

	string dataTypeDesc[] = { "DNA", "AA", "alphanumeric" };
	if (options->dataType < 0) {
		map<char, unsigned long> baseOccurences;
		for (unsigned int i = 0; i < _alignment.size(); i++) {
			string s = _alignment[i].getSequence();
			for (unsigned int j = 0; j < s.length(); j++)
				baseOccurences[s[j]]++;
		}

		string maps[] = { _DNA_MAP, _AA_MAP, _ALPHANUM_MAP };
		unsigned long counts[3];
		for (unsigned int i = 0; i < 3; i++) {
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

		cout << "Read " << getNumOfRows() << " sequences " << "which appear to be " << dataTypeDesc[_dataType] << "." << endl;
	} else {
		_dataType = options->dataType;
		cout << "Read " << getNumOfRows() << " sequences which have been defined to be " << dataTypeDesc[_dataType] << "." << endl;
	}

	if (_dataType == _DNA_DATA)
		_dim = 4;
	else if (_dataType == _AA_DATA)
		_dim = 20;
	else
		_dim = 36;

	for (unsigned int i = 0; i < _alignment.size(); i++) {
		_alignment[i].translateToNum(_dataType, options);
		if (verbose >= 2) {
			cout << i << ": ";
			for (unsigned int j = 0; j < _alignment[i].getLength() / options->groupLength; j++)
				cout << mapNumToChar(_alignment[i].getNumerical(j), _dataType, options->grouping.size()) << " ";
			cout << endl;
		}
	}
}

Alignment::~Alignment() {
}

void Alignment::testSymmetry(string prefix, int windowSize, int windowStep) {
	cout << endl << "Performing tests of pairwise symmetry" << endl;

	unsigned int len = _alignment.size();
	if (windowSize <= 0)
		windowSize = _cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	if (windowSize < (int) _cols)
		cout << "  WindowSize=" << windowSize << " StepWidth=" << windowStep << " _cols=" << _cols << endl;

	cout.precision(6);
	for (unsigned int windowStart = 0; windowStart + windowSize <= _cols; windowStart += windowStep) {
		vector<vector<double> > baseFrequencies;
		vector<unsigned int> dfList;
		vector<double> bowkerList;
		vector<double> pBowkerList;
		vector<double> dsList;
		vector<double> dmsList;
		vector<double> stuartList;
		vector<double> pStuartList;
		vector<double> pAbabnehList;
		vector<double> aitchisonList;

		for (unsigned int k = 0; k < len; k++) // 1st sequence
				{
			Sequence s1 = _alignment[k];
			vector<unsigned long> baseOccurences(_dim, 0);
			for (unsigned int m = windowStart; m < windowStart + windowSize; m++) {
				unsigned int c = s1.getNumerical(m);
				if (s1.charIsUnambiguous(c))
					baseOccurences[c]++;
			}

			if (verbose)
				cout << "baseFreq[" << k << "]:";
			vector<double> baseFreq(_dim);
			for (unsigned int c = 0; c < _dim; c++) {
				baseFreq[c] = ((double) baseOccurences[c] + (1.0 / _dim)) / (windowSize + 1);
				if (verbose)
					cout << baseFreq[c] << " ";
			}
			baseFrequencies.push_back(baseFreq);
			if (verbose)
				cout << endl;

			for (unsigned int l = k + 1; l < len; l++) // 2nd sequence
					{
				Sequence s2 = _alignment[l];
				unsigned int sum = 0;
				map<pair<unsigned int, unsigned int>, unsigned int> nmap;
				map<pair<unsigned int, unsigned int>, unsigned int>::iterator it;
				unsigned long id;
				for (unsigned int m = windowStart; m < windowStart + windowSize; m++) {
					unsigned int c1 = s1.getNumerical(m);
					unsigned int c2 = s2.getNumerical(m);
					if (s1.charIsUnambiguous(c1) && s2.charIsUnambiguous(c2)) {
						pair<unsigned int, unsigned int> p(c1, c2);
						if ((it = nmap.find(p)) != nmap.end())
							it->second = it->second + 1;
						else
							nmap.insert(pair<pair<unsigned int, unsigned int>, unsigned int>(p, 1));
						sum++;
					}
				}

				Matrix n(_dim);
				for (int i = 0; i < _dim; i++)
					for (int j = 0; j < _dim; j++) {
						pair<unsigned int, unsigned int> p(i, j);
						if ((it = nmap.find(p)) != nmap.end())
							n(i, j) = it->second;
						else
							n(i, j) = 0;
					}

				unsigned int df = 0;
				double bowker = .0;
				double deltaS = 0.0;
				for (int i = 0; i < _dim; i++) {
					for (int j = i + 1; j < _dim; j++) {
						if (n(i, j) + n(j, i) > 0) {
							df++;
							bowker += (double) ((n(i, j) - n(j, i)) * (n(i, j) - n(j, i))) / (n(i, j) + n(j, i));
						}

						double x = n(j, i) - n(i, j);
						deltaS += (x / sum) * (x / sum);
					}
				}
				dfList.push_back(df);
				bowkerList.push_back(bowker);
				dsList.push_back(sqrt(deltaS));

				double deltaMs = 0.0;
				for (int i = 0; i < _dim; i++) {
					double rowSum = n.getRowSum(i);
					double colSum = n.getColSum(i);
					deltaMs += ((rowSum - colSum) / sum) * ((rowSum - colSum) / sum);
				}
				dmsList.push_back(sqrt(deltaMs) / sqrt(2.0));

				Matrix V(_dim - 1);
				for (unsigned int i = 0; i < _dim - 1; i++)
					for (unsigned int j = 0; j < _dim - 1; j++)
						if (i == j)
							V(i, j) = n.getRowSum(i) + n.getColSum(i) - 2 * n(i, i);
						else
							V(i, j) = -(n(i, j) + n(j, i));
				double stuart = 0;
				try {
					V.inverse();
					for (unsigned int i = 0; i < _dim - 1; i++) {
						double d_i = n.getRowSum(i) - n.getColSum(i);
						for (unsigned int j = 0; j < _dim - 1; j++) {
							double d_j = n.getRowSum(j) - n.getColSum(j);
							stuart += V(i, j) * d_i * d_j;
						}
					}
				} catch (string& s) {
					if (verbose)
						cerr << "Error while inverting matrix for Stuart's test (" << k << "x" << l << "): " << s << endl;
					stuart = numeric_limits<double>::quiet_NaN();
				}

				stuartList.push_back(stuart);

				if (df > 0)
					pBowkerList.push_back(gammq(df / 2.0, (bowker / 2.0)));
				else
					pBowkerList.push_back(1.0);

				unsigned int dfS = _dim - 1;
				unsigned int dfA = df - dfS;
				if (stuart != stuart) { // NaN
					pStuartList.push_back(numeric_limits<double>::quiet_NaN());
					pAbabnehList.push_back(numeric_limits<double>::quiet_NaN());
				} else {
					if (dfS > 0)
						pStuartList.push_back(gammq(dfS / 2.0, (stuart / 2.0)));
					else
						pStuartList.push_back(1.0);

					if (dfA > 0)
						pAbabnehList.push_back(gammq(dfA / 2.0, ((bowker - stuart) / 2.0)));
					else
						pAbabnehList.push_back(1.0);
				}
			}
		}

		for (unsigned int k = 0; k < len; k++) // 1st sequence
			for (unsigned int l = k + 1; l < len; l++) // 2nd sequence
					{
				double aitchison = .0;
				for (unsigned int i = 0; i < _dim; i++) {
					double x_i = baseFrequencies[k][i];
					double y_i = baseFrequencies[l][i];
					for (unsigned int j = 0; j < _dim; j++) {
						double x_j = baseFrequencies[k][j];
						double y_j = baseFrequencies[l][j];
						double logDiff = log(x_i / x_j) - log(y_i / y_j);
						aitchison += (logDiff * logDiff) / (2 * _dim);
					}
				}
				aitchisonList.push_back(aitchison);
			}

		_df.push_back(dfList);
		_bowker.push_back(bowkerList);
		_pBowker.push_back(pBowkerList);
		_ds.push_back(dsList);
		_dms.push_back(dmsList);
		_stuart.push_back(stuartList);
		_pStuart.push_back(pStuartList);
		_pAbabneh.push_back(pAbabnehList);
		_aitchison.push_back(aitchisonList);
	}
}

void Alignment::writeSummary(string prefix, int windowSize, int windowStep) {
	string resultsFileName = prefix + ".summary.csv";
	ofstream resultsFile;
	resultsFile.open(resultsFileName.c_str(), ifstream::trunc);
	if (!resultsFile.is_open())
		throw("Error, cannot open file " + resultsFileName);

	cout << "Writing summary to: " << resultsFileName << endl;

	unsigned int len = _alignment.size();
	if (windowSize <= 0)
		windowSize = _cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	resultsFile << "Seq1\tSeq2";
	resultsFile << "\tBowker (B)\tdf_B\tp_B";
	resultsFile << "\tStuart (S)\tdf_S\tp_S";
	resultsFile << "\tAbabneh (A)\tdf_A\tp_A";
	resultsFile << "\tAitchison\tDelta_s\tDelta_ms";
	resultsFile << "\tStart\tEnd" << endl;

	unsigned int i = 0;
	for (unsigned int windowStart = 0; windowStart + windowSize <= _cols; windowStart += windowStep) {
		vector<unsigned int> count(10, 0);
		double minP = 1.0;
		unsigned int j = 0;
		for (unsigned int k = 0; k < len; k++)
			for (unsigned int l = k + 1; l < len; l++) {

				resultsFile << _alignment[k].getName() << "\t" << _alignment[l].getName() << scientific;

				unsigned int dfB = _df[i][j];
				resultsFile << "\t" << _bowker[i][j] << "\t" << dfB << "\t" << _pBowker[i][j];

				unsigned int dfS = _dim - 1;
				unsigned int dfA = dfB - dfS;
				if (_stuart[i][j] != _stuart[i][j]) { // NaN
					resultsFile << "\t" << "n/a" << "\t" << "n/a" << "\t" << "n/a";
					resultsFile << "\t" << "n/a" << "\t" << "n/a" << "\t" << "n/a";
				} else {
					resultsFile << "\t" << _stuart[i][j] << "\t" << dfS << "\t" << _pStuart[i][j];

					double ababneh = _bowker[i][j] - _stuart[i][j];
					resultsFile << "\t" << ababneh << "\t" << dfA << "\t" << _pAbabneh[i][j];
				}

				resultsFile << "\t" << _aitchison[i][j] << "\t" << _ds[i][j] << "\t" << _dms[i][j];

				resultsFile << "\t" << windowStart << "\t" << windowStart + windowSize - 1 << endl;

				if (_pBowker[i][j] < minP)
					minP = _pBowker[i][j];
				if (_pBowker[i][j] < 0.05)
					count[0]++;
				if (_pBowker[i][j] < 0.01)
					count[1]++;
				if (_pBowker[i][j] < 0.005)
					count[2]++;
				if (_pBowker[i][j] < 0.001)
					count[3]++;
				if (_pBowker[i][j] < 0.0005)
					count[4]++;
				if (_pBowker[i][j] < 0.0001)
					count[5]++;
				if (_pBowker[i][j] < 0.00005)
					count[6]++;
				if (_pBowker[i][j] < 0.00001)
					count[7]++;
				if (_pBowker[i][j] < 0.000005)
					count[8]++;
				if (_pBowker[i][j] < 0.000001)
					count[9]++;

				j++;
			}

		cout << endl << "Highlights from the analysis (window " << windowStart << "-" << windowStart + windowSize - 1 << "):" << endl;
		cout.precision(2);
		if (minP < 0.05)
			cout << "P-values < 0.05                 " << setw(8) << count[0] << " (" << fixed << (double) count[0] * 100 / j << "%)" << endl;
		if (minP < 0.01)
			cout << "P-values < 0.01                 " << setw(8) << count[1] << " (" << fixed << (double) count[1] * 100 / j << "%)" << endl;
		if (minP < 0.005)
			cout << "P-values < 0.005                " << setw(8) << count[2] << " (" << fixed << (double) count[2] * 100 / j << "%)" << endl;
		if (minP < 0.001)
			cout << "P-values < 0.001                " << setw(8) << count[3] << " (" << fixed << (double) count[3] * 100 / j << "%)" << endl;
		if (minP < 0.0005)
			cout << "P-values < 0.0005               " << setw(8) << count[4] << " (" << fixed << (double) count[4] * 100 / j << "%)" << endl;
		if (minP < 0.0001)
			cout << "P-values < 0.0001               " << setw(8) << count[5] << " (" << fixed << (double) count[5] * 100 / j << "%)" << endl;
		if (minP < 0.00005)
			cout << "P-values < 0.00005              " << setw(8) << count[6] << " (" << fixed << (double) count[6] * 100 / j << "%)" << endl;
		if (minP < 0.00001)
			cout << "P-values < 0.00001              " << setw(8) << count[7] << " (" << fixed << (double) count[7] * 100 / j << "%)" << endl;
		if (minP < 0.000005)
			cout << "P-values < 0.000005             " << setw(8) << count[8] << " (" << fixed << (double) count[8] * 100 / j << "%)" << endl;
		if (minP < 0.000001)
			cout << "P-values < 0.000001             " << setw(8) << count[9] << " (" << fixed << (double) count[9] * 100 / j << "%)" << endl;
		cout << "Number of tests                   " << setw(15) << j << endl;
		cout.precision(8);
		cout << "Smallest P-value                  " << setw(15) << scientific << minP << endl;

		i++;
	}

	resultsFile.close();
}

void Alignment::writeExtendedResults(string prefix, int windowSize, int windowStep) {
	string bowkerFileName = prefix + ".bowker.csv";
	string stuartFileName = prefix + ".stuart.csv";
	string ababnehFileName = prefix + ".ababneh.csv";
	string delta_sFileName = prefix + ".delta_s.dis";
	string delta_msFileName = prefix + ".delta_ms.dis";
	string aitchisonFileName = prefix + ".aitchison.dis";

	cout << endl;
	cout << "Writing extended results to:" << endl;
	cout << "  Bowker\'s test:               " << bowkerFileName << endl;
	cout << "  Stuart\'s test:               " << stuartFileName << endl;
	cout << "  Ababneh\'s test:              " << ababnehFileName << endl;
	cout << "  Aitchison\'s distance matrix: " << aitchisonFileName << endl;
	cout << "  delta_s distance matrix:     " << delta_sFileName << endl;
	cout << "  delta_ms distance matrix:    " << delta_msFileName << endl;

	unsigned int len = _alignment.size();
	if (windowSize <= 0)
		windowSize = _cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	ofstream bowkerFile, stuartFile, ababnehFile, aitchisonFile, delta_sFile, delta_msFile;
	bowkerFile.open(bowkerFileName.c_str(), ifstream::trunc);
	if (!bowkerFile.is_open())
		throw("Error, cannot open file " + bowkerFileName);

	stuartFile.open(stuartFileName.c_str(), ifstream::trunc);
	if (!stuartFile.is_open())
		throw("Error, cannot open file " + stuartFileName);

	ababnehFile.open(ababnehFileName.c_str(), ifstream::trunc);
	if (!ababnehFile.is_open())
		throw("Error, cannot open file " + ababnehFileName);

	aitchisonFile.open(aitchisonFileName.c_str(), ifstream::trunc);
	if (!aitchisonFile.is_open())
		throw("Error, cannot open file " + aitchisonFileName);

	delta_sFile.open(delta_sFileName.c_str(), ifstream::trunc);
	if (!delta_sFile.is_open())
		throw("Error, cannot open file " + delta_sFileName);

	delta_msFile.open(delta_msFileName.c_str(), ifstream::trunc);
	if (!delta_msFile.is_open())
		throw("Error, cannot open file " + delta_msFileName);

	bowkerFile.flags(ios::left);
	stuartFile.flags(ios::left);
	ababnehFile.flags(ios::left);
	aitchisonFile.flags(ios::left);
	delta_sFile.flags(ios::left);
	delta_msFile.flags(ios::left);

	unsigned int i = 0;
	for (unsigned int windowStart = 0; windowStart + windowSize <= _cols; windowStart += windowStep) {
		bowkerFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
		stuartFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
		ababnehFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
		aitchisonFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
		delta_sFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;
		delta_msFile << windowStart << "-" << setw(6) << windowStart + windowSize - 1;

		for (unsigned int l = 0; l < len; l++) {
			bowkerFile << "\t" << setw(12) << _alignment[l].getName();
			stuartFile << "\t" << setw(12) << _alignment[l].getName();
			ababnehFile << "\t" << setw(12) << _alignment[l].getName();
			aitchisonFile << "\t" << setw(12) << _alignment[l].getName();
			delta_sFile << "\t" << setw(12) << _alignment[l].getName();
			delta_msFile << "\t" << setw(12) << _alignment[l].getName();
		}
		bowkerFile << endl;
		stuartFile << endl;
		ababnehFile << endl;
		aitchisonFile << endl;
		delta_sFile << endl;
		delta_msFile << endl;

		for (unsigned int k = 0; k < len; k++) {
			bowkerFile.flags(ios::left);
			stuartFile.flags(ios::left);
			ababnehFile.flags(ios::left);
			aitchisonFile.flags(ios::left);
			delta_sFile.flags(ios::left);
			delta_msFile.flags(ios::left);
			bowkerFile << setw(12) << _alignment[k].getName();
			stuartFile << setw(12) << _alignment[k].getName();
			ababnehFile << setw(12) << _alignment[k].getName();
			aitchisonFile << setw(12) << _alignment[k].getName();
			delta_sFile << setw(12) << _alignment[k].getName();
			delta_msFile << setw(12) << _alignment[k].getName();
			bowkerFile.flags(ios::right);
			stuartFile.flags(ios::right);
			ababnehFile.flags(ios::right);
			aitchisonFile.flags(ios::right);
			delta_sFile.flags(ios::right);
			delta_msFile.flags(ios::right);
			for (unsigned int l = 0; l < len; l++) {
				if (k == l) {
					bowkerFile    << "\t0";
					stuartFile    << "\t0";
					ababnehFile   << "\t0";
					aitchisonFile << "\t0";
					delta_sFile   << "\t0";
					delta_msFile  << "\t0";
				} else {
					int m;
					if (k < l)
						m = k*(len-1) - (k-1)*k/2 + l - k - 1;
					else
						m = l*(len-1) - (l-1)*l/2 + k - l - 1;

					bowkerFile << "\t" << scientific << _pBowker[i][m];
					if (_pStuart[i][m] != _pStuart[i][m]) { // NaN
						stuartFile << "\t" << "n/a";
						ababnehFile << "\t" << "n/a";
					} else {
						stuartFile << "\t" << scientific << _pStuart[i][m];
						ababnehFile << "\t" << scientific << _pAbabneh[i][m];
					}
					aitchisonFile << "\t" << scientific << _aitchison[i][m];
					delta_sFile << "\t" << scientific << _ds[i][m];
					delta_msFile << "\t" << scientific << _dms[i][m];
				}
			}
			bowkerFile << endl;
			stuartFile << endl;
			ababnehFile << endl;
			aitchisonFile << endl;
			delta_sFile << endl;
			delta_msFile << endl;
		}
		i++;
	}

	bowkerFile.close();
	stuartFile.close();
	ababnehFile.close();
	aitchisonFile.close();
	delta_sFile.close();
	delta_msFile.close();
}
