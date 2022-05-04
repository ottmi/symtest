#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cctype>
#include <sstream>
#include <map>
#include <utility>
#include "AlignmentReader.h"
#include "Alignment.h"
#include "Matrix.h"
#include "chi.h"
#include "heatmap.h"
#include "helper.h"

Alignment::Alignment() {
	_cols = 0;
}

Alignment::Alignment(Options *options) {
	AlignmentReader alignmentReader(options->inputAlignment);
	_alignment = alignmentReader.getSequences(options->columnFrom, options->columnTo, options->listOfSequences);
	_cols = alignmentReader.getCols();
	if (_cols % options->groupLength != 0) {
		cout << "ERROR: Grouping is active, but the alignment length of " << _cols << " columns is not a multiple of the group length of " << options->groupLength << endl;
		exit(255);
	}
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

	if (_dataType != _DNA_DATA && options->grouping.size() != 1) {
		cout << "ERROR: The input data is " << dataTypeDesc[_dataType] << " and grouping is active. This makes no sense." << endl;
		exit(254);
	}

	_groupSize = options->grouping.size();

	int dim = 36;
	if (_dataType == _DNA_DATA)
		dim = 4;
	else if (_dataType == _AA_DATA)
		dim = 20;
	_dim = 1;
	for (int i = 0; i<_groupSize; i++) {
		_dim*= dim;
	}


	for (unsigned int i = 0; i < _alignment.size(); i++) {
		_alignment[i].translateToNum(_dataType, options);
		if (verbose >= 2) {
			cout << i << ": ";
			for (unsigned int j = 0; j < _alignment[i].getLength() / options->groupLength; j++)
				cout << mapNumToChar(_alignment[i].getNumerical(j), _dataType, _groupSize) << " ";
			cout << endl;
		}
	}

	cout << "Groupsize is " << _groupSize << ", matrix dimension will be " << _dim << endl;

	unsigned int n = 0;
	unsigned int i_max, j_max;
	if (_dataType == _DNA_DATA) {
		(_groupSize > 2) ? i_max = 4 : i_max = 0;
		(_groupSize > 1) ? j_max = 4 : j_max = 0;
		for (unsigned int i=0; i<=i_max; i++) {
			if (i == 4)
				n-= 16;
			for (unsigned int j=0; j<=j_max; j++) {
				unsigned int c = (i << 16) | (j << 8);
				if (j == 4)
					n-= 4;
				for (unsigned int k=0; k<=4; k++) {
					_mapEnumeration[c] = n;
					c++;
					if (k != 3)
						n++;
				}
			}
		}
	} else if (_dataType == _AA_DATA) {
		for (unsigned int c=0; c<20; c++) {
			_mapEnumeration[c] = c;
		}
	} else {
		for (unsigned int c=0; c<36; c++) {
			_mapEnumeration[c] = c;
		}
	}
}

Alignment::~Alignment() {
}

void Alignment::computeTestsOfSymmetry(int windowSize, int windowStep) {
	cout << endl << "Performing tests of pairwise symmetry" << endl;

	unsigned int len = _alignment.size();
	if (windowSize <= 0)
		windowSize = _cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	if (windowSize < (int) _cols)
		cout << "  WindowSize=" << windowSize << " StepWidth=" << windowStep << " _cols=" << _cols << endl;

	cout.precision(6);
	for (unsigned int windowStart = 0; windowStart < _cols; windowStart += windowStep) {
		vector<unsigned int> dfList;
		vector<double> bowkerList; // Bowker's test
		vector<double> pBowkerList; // Bowker's test probabilities
		vector<double> dsList;  // Euclidian Distance (full symmetry)
		vector<double> dmsList; // Euclidian Distance (marginal symmetry)
		vector<double> stuartList; // Stuart's test
		vector<double> pStuartList; // Stuart's test probabilities
		vector<double> pAbabnehList; // Ababneh's test probabilities
		vector<double> aitchisonMargList; // Aitchison's distance (marginal symmetry)
		vector<double> aitchisonFullList; // Aitchison's distance (full symmetry)

		for (unsigned int k = 0; k < len; k++) { // 1st sequence
			Sequence s1 = _alignment[k];
			for (unsigned int l = k + 1; l < len; l++) { // 2nd sequence
				Sequence s2 = _alignment[l];
				
				/* Counting the combinations of characters in the alignment */
				unsigned int sum = 0; // Sum of different combinations of characters
				map<pair<unsigned int, unsigned int>, unsigned int> nmap; // Map to count the different combinations of characters
				map<pair<unsigned int, unsigned int>, unsigned int>::iterator it;
				unsigned int windowEnd = windowStart + windowSize;
				if (windowEnd > _cols)
					windowEnd = _cols;
				for (unsigned int m = windowStart; m < windowEnd; m++) {
					unsigned int c1 = s1.getNumerical(m);
					unsigned int c2 = s2.getNumerical(m);
					if (s1.charIsUnambiguous(c1) && s2.charIsUnambiguous(c2)) {
						pair<unsigned int, unsigned int> p(_mapEnumeration[c1], _mapEnumeration[c2]);
						if ((it = nmap.find(p)) != nmap.end())
							it->second = it->second + 1;
						else
							nmap.insert(pair<pair<unsigned int, unsigned int>, unsigned int>(p, 1));
						sum++;
					}
				}

				/* Creating the divergence matrix that stores the combinations of characters */
				Matrix n(_dim);
				for (int i = 0; i < _dim; i++)
					for (int j = 0; j < _dim; j++) {
						pair<unsigned int, unsigned int> p(i, j);
						if ((it = nmap.find(p)) != nmap.end())
							n(i, j) = it->second;
						else
							n(i, j) = 0;
					}

				/* Compute the test statistic and degrees of freedom for Bowker's (1948) test */
				unsigned int df = 0;
				double bowker = .0;
				for (int i = 0; i < _dim; i++) {
					for (int j = i + 1; j < _dim; j++) {
						if (n(i, j) + n(j, i) > 0) {
							df++;
							bowker += (double) ((n(i, j) - n(j, i)) * (n(i, j) - n(j, i))) / (n(i, j) + n(j, i));
						}
					}
				}
				dfList.push_back(df);
				bowkerList.push_back(bowker);
				if (df > 0)
					pBowkerList.push_back(xChi_Square_Distribution_Tail(bowker, df));
				else
					pBowkerList.push_back(1.0);

				/* Compute Euclidean distance (full sym.) */
				double deltaS = 0.0;
				for (int i = 0; i < _dim; i++) {
					for (int j = i + 1; j < _dim; j++) {
						double x = n(j, i) - n(i, j);
						deltaS += (x / sum) * (x / sum);
					}
				}
				dsList.push_back(sqrt(deltaS));

				/* Compute Euclidean distance (mar. sym.) */
				double deltaMs = 0.0; // Euclidean distance (mar. sym.)
				for (int i = 0; i < _dim; i++) {
					double rowSum = n.getRowSum(i);
					double colSum = n.getColSum(i);
					deltaMs += ((rowSum - colSum) / sum) * ((rowSum - colSum) / sum);
				}
				dmsList.push_back(sqrt(deltaMs) / sqrt(2.0));

				/* Compute the test statistic for Stuart's test */
				Matrix V(_dim - 1);
				for (int i = 0; i < _dim - 1; i++)
					for (int j = 0; j < _dim - 1; j++)
						if (i == j)
							V(i, j) = n.getRowSum(i) + n.getColSum(i) - 2 * n(i, i);
						else
							V(i, j) = -(n(i, j) + n(j, i));
				double stuart = 0;
				try {
					V.inverse();
					for (int i = 0; i < _dim - 1; i++) {
						double d_i = n.getRowSum(i) - n.getColSum(i);
						for (int j = 0; j < _dim - 1; j++) {
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

				/* Compute probabilities for Stuart's and Ababneh's test */
				unsigned int dfS = _dim - 1;
				unsigned int dfA = df - dfS;
				if (stuart != stuart) { // NaN
					pStuartList.push_back(numeric_limits<double>::quiet_NaN());
					pAbabnehList.push_back(numeric_limits<double>::quiet_NaN());
				} else {
					if (dfS > 0)
						pStuartList.push_back(xChi_Square_Distribution_Tail(stuart, dfS));
					else
						pStuartList.push_back(1.0);

					if (dfA > 0)
						pAbabnehList.push_back(xChi_Square_Distribution_Tail(fabs(bowker - stuart), dfA));
					else
						pAbabnehList.push_back(1.0);
				}

				/* Compute Aitchison's distance (marg. sym.) */
				double prior = 1.0;
				double aitchisonMarg = 0.0;
				int c = 0;
				for (int i = 0; i < _dim; i++) {
					double x_i = n.getRowSum(i);
					double y_i = n.getColSum(i);
					if (verbose) {
						cout << "rowSum[" << i << "]=" << x_i << " colSum[" << i << "]=" << y_i << endl;
					}

					if (x_i > 0 && y_i > 0) {
						c++;
						x_i += prior;
						y_i += prior;
						for (int j = 0; j < _dim; j++) {
							double x_j = n.getRowSum(j) + prior;
							double y_j = n.getColSum(j) + prior;
							double logDiff = log(x_i / x_j) - log(y_i / y_j);
							aitchisonMarg += (logDiff * logDiff);
						}
					}
				}
				aitchisonMarg = sqrt(aitchisonMarg / (2 * c));
				aitchisonMargList.push_back(aitchisonMarg);

				/* Compute Aitchison's distance (full sym.) */
				double aitchisonFull = 0.0;
				c = 0;
				for (int g = 0; g < _dim; g++) {
					for (int h = g + 1; h < _dim; h++) {
						double y_i = n(g, h);
						double z_i = n(h, g);
						if (y_i > 0 && z_i > 0) {
							if (verbose) {
								cout << "y[" << c << "]=" << y_i << " z[" << c << "]=" << z_i << endl;
							}
							y_i += prior;
							z_i += prior;
							int d = 0;
							for (int i = 0; i < _dim; i++) {
								for (int j = i + 1; j < _dim; j++) {
									double y_j = n(i, j);
									double z_j = n(j, i);
									if (verbose) {
										cout << "  y[" << d << "]=" << y_j << " z[" << d << "]=" << z_j << endl;
									}
									y_j += prior;
									z_j += prior;
									double logDiff = log(y_i / y_j) - log(z_i / z_j);
									aitchisonFull += (logDiff * logDiff);
									d++;
								}
							}
							c++;
						}
					}
				}
				aitchisonFull = sqrt(aitchisonFull / (2 * c));
				aitchisonFullList.push_back(aitchisonFull);
			}
		}

		_df.push_back(dfList);
		_bowker.push_back(bowkerList);
		_pBowker.push_back(pBowkerList);
		_ds.push_back(dsList);
		_dms.push_back(dmsList);
		_stuart.push_back(stuartList);
		_pStuart.push_back(pStuartList);
		_pAbabneh.push_back(pAbabnehList);
		_aitchisonMarg.push_back(aitchisonMargList);
		_aitchisonFull.push_back(aitchisonFullList);
	}
}

void Alignment::printStatistics(vector<double> &pValues, char id) {
	double minP = 1.0;
	double sumP = 0;
	vector<unsigned int> count(10, 0);

	size_t total = pValues.size();
	for (unsigned int j = 0; j < total; j++) {
		if (pValues[j] < minP)
			minP = pValues[j];
		if (pValues[j] < 0.05)
			count[0]++;
		if (pValues[j] < 0.01)
			count[1]++;
		if (pValues[j] < 0.005)
			count[2]++;
		if (pValues[j] < 0.001)
			count[3]++;
		if (pValues[j] < 0.0005)
			count[4]++;
		if (pValues[j] < 0.0001)
			count[5]++;
		if (pValues[j] < 0.00005)
			count[6]++;
		if (pValues[j] < 0.00001)
			count[7]++;
		if (pValues[j] < 0.000005)
			count[8]++;
		if (pValues[j] < 0.000001)
			count[9]++;
		sumP += pValues[j];
	}

	vector<double> sortedP(pValues.begin(), pValues.end());
	sort(sortedP.begin(), sortedP.end());
	double medP = sortedP[sortedP.size() / 2];
	if (sortedP.size() % 2 == 0) {
		medP += sortedP[sortedP.size() / 2 - 1];
		medP /= 2;
	}

	cout.precision(2);
	if (minP < 0.05)
		cout << "  " << id << ": P-values < 0.05            " << setw(8) << count[0] << " (" << fixed << (double) count[0] * 100 / total << "%)" << endl;
	if (minP < 0.01)
		cout << "  " << id << ": P-values < 0.01            " << setw(8) << count[1] << " (" << fixed << (double) count[1] * 100 / total << "%)" << endl;
	if (minP < 0.005)
		cout << "  " << id << ": P-values < 0.005           " << setw(8) << count[2] << " (" << fixed << (double) count[2] * 100 / total << "%)" << endl;
	if (minP < 0.001)
		cout << "  " << id << ": P-values < 0.001           " << setw(8) << count[3] << " (" << fixed << (double) count[3] * 100 / total << "%)" << endl;
	if (minP < 0.0005)
		cout << "  " << id << ": P-values < 0.0005          " << setw(8) << count[4] << " (" << fixed << (double) count[4] * 100 / total << "%)" << endl;
	if (minP < 0.0001)
		cout << "  " << id << ": P-values < 0.0001          " << setw(8) << count[5] << " (" << fixed << (double) count[5] * 100 / total << "%)" << endl;
	if (minP < 0.00005)
		cout << "  " << id << ": P-values < 0.00005         " << setw(8) << count[6] << " (" << fixed << (double) count[6] * 100 / total << "%)" << endl;
	if (minP < 0.00001)
		cout << "  " << id << ": P-values < 0.00001         " << setw(8) << count[7] << " (" << fixed << (double) count[7] * 100 / total << "%)" << endl;
	if (minP < 0.000005)
		cout << "  " << id << ": P-values < 0.000005        " << setw(8) << count[8] << " (" << fixed << (double) count[8] * 100 / total << "%)" << endl;
	if (minP < 0.000001)
		cout << "  " << id << ": P-values < 0.000001        " << setw(8) << count[9] << " (" << fixed << (double) count[9] * 100 / total << "%)" << endl;
	cout << "  " << id << ": Number of tests              " << setw(15) << total << endl;
	cout.precision(8);
	cout << "  " << id << ": Median P-value               " << setw(15) << fixed << medP << endl;
	cout << "  " << id << ": Average P-value              " << setw(15) << fixed << sumP / total << endl;
	cout << "  " << id << ": Smallest P-value             " << setw(15) << fixed << minP << endl;
}

void Alignment::writeResults(Options* options) {
	string resultsFileName = options->prefix + ".summary.csv";
	ofstream resultsFile;
	resultsFile.open(resultsFileName.c_str(), ifstream::trunc);
	if (!resultsFile.is_open())
		throw("Error, cannot open file " + resultsFileName);

	cout << "Writing summary to: " << resultsFileName << endl;

	unsigned int len = _alignment.size();
	unsigned int windowSize = options->windowSize;
	unsigned int windowStep = options->windowStep;
	if (windowSize <= 0)
		windowSize = _cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	resultsFile << "Seq1" << CSV_SEPARATOR;
	resultsFile << "Seq2" << CSV_SEPARATOR;
	resultsFile << "Bowker (B)" << CSV_SEPARATOR << "df_B" << CSV_SEPARATOR << "p_B" << CSV_SEPARATOR;
	resultsFile << "Stuart (S)" << CSV_SEPARATOR << "df_S" << CSV_SEPARATOR << "p_S" << CSV_SEPARATOR;
	resultsFile << "Ababneh (A)" << CSV_SEPARATOR << "df_A" << CSV_SEPARATOR << "p_A" << CSV_SEPARATOR;
	resultsFile << "d_EMS" << CSV_SEPARATOR << "d_EFS" << CSV_SEPARATOR;
	resultsFile << "d_AMS" << CSV_SEPARATOR << "d_AFS" << CSV_SEPARATOR;
	resultsFile << "Start" << CSV_SEPARATOR << "End" << endl;

	unsigned int i = 0;
	for (unsigned int windowStart = 0; windowStart < _cols; windowStart += windowStep) {
		unsigned int windowEnd = windowStart + windowSize;
		if (windowEnd > _cols)
			windowEnd = _cols;

		unsigned int j = 0;
		for (unsigned int k = 0; k < len; k++) {
			for (unsigned int l = k + 1; l < len; l++) {
				resultsFile << _alignment[k].getName() << CSV_SEPARATOR << _alignment[l].getName() << fixed;

				unsigned int dfB = _df[i][j];
				resultsFile << CSV_SEPARATOR << _bowker[i][j] << CSV_SEPARATOR << dfB << CSV_SEPARATOR << _pBowker[i][j];

				unsigned int dfS = _dim - 1;
				unsigned int dfA = dfB - dfS;
				if (_stuart[i][j] != _stuart[i][j]) { // NaN
					resultsFile << CSV_SEPARATOR << "n/a" << CSV_SEPARATOR << "n/a" << CSV_SEPARATOR << "n/a";
					resultsFile << CSV_SEPARATOR << "n/a" << CSV_SEPARATOR << "n/a" << CSV_SEPARATOR << "n/a";
				} else {
					resultsFile << CSV_SEPARATOR << _stuart[i][j] << CSV_SEPARATOR << dfS << CSV_SEPARATOR << _pStuart[i][j];

					double ababneh = fabs(_bowker[i][j] - _stuart[i][j]);
					resultsFile << CSV_SEPARATOR << ababneh << CSV_SEPARATOR << dfA << CSV_SEPARATOR << _pAbabneh[i][j];
				}

				resultsFile << CSV_SEPARATOR << _dms[i][j] << CSV_SEPARATOR << _ds[i][j];
				resultsFile << CSV_SEPARATOR << _aitchisonMarg[i][j] << CSV_SEPARATOR << _aitchisonFull[i][j];
				resultsFile << CSV_SEPARATOR << windowStart << CSV_SEPARATOR << windowEnd - 1 << endl;

				j++;
			}
		}

		cout << endl << "Highlights from the analysis (window " << windowStart << "-" << windowEnd - 1 << "):" << endl;

		cout << "  Bowker's test:" << endl;
		printStatistics(_pBowker[i], 'B');
		cout << endl;

		cout << "  Stuarts's test:" << endl;
		printStatistics(_pStuart[i], 'S');
		cout << endl;

		cout << "  Ababneh et al.'s test:" << endl;
		printStatistics(_pAbabneh[i], 'A');

		i++;
	}

	resultsFile.close();

	if (options->writeExtendedTestResults) {
		cout << endl;
		cout << "Writing extended results/heatmaps to:" << endl;
		if (options->writeBowkerFile)
			writeExtendedResult("Bowker\'s test             ", options->prefix + ".bowker.", "csv", windowSize, windowStep, _pBowker);
		if (options->writeStuartFile)
			writeExtendedResult("Stuart\'s test             ", options->prefix + ".stuart.", "csv", windowSize, windowStep, _pStuart);
		if (options->writeAbabnehFile)
			writeExtendedResult("Ababneh\'s test            ", options->prefix + ".ababneh.", "csv", windowSize, windowStep, _pAbabneh);

		cout << "Writing extended distances to:" << endl;
		if (options->writeAmsFile)
			writeExtendedDistances("Aitchison\'s distances (marg)", options->prefix + ".AMS.", "dst", windowSize, windowStep, _aitchisonMarg);
		if (options->writeAfsFile)
			writeExtendedDistances("Aitchison\'s distances (full)", options->prefix + ".AFS.", "dst", windowSize, windowStep, _aitchisonFull);
		if (options->writeEmsFile)
			writeExtendedDistances("Euclidian distances (marg)   ", options->prefix + ".EMS.", "dst", windowSize, windowStep, _dms);
		if (options->writeEfsFile)
			writeExtendedDistances("Euclidian distances (full)   ", options->prefix + ".EFS.", "dst", windowSize, windowStep, _ds);
	}

}

void Alignment::writeExtendedResult(string title, string baseName, string ext, unsigned int windowSize, unsigned int windowStep, vector<vector<double> >& matrix) {
	ofstream outFile;
	string outFileName;
	cout.flags(ios::left);
	cout << "  " << title << endl;
	if (windowSize < _cols) {
		cout << "    Extended results           " << baseName << "<window>." << ext << endl;
		cout << "    Triangular matrix          " << baseName << "<window>.triangular.svg" << endl;
		cout << "    Full matrix                " << baseName << "<window>.full.svg" << endl;
	} else {
		outFileName = baseName + ext;
		cout << "    Extended results           " << outFileName << endl;
		cout << "    Triangular matrix          " << baseName << "triangular.svg" << endl;
		cout << "    Full matrix                " << baseName << "full.svg" << endl;
	}

	unsigned int len = _alignment.size();
	unsigned int i = 0;
	for (unsigned int windowStart = 0; windowStart < _cols; windowStart += windowStep) {
		if (windowSize < _cols) {
			unsigned int windowEnd = windowStart + windowSize - 1;
			if (windowEnd > _cols - 1)
				windowEnd = _cols - 1;

			stringstream ss;
			ss << windowStart << "-" << windowEnd << ".";
			outFileName = baseName + ss.str() + ext;
		}

		outFile.open(outFileName.c_str(), ifstream::trunc);
		if (!outFile.is_open())
			throw("Error, cannot open file " + outFileName);

		outFile.flags(ios::left);
		for (unsigned int l = 0; l < len; l++) {
			outFile << CSV_SEPARATOR << _alignment[l].getName();
		}
		outFile << endl;

		double *cij = new double[len*len];
		string *seqNames = new string[len];
		for (unsigned int k = 0; k < len; k++) {
			outFile.flags(ios::left);
			outFile << _alignment[k].getName();
			outFile.flags(ios::right);
			seqNames[k] = _alignment[k].getName();
			for (unsigned int l = 0; l < len; l++) {
				if (k == l) {
					outFile << CSV_SEPARATOR << "0";
				} else {
					int m;
					if (k < l)
						m = k * (len - 1) - (k - 1) * k / 2 + l - k - 1;
					else
						m = l * (len - 1) - (l - 1) * l / 2 + k - l - 1;
					cij[l*len + k] = matrix[i][m];
					if (matrix[i][m] != matrix[i][m]) { // NaN
						outFile << CSV_SEPARATOR << "n/a";
					} else {
						outFile << CSV_SEPARATOR << fixed << matrix[i][m];
					}
				}
			}
			outFile << endl;
		}
		outFile.close();
		outputFullHeatmap(outFileName.substr(0, outFileName.length()-ext.length()-1), seqNames, cij, len);
		outputTriHeatmap(outFileName.substr(0, outFileName.length()-ext.length()-1), seqNames, cij, len);
		i++;
		delete[] cij;
		delete[] seqNames;
	}
}

void Alignment::writeExtendedDistances(string title, string baseName, string ext, unsigned int windowSize, unsigned int windowStep, vector<vector<double> >& matrix) {
	ofstream outFile;
	string outFileName;
	cout.flags(ios::left);
	if (windowSize < _cols) {
		cout << "  " << setw(29) << title << baseName << "<window>." << ext << endl;
	} else {
		outFileName = baseName + ext;
		cout << "  " << setw(29) << title << outFileName << endl;
	}

	unsigned int len = _alignment.size();
	unsigned int i = 0;
	for (unsigned int windowStart = 0; windowStart < _cols; windowStart += windowStep) {
		if (windowSize < _cols) {
			unsigned int windowEnd = windowStart + windowSize - 1;
			if (windowEnd > _cols - 1)
				windowEnd = _cols - 1;

			stringstream ss;
			ss << windowStart << "-" << windowEnd << ".";
			outFileName = baseName + ss.str() + ext;
		}

		outFile.open(outFileName.c_str(), ifstream::trunc);
		if (!outFile.is_open())
			throw("Error, cannot open file " + outFileName);

		outFile.flags(ios::left);
		outFile << len << endl;

		for (unsigned int k = 0; k < len; k++) {
			string name = _alignment[k].getName();
			if (name.length() >= 10)
				outFile << name << " ";
			else
				outFile << setw(10) << left << name;

			outFile << right;
			for (unsigned int l = 0; l < len; l++) {
				if (k == l) {
					outFile << "0";
				} else {
					int m;
					if (k < l)
						m = k * (len - 1) - (k - 1) * k / 2 + l - k - 1;
					else
						m = l * (len - 1) - (l - 1) * l / 2 + k - l - 1;

					if (matrix[i][m] != matrix[i][m]) { // NaN
						outFile << "n/a";
					} else {
						outFile << fixed << matrix[i][m];
					}
				}
				if (l < len - 1)
					outFile << "\t";
			}
			outFile << endl;
		}
		outFile.close();
		i++;
	}
}

