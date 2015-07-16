#include<cstdlib>
#include<string>
#include<iostream>
#include<fstream>

using namespace std;

// Output the triangular heatmap
void outputTriHeatmap(string prefixOut, string* seqNames, double* cij, int seqNum);

// Output the full heatmap
void outputFullHeatmap(string prefixOut, string* seqNames, double* cij, int seqNum);
