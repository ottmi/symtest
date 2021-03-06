#ifndef HELPER_H_
#define HELPER_H_

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include <string>
#include <sstream>
using namespace std;

double factorial(int n);
string printTime(long t);
istream& safeGetline(istream& is, string& t);
string adjustString(string s, bool upercase=false);
unsigned int mapCharToNum(string s, int dataType);
string mapNumToChar(unsigned int n, int dataType, int groupSize);

template<class T> string str(const T& t)
{
	stringstream ss;
	ss << t;
	return ss.str();
}

#endif /* HELPER_H_ */
