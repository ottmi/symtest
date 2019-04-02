#include "helper.h"
#include "globals.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdio>
using namespace std;

string printTime(long t)
{
	stringstream s;
	if (t > 3600)
	{
		s << t / 3600 << ":" << setfill('0') << setw(2);
		t = t % 3600;
	}
	s << t / 60 << ":" << setfill('0') << setw(2) << t % 60;

	return s.str();
}

istream& safeGetline(istream& is, string& t)
{
	/* Courtesy of http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf */
	t.clear();
	istream::sentry se(is);
	streambuf* sb = is.rdbuf();

	for (;;)
	{
		int c = sb->sbumpc();
		switch (c)
		{
			case '\r':
				c = sb->sgetc();
				if (c == '\n') sb->sbumpc();
				return is;
			case '\n':
			case EOF:
				return is;
			default:
				t += (char) c;
		}
	}
}

string adjustString(string s, bool upercase)
{
	string r = "";

	for (unsigned int i = 0; i < s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != '\r' && c != ' ')
		{
			if (upercase)
				r += toupper(c);
			else
				r += c;
		}
	}

	return (r);
}

unsigned int mapCharToNum(string s, charMap_t& m)
{
	unsigned int d = 0;
	for (unsigned int i = 0; i < s.size(); i++) {
		d = d << 8;
		charMap_t::iterator it = m.find(toupper(s[i]));
		if (it == m.end()) {
			d+= m['?'];
		} else {
			d+= it->second;
		}
	}
	return d;
}

string mapNumToChar(unsigned int n, int groupSize, string m)
{
	string s;
	for (int i = 0; i < groupSize; i++)
	{
		s = m[n & 255] + s;
		n = n >> 8;
	}

	return s;
}

