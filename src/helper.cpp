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
				if (c == '\n')
					sb->sbumpc();
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


unsigned int mapCharToNum(string s, int dataType)
{
    unsigned int d = 0;
    if (dataType == _DNA_DATA)
    {
	for (unsigned int i = 0; i < s.size(); i++)
	{
	    char c = s[i];
	    d = d << 8;
	    switch (c)
	    {
		case 'A':
		case 'a':
		    d += 0x00;
		    break;

		case 'C':
		case 'c':
		    d += 0x01;
		    break;

		case 'G':
		case 'g':
		    d += 0x02;
		    break;

		case 'T':
		case 't':
		    d += 0x03;
		    break;

		case 'U':
		case 'u':
		    d += 0x04;
		    break;

		case 'R':
		case 'r':
		    d += 0x05;
		    break;

		case 'Y':
		case 'y':
		    d += 0x06;
		    break;

		case 'K':
		case 'k':
		    d += 0x07;
		    break;

		case 'M':
		case 'm':
		    d += 0x08;
		    break;

		case 'S':
		case 's':
		    d += 0x09;
		    break;

		case 'W':
		case 'w':
		    d += 0x0a;
		    break;

		case 'B':
		case 'b':
		    d += 0x0b;
		    break;

		case 'D':
		case 'd':
		    d += 0x0c;
		    break;

		case 'H':
		case 'h':
		    d += 0x0d;
		    break;

		case 'V':
		case 'v':
		    d += 0x0e;
		    break;

		case 'N':
		case 'n':
		    d += 0x0f;
		    break;

		case '?':
		    d += 0x10;
		    break;

		default:
		    d += 0x11;
		    break;
	    }
	}
    } else if (dataType == _AA_DATA)
    {
	for (unsigned int i = 0; i < s.size(); i++)
	{
	    char c = s[i];
	    d = d << 8;
	    switch (c)
	    {
		case 'A': // unambiguous
		case 'a':
		    d += 0x00;
		    break;

		case 'C':
		case 'c':
		    d += 0x01;
		    break;

		case 'D':
		case 'd':
		    d += 0x02;
		    break;

		case 'E':
		case 'e':
		    d += 0x03;
		    break;

		case 'F':
		case 'f':
		    d += 0x04;
		    break;

		case 'G':
		case 'g':
		    d += 0x05;
		    break;

		case 'H':
		case 'h':
		    d += 0x06;
		    break;

		case 'I':
		case 'i':
		    d += 0x07;
		    break;

		case 'K':
		case 'k':
		    d += 0x08;
		    break;

		case 'L':
		case 'l':
		    d += 0x09;
		    break;

		case 'M':
		case 'm':
		    d += 0x0A;
		    break;

		case 'N':
		case 'n':
		    d += 0x0B;
		    break;

		case 'P':
		case 'p':
		    d += 0x0C;
		    break;

		case 'Q':
		case 'q':
		    d += 0x0D;
		    break;

		case 'R':
		case 'r':
		    d += 0x0E;
		    break;

		case 'S':
		case 's':
		    d += 0x0F;
		    break;

		case 'T':
		case 't':
		    d += 0x10;
		    break;

		case 'V':
		case 'v':
		    d += 0x11;
		    break;

		case 'W':
		case 'w':
		    d += 0x12;
		    break;

		case 'Y':
		case 'y':
		    d += 0x13;
		    break;

		case 'B': // ambiguous
		case 'b':
		    d += 0x14;
		    break;

		case 'J':
		case 'j':
		    d += 0x15;
		    break;

		case 'X':
		case 'x':
		    d += 0x16;
		    break;

		case 'Z':
		case 'z':
		    d += 0x17;
		    break;

		case '?':
		    d += 0x18;
		    break;

		default: // missing
		    d += 0x19;
		    break;
	    }
	}
    } else
    {
	for (unsigned int i = 0; i < s.size(); i++)
	{
	    char c = s[i];
	    d = d << 8;

	    if (c >= 65 && c <= 90) // A-Z
		d += c - 65;
	    else if (c >= 97 && c <= 122) // a-z
		d += c - 97;
	    else if (c >= 48 && c <= 57) // 0-9
		d += c - 48 + 26;
	    else if (c == '?')
		d += 36;
	    else
		d += 37;
	}
    }

    return d;
}

string mapNumToChar(unsigned int n, int dataType, int groupSize)
{
    string map, s;
    if (dataType == _DNA_DATA
    )
	map = _DNA_MAP;
    else if (dataType == _AA_DATA
    )
	map = _AA_MAP;
    else
	map = _ALPHANUM_MAP;

    for (unsigned int i = 0; i < groupSize; i++)
    {
	s = map[n & 255] + s;
	n = n >> 8;
    }

    return s;
}

/*
 From Numerical Recipes
 */

#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define FPMIN 1.0e-30
#define ITMAX 100

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an, b, c, d, del, h;

	*gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i = 1; i <= ITMAX; i++)
	{
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = b + an / c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS)
			break;
	}
	if (i > ITMAX)
		cerr << "a too large, ITMAX too small in gcf" << endl;
	*gammcf = exp(-x + a * log(x) - (*gln)) * h;
}

/*
 From Numerical Recipes
 */

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum, del, ap;

	*gln = gammln(a);
	if (x <= 0.0)
	{
		if (x < 0.0)
			cerr << "x less than 0 in routine gser" << endl;
		*gamser = 0.0;
		return;
	} else
	{
		ap = a;
		del = sum = 1.0 / a;
		for (n = 1; n <= ITMAX; n++)
		{
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS)
			{
				*gamser = sum * exp(-x + a * log(x) - (*gln));
				return;
			}
		}
		cerr << "a too large, ITMAX too small in routine gser" << endl;
		return;
	}
}

/*
 From Numerical Recipes
 */

double gammln(double xx)
{
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.00000000019001;
	for (j = 0; j <= 5; j++)
		ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}

/*
 From Numerical Recipes
 */

double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0)
		cerr << "Invalid arguments in routine gammq" << endl;
	if (x < (a + 1.0))
	{
		gser(&gamser, a, x, &gln);
		return 1.0 - gamser;
	} else
	{
		gcf(&gammcf, a, x, &gln);
		return gammcf;
	}
}
