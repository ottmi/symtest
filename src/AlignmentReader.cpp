#include <sstream>
#include "globals.h"
#include "helper.h"
#include "AlignmentReader.h"

AlignmentReader::AlignmentReader(string fileName)
{
	_fileReader.open(fileName.c_str());
	if (! _fileReader.is_open())
		throw("Error, cannot open file " + fileName );

	cout << "Reading alignment file: " << fileName << endl;

	_rows = 0;
	_cols = 0;

	safeGetline(_fileReader, _lastLine);

	if (_lastLine[0] == '>')
	{
		cout << "The file appears to be in Fasta format." << endl;
		_format = _FASTA_FORMAT;
	} else
	{
		stringstream ss(_lastLine);
		ss >> _rows >> _cols;
		if (_rows && _cols)
		{
			cout << "The file appears to be in be in Phylip format (" << _rows << " rows, " << _cols << " columns)."<< endl;
			_format = _PHYLIP_FORMAT;
		} else
		{
			string ext = fileName.substr(fileName.find_last_of('.') + 1);
			if (!ext.compare("fsa") || !ext.compare("fst") || !ext.compare("fasta"))
			{
				cout << "According to its extension, this file should be in Fasta format." << endl;
				_format = _FASTA_FORMAT;
			} else if (!ext.compare("phy") || !ext.compare("phylip"))
			{
				cout << "According to its extension, this file should be in Phylip format." << endl;
				_format = _PHYLIP_FORMAT;
			} else
			{
				stringstream s;
				s << "Unable to detect alignment format.\n" << PROGNAME << " only supports the Fasta and sequential Phylip formats.";
				throw(s.str());
			}
		}
	}

	string whiteSpace = " \n\t";

	if (_format == _FASTA_FORMAT)
	{
		while ((!_fileReader.eof()) && _lastLine[0] != '>')
			safeGetline(_fileReader, _lastLine);

		while (!_fileReader.eof())
		{
			string header;
			string seq;
			header = _lastLine;
			_lastLine = "";
			while (!_fileReader.eof() && _lastLine[0] != '>')
			{
				safeGetline(_fileReader, _lastLine);
				if (_lastLine[0] != '>')
					seq += _lastLine;
			}
			seq = adjustString(seq, false);
			header = adjustString(header.substr(1), false);
			if (header.length() > 1 && seq.length())
			{
				Sequence s(header, seq);
				_sequences.push_back(s);
				_rows++;
				if (seq.length() > (unsigned int) _cols)
					_cols = seq.length();
			}
		}
	} else if (_format == _PHYLIP_FORMAT)
	{
		while (! _fileReader.eof())
	    {
	   		safeGetline(_fileReader, _lastLine);
	   		if (_lastLine.length())
	   		{
	   			int n = _lastLine.find_first_of(whiteSpace);
	   			string name, seq;
	   			if (n == -1) // there's no whitespace, so the sequence starts at pos 11
	   			{
	   	   			name = _lastLine.substr(0, 10);
	   	   			seq = _lastLine.substr(10);
	   			}
	   			else
	   			{
	   				name = _lastLine.substr(0, n);
	   				n = _lastLine.find_first_not_of(whiteSpace, n);
	   	   			seq = _lastLine.substr(n);
	   			}

	   			if ((int) seq.length() != _cols)
	   				cerr << "Sequence #" << _sequences.size() + 1 << " (" << name << ") consists of " << seq.length() << " characters when it should be " << _cols << "." << endl;
	   			seq = adjustString(seq, false);
	   			if (name.length() && seq.length())
	   			{
	   				Sequence s(name, seq);
	   				_sequences.push_back(s);
	   			}
	   		}
	    }
		if ((int) _sequences.size() < _rows)
			cerr << "The alignment contains only " << _sequences.size() << " rows, but it should be " << _rows << "."<< endl;
	}
}


AlignmentReader::~AlignmentReader()
{
	if (!_fileReader.is_open())
		_fileReader.close();
}


vector<Sequence> AlignmentReader::getSequences(int from, int to)
{
	if (from == -1 && to == -1)
		return _sequences;
	else
	{
		vector<Sequence> sequences;
		for (unsigned int i = 0; i < _sequences.size(); i++)
		{
			string name = _sequences[i].getName();
			string seq = _sequences[i].getSequence();
			if (to == -1)
			    sequences.push_back(Sequence(name, seq.substr(from-1)));
			else
			    sequences.push_back(Sequence(name, seq.substr(from-1, to-from+1)));
		}
		return sequences;
	}
}
