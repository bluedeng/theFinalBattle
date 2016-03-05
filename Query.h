/*
 * Query.h
 *
 *  Created on: Sep 14, 2010
 *      Author: xiaoliwang
 */

#ifndef _QUERY_H_
#define _QUERY_H_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>

#include "Gram.h"

using namespace std;

class CQuery 
{
public:
	const string sequence;
	const unsigned length;
	unsigned constraint;
	unsigned threshold;
	unsigned countBound;
	unsigned gramsize;
	//ngrams of the query, gramCodes and gramCount of the ngrams
	vector<string> ngrams;
	vector<unsigned> gramCodes;
	map<unsigned, int> gramCount;

	CGram* gramGenUpper;

	CQuery(string sequence, CGram* gramgenupper, unsigned constraint):
		sequence(sequence), 
		length(sequence.length()),
		gramGenUpper(gramgenupper),
		constraint(constraint)
    {
		setQueryGrams();
	}

private:
	void setQueryGrams() 
    {
		vector<string> stemp;
		vector<unsigned> temp;
		gramGenUpper->decompose(sequence, stemp, temp);

		//move the duplicate from temp and count it
		for (unsigned i = 0; i < temp.size(); i++)
			if (++gramCount[temp[i]] == 1) {
				gramCodes.push_back(temp[i]);
				ngrams.push_back(stemp[i]);
			}

		gramsize = stemp.size();
    }

public:
	void print(FILE* pf)
	{
		fprintf(pf, "The query string: %s\n", sequence.c_str());
		fprintf(pf, "The query gram codes: ");
		for(unsigned i = 0; i < ngrams.size(); i++) {
			fprintf(pf, " (%s, %u)", ngrams[i].c_str(), gramCodes[i]);
		}
		fprintf(pf, "\n");
	}
};

#endif
