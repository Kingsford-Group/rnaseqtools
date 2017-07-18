#ifndef __BAMKIT_H__
#define __BAMKIT_H__

#include "hit.h"

using namespace std;

class bamkit
{
public:
	bamkit(const string &bamfile);
	~bamkit();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

	int qcnt;			// single reads
	double qlen;		// single reads
	vector<int> ivec;	// insert size

public:
	int solve_count();
	int solve_strand();

};

#endif
