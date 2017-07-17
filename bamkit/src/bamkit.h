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

	int qcnt;		// single reads
	int pcnt;		// paired reads
	double qlen;	// single reads
	double plen;	// paired reads

public:
	int solve();

};

#endif
