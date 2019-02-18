#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"
#include "genome1.h"

using namespace std;

int main(int argc, const char **argv)
{
 	if(argc == 1)
	{
		cout<<"usage: " << endl;
		cout<<"       " << argv[0] << " RPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " FPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " shrink <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " filter <min-transcript-coverage> <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " top <integer-n> <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " stats-exons <in-gtf-file> <exons-bins>"<<endl;
		cout<<"       " << argv[0] << " stats-length <in-gtf-file>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "FPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_FPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "RPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_RPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "shrink")
	{
		genome1 gm;
		gm.shrink(argv[2]);
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "filter")
	{
		genome1 gm;
		gm.filter(argv[3], atof(argv[2]));
		gm.write(argv[4]);
	}

	if(string(argv[1]) == "top")
	{
		genome1 gm;
		gm.top(argv[3], atoi(argv[2]));
		gm.write(argv[4]);
	}

	if(string(argv[1]) == "stats-exons")
	{
		genome1 gm;
		gm.stats_exons(argv[2], atoi(argv[3]));
	}

	if(string(argv[1]) == "stats-length")
	{
		genome1 gm;
		gm.stats_length(argv[2]);
	}

    return 0;
}
