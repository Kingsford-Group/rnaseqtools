#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;
typedef pair<int32_t, double> PID32;
typedef pair<PI32, double> PPD32;
typedef map<int32_t, double> MID32;
typedef map<PI32, double> MPD32;

typedef pair<string, int> PSI;
typedef map<string, int> MSI;
typedef map<string, MID32> MSID;
typedef map<string, MPD32> MSPD;

class genome1
{
public:
	genome1();

public:
	vector<transcript> transcripts;
	MSI chain_hashing;
	MSID mboundary;
	MSPD mjunction;

public:
	int compare_boundary(const genome1 &gm);
	int compare_junction(const genome1 &gm);

public:
	int add_transcript(const transcript &t);
	int build(const string &file);
	int clear();
	int write(const string &file);
	int add_suffix(const string &p);
	int print(int index);
	int print_hashing();

private:
	int add_boundary(int32_t p, string chrm, double coverage);
	int add_junction(PI32 p, string chrm, double coverage);
	int compare_boundary(string, const MID32&, const MID32&);
	int compare_junction(string, const MPD32&, const MPD32&);
};

string tostring(int32_t p);
string compute_chain_hashing(const transcript &t);

#endif
