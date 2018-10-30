#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<string, int> PSI;
typedef map<string, int> MSI;

class genome1
{
public:
	genome1();
	genome1(const string &file);

public:
	vector<transcript> transcripts;
	MSI intron_hashing;

public:
	int add_transcript(const transcript &t);
	int build(const string &file);
	int build(const vector<transcript> &v);
	int build_union(const genome1 &gm);
	int build_intersection(const genome1 &gm, genome1 &out);
	int clear();
	int write(const string &file);
	int add_suffix(const string &p);
	int remove_redundancy();
	int print(int index);
	int print_hashing();

private:
	int build_multiexon_transcripts(const string &file);
};

string tostring(int p);
string compute_intron_hashing(const transcript &t);

#endif
