#ifndef __GENOME1_H__
#define __GENOME1_H__

#include "genome.h"

using namespace std;

typedef pair<PI32, set<int> > PPIS;
typedef map<PI32, set<int> > MPIS;
typedef pair<int, int> PII;
typedef map<int, int> MII;

class genome1
{
public:
	genome1();

public:
	vector<transcript> transcripts;

private:
	MPIS intron_index;	

public:
	int build(const string &file);
	int clear();
	int print(int index);
	int write(const string &file);
	int compare(const genome1 &gy);
	int shrink(const string &file);
	int filter(const string &file, double c);
	int select_by_length(const string &file, int min_length, int max_length);
	int top(const string &file, int n);
	int stats_exons(const string &file, int n);
	int stats_length(const string &file);
    int write_all(const string &input, const string &output, set<string> expressedTr);

private:
	int build_multiexon_transcripts(const string &file);
	int build_all_transcripts(const string &file);
	int build_intron_index();
	int query(const transcript &t, const set<int> &fb);
	int query(const transcript &t, int min_index);
};

bool transcript_cmp_coverage(const transcript &x, const transcript &y);
bool transcript_cmp_intron_chain(const transcript &x, const transcript &y);
string compute_intron_hashing(const transcript &t);
string tostring(int32_t p);

#endif
