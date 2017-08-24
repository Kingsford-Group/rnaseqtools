#ifndef __GTFMERGE_H__
#define __GTFMERGE_H__

#include "genome1.h"

using namespace std;

//int do_union(const string &file, genome1 &gm);
int do_union(const vector<string> &v, genome1 *gm);

class gtfmerge
{
public:
	genome1 gm;

public:
	int build_union(const string &file);
};

#endif
