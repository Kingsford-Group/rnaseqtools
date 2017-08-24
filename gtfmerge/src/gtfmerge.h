#ifndef __GTFMERGE_H__
#define __GTFMERGE_H__

#include "genome1.h"

using namespace std;

class gtfmerge
{
public:
	genome1 gm;

public:
	int build_union(const string &file);
};

#endif
