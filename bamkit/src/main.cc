/*
Part of bamkit 
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "bamkit.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc != 3)
	{
		printf("usage: \n");
		printf(" %s count <bam-file>\n", argv[0]);
		printf(" %s strand <bam-file>\n", argv[0]);
		return 0;
	}

	if(string(argv[1]) == "count")
	{
		bamkit bk(argv[1]);
		bk.solve_count();
	}

	if(string(argv[1]) == "strand")
	{
		bamkit bk(argv[2]);
		bk.solve_strand();
	}

	return 0;
}
