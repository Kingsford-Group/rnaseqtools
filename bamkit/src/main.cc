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

	if(argc != 2)
	{
		printf("usage: %s <bam-file>\n", argv[0]);
		return 0;
	}

	bamkit bk(argv[1]);
	bk.solve();

	return 0;
}
