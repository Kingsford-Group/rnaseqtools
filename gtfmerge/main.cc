#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "gtfmerge.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		print_help();
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "union")
	{
		assert(argc >= 4);
		gtfmerge gm;
		gm.build_union(argv[2]);
		gm.gm.write(argv[3]);
	}

    return 0;
}
