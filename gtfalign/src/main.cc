#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "genome1.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		cout<<"usage: "<<argv[0]<< " intron <gtf-file1> <gtf-file2> [options]"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "intron")
	{
		genome1 g1;
		genome1 g2;
	}

    return 0;
}
