#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome1.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc < 3)
	{
		cout<<"usage: "<<argv[0]<< " <gtf-file-1> <gtf-file-2> [-l length -a algo -s]"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	genome1 gm1(argv[1]);
	genome1 gm2(argv[2]);

	//gm1.remove_redundancy();
	//gm2.remove_redundancy();
	gm1.compare(gm2);

    return 0;
}
