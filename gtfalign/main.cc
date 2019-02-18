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
		cout<<"usage: "<<argv[0]<< " boundary <gtf-file1> <gtf-file2> [options]"<<endl;
		cout<<"usage: "<<argv[0]<< " junction <gtf-file1> <gtf-file2> [options]"<<endl;
		cout<<"usage: "<<argv[0]<< " chain <gtf-file1> <gtf-file2> <novel-chain.gtf>"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "boundary")
	{
		genome1 g1;
		genome1 g2;

		g1.build_chain_hashing(argv[2]);
		g2.build_chain_hashing(argv[3]);

		g1.compare_boundary(g2);
	}

	if(string(argv[1]) == "junction")
	{
		genome1 g1;
		genome1 g2;

		g1.build_chain_hashing(argv[2]);
		g2.build_chain_hashing(argv[3]);

		g1.compare_junction(g2);
	}

	if(string(argv[1]) == "chain")
	{
		genome1 g1;
		genome1 g2;

		g1.build_boundary_hashing(argv[2]);
		g2.build_boundary_hashing(argv[3]);

		g1.compare_chain(g2, argv[4]);
	}

    return 0;
}
