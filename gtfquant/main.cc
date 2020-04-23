#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "genome.h"
#include "gtfquant.h"

using namespace std;

void print_usage(const char* name) {
	cout << "usage: "<< name << " filter|compare <salmon.quan> <gtffile> [min-tpm] [min-numreads] [TPM/RPKM]"<<endl;
}

int main(int argc, const char **argv)
{
	if(argc != 7 && argc != 6)
	{
		print_usage(argv[0]);
		return 1;
	}

	if(string(argv[1]) == "filter") {
		if(argc != 6) {
			cout << "Error: expected 4 argument to subcommand 'filter'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfquant roc(argv[2], argv[3], atof(argv[4]), atof(argv[5]));
		roc.filter();
		return 0;
        }
        else if(string(argv[1]) == "compare") {
		if(argc != 7) {
			cout << "Error: expected 5 argument to subcommand 'compare'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfquant roc(argv[2], argv[3], atof(argv[4]), atof(argv[5]));
		if(string(argv[6]) == "RPKM") roc.gm.assign_TPM_by_RPKM();
		roc.compare();
	} else {
		std::cout << "Error: unknown subcommand '" << argv[1] << "'\n";
		print_usage(argv[0]);
		return 1;
	}

    return 0;
}
