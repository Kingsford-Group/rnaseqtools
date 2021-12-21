#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "item.h"
#include "gtfcuff.h"

using namespace std;

void print_usage(const char* name) {
	cout <<"usage: \n"
	     << "       " << name << " roc <cuff.tmap> <ref-size> <cov|TPM|FPKM>\n"
	     << "       " << name << " roc-salmon <cuff.tmap> <ref-size> <salmon-quant-file>\n"
	     << "       " << name << " auc <cuff.tmap> <ref-size>\n"
	     << "       " << name << " acc <cuff.tmap> <ref-size>\n"
	     << "       " << name << " acc-single <cuff.tmap>\n"
	     << "       " << name << " roc-trunc <cuff.tmap> <ref-size> <min-coverage> <max-coverage>\n"
	     << "       " << name << " roc-quant <cuff.tmap> <quant-file> <min-tpm> <max-tpm>\n"
	     << "       " << name << " acc-quant <cuff.tmap> <quant-file> <tpm-threshold>\n"
	     << "       " << name << " classify <cuff.tmap> <pred-gtf-file>\n"
	     << "       " << name << " quant <cuff.tmap> <pred-gtf-file> <ref-gtf-file>\n"
	     << "       " << name << " split <cuff.tmap> <pred-gtf-file> <true-file> <false-file>\n"
	     << "       " << name << " union <cuff.tmap> <pred-gtf-file> <ref-gtf-file> <merged-file>\n"
	     << "       " << name << " puniq <cuff.tmap> <pred-gtf-file> <ref-gtf-file> <puniq-file>\n"
	     << "       " << name << " match-precision <cuff.tmap> <ref-size> <balanced-precision>\n"
	     << "       " << name << " match-sensitivity <cuff.tmap> <ref-size> <balanced-sensitivity>\n"
	     << "       " << name << " match-correct <cuff.tmap> <ref-size> <balanced-correct>\n"
	     << "       " << name << " split-class <cuff.tmap> <gtf-file> <prefix>"<<std::endl;
}

int main(int argc, const char **argv)
{
	if(argc == 1) {
		print_usage(argv[0]);
		return 1;
	}
	else if(string(argv[1]) == "auc") {
		if(argc != 4) {
			cout << "Error: expected 2 arguments for subcommand 'auc'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.auc(atoi(argv[3]));
	}
	else if(string(argv[1]) == "roc") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'roc'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2], argv[4]);
		cuff.roc(atoi(argv[3]));
	}
	else if(string(argv[1]) == "roc_salmon") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'roc_salmon'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.roc_salmon(atoi(argv[3]), argv[4]);
	}
	else if(string(argv[1]) == "acc") {
		if(argc != 4) {
			cout << "Error: expected 2 arguments for subcommand 'acc'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.acc(atoi(argv[3]));
	}
	else if(string(argv[1]) == "acc-single") {
		if(argc != 3) {
			cout << "Error: expected 1 arguments for subcommand 'acc-single'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.compute_single_accuracy();
	}
	else if(string(argv[1]) == "match-precision") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'match-precision'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.match_precision(atoi(argv[3]), atof(argv[4]));
	}
	else if(string(argv[1]) == "match-sensitivity") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'match-sensitivity'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.match_sensitivity(atoi(argv[3]), atof(argv[4]));
	}
	else if(string(argv[1]) == "match-correct") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'match-correct'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.match_correct(atoi(argv[3]), atoi(argv[4]));
	}
	else if(string(argv[1]) == "roc-trunc") {
		if(argc != 6) {
			cout << "Error: expected 4 arguments for subcommand 'roc-trunc'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.roc_trunc(atoi(argv[3]), atof(argv[4]), atof(argv[5]), false);
	}
	else if(string(argv[1]) == "roc-quant") {
		if(argc != 6) {
			cout << "Error: expected 4 arguments for subcommand 'roc-quant'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.roc_quant(argv[3], atof(argv[4]), atof(argv[5]));
	}
	else if(string(argv[1]) == "acc-quant") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'acc-quant'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.acc_quant(argv[3], atof(argv[4]));
	}
	else if(string(argv[1]) == "split") {
		if(argc != 6) {
			cout << "Error: expected 4 arguments for subcommand 'split'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.split(argv[4], argv[5]);
	}
	else if(string(argv[1]) == "union") {
		if(argc != 6) {
			cout << "Error: expected 4 arguments for subcommand 'union'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.build_union(argv[5]);
	}
	else if(string(argv[1]) == "puniq") {
		if(argc != 6) {
			cout << "Error: expected 4 arguments for subcommand 'puniq'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.build_pred_unique(argv[5]);
	}
	else if(string(argv[1]) == "classify") {
		if(argc != 4) {
			cout << "Error: expected 2 arguments for subcommand 'classify'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.classify();
	}
	else if(string(argv[1]) == "quant") {
		if(argc != 5) {
			cout << "Error: expected 3 arguments for subcommand 'quant'\n";
			print_usage(argv[0]);
			return 1;
		}
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.quant();
	}
        else if(string(argv[1]) == "split-class") {
                if(argc != 5) {
                        cout << "Error: expected 3 arguments for subcommand 'split-class'\n";
                        print_usage(argv[0]);
                        return 1;
                }
                gtfcuff cuff(argv[2]);
                cuff.assign_ref(argv[3]);
                cuff.split_class(argv[4]);
        }
	else {
		cout << "Error: unknown subcommand '" << argv[1] << "'\n";
		print_usage(argv[0]);
		return 1;
	}

	return 0;
}
