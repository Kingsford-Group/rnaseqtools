#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "item.h"
#include "gtfcuff.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		cout<<"usage: " <<endl;
		cout<<"       " <<argv[0] << " roc <cuff.tmap> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " acc <cuff.tmap> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " acc-single <cuff.tmap>"<<endl;
		cout<<"       " <<argv[0] << " roc-trunc <cuff.tmap> <ref-size> <min-coverage> <max-coverage>"<<endl;
		cout<<"       " <<argv[0] << " roc-quant <cuff.tmap> <quant-file> <min-tpm> <max-tpm>"<<endl;
		cout<<"       " <<argv[0] << " acc-quant <cuff.tmap> <quant-file> <tpm-threshold>"<<endl;
		cout<<"       " <<argv[0] << " classify <cuff.tmap> <pred-gtf-file>"<<endl;
		cout<<"       " <<argv[0] << " quant <cuff.tmap> <pred-gtf-file> <ref-gtf-file>"<<endl;
		cout<<"       " <<argv[0] << " split <cuff.tmap> <pred-gtf-file> <true-file> <false-file>"<<endl;
		cout<<"       " <<argv[0] << " union <cuff.tmap> <pred-gtf-file> <ref-gtf-file> <merged-file>"<<endl;
		cout<<"       " <<argv[0] << " puniq <cuff.tmap> <pred-gtf-file> <ref-gtf-file> <puniq-file>"<<endl;
		cout<<"       " <<argv[0] << " match-precision <cuff.tmap> <ref-size> <balanced-precision>"<<endl;
		cout<<"       " <<argv[0] << " match-sensitivity <cuff.tmap> <ref-size> <balanced-sensitivity>"<<endl;
		cout<<"       " <<argv[0] << " match-correct <cuff.tmap> <ref-size> <balanced-correct>"<<endl;
		return 0;
	}

	if(string(argv[1]) == "roc")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.roc(atoi(argv[3]));
	}

	if(string(argv[1]) == "acc")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.acc(atoi(argv[3]));
	}

	if(string(argv[1]) == "acc-single")
	{
		assert(argc == 3);
		gtfcuff cuff(argv[2]);
		cuff.compute_single_accuracy();
	}

	if(string(argv[1]) == "match-precision")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.match_precision(atoi(argv[3]), atof(argv[4]));
	}

	if(string(argv[1]) == "match-sensitivity")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.match_sensitivity(atoi(argv[3]), atof(argv[4]));
	}

	if(string(argv[1]) == "match-correct")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.match_correct(atoi(argv[3]), atoi(argv[4]));
	}

	if(string(argv[1]) == "roc-trunc")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.roc_trunc(atoi(argv[3]), atof(argv[4]), atof(argv[5]));
	}

	if(string(argv[1]) == "roc-quant")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.roc_quant(argv[3], atof(argv[4]), atof(argv[5]));
	}

	if(string(argv[1]) == "acc-quant")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.acc_quant(argv[3], atof(argv[4]));
	}

	if(string(argv[1]) == "split")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.split(argv[4], argv[5]);
	}

	if(string(argv[1]) == "union")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.build_union(argv[5]);
	}

	if(string(argv[1]) == "puniq")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.build_pred_unique(argv[5]);
	}

	if(string(argv[1]) == "classify")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.classify();
	}

	if(string(argv[1]) == "quant")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.assign_ref(argv[4]);
		cuff.quant();
	}

    return 0;
}
