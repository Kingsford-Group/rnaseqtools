#include "gtfmerge.h"
#include "config.h"

#include <string>
#include <thread>

int do_union(const vector<string> &v, genome1 *gm)
{
	gm->clear();
	int cnt = 0;	
	for(int k = 0; k < v.size(); k++)
	{
		string s = v[k];
		genome1 g(s);
		gm->build_union(g);
		cout << "union genome " << ++cnt << " " << s.c_str() << endl << flush;
	}
	return 0;
}

int gtfmerge::build_union(const string &file)
{
	vector<genome1> gv(num_threads);
	vector< vector<string> > sv(num_threads);

	ifstream fin(file.c_str());
	if(fin.fail()) return 0;
	string line;
	int cnt = 0;
	while(getline(fin, line))
	{
		if(line == "") continue;
		int k = cnt % num_threads;
		sv[k].push_back(line);
		cnt++;
	}
	fin.close();

	vector<thread> threads;
	for(int k = 0; k < sv.size(); k++)
	{
		printf("thread %d processes %lu genomes\n", k, sv[k].size());
		threads.emplace_back(do_union, sv[k], &(gv[k]));
	}
	for(int k = 0; k < threads.size(); k++)
	{
		threads[k].join();
	}

	gm.clear();
	for(int k = 0; k < gv.size(); k++)
	{
		printf("final union genome %d\n", k);
		gm.build_union(gv[k]);
	}

	return 0;
}
