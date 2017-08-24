#include "gtfmerge.h"
#include <string>

int gtfmerge::build_union(const string &file)
{
	gm.clear();
	ifstream fin(file.c_str());
	if(fin.fail()) return 0;

	string line;
	int cnt = 0;
	while(getline(fin, line))
	{
		if(line == "") continue;
		genome1 g(line);
		gm.build_union(g);
		cout << "union genome " << ++cnt << " " << line.c_str() << endl << flush;
	}

	fin.close();
	return 0;
}
