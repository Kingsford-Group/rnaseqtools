#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>
#include <algorithm>

genome1::genome1()
{}

genome1::genome1(const string &file)
{
	build(file);
}

int genome1::build(const string &file)
{
	build_multiexon_transcripts(file);
	sort(transcripts.begin(), transcripts.end(), transcript_cmp);
	build_intron_index();
	return 0;
}

int genome1::clear()
{
	transcripts.clear();
	intron_index.clear();
	return 0;
}

int genome1::build_multiexon_transcripts(const string &file)
{
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			transcripts.push_back(t);
		}
	}
	return 0;
}

int genome1::build_intron_index()
{
	intron_index.clear();
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		assert(t.exons.size() >= 2);
		PI32 p = t.get_first_intron();
		if(intron_index.find(p) == intron_index.end()) 
		{
			set<int> s;
			s.insert(i);
			intron_index.insert(PPIS(p, s));
		}
		else
		{
			intron_index[p].insert(i);
		}
	}
	return 0;
}

int genome1::query(const transcript &t, const set<int> &fb)
{
	if(t.exons.size() <= 1) return -1;
	PI32 p = t.get_first_intron();
	if(intron_index.find(p) == intron_index.end()) return -1;
	set<int> s = intron_index[p];
	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	for(int i = v.size() - 1; i >= 0; i--)
	{
		int k = v[i];
		transcript &x = transcripts[k];
		if(x.strand != t.strand) continue;
		if(x.seqname != t.seqname) continue;
		if(x.exons.size() != t.exons.size()) continue;
		if(x.intron_chain_match(t) == false) continue;
		if(fb.find(k) != fb.end()) continue;
		return k;
	}
	return -1;
}

int genome1::query(const transcript &t, int min_index)
{
	if(t.exons.size() <= 1) return -1;
	PI32 p = t.get_first_intron();
	if(intron_index.find(p) == intron_index.end()) return -1;
	set<int> s = intron_index[p];
	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	for(int i = v.size() - 1; i >= 0; i--)
	{
		int k = v[i];
		if(k < min_index) continue;
		transcript &x = transcripts[k];
		if(x.strand != t.strand) continue;
		if(x.seqname != t.seqname) continue;
		if(x.exons.size() != t.exons.size()) continue;
		if(x.intron_chain_match(t) == false) continue;
		return k;
	}
	return -1;
}

int genome1::remove_redundancy()
{
	set<int> rd;
	for(int i = transcripts.size() - 1; i >= 0; i--)
	{
		transcript &t = transcripts[i];	
		int k = query(t, i + 1);
		if(k == -1) continue;
		rd.insert(i);
	}

	vector<transcript> v;
	for(int i = 1; i < transcripts.size(); i++)
	{
		if(rd.find(i) != rd.end()) continue;
		transcript &t = transcripts[i];	
		v.push_back(t);
	}
	transcripts = v;
	build_intron_index();
	return 0;
}

int genome1::compare(const genome1 &gy)
{
	MII x2y;
	MII y2x;
	set<int> fb;
	for(int i = gy.transcripts.size() - 1; i >= 0; i--)
	{
		const transcript &t = gy.transcripts[i];
		int k = query(t, fb);
		if(k == -1) continue;
		x2y.insert(PII(k, i));
		y2x.insert(PII(i, k));
		fb.insert(k);
	}

	int correct = x2y.size();
	int refsize = transcripts.size();
	int prdsize = gy.transcripts.size();

	double sen0 = correct * 100.0 / refsize;
	for(int i = 0; i < gy.transcripts.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (prdsize - i);

		if(sen * 2.0 < sen0) break;

		if(i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf\n",
				refsize, prdsize - i, correct, sen, pre, gy.transcripts[i].coverage);
		}

		if(y2x.find(i) != y2x.end()) correct--;
	}
	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_index.size());
	return 0;
}

bool transcript_cmp(const transcript &x, const transcript &y)
{
	if(x.coverage < y.coverage) return true;
	else return false;
}
