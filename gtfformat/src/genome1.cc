#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>
#include <algorithm>

genome1::genome1()
{}

int genome1::clear()
{
	transcripts.clear();
	intron_index.clear();
	return 0;
}

int genome1::build_all_transcripts(const string &file)
{
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			transcripts.push_back(t);
		}
	}
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

int genome1::shrink(const string &file)
{
	map<int, int> mm;
	vector<transcript> vv;
	build_multiexon_transcripts(file);
	int n = transcripts.size();
	if(n <= 0) return 0;

	sort(transcripts.begin(), transcripts.end(), transcript_cmp_intron_chain);

	vv.push_back(transcripts[0]);

	int c = 1;
	for(int k = 1; k < n; k++)
	{
		bool b = transcripts[k].intron_chain_match(transcripts[k - 1]);
		if(b == true) c++;
		if(b == true) continue;

		vv.push_back(transcripts[k]);
		if(mm.find(c) == mm.end()) mm.insert(PII(c, 1));
		else mm[c]++;
		c = 1;
	}

	for(MII::iterator it = mm.begin(); it != mm.end(); it++)
	{
		printf("size %d : %d groups (by identifical intron-chain)\n", it->first, it->second);
	}

	transcripts = vv;
	return 0;
}

int genome1::filter(const string &file, double c)
{
	build_multiexon_transcripts(file);
	vector<transcript> vv;
	for(int i = 0; i < transcripts.size(); i++)
	{
		if(transcripts[i].coverage < c) continue;
		vv.push_back(transcripts[i]);
	}
	transcripts = vv;
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

int genome1::stats(const string &file, int n)
{
	build_all_transcripts(file);
	vector<int> counts;
	vector<int> length;
	counts.assign(n, 0);
	length.assign(n, 0);
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		int k = t.exons.size();
		int l = t.length();
		if(k >= n) k = n;
		counts[k - 1] ++;
		length[k - 1] += l;
	}

	for(int k = 0; k < n; k++)
	{
		printf("transcripts with %d exons: %d\n", k + 1, counts[k]);
	}

	/*
	for(int k = 0; k < n; k++)
	{
		double ave = -1;
		if(counts[k] >= 1) ave = length[k] * 1.0 / counts[k];
	}

	for(int k = 0; k < n; k++)
	{
		if(counts[k] == 0) printf("0.0 ");
		else printf("%.2lf ", length[k] * 1.0 / counts[k]);
	}
	*/

	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_index.size());
	return 0;
}

int genome1::write(const string &file)
{
	ofstream fout(file.c_str());
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].write(fout);
	}
	fout.close();
	return 0;
}

bool transcript_cmp_coverage(const transcript &x, const transcript &y)
{
	if(x.coverage < y.coverage) return true;
	else return false;
}

bool transcript_cmp_intron_chain(const transcript &x, const transcript &y)
{
	if(x.strand < y.strand) return true;
	if(x.strand > y.strand) return false;
	if(x.seqname < y.seqname) return true;
	if(x.seqname > y.seqname) return false;
	if(x.exons.size() < y.exons.size()) return true;
	if(x.exons.size() > y.exons.size()) return false;

	int n = x.exons.size() - 1;
	if(n >= 0 && x.exons[0].second < y.exons[0].second) return true;
	if(n >= 0 && x.exons[0].second > y.exons[0].second) return false;
	if(n >= 0 && x.exons[n].first < y.exons[n].first) return true;
	if(n >= 0 && x.exons[n].first > y.exons[n].first) return false;

	for(int k = 1; k < n - 1; k++)
	{
		int32_t x1 = x.exons[k].first;
		int32_t y1 = y.exons[k].first;
		int32_t x2 = x.exons[k].second;
		int32_t y2 = y.exons[k].second;
		if(x1 < y1) return true;
		if(x1 > y1) return false;
		if(x2 < y2) return true;
		if(x2 > y2) return false;
	}
	if(x.length() > y.length()) return true;
	if(x.length() < y.length()) return false;
	if(x.get_bounds().first < y.get_bounds().first) return true;
	if(x.get_bounds().first > y.get_bounds().first) return false;
	if(x.transcript_id.compare(y.transcript_id) < 0) return true;
	else return false;
}
