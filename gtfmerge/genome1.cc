#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>

genome1::genome1()
{}

genome1::genome1(const string &file)
{
	build(file);
}

int genome1::build(const string &file)
{
	clear();
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			transcript t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			if(merge_coverage_as_counts == true) t.coverage = 1.0;
			add_transcript(t);
		}
	}
	return 0;
}

int genome1::build(const vector<transcript> &v)
{
	clear();
	for(int i = 0; i < v.size(); i++)
	{
		add_transcript(v[i]);
	}
	return 0;
}

int genome1::clear()
{
	transcripts.clear();
	intron_hashing.clear();
	return 0;
}

int genome1::add_transcript(const transcript &t)
{
	string s = compute_intron_hashing(t);
	if(intron_hashing.find(s) == intron_hashing.end())
	{
		intron_hashing.insert(PSI(s, transcripts.size()));
		transcripts.push_back(t);
	}
	else
	{
		int k = intron_hashing[s];
		assert(k >= 0 && k < transcripts.size());
		transcript tt = transcripts[k];
		if(transcripts[k].length() > t.length()) 
		{
			transcripts[k].coverage += t.coverage;
		}
		else
		{
			double c = transcripts[k].coverage + t.coverage;
			transcripts[k] = t;
			transcripts[k].coverage = c;
		}
	}
	return 0;
}

int genome1::build_intersection(const genome1 &gm, genome1 &out)
{
	out.clear();
	for(MSI::iterator it = intron_hashing.begin(); it != intron_hashing.end(); it++)
	{
		string s = it->first;
		int k1 = it->second;
		transcript t = transcripts[k1];
		MSI::const_iterator x = gm.intron_hashing.find(s);
		if(x == gm.intron_hashing.end()) continue;
		int k2 = x->second;
		t.coverage += gm.transcripts[k2].coverage;
		out.add_transcript(t);
	}
	return 0;
}

int genome1::build_union(const genome1 &gm)
{
	for(int k = 0; k < gm.transcripts.size(); k++)
	{
		const transcript &t = gm.transcripts[k];
		add_transcript(t);
	}
	return 0;
}

int genome1::add_suffix(const string &p)
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		t.transcript_id.append("-").append(p);
		t.gene_id.append("-").append(p);
	}
	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_hashing.size());
	return 0;
}

int genome1::print_hashing()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		string s = compute_intron_hashing(transcripts[i]);
		printf("hash = %s\n", s.c_str());
	}
	return 0;
}

int genome1::write(const string &file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return 0;
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		if(t.coverage < min_transcript_coverage) continue;
		t.write(fout);
	}
	fout.close();
}

string compute_intron_hashing(const transcript &t)
{
	string h = t.seqname;
	
	if(t.strand == '.') h.append("0");
	if(t.strand == '+') h.append("1");
	if(t.strand == '-') h.append("2");

	if(t.exons.size() <= 1) return h;
	int32_t p = t.exons[0].second;
	h.append(to_string(p));

	for(int k = 1; k < t.exons.size(); k++)
	{
		int32_t q1 = t.exons[k].first;
		int32_t q2 = t.exons[k].second;
		h.append(to_string(q1 - p));
		if(k == t.exons.size() - 1) break;
		h.append(to_string(q2 - q1));
		p = q2;
	}
	return h;
}
