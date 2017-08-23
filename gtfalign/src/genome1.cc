#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>

genome1::genome1()
{}

int genome1::build(const string &file)
{
	clear();
	genome gm(file);
	for(int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for(int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			if(t.exons.size() <= 1) continue;
			add_transcript(t);
		}
	}
	return 0;
}

int genome1::clear()
{
	transcripts.clear();
	chain_hashing.clear();
	return 0;
}

int genome1::compare_boundary(const genome1 &gm)
{
	for(MSID::iterator it = mboundary.begin(); it != mboundary.end(); it++)
	{
		string chrm = it->first;
		MSID::const_iterator it2 = gm.mboundary.find(chrm);

		if(it2 == gm.mboundary.end()) continue;

		compare_boundary(chrm, it->second, it2->second);
	}
	return 0;
}

int genome1::compare_boundary(string chrm, const MID32 &mx, const MID32 &my)
{
	return 0;
}

int genome1::add_transcript(const transcript &t)
{
	if(t.exons.size() <= 1) return 0;

	// add boundaries and junctions
	vector<PI32> chain = t.get_intron_chain();
	for(int k = 0; k < chain.size(); k++)
	{
		add_boundary(chain[k].first, t.seqname, t.coverage);
		add_boundary(chain[k].second, t.seqname, t.coverage);
		add_junction(chain[k], t.seqname, t.coverage);
	}
	return 0;

	// add transcript and intron hashing
	string s = compute_chain_hashing(t);
	if(chain_hashing.find(s) == chain_hashing.end())
	{
		chain_hashing.insert(PSI(s, transcripts.size()));
		transcripts.push_back(t);
	}
	else
	{
		int k = chain_hashing[s];
		assert(k >= 0 && k < transcripts.size());
		transcripts[k].coverage += t.coverage;
	}
	return 0;
}

int genome1::add_boundary(int32_t p, string chrm, double coverage)
{
	if(mboundary.find(chrm) == mboundary.end())
	{
		MID32 m;
		m.insert(PID32(p, coverage));
		mboundary.insert(pair<string, MID32>(chrm, m));
	}
	else
	{
		MID32 &m = mboundary[chrm];
		if(m.find(p) == m.end()) m.insert(PID32(p, coverage));
		else m[p] += coverage;
	}
	return 0;
}

int genome1::add_junction(PI32 p, string chrm, double coverage)
{
	if(mjunction.find(chrm) == mjunction.end())
	{
		MPD32 m;
		m.insert(PPD32(p, coverage));
		mjunction.insert(pair<string, MPD32>(chrm, m));
	}
	else
	{
		MPD32 &m = mjunction[chrm];
		if(m.find(p) == m.end()) m.insert(PPD32(p, coverage));
		else m[p] += coverage;
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
	printf("genome %d: %lu transcripts\n", index, transcripts.size());
	return 0;
}

int genome1::print_hashing()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		string s = compute_chain_hashing(transcripts[i]);
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

string compute_chain_hashing(const transcript &t)
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
