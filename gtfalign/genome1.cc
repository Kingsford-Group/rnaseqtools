#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>

genome1::genome1()
{}

int genome1::build_chain_hashing(const string &file)
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
			add_transcript1(t);
		}
	}
	return 0;
}

int genome1::build_boundary_hashing(const string &file)
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
			add_transcript2(t);
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

int genome1::compare_chain(genome1 &gm)
{
	MSPIM64 &m = gm.boundary_hashing;

	for(int k = 0; k < transcripts.size(); k++)
	{
		transcript &t = transcripts[k];
		vector<PI32> chain = t.get_intron_chain();
		assert(chain.size() >= 1);
		int32_t p1 = chain.front().first;
		int32_t p2 = chain.back().second;
		int64_t pp = pack(p1, p2);

		string type;
		string match = "N/A";

		if(m.find(t.seqname) == m.end())
		{
			type = "chrm_mismatch";
			continue;
		}

		PIM64 &pim = m[t.seqname];
		if(pim.find(pp) == pim.end())
		{
			type = "boundary_mismatch";
			continue;
		}

		type = "chain_mismatch";

		set<int> &s = pim[pp];
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			transcript &tt = gm.transcripts[*it];
			if(tt.intron_chain_match(t) == false) continue;

			match = tt.transcript_id;
			type = "strand_mismatch";
			
			if(tt.strand != t.strand) continue;
			type = "identical";

			break;
		}

		printf("%s %s %s\n", t.transcript_id.c_str(), type.c_str(), match.c_str());
	}
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
	for(MID32::const_iterator it = mx.begin(); it != mx.end(); it++)
	{
		int32_t p1 = it->first;
		double c1 = it->second;

		MID32::const_iterator it2 = my.find(p1);

		if(it2 == my.end())
		{
			printf("type1 %s:%d %.3lf %.3lf\n", chrm.c_str(), p1, c1, 0.0);
		}
		else
		{
			double c2 = it2->second;
			printf("type2 %s:%d %.3lf %.3lf\n", chrm.c_str(), p1, c1, c2);
		}
	}

	for(MID32::const_iterator it2 = my.begin(); it2 != my.end(); it2++)
	{
		int32_t p2 = it2->first;
		double c2 = it2->second;

		MID32::const_iterator it = mx.find(p2);

		if(it == mx.end())
		{
			printf("type3 %s:%d %.3lf %.3lf\n", chrm.c_str(), p2, 0.0, c2);
		}
	}

	return 0;
}

int genome1::compare_junction(const genome1 &gm)
{
	for(MSPD::iterator it = mjunction.begin(); it != mjunction.end(); it++)
	{
		string chrm = it->first;
		MSPD::const_iterator it2 = gm.mjunction.find(chrm);

		if(it2 == gm.mjunction.end()) continue;

		compare_junction(chrm, it->second, it2->second);
	}
	return 0;
}

int genome1::compare_junction(string chrm, const MPD32 &mx, const MPD32 &my)
{
	for(MPD32::const_iterator it = mx.begin(); it != mx.end(); it++)
	{
		PI32 p1 = it->first;
		double c1 = it->second;

		MPD32::const_iterator it2 = my.find(p1);

		if(it2 == my.end())
		{
			printf("type1 %s:%d-%d %.3lf %.3lf\n", chrm.c_str(), p1.first, p1.second, c1, 0.0);
		}
		else
		{
			double c2 = it2->second;
			printf("type2 %s:%d-%d %.3lf %.3lf\n", chrm.c_str(), p1.first, p1.second, c1, c2);
		}
	}

	for(MPD32::const_iterator it2 = my.begin(); it2 != my.end(); it2++)
	{
		PI32 p2 = it2->first;
		double c2 = it2->second;

		MPD32::const_iterator it = mx.find(p2);

		if(it == mx.end())
		{
			printf("type3 %s:%d-%d %.3lf %.3lf\n", chrm.c_str(), p2.first, p2.second, 0.0, c2);
		}
	}

	return 0;
}

int genome1::add_transcript1(const transcript &t)
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

int genome1::add_transcript2(const transcript &t)
{
	if(t.exons.size() <= 1) return 0;
	int index = transcripts.size();

	// add boundaries and junctions
	vector<PI32> chain = t.get_intron_chain();
	assert(chain.size() >= 1);

	int32_t p1 = chain.front().first;
	int32_t p2 = chain.back().second;
	int64_t pp = pack(p1, p2);

	string h = t.seqname;
	
	/*
	if(t.strand == '.') h.append("0");
	if(t.strand == '+') h.append("1");
	if(t.strand == '-') h.append("2");
	*/

	if(boundary_hashing.find(h) == boundary_hashing.end())
	{
		PIM64 pim;
		set<int> s;
		s.insert(index);
		pim.insert(PIS64(pp, s));
		boundary_hashing.insert(PSPIM64(h, pim));
	}
	else if(boundary_hashing[h].find(pp) == boundary_hashing[h].end())
	{
		set<int> s;
		s.insert(index);
		boundary_hashing[h].insert(PIS64(pp, s));
	}
	else
	{
		boundary_hashing[h][pp].insert(index);
	}

	transcripts.push_back(t);
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

string tostring(int32_t p)
{
	char buf[10240];
	sprintf(buf, "%d", p);
	return string(buf);
}

string compute_chain_hashing(const transcript &t)
{
	string h = t.seqname;
	
	if(t.strand == '.') h.append("0");
	if(t.strand == '+') h.append("1");
	if(t.strand == '-') h.append("2");

	if(t.exons.size() <= 1) return h;
	int32_t p = t.exons[0].second;
	h.append(tostring(p));

	for(int k = 1; k < t.exons.size(); k++)
	{
		int32_t q1 = t.exons[k].first;
		int32_t q2 = t.exons[k].second;
		h.append(tostring(q1 - p));
		if(k == t.exons.size() - 1) break;
		h.append(tostring(q2 - q1));
		p = q2;
	}
	return h;
}
