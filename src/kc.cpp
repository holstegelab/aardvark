#include "arktools.hpp"
#include "kc.hpp"
#include "kseq.h"
#include <zlib.h>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <utility>


KSEQ_INIT(gzFile, gzread);

class CKmer {
	public:
		uint64_t o;
		uint64_t m;
		uint32_t n;
		CKmer(uint64_t _o, uint64_t _m, uint64_t _n): o(_o), m(_m), n(_n){};
};


void count_pairs(const std::vector<std::pair<uint64_t, uint64_t>>& values, std::vector<CKmer>& output)
{
	auto it = values.begin();
	uint64_t last_o = it->first;
	uint64_t last_g = it->second;
	uint32_t count = 1;
	++it;
	while(it != values.end()){
		if(it->first == last_o && it->second == last_g) ++count;
		else{
			output.emplace_back(last_o, last_g, count);
			count = 1;
			last_o = it->first;
			last_g = it->second;
		}
		++it;
	}
	output.emplace_back(last_o, last_g, count);	
}

void count_cks(std::vector<CKmer>& values)
{
	auto it = values.begin();
	auto last = it;
	++it;
	while(it != values.end()){
		if(it->m == last->m) {
			if(it->n > last->n) last->o = it->o;
			last->n += it->n;
			it = values.erase(it);
		}
		else{
			last = it;
			++it;	
		}
	}
}

void kc(const std::vector<std::string>& fastas, const uint32_t k, const uint32_t m)
{
	init_encoding();
	KS_FULL_COMMENT = true; 
	uint64_t seed = 1;
	int seq_l;

	for(uint32_t k_l = k; k_l < m; ++k_l){
		//std::vector<uint64_t> kmer_hashes;
		std::vector<std::pair<uint64_t, uint64_t>> kmer_indeces;

		for(const auto& fasta: fastas){
			gzFile fp = gzopen(fasta.c_str(), "r");
			kseq_t *seq = kseq_init(fp);
			while ((seq_l = kseq_read(seq)) >= (int)k_l ) {
				std::string name = seq->name.s;
				std::string local_sequence = seq->seq.s;
				get_kmer_min_indeces(local_sequence, k_l, kmer_indeces);
			}
			
			kseq_destroy(seq);
			gzclose(fp);
		}

		sort(kmer_indeces.begin(), kmer_indeces.end(),[](const auto& a, const auto& b) -> bool{return a.first < b.first;});

		std::vector<CKmer> ck_counts;
		count_pairs(kmer_indeces, ck_counts);

		sort(ck_counts.begin(), ck_counts.end(),[](const auto& a, const auto& b) -> bool{
			if(a.m == b.m) return a.o < b.o;
			else return a.m < b.m;
		});

		count_cks(ck_counts);
		for(const auto& ck : ck_counts) {
			std::string seq(k_l, 'N');
			index2kmer(seq, ck.o, k_l);
			std::cout << seq << '\t' << ck.n << '\n';
		}
	}
}