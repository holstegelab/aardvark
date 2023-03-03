#include "arktools.hpp"
#include "mk.hpp"
#include "kseq.h"
#include <zlib.h>
#include <vector>
#include <string>
#include <iostream>

KSEQ_INIT(gzFile, gzread);


void mk(const std::vector<std::string>& fastas, const std::string& mask_fasta)
{
	KS_FULL_COMMENT = true; 
	gzFile fp;
	kseq_t* seq;
	init_encoding();

	Trie mask_trie;

	int seq_l;
	fp = gzopen(mask_fasta.c_str(), "r");
	seq = kseq_init(fp);
	while ((seq_l = kseq_read(seq)) > 0) {
		std::string name = seq->name.s;
		std::string local_sequence = seq->seq.s;
		mask_trie.set_seq(local_sequence);
	}
	
	kseq_destroy(seq);
	gzclose(fp);

	//std::vector<std::string>masks;
	//mask_trie.collect_seqs(masks);

	for(const auto& fasta : fastas){
		fp = gzopen(fasta.c_str(), "r");
		seq = kseq_init(fp);
		while ((seq_l = kseq_read(seq)) > 0) {
			std::string name = seq->name.s;
			std::string local_sequence = seq->seq.s;
			for(int i = 0; i < (int)local_sequence.size();){
				int max_i = mask_trie.contains(local_sequence, local_sequence.size(), i);
				if(max_i <= i) ++i;
				else{
					for(int j = i; j < max_i; ++j) local_sequence[j] = 'N';
					i = max_i;
				}
			}
			std::cout << '>' << name << '\n' << local_sequence << '\n';
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

}