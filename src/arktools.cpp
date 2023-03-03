#include "arktools.hpp"
#include <string>
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <utility>

void init_encoding()
{
	for(int i = 0; i < 256; ++i) nt2encoding[i] = 5;
	nt2encoding['A'] = nt2encoding['a'] = 0; 
	nt2encoding['C'] = nt2encoding['c'] = 1;
	nt2encoding['G'] = nt2encoding['g'] = 2; 
	nt2encoding['T'] = nt2encoding['t'] = 3;
	nt2encoding['N'] = nt2encoding['n'] = 4;
	for(int i = 0; i < 4; ++i) encoding2nt[i] = 'N';
	encoding2nt[0] = 'A';
	encoding2nt[1] = 'C';
	encoding2nt[2] = 'G';
	encoding2nt[3] = 'T';
	for(int i = 0; i < 4; ++i) encoding2rcnt[i] = 'N';
	encoding2rcnt[0] = 'T';
	encoding2rcnt[1] = 'G';
	encoding2rcnt[2] = 'C';
	encoding2rcnt[3] = 'A';
}

bool is_ambi(const char* nt, const uint32_t& k)
{
	for(uint32_t i = 0; i < k; ++i) {
		uint8_t encoding = nt2encoding[(uint8_t)*(nt + i)];
		if(encoding > 3) return true;
	}
	return false;
}

uint64_t _kmer2index(const char* nt, uint32_t k){
  if(k == 0) return 0;
  uint64_t index = (uint64_t)nt2encoding[(uint8_t)(*nt)];
  return 4 * _kmer2index(nt - 1, k - 1) + index;
}

uint64_t kmer2index(const char* nt, uint32_t k){
  return _kmer2index(nt + k - 1, k);
}

void index2kmer(std::string& seq, uint64_t index, uint32_t k){
  if(k == 1) {
  	uint64_t nt_index = encoding2nt[(uint8_t)index];
  	seq[k-1] = nt_index;
  }
  else {
    uint64_t prefixindex = index / 4;
    uint8_t r = (uint8_t)(index % 4);
    uint64_t nt_index = encoding2nt[(uint8_t)r];
    seq[k-1] = nt_index;
    index2kmer(seq, prefixindex, k - 1);
  }
}

std::pair<uint64_t, uint64_t> min_permute_index(const char* nt, const uint32_t& k)
{
	std::string minseq(k, 'N');
	for(int i = 0; i < (int)k; ++i) minseq[i] = nt[i];

	uint64_t min_index = kmer2index(nt, k);
	uint64_t og_index = min_index;
	for(uint32_t i = 1; i < k; ++i){
		std::string seq(k, 'N');
		int seq_i = 0;
		for(uint32_t j = i; j < k; ++j){
			seq[seq_i] = *(nt + j);
			++seq_i;
			//seq.push_back(*(nt + j));
		}
		for(uint32_t j = 0; j < i; ++j){
			seq[seq_i] = *(nt + j);
			++seq_i;
			//seq.push_back(*(nt + j));
		}
		uint64_t index = kmer2index(seq.c_str(), k);
		if(index < min_index) {
			min_index = index;
			minseq = seq;
		}
	}
	return std::make_pair(og_index, min_index);
}

void get_kmer_min_indeces(const std::string& seq, const uint32_t& k, std::vector<std::pair<uint64_t,uint64_t>>& output)
{
	const char* nt = seq.c_str();
	for(uint32_t i = 0; i < seq.size() - k + 1; ++i) if(!is_ambi(nt + i, k)) output.emplace_back(min_permute_index(nt + i, k));
}

Trie::Trie():is_end(false){};

int Trie::contains(const std::string& seq, const int& l, const int i) const
{
	return Trie::_contains(seq, l, i, i);
}

int Trie::_contains(const std::string& seq, const int& l, const int max_i, const int i) const
{
	if(i >= l) if(is_end) return i; else return max_i;
	else{
		uint8_t index = nt2encoding[seq[i]];
		int updated_max_i = max_i;
		if(is_end) updated_max_i = i;
		if(index == 4) return updated_max_i;
		if(!nodes[index]) return updated_max_i;
		return nodes[index]->_contains(seq, l, updated_max_i, i + 1);
	}
}

void Trie::set_seq(const std::string& seq)
{
	int i = 0;
	int l = (int)seq.size();
	_set_seq(seq, l, i);
}

void Trie::_set_seq(const std::string& seq, const int& l, int i)
{
	if(i >= l) is_end = true;
	else {
		uint8_t index = nt2encoding[seq[i]];
		if(index == 4) exit(1);
		if(!nodes[index]) nodes[index].reset(new Trie());
		nodes[index]->_set_seq(seq, l, i + 1);
	}
}

void Trie::collect_seqs(std::vector<std::string>& output) const
{
	std::string acc;
	_collect_seqs(acc, output);
}

void Trie::_collect_seqs(std::string acc_seq, std::vector<std::string>& output) const
{
	if(is_end) output.emplace_back(acc_seq);
	for(uint8_t i = 0; i < (uint8_t)nodes.size(); ++i) {
		if(nodes[i]) nodes[i]->_collect_seqs(acc_seq + encoding2nt[i], output);
	}
}

