#include <string>
#include <vector>
#include <array>
#include <memory>
#include <utility>

inline std::array<uint8_t, 256> nt2encoding;
inline std::array<char, 5> encoding2nt;
inline std::array<char, 5> encoding2rcnt;

void init_encoding();

void get_kmer_min_indeces(
	const std::string&,
	const uint32_t&, 
	std::vector<std::pair<uint64_t,uint64_t>>&
);

void index2kmer(
	std::string&, 
	uint64_t, 
	uint32_t
);

class Trie{
	public:
		std::array<std::unique_ptr<Trie>, 4> nodes;
		bool is_end;
		Trie();
		int contains(const std::string&, const int&, const int) const;
		void set_seq(const std::string&);
		void collect_seqs(std::vector<std::string>&) const;

	private:
		int _contains(const std::string&, const int&, const int, const int) const;
		void _set_seq(const std::string&, const int&, int);
		void _collect_seqs(std::string, std::vector<std::string>&) const;
};

