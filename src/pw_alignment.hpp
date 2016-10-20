#ifndef PW_ALIGNMENT_HPP
#define PW_ALIGNMENT_HPP

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <math.h> 
#include <limits>

#include <cassert>
#include <cstdlib>



#define ALIGNMENT_BLOCK_SIZE 1000
// longer piece with gaps
struct alblock {
	size_t cfrom;
	size_t cto;
	double c1;
	double c2;
	double m1;
	double m2;
};



class pw_alignment {
	public:
	
	pw_alignment(std::string sample1, std::string sample2, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end,size_t sample1reference, size_t sample2reference);
	pw_alignment(const size_t & begin1, const size_t & begin2, const size_t & end1, const size_t & end2, const size_t & reference1, const size_t & reference2);
	pw_alignment();
	virtual ~pw_alignment();
	pw_alignment(const pw_alignment & p);


	virtual size_t alignment_length() const;
	virtual void alignment_col(size_t c, char & s1, char & s2) const;

	// TODO remove this function (slow copy, bad oo design)
	std::vector<bool> getsample1()const;
	// TODO remove this function 
	std::vector<bool> getsample2()const;

	virtual size_t getbegin1()const;
	virtual size_t getbegin2() const;
	virtual size_t getend1()const;
	virtual size_t getend2()const;
	virtual size_t getreference1() const;
	virtual size_t getreference2() const;
	void split(bool sample, size_t position, pw_alignment & first_part, pw_alignment & second_part) const;
	void remove_end_gaps(pw_alignment & res) const;

	// TODO remove
	void set_alignment_bits(std::vector<bool> s1, std::vector<bool> s2);
	
	// TODO slow, bad design, remove
	std::vector<bool> getsample(size_t id)const;


	virtual size_t getbegin(size_t id)const;
	virtual size_t getend(size_t id)const;
	virtual size_t getreference(size_t id)const;

	static void get_bits(char&, std::vector<bool> &);
	static void get_base(char&, bool&, bool&, bool&);

	virtual std::string get_al_ref1() const;
	virtual std::string get_al_ref2() const;
	
	virtual void setbegin1(size_t begin1);
	virtual void setbegin2(size_t begin2);
	
	virtual void setend1(size_t end1);
	virtual void setend2(size_t end2);

	virtual void setreference1(size_t ref1);
	virtual void setreference2(size_t ref2);

	virtual void get_lr_on_reference(size_t sequence, size_t & left, size_t & right) const;
	virtual void get_lr1(size_t & left, size_t & right) const;
	virtual void get_lr2(size_t & left, size_t & right) const;
	virtual void split_on_reference(size_t sequence, size_t pos, pw_alignment & fp, pw_alignment & sp) const;
	virtual void single_ref_intersect(size_t alref, const pw_alignment & other, size_t otherref, pw_alignment & result, pw_alignment & otherresult) const;
	virtual bool is_same_direction() const;

	virtual void print() const;

	virtual void set_cost(const std::vector<double> & create, const std::vector<double> & modify) const;
	virtual bool is_cost_cached() const;
	virtual double get_create1() const;
	virtual double get_create2() const;
	virtual double get_modify1() const;
	virtual double get_modify2() const;


	virtual void is_block_cached(std::vector<bool> & isbc) const;
	virtual void set_block_cache(const size_t & bnum, const double & c1, const double & c2, const double & m1, const double & m2) const;
	virtual void set_block_cached(const size_t & bl) const;
	
	virtual bool equals(const pw_alignment & al) const;
	virtual std::vector<alblock> & getblocks() const;

	// inefficient in pw_aligment, we assert that it is never used
	virtual void alignment_position(size_t col, size_t & r1pos, size_t & r2pos) const;

	virtual pw_alignment & operator=(const pw_alignment & al);

	// TODO remove this
	std::vector<bool> reverse_complement(std::vector<bool> & ); 

	virtual void get_reverse(pw_alignment &);
	virtual void get_reverse_complement_sample(std::vector<std::vector<bool> > &);
	virtual void simulate_split(const size_t & ref, const size_t & pos, size_t & sp1, size_t & sp2) const ;
	virtual bool onlyGapSample() const;
	protected:
	std::vector<size_t> begins;
	std::vector<size_t> ends;
	std::vector<size_t> references;

	static inline void base_translate(char base, bool & bit1, bool & bit2, bool & bit3);
	static inline char base_translate_back(bool bit1, bool bit2, bool bit3);
	static inline char complement(char &);


	// cache costs here:
	mutable bool costs_cached;
	mutable std::vector<double> create_costs;
	mutable std::vector<double> modify_costs;
	mutable std::vector<alblock> blocks;

	private:

	
	std::vector<std::vector<bool> > samples;
};

// TODO if pw_alignment is not always fixed on a reference sequence, we should not use these comparators on it
class compare_pointer_pw_alignment {
	public:
	bool operator()(const pw_alignment *ar, const pw_alignment *br) const;
};

class compare_pw_alignment {
	public:
	bool operator()(const pw_alignment &a, const pw_alignment &b) const ;
};
class sort_pw_alignment{
	public:
	bool operator () (const pw_alignment &p1, const pw_alignment &p2 )const;
};
class sort_right_pw_alignment{
	public:
	bool operator () (const pw_alignment &p1, const pw_alignment &p2 )const;
};
class wrapper{
	public:
		wrapper();
		~wrapper();
		void encode(unsigned int&, unsigned int &, unsigned int&);
		void context(int &);
				
	private:
		std::ofstream encodeout;
		std::ofstream al_encode;
	};

class decoding_wrapper{
	public:
		decoding_wrapper();
		~decoding_wrapper();
		void decode(unsigned int &, unsigned int &, unsigned int&);
		void decodeContext(int &);
	private:
		std::ofstream decodeout;
		std::ofstream al_decode;
};


class dnastring;
class all_data;


class located_alignment : public pw_alignment {
public:
//	virtual located_alignment(const all_data * data, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end,size_t sample1reference, size_t sample2reference);
	virtual ~located_alignment();
	located_alignment(const all_data * data, const pw_alignment & p);
	located_alignment(const located_alignment & p);


	virtual size_t alignment_length() const;
	virtual void alignment_col(size_t c, char & s1, char & s2) const;
	virtual void alignment_col(size_t c, char & c1, char & c2, size_t & r1at, size_t & r2at) const;

	virtual void alignment_position(size_t col, size_t & r1pos, size_t & r2pos) const;

	virtual void is_block_cached(std::vector<bool> & isbc) const;
	virtual void set_block_cached(const size_t & bl) const;

	virtual void get_column(const size_t & ref, const size_t & refpos, size_t & alcol) const;
	virtual void simulate_split(const size_t & ref, const size_t & pos, size_t & sp1, size_t & sp2) const ;
	virtual void split(bool sample, size_t position, located_alignment & first_part, located_alignment & second_part) const;
	virtual void printseg() const;
private:
	// short piece without gaps
	struct segment {
		size_t begins[3]; // begin in r1, r2, and alignment col
		size_t length;
	};

	const all_data * data;
	size_t length;
	mutable size_t last_segment;
	std::vector<segment> segments; 
// we cache costs per block (BLOCK_SIZE alignment columns) to be much faster in splitting huge alignments and in recomputing cost of split parts
// edge effects between blocks can be ignored as BLOCK SIZE is large
	
	mutable std::vector<bool> block_cached;




	virtual void make_blocks() const;
	void bsearch_col(const size_t & col, const size_t & from, const size_t & to, size_t & seg, size_t & r1at, size_t & r2at, char & c1, char & c2) const;
	void bsearch_ref(const size_t & ref, const size_t & refpos, const size_t & from, const size_t & to, size_t & seg, size_t & alcolumn) const;


};
		
#include "data.hpp" // need data classes, that were forward declared (cross-dependency)
#endif 



