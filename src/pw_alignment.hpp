#ifndef PW_ALIGNMENT_HPP
#define PW_ALIGNMENT_HPP

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <math.h> 

#include <cassert>
#include <cstdlib>

using namespace std;


class pw_alignment {
	public:
	
	pw_alignment(string sample1, string sample2, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end,size_t sample1reference, size_t sample2reference);
	pw_alignment();
	~pw_alignment();
	pw_alignment(const pw_alignment & p);


	size_t alignment_length() const;
	void alignment_col(size_t c, char & s1, char & s2) const;
	vector<bool> getsample1()const;
	vector<bool> getsample2()const;
	size_t getbegin1()const;
	size_t getbegin2() const;
	size_t getend1()const;
	size_t getend2()const;
	size_t getreference1() const;
	size_t getreference2() const;
	void split(bool sample, size_t position, pw_alignment & first_part, pw_alignment & second_part) const;
	void remove_end_gaps(pw_alignment & res) const;
	void set_alignment_bits(vector<bool> s1, vector<bool> s2);
	vector<bool> getsample(size_t id)const;
	size_t getbegin(size_t id)const;
	size_t getend(size_t id)const;
	size_t getreference(size_t id)const;

	string get_al_ref1() const;
	string get_al_ref2() const;
	
	void setbegin1(size_t begin1);
	void setbegin2(size_t begin2);
	
	void setend1(size_t end1);
	void setend2(size_t end2);

	void setreference1(size_t ref1);
	void setreference2(size_t ref2);

	void get_lr_on_reference(size_t sequence, size_t & left, size_t & right) const;
	void get_lr1(size_t & left, size_t & right) const;
	void get_lr2(size_t & left, size_t & right) const;
	void split_on_reference(size_t sequence, size_t pos, pw_alignment & fp, pw_alignment & sp) const;

	void print() const;

	void set_cost(vector<double> create, vector<double> modify);

	double get_create1() const;
	double get_create2() const;
	double get_modify1() const;
	double get_modify2() const;

	private:
	vector<vector<bool> > samples;
	vector<size_t> begins;
	vector<size_t> ends;
	vector<size_t> references;
	vector<double> create_costs;
	vector<double> modify_costs;
	
	static inline void base_translate(char base, bool & bit1, bool & bit2, bool & bit3);
	static inline char base_translate_back(bool bit1, bool bit2, bool bit3);
};


class compare_pw_alignment {
	public:
	bool operator()(const pw_alignment *const &a, const pw_alignment *const &b) const ;
};

class wrapper{
	public:
		wrapper();
		~wrapper();
		void encode(unsigned int&, unsigned int &, unsigned int&);
		void decode(unsigned int &, unsigned int &, unsigned int&);

	private:
		ofstream encodeout;
		ofstream decodeout;
};
		
#endif 



