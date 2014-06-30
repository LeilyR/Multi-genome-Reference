#ifndef PW_ALIGNMENT_HPP
#define PW_ALIGNMENT_HPP

#include <vector>
#include <iostream>
#include <string>
#include <map>

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
	void set_alignment_bits(vector<bool> s1, vector<bool> s2);
	private:
	vector<bool> sample1;
	vector<bool> sample2;
	size_t sample1_begin;
	size_t sample2_begin;
	size_t sample1_end;
	size_t sample2_end;
	size_t sample1reference;
	size_t sample2reference;
	static inline void base_translate(char base, bool & bit1, bool & bit2, bool & bit3);
	static inline char base_translate_back(bool bit1, bool bit2, bool bit3);
};

class overlap{
public:
	overlap();
	~overlap();
	void add_alignment(pw_alignment & new_alignment);


	typedef multimap <size_t, pw_alignment> map;


	pw_alignment new_alignment;
	
};

		
#endif 



