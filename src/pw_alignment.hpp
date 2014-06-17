#ifndef PW_ALIGNMENT_HPP
#define PW_ALIGNMENT_HPP

#include <vector>
#include <iostream>
#include <string>

using namespace std;

class pw_alignment {
	public:
	
	pw_alignment(string sample1, string sample2, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end);
	~pw_alignment();
	pw_alignment(const pw_alignment & p);


	size_t alignment_length() const;
	void alignment_col(size_t c, char & s1, char & s2) const;



	private:
	vector<bool> sample1;
	vector<bool> sample2;
	size_t sample1_begin;
	size_t sample2_begin;
	size_t sample1_end;
	size_t sample2_end;



	static inline void base_translate(char base, bool & bit1, bool & bit2, bool & bit3);
	static inline char base_translate_back(bool bit1, bool bit2, bool bit3);

};





#endif 



