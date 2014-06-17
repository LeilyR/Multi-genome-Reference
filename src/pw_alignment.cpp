#include "pw_alignment.hpp"
#include <cassert>
#include <cstdlib>

pw_alignment::pw_alignment(string sample1str, string sample2str, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end) {
	assert(sample1str.length() == sample2str.length());
	size_t bitlength = 3*sample1str.size();
	sample1.resize(bitlength);
	sample2.resize(bitlength);

	for(size_t i=0; i<sample1str.length(); ++i) {
		char base1 = sample1str.at(i);
		char base2 = sample2str.at(i);
		bool bit1;
		bool bit2;
		bool bit3;
		base_translate(base1, bit1, bit2, bit3);
		sample1.at(3*i) = bit1;
		sample1.at(3*i+1) = bit2;
		sample1.at(3*i+2) = bit3;
		base_translate(base2, bit1, bit2, bit3);
		sample2.at(3*i) = bit1;
		sample2.at(3*i+1) = bit2;
		sample2.at(3*i+2) = bit3;
	}
}

pw_alignment::~pw_alignment() {

}


pw_alignment::pw_alignment(const pw_alignment & p) {

}

size_t pw_alignment::alignment_length() const {
	return sample1.size() / 3;
}


void pw_alignment::alignment_col(size_t c, char & s1, char & s2) const {
	bool bit1;
	bool bit2;
	bool bit3;
	bit1 = sample1.at(3*c);
	bit2 = sample1.at(3*c+1);
	bit3 = sample1.at(3*c+2);
	s1 = base_translate_back(bit1, bit2, bit3);
	bit1 = sample2.at(3*c);
	bit2 = sample2.at(3*c+1);
	bit3 = sample2.at(3*c+2);
	s2 = base_translate_back(bit1, bit2, bit3);
}




void pw_alignment::base_translate(char base, bool  & bit1, bool & bit2, bool & bit3) {
	switch(base) {
			case 'A':
			case 'a':
				bit1 = false;
				bit2 = false;
				bit3 = false;

			break;
			case 'C':
			case 'c':
				bit1 = true;
				bit2 = false;
				bit3 = false;

			break;
			case 'T':
			case 't':
				bit1 = false;
				bit2 = true;
				bit3 = false;

			break;
			case 'G':
			case 'g':
				bit1 = true;
				bit2 = true;
				bit3 = false;

			break;
			case 'N':
			case 'n':
				bit1 = false;
				bit2 = false;
				bit3 = true;

			break;
			case '-':
			case '.':
				bit1 = true;
				bit2 = false;
				bit3 = true;

			break;
			default:
			cerr << "Error: unknown character: " << base << endl;
			exit(1);
			break;
		}
}

char pw_alignment::base_translate_back(bool bit1, bool bit2, bool bit3) {

	if(bit1) {
		if(bit2) {
			if(bit3) {
				cerr << "Error: wrong bit vector 111" << endl;
				exit(1);
			} else {
				return 'G';
			}
		
		} else {
			if(bit3) {
				return '-';
			} else {
				return 'C';
			}

		}	
	
	} else {
		if(bit2) {
			if(bit3) {
				cerr << "Error: wrong bit vector 011" << endl;
				exit(1);
			} else {
				return 'T';
			}

		} else {
			if(bit3) {
				return 'N';
			} else {
				return 'A';
			}

		}
	}
}



