#include <iostream>
#include <string>

#include "pw_alignment.hpp"

using namespace std;



int main(int argc, char * argv[]) {
	cout << " hello " << endl;

	pw_alignment p(string("ATT----TTCTT"), string("AGTGATAT---"), 12, 15, 23, 26);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	  
	}

	pw_alignment p1;
	pw_alignment p2;
	pw_alignment p3;
	pw_alignment p3;
	p.split(4, p1, p2, p3, p4);
        cout << "first part_s1" << pw_alignment p1 "second part_s1" << pw_alignment p2 << endl;
        cout << "first part_s2" << pw_alignment p3 "second part_s2" << pw_alignment p4 << endl;

	return 0;

}
