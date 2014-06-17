#include <iostream>
#include <string>

#include "pw_alignment.hpp"

using namespace std;



int main(int argc, char * argv[]) {
	cout << " hello " << endl;

	pw_alignment p(string("ACT"), string("AGT"), 12, 15, 14, 17);

	for(size_t i=0; i<p.alignment_length(); ++i) {
		char s1;
		char s2;
		p.alignment_col(i, s1, s2);
		cout << "pos " << i << " s1 " << s1 << " s2 " << s2 << endl;
	  
	}

	return 0;

}
