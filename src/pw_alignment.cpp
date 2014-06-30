#include "pw_alignment.hpp"
#include <cassert>
#include <cstdlib>

pw_alignment::pw_alignment(string sample1str, string sample2str, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end, size_t sample1reference, size_t sample2reference ) {
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

pw_alignment::pw_alignment() {}

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
void  pw_alignment::set_alignment_bits(vector<bool> s1, vector<bool> s2){
	sample1 = s1;
	sample2 = s2;
}

void pw_alignment::split(bool sample, size_t position, pw_alignment & first_part, pw_alignment & second_part ) const{
	vector<bool> fps1(sample1.size());
	vector<bool> fps2(sample1.size());
	vector<bool> sps1(sample1.size());
	vector<bool> sps2(sample1.size());

cout << "size " << sample1.size() << endl;


	size_t s = 0;
	if (sample == true) {
	for (size_t i=0; i<sample1.size()/3; ++i){
cout << "bits " << sample1.at(i*3) << sample1.at(i*3+1) << sample1.at(i*3+2)<< endl;
 		if (sample1.at(0+i*3)== true && sample1.at(1+i*3)==false && sample1.at(2+i*3)==true) {
cout << "gap at " << i << endl;
			s=s+1;
		}
	fps1.at(3*i) = sample1.at(3*i);
	fps2.at(3*i) = sample2.at(3*i);
	fps1.at(3*i+1) = sample1.at(3*i+1);
	fps2.at(3*i+1) = sample2.at(3*i+1);
	fps1.at(3*i+2) = sample1.at(3*i+2);
	fps2.at(3*i+2) = sample2.at(3*i+2);
	if (i == position-1+s)
		break;
	
}

	fps1.resize(3*(position+s));
	fps2.resize(3*(position+s));
	first_part.set_alignment_bits(fps1,fps2); 

for (size_t i=position+s;i<sample1.size()/3;++i){
      	sps1.at(3*(i-position-s)) = sample1.at(3*i);
	sps2.at(3*(i-position-s)) = sample2.at(3*i);
	sps1.at(3*(i-position-s)+1) = sample1.at(3*i+1);
	sps2.at(3*(i-position-s)+1) = sample2.at(3*i+1);
	sps1.at(3*(i-position-s)+2) = sample1.at(3*i+2);
	sps2.at(3*(i-position-s)+2) = sample2.at(3*i+2);

}

	sps1.resize(sample1.size()-3*(position+s));
	sps2.resize(sample1.size()-3*(position+s));
	cout << "size " << sps1.size() << endl;

	second_part.set_alignment_bits(sps1,sps2);


}

	else{
    for (size_t i=0; i<sample2.size()/3; ++i){
cout << "bits " << sample2.at(i*3) << sample2.at(i*3+1) << sample2.at(i*3+2)<< endl;
 		if (sample2.at(0+i*3)== true && sample2.at(1+i*3)==false && sample2.at(2+i*3)==true) {
cout << "gap at " << i << endl;
			s=s+1;
		}
	fps1.at(3*i) = sample1.at(3*i);
	fps2.at(3*i) = sample2.at(3*i);
	fps1.at(3*i+1) = sample1.at(3*i+1);
	fps2.at(3*i+1) = sample2.at(3*i+1);
	fps1.at(3*i+2) = sample1.at(3*i+2);
	fps2.at(3*i+2) = sample2.at(3*i+2);
	if (i == position-1+s)
		break;
	
}

	fps1.resize(3*(position+s));
	fps2.resize(3*(position+s));
	first_part.set_alignment_bits(fps1,fps2); 

for (size_t i=position+s;i<sample2.size()/3;++i){
      	sps1.at(3*(i-position-s)) = sample1.at(3*i);
	sps2.at(3*(i-position-s)) = sample2.at(3*i);
	sps1.at(3*(i-position-s)+1) = sample1.at(3*i+1);
	sps2.at(3*(i-position-s)+1) = sample2.at(3*i+1);
	sps1.at(3*(i-position-s)+2) = sample1.at(3*i+2);
	sps2.at(3*(i-position-s)+2) = sample2.at(3*i+2);

}

	sps1.resize(sample2.size()-3*(position+s));
	sps2.resize(sample2.size()-3*(position+s));
	cout << "size " << sps1.size() << endl;

	second_part.set_alignment_bits(sps1,sps2);


}
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
	vector<bool> pw_alignment::getsample1()const{
	return sample1;
}
	vector<bool> pw_alignment::getsample2()const{
	return sample2;
}
	size_t pw_alignment::getbegin1()const{
	return sample1_begin;
}
	size_t pw_alignment::getbegin2()const{
	return sample2_begin;
}
	size_t pw_alignment::getend1()const{
	return sample1_end;
}
	size_t pw_alignment::getend2()const{
	return sample2_end;
}


	overlap::overlap(){

	}

	void overlap::add_alignment(pw_alignment & new_alignment){

	pw_alignment p1;
	pw_alignment p2;
	pw_alignment p3;
	pw_alignment p4;
	map alignments_on_reference1; 
	map alignments_on_reference2;
	//if( allalignit->second.getreference1()){

	size_t leftinsert = new_alignment.getbegin1();
	size_t rightinsert = new_alignment.getend1();
	if(new_alignment.getend1() < leftinsert) {
		leftinsert = new_alignment.getend1();
		rightinsert = new_alignment.getbegin1();
}
	map::const_iterator allalignit = alignments_on_reference1.lower_bound(leftinsert);


	if(allalignit == alignments_on_reference1.end()) allalignit = alignments_on_reference1.begin();

	for(; allalignit != alignments_on_reference1.end(); ++allalignit) {
		const pw_alignment & al = allalignit->second;
		size_t alleft;
		size_t alright = 0;
		alleft = al.getbegin1();
		
		if(alleft > al.getend1()) alleft = al.getend1();
			
		
		if(leftinsert <= alleft && rightinsert >= alright) {

		al.split (true,rightinsert,p1,p2);
		p2.split(true,rightinsert,p3,p4);

			}

			

		}


		// break loop if too far to the  right (??)
	
}



//}
	overlap::~overlap(){
	
}
	
