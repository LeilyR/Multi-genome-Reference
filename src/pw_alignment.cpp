#include "pw_alignment.hpp"

pw_alignment::pw_alignment(string sample1str, string sample2str, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end, size_t sample1reference, size_t sample2reference ): samples(2), begins(2), ends(2), references(2) {
	assert(sample1str.length() == sample2str.length());


	size_t bitlength = 3*sample1str.size();
	samples.at(0).resize(bitlength);
	samples.at(1).resize(bitlength);

	for(size_t i=0; i<sample1str.length(); ++i) {
		char base1 = sample1str.at(i);
		char base2 = sample2str.at(i);
		bool bit1;
		bool bit2;
		bool bit3;
		base_translate(base1, bit1, bit2, bit3);
		samples.at(0).at(3*i) = bit1;
		samples.at(0).at(3*i+1) = bit2;
		samples.at(0).at(3*i+2) = bit3;
		base_translate(base2, bit1, bit2, bit3);
		samples.at(1).at(3*i) = bit1;
		samples.at(1).at(3*i+1) = bit2;
		samples.at(1).at(3*i+2) = bit3;
	}
	begins.at(0) = sample1_begin;
	begins.at(1) = sample2_begin;
	ends.at(0) = sample1_end;
	ends.at(1) = sample2_end;
	references.at(0) = sample1reference;
	references.at(1) = sample2reference;
}

pw_alignment::pw_alignment(): samples(2), begins(2), ends(2), references(2) {}

pw_alignment::~pw_alignment() {

}


pw_alignment::pw_alignment(const pw_alignment & p) {
	begins = p.begins;
	ends = p.ends;
	samples = p.samples;
	references = p.references;

}

size_t pw_alignment::alignment_length() const {
	return samples.at(0).size() / 3;
}


void pw_alignment::alignment_col(size_t c, char & s1, char & s2) const {
	bool bit1;
	bool bit2;
	bool bit3;
	bit1 = samples.at(0).at(3*c);
	bit2 = samples.at(0).at(3*c+1);
	bit3 = samples.at(0).at(3*c+2);
	s1 = base_translate_back(bit1, bit2, bit3);
	bit1 = samples.at(1).at(3*c);
	bit2 = samples.at(1).at(3*c+1);
	bit3 = samples.at(1).at(3*c+2);
	s2 = base_translate_back(bit1, bit2, bit3);
}
void  pw_alignment::set_alignment_bits(vector<bool> s1, vector<bool> s2){
	samples.at(0) = s1;
	samples.at(1) = s2;
}


/* cut before position
*/
void pw_alignment::split(bool sample, size_t position, pw_alignment & first_part, pw_alignment & second_part ) const{
	vector<bool> fps1(samples.at(0).size());
	vector<bool> fps2(samples.at(0).size());
	vector<bool> sps1(samples.at(0).size());
	vector<bool> sps2(samples.at(0).size());
	size_t s = 0;
	size_t s2 = 0;
	if (sample == true) {
		if (getbegin1()<getend1()){
			assert(position > getbegin1());
				for (size_t i=0; i<samples.at(0).size()/3; ++i){
 					if (samples.at(0).at(0+i*3)== true && samples.at(0).at(1+i*3)==false && samples.at(0).at(2+i*3)==true) {
						s=s+1;
					}
					if (samples.at(1).at(0+i*3)== true && samples.at(1).at(1+i*3)==false && samples.at(1).at(2+i*3)==true) {
						s2=s2+1;
 					}
			
					fps1.at(3*i) = samples.at(0).at(3*i);
					fps2.at(3*i) = samples.at(1).at(3*i);
					fps1.at(3*i+1) = samples.at(0).at(3*i+1);
					fps2.at(3*i+1) = samples.at(1).at(3*i+1);
					fps1.at(3*i+2) = samples.at(0).at(3*i+2);
					fps2.at(3*i+2) = samples.at(1).at(3*i+2);
					if (i == position+s-getbegin1()-1)
						break;
				}
				fps1.resize(3*(position+s-getbegin1()));
				fps2.resize(3*(position+s-getbegin1()));
				first_part.set_alignment_bits(fps1,fps2); 
				first_part.setbegin1(getbegin1());
				first_part.setend1(position-1);
				first_part.setbegin2(getbegin2());
				if(getbegin2()<getend2()){
					first_part.setend2(getbegin2()+position+s-s2-getbegin1()-1);
				}
				else{
		//			cout<<"s: "<<s<<endl;
		//			cout<<"s2: "<<s2<<endl;
					first_part.setend2(getbegin2()-position-s+s2+getbegin1()+1);
				}
				
				first_part.setreference1(getreference1());
				first_part.setreference2(getreference2());
				
				for (size_t i=position+s-getbegin1();i<samples.at(0).size()/3;++i){
      					sps1.at(3*(i-position-s+getbegin1())) = samples.at(0).at(3*i);
					sps2.at(3*(i-position-s+getbegin1())) = samples.at(1).at(3*i);
					sps1.at(3*(i-position-s+getbegin1())+1) = samples.at(0).at(3*i+1);
					sps2.at(3*(i-position-s+getbegin1())+1) = samples.at(1).at(3*i+1);
					sps1.at(3*(i-position-s+getbegin1())+2) = samples.at(0).at(3*i+2);
					sps2.at(3*(i-position-s+getbegin1())+2) = samples.at(1).at(3*i+2);
				}
				sps1.resize(samples.at(0).size()-fps1.size());
				sps2.resize(samples.at(0).size()-fps1.size());
				second_part.set_alignment_bits(sps1,sps2);
				second_part.setbegin1(position);
				if(getbegin2()<getend2()){
					second_part.setbegin2(getbegin2()+first_part.alignment_length()-s2);
				}
				else{// cout<<"hey"<<endl;							
					second_part.setbegin2(getbegin2()-first_part.alignment_length()+s2);
				}
				second_part.setend1(getend1());
				second_part.setend2(getend2());
				second_part.setreference1(getreference1());
				second_part.setreference2(getreference2());
			
		}
		else{
			assert(position <= getbegin1());
			assert(position > getend1());
				for (size_t i=0; i<samples.at(0).size()/3; ++i){
 					if (samples.at(0).at(0+i*3)== true && samples.at(0).at(1+i*3)==false && samples.at(0).at(2+i*3)==true) {
						s=s+1;
					}
					if (samples.at(1).at(0+i*3)== true && samples.at(1).at(1+i*3)==false && samples.at(1).at(2+i*3)==true) {
						s2=s2+1;
					}
					fps1.at(3*i) = samples.at(0).at(3*i);
					fps2.at(3*i) = samples.at(1).at(3*i);
					fps1.at(3*i+1) = samples.at(0).at(3*i+1);
					fps2.at(3*i+1) = samples.at(1).at(3*i+1);
					fps1.at(3*i+2) = samples.at(0).at(3*i+2);
					fps2.at(3*i+2) = samples.at(1).at(3*i+2);
		
					if (i == getbegin1()-position+s)
						break;
				}
				fps1.resize(3*(getbegin1()-position+s+1));
				fps2.resize(3*(getbegin1()-position+s+1));
				first_part.set_alignment_bits(fps1,fps2); 
				first_part.setbegin1(getbegin1());
				first_part.setbegin2(getbegin2());
				first_part.setend1(position);
				if(getbegin2()<getend2()){
					first_part.setend2(getbegin2()+getbegin1()-position+s-s2);
				}else 
					first_part.setend2(getbegin2()-getbegin1()+position-s+s2);
				first_part.setreference1(getreference1());
				first_part.setreference2(getreference2());

				for (size_t i=getbegin1()-position+1+s; i < samples.at(0).size();++i){
      					sps1.at(3*(i-getbegin1()+position-1-s)) = samples.at(0).at(3*i);
					sps2.at(3*(i-getbegin1()+position-1-s)) = samples.at(1).at(3*i);
					sps1.at(3*(i-getbegin1()+position-1-s)+1) = samples.at(0).at(3*i+1);
					sps2.at(3*(i-getbegin1()+position-1-s)+1) = samples.at(1).at(3*i+1);
					sps1.at(3*(i-getbegin1()+position-1-s)+2) = samples.at(0).at(3*i+2);
					sps2.at(3*(i-getbegin1()+position-1-s)+2) = samples.at(1).at(3*i+2);
				}

				sps1.resize(samples.at(0).size()-fps1.size());
				sps2.resize(samples.at(0).size()-fps1.size());

				second_part.set_alignment_bits(sps1,sps2);
				second_part.setbegin1(position-1);
				if(getbegin2()<getend2()){
					second_part.setbegin2(getbegin2()+first_part.alignment_length()-s2);
				}else
					second_part.setbegin2(getbegin2()-first_part.alignment_length()+s2);
				second_part.setend1(getend1());
				second_part.setend2(getend2());
				second_part.setreference1(getreference1());
				second_part.setreference2(getreference2());
			
		}

	}

	else{
		if(getbegin2()<getend2()){
		assert(position > getbegin2());	
   				 for (size_t i=0; i<samples.at(1).size()/3; ++i){
 					if (samples.at(1).at(0+i*3)== true && samples.at(1).at(1+i*3)==false && samples.at(1).at(2+i*3)==true) {
						s=s+1;
					}
					if (samples.at(0).at(0+i*3)== true && samples.at(0).at(1+i*3)==false && samples.at(0).at(2+i*3)==true) {
						s2=s2+1;
					}
					fps1.at(3*i) = samples.at(0).at(3*i);
					fps2.at(3*i) = samples.at(1).at(3*i);
					fps1.at(3*i+1) = samples.at(0).at(3*i+1);
					fps2.at(3*i+1) = samples.at(1).at(3*i+1);
					fps1.at(3*i+2) = samples.at(0).at(3*i+2);
					fps2.at(3*i+2) = samples.at(1).at(3*i+2);
					if (i == position+s-getbegin2()-1)
						break;
				}

				fps1.resize(3*(position+s-getbegin2()));
				fps2.resize(3*(position+s-getbegin2()));
				first_part.set_alignment_bits(fps1,fps2);
				first_part.setbegin2(getbegin2());
				first_part.setbegin1(getbegin1());
				if(getbegin1()<getend1()){
					first_part.setend1(getbegin1()+(position+s-s2-getbegin2())-1);
				}else{
					first_part.setend1(getbegin1()-position-s+s2+getbegin2()+1);
				}
			
				first_part.setend2(position-1);
				first_part.setreference1(getreference1());
				first_part.setreference2(getreference2());
	
				for (size_t i=position+s-getbegin2();i<samples.at(1).size()/3;++i){
					sps1.at(3*(i-position-s+getbegin2())) = samples.at(0).at(3*i);
					sps2.at(3*(i-position-s+getbegin2())) = samples.at(1).at(3*i);
					sps1.at(3*(i-position-s+getbegin2())+1) = samples.at(0).at(3*i+1);
					sps2.at(3*(i-position-s+getbegin2())+1) = samples.at(1).at(3*i+1);
					sps1.at(3*(i-position-s+getbegin2())+2) = samples.at(0).at(3*i+2);
					sps2.at(3*(i-position-s+getbegin2())+2) = samples.at(1).at(3*i+2);
				}

				sps1.resize(samples.at(1).size()-fps2.size());
				sps2.resize(samples.at(1).size()-fps2.size());
				second_part.set_alignment_bits(sps1,sps2);
				if(getbegin1()<getend1()){
					second_part.setbegin1(getbegin1()+first_part.alignment_length()-s2);
				}else
					second_part.setbegin1(getbegin1()-first_part.alignment_length()+s2);					
				second_part.setbegin2(position);
				second_part.setend1(getend1());
				second_part.setend2(getend2());
				second_part.setreference1(getreference1());
				second_part.setreference2(getreference2());
		
	
		}
		else{
			assert(position <= getbegin2());
			assert(position > getend2());

			for (size_t i=0; i<samples.at(1).size()/3; ++i){
 			
				if (samples.at(1).at(0+i*3)== true && samples.at(1).at(1+i*3)==false && samples.at(1).at(2+i*3)==true) {
				s=s+1;
				}
				if (samples.at(0).at(0+i*3)== true && samples.at(0).at(1+i*3)==false && samples.at(0).at(2+i*3)==true) {
				s2=s2+1;
				}
				

				fps1.at(3*i) = samples.at(0).at(3*i);
				fps2.at(3*i) = samples.at(1).at(3*i);
				fps1.at(3*i+1) = samples.at(0).at(3*i+1);
				fps2.at(3*i+1) = samples.at(1).at(3*i+1);
				fps1.at(3*i+2) = samples.at(0).at(3*i+2);
				fps2.at(3*i+2) = samples.at(1).at(3*i+2);
		
				if (i == getbegin2()-position+s)
				break;
			}

			fps1.resize(3*(getbegin2()-position+s+1));
			fps2.resize(3*(getbegin2()-position+s+1));
			first_part.set_alignment_bits(fps1,fps2); 
			first_part.setbegin1(getbegin1());
			first_part.setbegin2(getbegin2());
			if(getbegin1()<getend1()){
				first_part.setend1(getbegin1()+getbegin2()-position+s-s2);			
			}else 
				first_part.setend1(getbegin1()-getbegin2()+position-s+s2);
			first_part.setend2(position);
			first_part.setreference1(getreference1());
			first_part.setreference2(getreference2());

			for (size_t i=getbegin2()-position+s+1; i < samples.at(1).size()/3;++i){

     		 		sps1.at(3*(i-getbegin2()+position-s-1)) = samples.at(0).at(3*i);
				sps2.at(3*(i-getbegin2()+position-s-1)) = samples.at(1).at(3*i);
				sps1.at(3*(i-getbegin2()+position-s-1)+1) = samples.at(0).at(3*i+1);
				sps2.at(3*(i-getbegin2()+position-s-1)+1) = samples.at(1).at(3*i+1);
				sps1.at(3*(i-getbegin2()+position-s-1)+2) = samples.at(0).at(3*i+2);
				sps2.at(3*(i-getbegin2()+position-s-1)+2) = samples.at(1).at(3*i+2);
			}
	//		cout<<"sample size: "<<samples.at(1).size()<< "first part size: "<< fps2.size()<<endl;
			sps1.resize(samples.at(1).size()-fps2.size());
			sps2.resize(samples.at(1).size()-fps2.size());
			second_part.set_alignment_bits(sps1,sps2);
			if(getbegin1()<getend1()){
				second_part.setbegin1(getbegin1()+first_part.alignment_length()-s2);
			}else
				second_part.setbegin1(getbegin1()-first_part.alignment_length()+s2);
			second_part.setbegin2(position-1);
			second_part.setend1(getend1());
			second_part.setend2(getend2());
			second_part.setreference1(getreference1());
			second_part.setreference2(getreference2());

			}
	
	
	}
}


void pw_alignment::get_lr_on_reference(size_t sequence, size_t & left, size_t & right) const {
	size_t ref_idx = 2;
	assert(references.at(0)!=references.at(1));
	if(sequence == references.at(0)) ref_idx= 0;
	if(sequence == references.at(1)) ref_idx= 1;
	assert(ref_idx < 2);
	left = begins.at(ref_idx);
	right = ends.at(ref_idx);
	if(left > right) {
		left = ends.at(ref_idx);
		right = begins.at(ref_idx);
	}
}


void pw_alignment::get_lr1(size_t & left, size_t & right) const {
	left = begins.at(0);
	right = ends.at(0);
	if(left > right) {
		left = ends.at(0);
		right = begins.at(0);
	}
}

void pw_alignment::get_lr2(size_t & left, size_t & right) const {
	left = begins.at(1);
	right = ends.at(1);
	if(left > right) {
		left = ends.at(1);
		right = begins.at(1);
	}

}


void pw_alignment::split_on_reference(size_t sequence, size_t pos, pw_alignment & fp, pw_alignment & sp) const {
	size_t ref_idx = 2;
	if(sequence == references.at(0)) ref_idx= 0;
	if(sequence == references.at(1)) ref_idx= 1;
	assert(ref_idx < 2);
	split(1-ref_idx, pos, fp, sp);
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
			case 'R':
			case 'r':
			case 'Y':
			case 'y':
			case 'M':
			case 'm':
			case 'K':
			case 'k':
			case 'W':
			case 'w':
			case 'S':
			case 's':
			case 'B':
			case 'b':
			case 'D':
			case 'd':
			case 'H':
			case 'h':
			case 'V':
			case 'v':
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
			cerr << "Error: Illegal character in DNA sequence: " << base << endl;
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
	return samples.at(0);
}
	vector<bool> pw_alignment::getsample2()const{
	return samples.at(1);
}
	size_t pw_alignment::getbegin1()const{
	return begins.at(0);
}
	size_t pw_alignment::getbegin2()const{
	return begins.at(1);
}
	size_t pw_alignment::getend1()const{
	return ends.at(0);
}
	size_t pw_alignment::getend2()const{
	return ends.at(1);
}
	size_t pw_alignment::getreference1() const{
	return references.at(0);
}
	size_t pw_alignment::getreference2() const{
	return references.at(1);
}
	vector<bool> pw_alignment::getsample(size_t id)const{
	return samples.at(id);
}
	size_t pw_alignment::getbegin(size_t id)const{
	return begins.at(id);
}
	size_t pw_alignment::getend(size_t id)const{
	return ends.at(id);
}
	size_t pw_alignment::getreference(size_t id)const{
	return references.at(id);
}
	void pw_alignment::setbegin1(size_t begin1){
		begins.at(0) = begin1;
}
	void pw_alignment::setbegin2(size_t begin2){
		begins.at(1) = begin2;
}

	void pw_alignment::setend1(size_t end1){
		ends.at(0) = end1;
}

	void pw_alignment::setend2(size_t end2){
		ends.at(1) = end2;
}
	void pw_alignment::setreference1(size_t ref1){
		references.at(0)=ref1;
}
	void pw_alignment::setreference2(size_t ref2){
		references.at(1)=ref2;
}

	void pw_alignment::print()const{
 
//	cout << "al1 seq "<< getreference1() << " b " << getbegin1() << " e " << getend1() << endl;
//	cout << "al2 seq "<< getreference2() << " b " << getbegin2() << " e " << getend2() << endl;

	for(size_t col = 0; col < alignment_length(); col++) {
		char c1;
		char c2;
		alignment_col(col, c1, c2);
		cout <<col <<"/t"<< c1<<"\t"<<c2<<endl;
	if(col>50) break;	
	}
	} 





bool compare_pw_alignment::operator()(const pw_alignment *const &a, const pw_alignment *const &b) const {
	size_t asmaller = 0;
	size_t abigger = 1;
	if(a->getreference(0) > a->getreference(1)) {
		asmaller = 1;
		abigger = 0;
	} else if(a->getreference(0) == a->getreference(1)) {
		if(a->getbegin(0) > a->getbegin(1)) {
			asmaller =1;
			abigger = 0;
		} else if(a->getbegin(0) == a->getbegin(1)) {
			 if(a->getend(0) > a->getend(1)) {
				asmaller =1;
				abigger = 0;
			}
			assert(a->getend(0)!=a->getend(1));
		}
	}

	
	size_t bsmaller = 0;
	size_t bbigger = 1;
	if (b->getreference(0) > b->getreference(1)){
		bsmaller = 1;
		bbigger = 0;
} 
	else if (b->getreference(0) == b->getreference(1)){
		if(b->getbegin(0) > b->getbegin(1)) {
			bsmaller =1;
			bbigger = 0;
		}else if(b->getbegin(0) == b->getbegin(1)) {
			 if(b->getend(0) > b->getend(1)) {
				bsmaller =1;
				bbigger = 0;
			}
			assert(b->getend(0)!=b->getend(1));
		}	
	}


	if (a->getbegin(asmaller) < b->getbegin(bsmaller))return true;
	if (a->getbegin(asmaller) > b->getbegin(bsmaller))return false;

	if (a->getend(asmaller) < b->getend(bsmaller))return true;
	if (a->getend(asmaller) > b->getend(bsmaller))return false;

	if ( a->getbegin(abigger) < b->getbegin(bbigger)) return true;
	if ( a->getbegin(abigger) > b->getbegin(bbigger)) return false;
	
	if ( a->getend(abigger) < b->getend(bbigger)) return true;
	if ( a->getend(abigger) > b->getend(bbigger)) return false;
	return false;
}

	void pw_alignment::set_cost(vector<double> create, vector<double> modify){
		create_costs = create;
		modify_costs = modify; 
	}


	double pw_alignment::get_create1() const {
		return create_costs.at(0);
	}
	double pw_alignment::get_create2() const {
		return create_costs.at(1);
	}
	double pw_alignment::get_modify1() const {
		return modify_costs.at(0);
	}
	double pw_alignment::get_modify2() const {
		return modify_costs.at(1);
	}
	
