#include "pw_alignment.hpp"

pw_alignment::pw_alignment(std::string sample1str, std::string sample2str, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end, size_t sample1reference, size_t sample2reference ): samples(2), begins(2), ends(2), references(2), costs_cached(false), create_costs(0), modify_costs(0){
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

pw_alignment::pw_alignment(): samples(2), begins(2), ends(2), references(2), costs_cached(false), create_costs(0), modify_costs(0) {}

pw_alignment::~pw_alignment() {

}


pw_alignment::pw_alignment(const pw_alignment & p) {
	begins = p.begins;
	ends = p.ends;
	samples = p.samples;
	references = p.references;
	costs_cached = p.costs_cached;
	create_costs = p.create_costs;
	modify_costs = p.modify_costs;

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
void  pw_alignment::set_alignment_bits(std::vector<bool> s1, std::vector<bool> s2){
	samples.at(0) = s1;
	samples.at(1) = s2;
	costs_cached = false;
}


/* cut before position
*/
void pw_alignment::split(bool sample, size_t position, pw_alignment & first_part, pw_alignment & second_part ) const{
	std::vector<bool> fps1(samples.at(0).size());
	std::vector<bool> fps2(samples.at(0).size());
	std::vector<bool> sps1(samples.at(0).size());
	std::vector<bool> sps2(samples.at(0).size());
	size_t s = 0;
	size_t s2 = 0;
	if (sample == true) {
		if (getbegin1()<getend1()){
			assert(position > getbegin1());
			assert(position <= getend1());
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
		//			std::cout<<"s: "<<s<<std::endl;
		//			std::cout<<"s2: "<<s2<<std::endl;
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
				else{// std::cout<<"hey"<<std::endl;							
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
		assert(position <= getend2());	
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
	//		std::cout<<"sample size: "<<samples.at(1).size()<< "first part size: "<< fps2.size()<<std::endl;
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

/*

   	remove gap containing columns at the ends of an alignment

	ATCTCTAAT--
	-TCTTTAATTT
	to
	 TCTCTAAT
	 TCTTTAAT
*/   
void pw_alignment::remove_end_gaps(pw_alignment & res) const {

//	std::cout << "reg: " << std::endl;
//	print();
//	std::cout << std::endl;

	// get first/last alignment column without gaps
	size_t first_col = 0;
	size_t last_col = alignment_length() -1;
	size_t start_gaps1 = 0;
	size_t start_gaps2 = 0;
	size_t end_gaps1 = 0;
	size_t end_gaps2 = 0;
	for(size_t i=0; i<=last_col; ++i) {
		char c1;
		char c2;
		alignment_col(i, c1, c2);
		if(c1=='-') start_gaps1++;
		if(c2=='-') start_gaps2++;
		if(c1!='-' && c2!='-') {
			first_col = i;
			break;
		}	
	}
	for(size_t i=last_col; i>=0; i--) {
		char c1;
		char c2;
		alignment_col(i, c1, c2);
		if(c1=='-') end_gaps1++;
		if(c2=='-') end_gaps2++;
		if(c1!='-' && c2!='-') {
			last_col = i;
			break;
		}
	}

//	std::cout << " start gaps " << start_gaps1 << " " << start_gaps2 << std::endl;
//	std::cout << " end gaps " << end_gaps1 << " " << end_gaps2 << std::endl;
//	std::cout << " first col " << first_col << " last col " << last_col << std::endl;


	assert(start_gaps1==0 || start_gaps2==0);
	assert(end_gaps1==0 || end_gaps2==0);

	res = *this;
	// cut alignment first part
	if(first_col>0) {
		pw_alignment gaps;
		if(start_gaps1>0) {
			if(getbegin2()<getend2()) {
//				std::cout << " split 1 " << std::endl;
				split(false, getbegin2()+start_gaps1, gaps, res);	
			} else {
//				std::cout << " split 2 " << std::endl;
				split(false, getbegin2()-start_gaps1+1, gaps, res); // falsch?
			}
		} else { // gaps on second ref
			if(getbegin1()<getend1()) {
//				std::cout << " split 3 " << std::endl;
				split(true, getbegin1()+start_gaps2, gaps, res);
			} else {
//				std::cout << " split 4 " << std::endl;
				split(true, getbegin1()-start_gaps2+1, gaps, res);
			}
		}
	}
	// cut alignment last part
	if(last_col + 1 < alignment_length()) {
		pw_alignment gaps;
		if(end_gaps1>0) {
			if(getbegin2()<getend2()) {
//				std::cout << " split 5 " << std::endl;
				res.split(false, getend2()-end_gaps1+1, res, gaps);
			} else {
//				std::cout << " split 6 " << std::endl;
				res.split(false, getend2()+end_gaps1, res, gaps);
			}
		} else { // gaps on first ref
			if(getbegin1()<getend1()) {
//				std::cout << " split 7 " << std::endl;
				res.split(true, getend1()-end_gaps2+1, res, gaps);
			} else {
//				std::cout << " split 8 " << std::endl;
				res.split(true, getend2()+end_gaps2, res, gaps);
			}
		}
	}

	res.costs_cached = false;

//	std::cout << " regres " << std::endl;
//	res.print();
//	std::cout << std::endl;
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
			std::cerr << "Error: Illegal character in DNA sequence: " << base << std::endl;
			exit(1);
			break;
		}
}

char pw_alignment::base_translate_back(bool bit1, bool bit2, bool bit3) {

	if(bit1) {
		if(bit2) {
			if(bit3) {
				std::cerr << "Error: wrong bit std::vector 111" << std::endl;
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
				std::cerr << "Error: wrong bit std::vector 011" << std::endl;
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
	std::vector<bool> pw_alignment::getsample1()const{
	return samples.at(0);
}
	std::vector<bool> pw_alignment::getsample2()const{
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
	std::vector<bool> pw_alignment::getsample(size_t id)const{
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

std::string pw_alignment::get_al_ref1() const {
	std::stringstream res;
	for(size_t i=0; i<alignment_length(); ++i) {
		char c1, c2;
		alignment_col(i, c1, c2);
		res << c1;
	}
	return res.str();
}
	
std::string pw_alignment::get_al_ref2() const {
	std::stringstream res;
	for(size_t i=0; i<alignment_length(); ++i) {
		char c1, c2;
		alignment_col(i, c1, c2);
		res << c2;
	}
	return res.str();
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
	std::cout << "al1 seq "<< getreference1() << " b " << getbegin1() << " e " << getend1() <<  " l "  << alignment_length() << std::endl;
	std::cout << "al2 seq "<< getreference2() << " b " << getbegin2() << " e " << getend2() << std::endl;
	/*
	for(size_t col = 0; col < alignment_length(); col++) {
	//	if(col < 3 || (alignment_length() - col < 3)) {
			char c1;
			char c2;
			alignment_col(col, c1, c2);
			std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
	//	}
	}
	*/
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

void pw_alignment::set_cost(const std::vector<double> & create, const std::vector<double> & modify) const{
	assert(create.size()==2 && modify.size()==2);
	create_costs = create;
	modify_costs = modify; 
	costs_cached = true;
}

bool pw_alignment::is_cost_cached() const {
	return costs_cached;
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




        wrapper::wrapper():encodeout("enc.txt"),decodeout("dec.txt"),al_encode("al_encode.txt"),al_decode("al_decode.txt"){
        }
        wrapper::~wrapper(){}
        void wrapper::encode(unsigned int& low, unsigned int& high, unsigned int & total){
                encodeout << "l"<< low << "h"<< high;
                al_encode << " low: "<< low << " high: "<< high <<std::endl;

        }
        void wrapper::decode(unsigned int& low, unsigned int& high, unsigned int & total){
                decodeout << "l"<< low << "h"<< high;
                al_decode << " low: "<< low << " high: "<< high <<std::endl;

        }
        void wrapper::context(size_t & pos , int & context){
                al_encode << "position: " << pos << " context: " <<  context <<std::endl;

        }
        void wrapper::decodeContext(size_t & pos, int & context){
                al_decode << "position: " << pos << " context: " <<  context <<std::endl;
        }



	
