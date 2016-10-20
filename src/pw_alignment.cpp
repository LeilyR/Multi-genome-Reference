#include "pw_alignment.hpp"

pw_alignment::pw_alignment(std::string sample1str, std::string sample2str, size_t sample1_begin, size_t sample2_begin, size_t sample1_end, size_t sample2_end, size_t sample1reference, size_t sample2reference ):  begins(2), ends(2), references(2), create_costs(0), modify_costs(0), costs_cached(false), samples(2) {
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

pw_alignment::pw_alignment(const size_t & begin1, const size_t & begin2, const size_t & end1, const size_t & end2, const size_t & reference1, const size_t & reference2):begins(2), ends(2), references(2), costs_cached(false), create_costs(0), modify_costs(0), samples(0)  {
	begins.at(0) = begin1;
	begins.at(1) = begin2;
	ends.at(0) = end1;
	ends.at(1) = end2;
	references.at(0) = reference1;
	references.at(1) = reference2;
	
}
pw_alignment::pw_alignment(): begins(2), ends(2), references(2), costs_cached(false), create_costs(0), modify_costs(0), samples(2) {}

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
	requires: overlap between *this (ref alref) and other (ref alref) 
	result: is a maximal part of *this which is contained one of the references of other

*/
void pw_alignment::single_ref_intersect(size_t alref, const pw_alignment & other, size_t otherref, pw_alignment & alresult, pw_alignment & otherresult) const {
	size_t intersect_ref;
	size_t all, alr, otl, otr;
	
	if(alref) {
		get_lr2(all, alr);
	} else {
		get_lr1(all, alr);
	}
	if(otherref) {
		other.get_lr2(otl, otr);
	} else {
		other.get_lr1(otl, otr);
	}
	assert( all <= otr && otl <= alr);
	bool alforward = true;
	bool otforward = true;
	if(getbegin(alref) > getend(alref)) {
		alforward = false;
	}
	if(other.getbegin(otherref) > other.getend(otherref)) {
		otforward = false;
	}
	if( (all==otl) && (alr==otr) ) {
		alresult = *this;
		otherresult = other;
		return;
	}

	pw_alignment trash1;
	pw_alignment trash2;

	if(all < otl) {
		if(otr > alr) {
		// alalalalalal
		//     otototototot
		// intersection: otl to alr	
			if(alforward)
				this->split(1 - (bool) alref, otl, trash1, alresult);
			else 
				this->split(1 - (bool) alref, otl, alresult, trash1);
			if(otforward)
				other.split(1 - (bool) otherref, alr + 1, otherresult, trash2);
			else 
				other.split(1 - (bool) otherref, alr + 1, trash2, otherresult);
		} else if (otr < alr) {
		// alalalalalalalalalal
		//     otototototot
		// intersection: otl to otr
			pw_alignment altmp;
			if(alforward) {
				this->split(1 - (bool) alref,  otl, trash1, altmp);
				altmp.split(1 - (bool) alref,  otr + 1, alresult, trash1);
			} else {
				this->split(1 - (bool) alref,  otl, altmp, trash1);
				altmp.split(1 - (bool) alref,  otr + 1, trash1, alresult);
			}	
			otherresult = other;
		} else {
		// alalalalalalal
		//       otototot
		// intersection: otl to alr/otr
			if(alforward)
				this->split(1 - (bool) alref, otl, trash1, alresult);
			else
				this->split(1 - (bool) alref, otl, alresult, trash1);
			otherresult = other;
		}
	} else if (otl < all) {
		if( otr > alr) {
		//       alalala
		// ototototototototot
		// intersection: all to alr
			alresult = *this; 
			pw_alignment ottmp;
			if(otforward) {
				other.split(1 - (bool) otherref, all, trash1, ottmp);
				ottmp.split(1 - (bool) otherref, alr + 1, otherresult, trash1); 
			} else {
				other.split(1 - (bool) otherref, all, ottmp, trash1);
				ottmp.split(1 - (bool) otherref, alr + 1, trash1, otherresult); 
			}
		} else if (otr < alr) {
		//     alalalalalalalal
		// otototototototot
		// intersection: all to alr
			if(alforward)
				this->split(1 - (bool) alref, otr + 1, alresult, trash1);
			else 
				this->split(1 - (bool) alref, otr + 1, trash1, alresult);
			if(otforward)
				other.split(1 - (bool) otherref, all, trash2, otherresult);
			else 
				other.split(1 - (bool) otherref, all, otherresult, trash2);
		} else {
		//   alalalalalalalal
		// ototototototototot
		// intersection: all to alr/otr
			alresult = *this;
			if(otforward)
				other.split(1 - (bool) otherref, all, trash2, otherresult);
			else
				other.split(1 - (bool) otherref, all, otherresult, trash2);
		}
	} else {
		if( otr > alr) {
		// alalalalalalal
		// otototototototototot
		// intersection: all/otl to alr
			alresult = *this;
			if(otforward)
				other.split(1 - (bool) otherref, alr + 1, otherresult, trash2);
			else 
				other.split(1 - (bool) otherref, alr + 1, trash2, otherresult);
		} else {
		// alalalalalalalalalalal
		// otototototototot
		// intersection: all/otl to otr
			if(alforward)
				this->split(1 - (bool) alref, otr + 1, alresult, trash1);
			else 
				this->split(1 - (bool) alref, otr + 1, trash1, alresult);
			otherresult = other;
		}
	}
	
	

}

/*
	returns true if both references go in the same direction
*/
bool pw_alignment::is_same_direction() const {
	if(begins.at(0) < ends.at(0) && begins.at(1) < ends.at(1)) return true;
	if(begins.at(0) > ends.at(0) && begins.at(1) > ends.at(1)) return true;
	return false;
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

void pw_alignment::get_bits(char& base, std::vector<bool> & three_bits){
	bool b1,b2,b3;
	base_translate(base,b1,b2,b3);
	three_bits.push_back(b1);
	three_bits.push_back(b2);
	three_bits.push_back(b3);
//	std::cout<< "b1,b2,b3 "<< b1 << " "<<b2 << " "<< b3 << std::endl;

}
void pw_alignment::get_base(char & base, bool& bit1, bool& bit2, bool& bit3){

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
//	if(alignment_length() > 5){
	size_t counter = 0;
	size_t counter1 = 0;
//	for(size_t col = 0; col < 5; col++) {	
/*	for(size_t col = 0; col < alignment_length(); col++) {
//	//	if(col < 3 || (alignment_length() - col < 3)) {
			char c1;
			char c2;
			alignment_col(col, c1, c2);
			std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
			if(c2 == '-'){
				counter ++;
			}
			if(c1 == '-'){
				counter1 ++;
			}
	//	}
	}*/
//	std::cout << "counter "<<counter << " counter1 "<< counter1 <<std::endl;
//}
	
}

bool compare_pointer_pw_alignment::operator()(const pw_alignment * a, const pw_alignment *b) const {
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


	
	if(a->getreference(asmaller) < b->getreference(bsmaller)) {
		return true;
	} else if(a->getreference(asmaller) > b->getreference(bsmaller)) {
		return false;	
	}
	if(a->getreference(abigger) < b->getreference(bbigger)) {
		return true;
	} else if(a->getreference(abigger) > b->getreference(bbigger)) {
		return false;	
	}
//This is the added part:
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

bool compare_pw_alignment::operator()(const pw_alignment &ar, const pw_alignment &br) const {
// if(ar.getbegin(0)
	const pw_alignment * a = &ar;
	const pw_alignment * b = &br;
	size_t asmaller = 0;//TODO replace it with calling the new one
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


	
	if(a->getreference(asmaller) < b->getreference(bsmaller)) {
		return true;
	} else if(a->getreference(asmaller) > b->getreference(bsmaller)) {
		return false;	
	}
	if(a->getreference(abigger) < b->getreference(bbigger)) {
		return true;
	} else if(a->getreference(abigger) > b->getreference(bbigger)) {
		return false;	
	}
//This is the added part:
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

bool sort_pw_alignment::operator() (const pw_alignment &p1, const pw_alignment &p2)const{
	size_t l1,r1,l2,r2;
	p1.get_lr2(l1,r1);
	p2.get_lr2(l2,r2);
	return (l1 <= l2);
}

bool sort_right_pw_alignment::operator() (const pw_alignment &p1, const pw_alignment &p2)const{
	size_t l1,r1,l2,r2;
	p1.get_lr2(l1,r1);
	p2.get_lr2(l2,r2);
	return (r1 <= r2);
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


bool pw_alignment::equals(const pw_alignment & al) const {
	if(al.getreference1()==references.at(0) && al.getreference2()==references.at(1)) {
		if(
			al.getbegin1() == begins.at(0) &&
			al.getend1() == ends.at(0) &&
			al.getbegin2() == begins.at(1) &&
			al.getend2() == ends.at(1)		
		) {
			return true;
		}
	
	} else if(al.getreference2()==references.at(0) && al.getreference1()==references.at(1)) {
		if(
			al.getbegin1() == begins.at(1) &&
			al.getend1() == ends.at(1) &&
			al.getbegin2() == begins.at(0) &&
			al.getend2() == ends.at(0)	
		) {
			return true;
		} 	
	} 


// when both refs are equal we can flip content (nonflipped case already handled in first condition)
       if( (al.getreference1() == al.getreference2()) &&
           (references.at(0) == references.at(1)) &&
           (references.at(0) == al.getreference1()) ) {
               if(
                       al.getbegin1() == begins.at(0) &&
                       al.getend1() == ends.at(0) &&
                       al.getbegin2() == begins.at(1) &&
                       al.getend2() == ends.at(0) 


               ) return true;


       }


	return false;
}


pw_alignment & pw_alignment::operator=(const pw_alignment & al) {
	if(&al!=this) {
		samples = al.samples;
		begins = al.begins;
		ends = al.ends;
		references = al.references;
		blocks=al.blocks;

		create_costs = al.create_costs;
		modify_costs = al.modify_costs;
		costs_cached = al.costs_cached;
// TODO should this also copy located_alignment s
	}

	return *this;
}
char pw_alignment::complement(char & c){
	switch(c) {
		case 'a':
		case 'A':
		return 'T';
		break;
		case 't':
		case 'T':
		return 'A';
		break;
		case 'c':
		case 'C':
		return 'G';
		break;
		case 'g':
		case 'G':
		return 'C';
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
		return 'N';
		break;
		case '.':
		case '-':
		return '-';
		break;

	}
	return 'X';
}
std::vector<bool> pw_alignment::reverse_complement(std::vector<bool> & sample){
	std::vector<bool> temp;
	for(int i =sample.size()-1; i > 1;i--){
		char s = base_translate_back(sample.at(i-2),sample.at(i-1),sample.at(i));
		char rev = complement(s);
		bool bit1;
		bool bit2;
		bool bit3;
		base_translate(rev, bit1, bit2, bit3);
		temp.push_back(bit1);
		temp.push_back(bit2);
		temp.push_back(bit3);
		i = i-2;
	}
	return temp;
}
void pw_alignment::get_reverse(pw_alignment & al){
	al.references = references;
	al.ends = begins;
	al.begins = ends;
	std::vector< std::vector<bool> >temp = samples;
	al.samples.at(0) = reverse_complement(temp.at(0)); 
	al.samples.at(1) = reverse_complement(temp.at(1)); 
}
void pw_alignment::get_reverse_complement_sample(std::vector<std::vector<bool> > & sample){
	sample.push_back(reverse_complement(samples.at(0)));
	sample.push_back(reverse_complement(samples.at(1)));
}

        wrapper::wrapper():encodeout("enc1.txt"),al_encode("al_encode1.txt"){
        }
        wrapper::~wrapper(){}
        void wrapper::encode(unsigned int& low, unsigned int& high, unsigned int & total){
                encodeout << "l"<< low << "h"<< high;
//              al_encode << " low: "<< low << " high: "<< high <<std::endl;

        }
        void wrapper::context(int & context){
		al_encode<<"c"<<context;
             //   al_encode << "position: " << pos << " context: " <<  context <<std::endl;

        }
       	decoding_wrapper::decoding_wrapper():decodeout("dec1.txt"),al_decode("al_decode1.txt"){}
	decoding_wrapper::~decoding_wrapper(){}
	void decoding_wrapper::decode(unsigned int& low, unsigned int& high, unsigned int & total){
                decodeout << "l"<< low << "h"<< high;
//		std::cout << "high in wrapper " << high <<std::endl;
//              al_decode << " low: "<< low << " high: "<< high <<std::endl;
        }
        void decoding_wrapper::decodeContext(int & context){
		al_decode<<"c"<<context;
	//	std::cout << "context "<< context << std::endl;
            //    al_decode << "position: " << pos << " context: " <<  context <<std::endl;
        }

void pw_alignment::is_block_cached(std::vector<bool> & isbc) const {
	isbc.clear();
}


void located_alignment::is_block_cached(std::vector<bool> & isbc) const {
	isbc = block_cached;
}

void pw_alignment::set_block_cache(const size_t & bnum, const double & c1, const double & c2, const double & m1, const double & m2) const {
	assert(0);	
}



void located_alignment::make_blocks() const {
	size_t len = alignment_length();
	size_t num_blocks = (len/ALIGNMENT_BLOCK_SIZE) + 1;

	for(size_t i=0; i<num_blocks; ++i) {
		alblock bl;
		bl.cfrom = i*ALIGNMENT_BLOCK_SIZE;
		bl.cto = bl.cfrom + ALIGNMENT_BLOCK_SIZE - 1;
		if(bl.cto >= alignment_length()) {
			bl.cto = alignment_length() - 1;
		} 
		blocks.push_back(bl);
	}



	block_cached = std::vector<bool>(num_blocks, false);

}
located_alignment::~located_alignment() {}

located_alignment::located_alignment(const located_alignment & p) {
	data = p.data;
	length = p.length;
	last_segment = p.last_segment;
	segments = p.segments;
	blocks = p.blocks;
	block_cached = p.block_cached;

	begins = p.begins;
	ends = p.ends;
	references = p.references;
	costs_cached = p.costs_cached;
       	create_costs = p.create_costs;
       	modify_costs = p.modify_costs;



}



/*
	create from pairwise alignment. la is always without gaps at ends
*/
located_alignment::located_alignment(const all_data * data, const pw_alignment & p): data(data), length(p.alignment_length()), last_segment(0), pw_alignment(p.getbegin1(), p.getbegin2(), p.getend1(), p.getend2(), p.getreference1(), p.getreference2()) {
// separate content of p into segments without gap
	size_t last_begin1 = std::numeric_limits<size_t>::max();
	size_t last_begin2 = std::numeric_limits<size_t>::max();
	size_t last_begin3 = std::numeric_limits<size_t>::max();
	size_t r1_at = begins.at(0); // position of the current character on its reference (or previous non-gap character)
	size_t r2_at = begins.at(1);
	bool r1forward = true;
	bool r2forward = true;
	if(begins.at(0) > ends.at(0)) r1forward = false;
	if(begins.at(1) > ends.at(1)) r2forward = false;
	size_t in_seg_l = 0;
	size_t r1_gap_l = 0;
	size_t r2_gap_l = 0;

//	std::cout << " r fr al  len " << length << std::endl;
	for(size_t i=0; i<p.alignment_length(); ++i) {
		char c1;
		char c2;
		p.alignment_col(i, c1, c2);
//		std::cout << " i " << i << " al " << c1 << " " << c2 <<" at " << r1_at << " " << r2_at<< " l " << in_seg_l << " g1 " << r1_gap_l << " g2 " << r2_gap_l << std::endl;


// within loop order: first all things that end (have ended at end of previous loop. c1,c2 is first column after)
// either end a segment or continue it

		if(r1_gap_l != 0) {
			assert(r2_gap_l==0 && in_seg_l == 0);
			if(c1!='-') {
// end gap r1
				r1_gap_l = 0;
			} else {
				r1_gap_l++;
			}


		}
		if(r2_gap_l != 0) {
			assert(r1_gap_l==0 && in_seg_l ==0);
			if(c2!='-') {
// end gap r2
				r2_gap_l = 0;
			} else {
				r2_gap_l++;
			}


		}
		if(in_seg_l !=0) {
			assert(r1_gap_l == 0 && r2_gap_l == 0);
			if(c1=='-' || c2=='-') {
// end segment:		
				assert(last_begin1 < std::numeric_limits<size_t>::max() && last_begin2 < std::numeric_limits<size_t>::max());
				segment s;
				s.begins[0] = last_begin1;
				s.begins[1] = last_begin2;
				s.begins[2] = last_begin3;
				s.length = in_seg_l;
//				std::cout << " end seg " << " c " << s.begins[2] << " l " << s.length << " to " << s.begins[2] + s.length -1 << " r1 " << s.begins[0] << " r2 " << s.begins[1] << std::endl;
				segments.push_back(s);
				in_seg_l = 0;
			} else {
				in_seg_l++;
			}
		} 


// then all segments that begin with i
		if(in_seg_l==0) {		
			if(c1!='-' && c2!='-') {
				assert(r1_gap_l==0 || r2_gap_l==0);
//				std::cout << " begin seg at " << i << " r1 " << r1_at << " r2 " << r2_at << std::endl;
// begin segment 
				last_begin1 = r1_at;
				last_begin2 = r2_at;
				last_begin3 = i;
				in_seg_l = 1;
				r1_gap_l = 0;		
				r2_gap_l = 0;
			}
		}

		if(r1_gap_l==0) {
			if(c1=='-') {
// begin r1 gap
				last_begin1 = std::numeric_limits<size_t>::max();
				last_begin2 = std::numeric_limits<size_t>::max();
				last_begin3 = std::numeric_limits<size_t>::max();
				r1_gap_l = 1;

			}
		}

		if(r2_gap_l==0) {
			if(c2=='-') {
// begin r2 gap
				last_begin1 = std::numeric_limits<size_t>::max();
				last_begin2 = std::numeric_limits<size_t>::max();
				last_begin3 = std::numeric_limits<size_t>::max();
				r2_gap_l = 1;

			}
		}

// continue to position in next loop
		if(c1!='-') {
			if(r1forward) r1_at++;
			else r1_at--;
		}
		if(c2!='-') {
			if(r2forward) r2_at++;
			else r2_at--;
		}
	} // for i


	if(in_seg_l > 0) {
		assert(last_begin1 < std::numeric_limits<size_t>::max() && last_begin2 < std::numeric_limits<size_t>::max());
		segment s;
		s.begins[0] = last_begin1;
		s.begins[1] = last_begin2;
		s.begins[2] = last_begin3;
		s.length = in_seg_l;
		segments.push_back(s);
		in_seg_l = 0;
	}


	for(size_t i=0; i<segments.size(); ++i) {
		const segment & s = segments.at(i);
//		std::cout << " s " << i << " c " << s.begins[2] << " l " << s.length << " to " << s.begins[2] + s.length -1 << " r1 " << s.begins[0] << " r2 " << s.begins[1] << std::endl;

	}

	make_blocks();
}


size_t located_alignment::alignment_length() const {
	return length;
}
void located_alignment::alignment_col(size_t c, char & s1, char & s2) const {
// TODO 
	size_t r1, r2;
	alignment_col(c, s1, s2, r1, r2);
}


void located_alignment::alignment_col(size_t c, char & c1, char & c2, size_t & r1at, size_t & r2at) const {
	size_t fr = 0;
	size_t to = segments.size() - 1;

	size_t seg;
	c1 = 'x'; // bsearch never sets to x
	// to speed up sequential access, we try first on the segment that was previously used
//	std::cout << " al col initial on " << last_segment << std::endl;
	bsearch_col(c, last_segment, last_segment, seg, r1at, r2at, c1, c2);
	if(c1=='x') {
//		std::cout << " al col binary search " << std::endl;
	// then we do full binary search
		bsearch_col(c, fr, to, seg, r1at, r2at, c1, c2);
	}
	last_segment = seg;

}


/*
	search for alignment column

	from segment from to segment to

	return col is contained in segments.at(seg) or the gap afterwards

*/
void located_alignment::bsearch_col(const size_t & col, const size_t & from, const size_t & to, size_t & seg, size_t & r1at, size_t & r2at, char & c1, char & c2) const{
	size_t dist = to - from;
	size_t midpos = from + (dist/2);


//	std::cout << " al column binary search from seg " << from << " to " << to << " at " << midpos << " for " << col <<std::endl;

	const segment & mseg = segments.at(midpos);
	size_t mpcol_fr = mseg.begins[2];
	size_t mpcol_to = mpcol_fr + mseg.length - 1;

	if(col >= mpcol_fr && col <= mpcol_to) {
// col is in gapless segment
		
//		std::cout << " in segment " << std::endl;
		
		seg = midpos;
		size_t offs = col - mpcol_fr;
		bool r1forward = true;
		bool r2forward = true;
		if(begins.at(0) > ends.at(0)) r1forward = false;
		if(begins.at(1) > ends.at(1)) r2forward = false;
		const dnastring & s1 = data->getSequence(references.at(0));
		const dnastring & s2 = data->getSequence(references.at(1));
		if(r1forward) {
			r1at = mseg.begins[0] + offs;
			c1 = s1.at(r1at);
		} else {
			r1at = mseg.begins[0] - offs;
			c1 = dnastring::complement(s1.at(r1at));
		}
		if(r2forward) {
			r2at = mseg.begins[1] + offs;
			c2 = s2.at(r2at);
		} else {
			r2at = mseg.begins[1] - offs;
			c2 = dnastring::complement(s2.at(r2at));
		}



	} else if(col > mpcol_to) {
// now find out if col is in the gaps after the midpos segment
		if(midpos + 1 < segments.size()) { // are there any gaps after midpos


//			std::cout << " check gaps after segment " << std::endl;

			bool r1forward = true;
			bool r2forward = true;
			if(begins.at(0) > ends.at(0)) r1forward = false;
			if(begins.at(1) > ends.at(1)) r2forward = false;
			const segment & nseg = segments.at(midpos + 1);
//			std::cout << " between " << midpos << " and " << midpos+1<< std::endl;
//			std::cout << " col " << col << " seg " << mpcol_fr << " " << mpcol_to << " next " << nseg.begins[2] << " l " << nseg.length << std::endl;
//			std::cout << " r1 " << r1forward << " seg beg " << mseg.begins[0] << " l " << mseg.length << " "<< mseg.begins[0] + mseg.length -1  << " n " << nseg.begins[0] << std::endl;
//			std::cout << " r2 " << r2forward << " seg beg " << mseg.begins[1] << " l " << mseg.length << " "<< mseg.begins[1] + mseg.length -1 << " n " << nseg.begins[1] << std::endl;
			size_t r1gaps, r2gaps; // gapped bases on r1/r2, the gaps are on the other reference
			if(r1forward) {
				r1gaps = nseg.begins[0] - mseg.length - mseg.begins[0];
			} else {
				r1gaps = mseg.begins[0] - mseg.length - nseg.begins[0];
			}
			if(r2forward) {
				r2gaps = nseg.begins[1] - mseg.length - mseg.begins[1];
			} else {
				r2gaps = mseg.begins[1] - mseg.length - nseg.begins[1];
			}
			assert(r1gaps == 0 || r2gaps==0);

//			std::cout << " gaps r1 " << r1gaps << " r2 " << r2gaps << std::endl;
			size_t offs = col - mpcol_to;
			if(r1gaps > 0) {
				if(col <= mpcol_to + r1gaps) { 
					seg = midpos;
					const dnastring & s1 = data->getSequence(references.at(0));


					if(r1forward) {
						r1at = mseg.begins[0] + mseg.length - 1 + offs;
						c1 = s1.at(r1at);
					} else {
						r1at = mseg.begins[0] + 1 - mseg.length - offs; 
						c1 = dnastring::complement(s1.at(r1at));
						
					}
					c2 = '-';
					if(r2forward) {
						r2at = mseg.begins[1] + mseg.length - 1;
					} else {
						r2at = mseg.begins[1] + 1 - mseg.length - r1gaps; 
					}


//					std::cout << " found (in r1) " << r1at << " " << r2at << " chars " << c1 << " " << c2 << std::endl;
					return;
				}	
			}
			if(r2gaps > 0) {
				if(col <= mpcol_to + r2gaps) { 
					seg = midpos;
					const dnastring & s2 = data->getSequence(references.at(1));
					if(r1forward) {
						r1at = mseg.begins[0] + mseg.length - 1;
					} else {
						r1at = mseg.begins[0] + 1 - mseg.length; 
					}
					c1 = '-';
					if(r2forward) {
						r2at = mseg.begins[1] + mseg.length - 1 + offs;
						c2 = s2.at(r2at);
					} else {
						r2at = mseg.begins[1] + 1 - mseg.length - offs; 
						c2 = dnastring::complement(s2.at(r2at));
					}
//					std::cout << " found (in r2) " << r1at << " " << r2at << " chars " << c1 << " " << c2 << std::endl;
					return;
				}	
			}	
			// if col was not found in gaps after mseg, recursion goes on
		}
		if(midpos < to) {

//			std::cout << " continue recursion after " << std::endl;
			bsearch_col(col , midpos + 1, to, seg, r1at, r2at, c1, c2); 
		}
	} else if(col < mpcol_fr) {
		if(from < midpos) {
//			std::cout << " continue recursion before " << std::endl;
			bsearch_col(col, from, midpos - 1, seg, r1at, r2at, c1, c2);
		}
	}	
}
/*
	search for alignment position on reference ref at position refpos

	from segment from to segment to

	return col is contained in segments.at(seg) or the gap afterwards

*/

void located_alignment::bsearch_ref(const size_t & ref, const size_t & refpos, const size_t & from, const size_t & to, size_t & seg, size_t & alcolumn) const {
	assert(ref < 2);

	// look at segment in the center of search area
	size_t dist = to - from;
	size_t midpos = from + (dist/2);


	bool refforward = true;
	if(begins.at(ref) > ends.at(ref)) refforward = false;

	const segment & mseg = segments.at(midpos);
	size_t mref_l, mref_r;
	if(refforward) {
		mref_l = mseg.begins[ref];
		mref_r = mref_l + mseg.length -1;
	} else {
		mref_r = mseg.begins[ref];
		mref_l = mref_r - mseg.length +1;
	}

// std::cout << " search for ref " << ref << " pos " << refpos << " on " << midpos << " from " << mref_l << " to " << mref_r << " forward " << refforward << std::endl;

	if(mref_l <= refpos && refpos <= mref_r) {
// refpos on the segment
		seg = midpos;
		size_t offs = refpos - mref_l;
		if(refforward) {
			alcolumn = mseg.begins[2] + offs;
		} else {
			alcolumn = mseg.begins[2] + mseg.length - 1 - offs;
		}
	} else if( refpos > mref_r) {

		size_t bases_after_mr = refpos - mref_r; // how far to look after mseg (in bases on ref) 
//	std::cout << " pos is " << bases_after_mr << " bases after seg " << std::endl;
// is it in the gap after mseg?
		if(refforward) {
			if(midpos + 1 < segments.size()) {
				const segment & nseg = segments.at(midpos + 1);
// are there bases on ref between mseg and nseg?
				size_t num_gapb = nseg.begins[ref] - mref_r -1;
//	std::cout << " there are " <<num_gapb << " forward gapped bases " << std::endl;
	
				if(bases_after_mr <= num_gapb) {
					seg = midpos;
					alcolumn = mseg.begins[2] + mseg.length - 1 + bases_after_mr;
					return;
				} 
			}
		} else {
// if ref is backwards, positions after refpos, come before it on the alignment
			if(midpos > 0) {
				const segment & nseg = segments.at(midpos - 1);	
				size_t nseg_left = nseg.begins[ref] - nseg.length + 1;
				size_t num_gapb = nseg_left - mref_r - 1;
				if(bases_after_mr <= num_gapb) {
					seg = midpos - 1; // now on alignment before mseg, after nseg
					alcolumn = mseg.begins[2] - bases_after_mr;
					return;
				}
			}

		}
// if it was not in gaps after mseg (on ref) we recurse	
		if(refforward) {
			if(midpos < to) {
				bsearch_ref(ref, refpos, midpos + 1, to, seg, alcolumn);
			}	
		} else {
			if(from < midpos) {
				bsearch_ref(ref, refpos, from, midpos - 1, seg, alcolumn);

			}
		}
	} else {
		if(refforward) {
			if(from < midpos) {
				bsearch_ref(ref, refpos, from, midpos - 1, seg, alcolumn);
			}
		} else {
			if(midpos < to) {
				bsearch_ref(ref, refpos, midpos +1, to, seg, alcolumn);
			}
		}
	}
}

void located_alignment::get_column(const size_t & ref, const size_t & refpos, size_t & alcol) const {
	assert(ref < 2);
	alcol = std::numeric_limits<size_t>::max();
	size_t segment;
//	std::cout << " initial search at " << last_segment << std::endl;
	bsearch_ref(ref, refpos, last_segment, last_segment, segment, alcol);
	if(alcol == std::numeric_limits<size_t>::max()) {
//		std::cout << " binary search " << std::endl;
		bsearch_ref(ref, refpos, 0, segments.size(), segment, alcol);
	}
	last_segment = segment;
}



void located_alignment::alignment_position(size_t col, size_t & r1pos, size_t & r2pos) const {
	char c1, c2;
	size_t seg;
	size_t from = 0;
	size_t to = segments.size() - 1;
	bsearch_col(col, from, to, seg, r1pos, r2pos, c1, c2);



}

void pw_alignment::alignment_position(size_t col, size_t & r1pos, size_t & r2pos) const {
	assert(0);
}

std::vector<alblock> & pw_alignment::getblocks() const {
	return blocks;
}

void located_alignment::set_block_cached(const size_t & bl) const {
	block_cached.at(bl) = true;
}
void pw_alignment::set_block_cached(const size_t & bl) const {
	assert(0);
}
	
void located_alignment::printseg() const {
	bool r1forward = true;
	bool r2forward = true;
	if(begins.at(0) > ends.at(0)) r1forward = false;
	if(begins.at(1) > ends.at(1)) r2forward = false;

	for(size_t i=0; i<segments.size(); ++i) {
		segment s = segments.at(i);
		size_t r1to, r2to;
		if(r1forward) {
			r1to = s.begins[0] + s.length - 1;
		} else {
			r1to = s.begins[0] - s.length + 1;
		}
		if(r2forward) {
			r2to = s.begins[1] + s.length - 1;
		} else {
			r2to = s.begins[1] - s.length + 1;
		}
		
		
//		std::cout << " seg " << i << " r1 " << s.begins[0] << " to " << r1to << " r2 " << s.begins[1] << " to " << r2to << " col " << s.begins[2] << " to " << s.begins[2] + s.length -1 << " len " << s.length << std::endl;
		if(i+1 < segments.size()) {
			segment n = segments.at(i+1);
			size_t gc = n.begins[2] - s.begins[2] - s.length;
			size_t g1, g2;
			if(r1forward) {
				g1 = n.begins[0] - s.begins[0] - s.length;
			} else {
				g1 = s.begins[0] - n.begins[0] - n.length;
			}
			if(r2forward) {
				g2 = n.begins[1] - s.begins[1] - s.length;
			} else {
				g2 = s.begins[1] - n.begins[1] - n.length;
			}
//		std::cout << " gaps   r1 " << g1 << " r2 " << g2 << " col " << gc << std::endl; 
		}  

	}

}

bool pw_alignment::onlyGapSample() const {
               bool gapOnSample1 = true;
               bool gapOnSample2 = true;
               for(size_t i = 0; i< alignment_length();i++){
                       char p1char = 'X';
                       char p2char = 'X';
                       alignment_col(i, p1char, p2char);       
                       if(p1char !='-' ) gapOnSample1 = false;
                       if(p2char !='-' ) gapOnSample2 = false;                 
                       if(! gapOnSample1 && !gapOnSample2) {
                               return false;
                       }               
               }
               return gapOnSample1 || gapOnSample2;
}

/* old code compatibility function, this is slow
-       for fast splitting/simulations use located_alignment
-*/
void pw_alignment::simulate_split(const size_t & al_ref, const size_t & position, size_t & sp1, size_t & sp2) const {
       pw_alignment fp;
       pw_alignment sp;
       bool alref_for_split = 1 - (bool) al_ref; // even more confusing: split splits ref 2 if the flag is false
       split(alref_for_split, position, fp, sp);
       pw_alignment fpe;
       pw_alignment spe;


       size_t l, r; // positions on the reference with the split point
       bool alforward = true; // is al in forward direction on al_other_seq?
       size_t al_ref_seq, al_other_seq;
       if(al_ref) {
               get_lr2(l,r);
               if(getbegin1() > getend1()) alforward = false;
               al_ref_seq = getreference2();
               al_other_seq = getreference1();
               
       } else {
               get_lr1(l,r);
               if(getbegin2() > getend2()) alforward = false;
               al_ref_seq = getreference1();
               al_other_seq = getreference2();
       }

       // do splitted parts have only gaps on one of their reference sequences
       bool fgaps = fp.onlyGapSample();
       bool sgaps = sp.onlyGapSample();
                                       
       // find split part alignment ends on the other reference after removing end gaps
       size_t fpeleft= (size_t) -1;
       size_t fperight = (size_t) -1;
       size_t speleft = (size_t) -1;
       size_t speright = (size_t) -1;
       if(!fgaps) {
               fp.remove_end_gaps(fpe);
               if(al_ref)
                       fpe.get_lr1(fpeleft, fperight);
               else
                       fpe.get_lr2(fpeleft, fperight);
#if SPLITPRINT
                       std::cout << "newal fpe " << std::endl;
                       fpe.print();
                       std::cout  << std::endl;
#endif
               }
               if(!sgaps) {
                       sp.remove_end_gaps(spe);
                       if(al_ref)
                               spe.get_lr1(speleft, speright);
                       else 
                               spe.get_lr2(speleft, speright);
#if SPLITPRINT
                       std::cout << "newal spe " << std::endl;
                       spe.print();
                       std::cout << std::endl;
#endif
               }

       



}


/*
	ref:0/1

	pos is on ref 0 or ref 1 of this alignment

	splitting the alignment before column pos, will cause split(s) on the other reference of this alignment. 
	
	if it causes splits, we will set sp1/sp2
	
	we set sp1/sp2 to numeric_limits<size_t>::max if fewer split points were caused
	

*/
void located_alignment::simulate_split(const size_t & ref, const size_t & pos, size_t & sp1, size_t & sp2) const {
// splitpoint at pos means we split left of pos on the sequence
	assert(ref < 2);
	size_t l, r;
	if(ref == 0) {
		get_lr1(l,r);
	}
	if(ref == 1) {
		get_lr2(l,r);
	}

	assert(l<=pos);
	assert(pos<=r);

	size_t split_col;
	get_column(ref, pos, split_col);
	size_t r1pos, r2pos;
	char c1, c2;
	alignment_col(split_col, c1, c2, r1pos, r2pos);
	size_t segid = last_segment;
	const segment & seg = segments.at(segid);

//     std::cout << " simulate split at " << split_col << " r1 " << r1pos << " r2 " << r2pos << " content " << c1 << " " << c2 << std::endl;
//     printseg();

	
	bool rforward = true;
	bool oforward = true;
	if(begins.at(ref) > ends.at(ref)) rforward = false;
	if(begins.at(1-ref) > ends.at(1-ref)) oforward = false;
	
	size_t seg_l, seg_r;
	if(rforward) {
		seg_l = seg.begins[ref];
		seg_r = seg_l + seg.length - 1;
	} else {
		seg_r = seg.begins[ref];
		seg_l = seg_r - seg.length + 1;
	}
//  	std::cout << " split point on ref " << pos << " close to segment l " << seg_l << " r " << seg_r << std::endl;	

	if(seg_l < pos && pos <= seg_r) { // we would cut within seg
		// cut before this pos within the segment (relative to segment alignment)
		size_t cut_col;
		if(rforward) {
			cut_col = pos - seg.begins[ref];	
		} else {
			cut_col = seg.begins[ref] - pos + 1;
		}
//		std::cout << " cut in that segment before column  " << cut_col<< std::endl;
		assert(cut_col > 0 && cut_col < seg.length);
		// transfer cut point to other ref
		if(oforward) {
			sp1 = seg.begins[1-ref] + cut_col;
			sp2 = std::numeric_limits<size_t>::max();
		} else {
			sp1 = seg.begins[1-ref] - cut_col + 1;
//			sp1 = seg.begins[1-ref] + seg.length - cut_col;
			sp2 = std::numeric_limits<size_t>::max();
		}
		return;
	} 
	// cut in gaps between segf and segs (first and second in alignment space)
	segment segf;
	segment segs;
	if(rforward) {
		if(pos <= seg_l) {
			assert(segid > 0);
			segs = seg;
			segf = segments.at(segid - 1);
		} else {
			assert(segid + 1 < segments.size());
			segf = seg;
			segs = segments.at(segid + 1);
		}
	} else {
		if(pos <= seg_l) {	
			assert(segid + 1 < segments.size());
			segf = seg;
			segs = segments.at(segid + 1);
		} else {
			assert(segid > 0);
			segs = seg;
			segf = segments.at(segid - 1);
		}
	}
// 	std::cout << " cut between " << segf.begins[2] << " and " << segs.begins[2] << std::endl;
	if(oforward) {
		sp1 = segf.begins[1-ref] + segf.length; // just after segf
		sp2 = segs.begins[1-ref]; // before begin of segs
	} else {
		sp1 = segf.begins[1-ref] - segf.length + 1; // just left of f
                sp2 = segs.begins[1-ref] + 1; // just right of s

	}


	


}

void located_alignment::split(bool sample, size_t pos, located_alignment & first_part, located_alignment & second_part) const {
       size_t ref = 1 - (size_t) sample;
       assert(ref < 2);
	assert(0);
}
