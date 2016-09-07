
#include "overlap.hpp"


overlap_interval_tree::overlap_interval_tree(const all_data & d): data(d), als_on_reference(d.numSequences()),alind(d.numSequences()){

	}

overlap_interval_tree::overlap_interval_tree(const overlap_interval_tree & o): data(o.data),als_on_reference(o.data.numSequences()),alind(o.data.numSequences()){}


	
overlap_interval_tree::~overlap_interval_tree(){
}

void overlap_interval_tree::remove_alignment(const pw_alignment & removeR){
//	std::cout << "removeR"<<std::endl;
	std::set<pw_alignment, compare_pw_alignment>::iterator findr = alignments.find(removeR);
//	assert(findr != alignments.end()); XXX
	if(findr != alignments.end()){
		const pw_alignment & al = *findr;
		alind.erase(&al);
		alignments.erase(removeR);
	}
}


void overlap_interval_tree::test_all() const {

	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment & al = *it;
		bool alok = data.alignment_fits_ref(al);
		if(!alok) {
			std::cout << " test all failed " << std::endl;
			al.print();
			exit(1);
		}
		
	
	}
}




void overlap_interval_tree::print_all_alignment() const {
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
		const pw_alignment * al = &(*it);
		//std::cout << " x " << al << std::endl;
		al->print();
		std::cout << std::endl;
	}

}

void overlap_interval_tree::test_multimaps() {

	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it ) {
		const pw_alignment * remove = &(*it);
		std::pair< std::multimap<size_t, const pw_alignment &>::iterator, std::multimap<size_t, const pw_alignment &>::iterator > eqrb1 =
		als_on_reference.at(remove->getreference1()).equal_range(remove->getbegin1());
		size_t found =0;
		for(std::multimap<size_t, const pw_alignment &>::iterator it = eqrb1.first; it!=eqrb1.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}
		std::pair< std::multimap<size_t,const pw_alignment &>::iterator, std::multimap<size_t, const pw_alignment &>::iterator > eqre1 =
		als_on_reference.at(remove->getreference1()).equal_range(remove->getend1());
		for(std::multimap<size_t, const pw_alignment &>::iterator it = eqre1.first; it!=eqre1.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}

		std::pair< std::multimap<size_t, const pw_alignment &>::iterator, std::multimap<size_t, const pw_alignment &>::iterator > eqrb2 =
		als_on_reference.at(remove->getreference2()).equal_range(remove->getbegin2());
		for(std::multimap<size_t, const pw_alignment &>::iterator it = eqrb2.first; it!=eqrb2.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}

		std::pair< std::multimap<size_t, const pw_alignment &>::iterator, std::multimap<size_t, const pw_alignment &>::iterator > eqre2 =
		als_on_reference.at(remove->getreference2()).equal_range(remove->getend2());
		for(std::multimap<size_t, const pw_alignment &>::iterator it = eqre2.first; it!=eqre2.second; ++it) {
			if ( remove->equals(it->second) ) {
				found++;
				break;
			}		
		}
		assert(found==4);

	}
	for(size_t s=0; s<data.numSequences(); s++) {
		for(std::multimap< size_t, const pw_alignment & >::const_iterator it = als_on_reference.at(s).begin(); it!=als_on_reference.at(s).end(); ++it) {
			const pw_alignment & p = it->second;
	
			//std::cout << p << std::endl;
			std::set<pw_alignment, compare_pw_alignment>::iterator findp = alignments.find(p);
			if(findp==alignments.end()){
				std::cout<<"wrong one: "<<std::endl;
				p.print();
			}
			assert(findp!=alignments.end());
		}

	}
		
}



void overlap_interval_tree::insert_without_partial_overlap(const pw_alignment & p){
//	std::cout << "insert_without_partial_overlap"<<std::endl;
//	alind.insert(&p);
//	std::cout << &p << std::endl;
//	p.print();
	assert(alignments.find(p) == alignments.end());
	std::pair<std::set<pw_alignment, compare_pw_alignment>::iterator, bool > npp = alignments.insert(p);
	std::set<pw_alignment, compare_pw_alignment>::iterator it = alignments.find(p);
	assert(it != alignments.end());
	const pw_alignment & al = *it;
//	std::cout << "alignment size is "<< alignments.size() << std::endl;
//	std::cout << &al <<std::endl;
	alind.insert(&al);
}

const pw_alignment * overlap_interval_tree::get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const {
/*		const std::multimap<size_t, const pw_alignment &> & r1map = als_on_reference.at(ref1);
		std::pair<std::multimap<size_t, const pw_alignment &>::const_iterator,std::multimap<size_t, const pw_alignment &>::const_iterator> eqr = r1map.equal_range(left1);
		for( std::multimap<size_t, const pw_alignment &>::const_iterator it = eqr.first; it!= eqr.second; ++it) {
			const pw_alignment & calr = it->second;
			const pw_alignment * cal = & calr;
			if(cal->getreference1()==ref1) {
				size_t cal_left1 = cal->getbegin1();
				if(cal->getend1() < cal_left1) cal_left1 = cal->getend1();
				if(cal_left1!=left1) continue;
				if(ref2!=cal->getreference2()) continue;
				size_t cal_left2 = cal->getbegin2();
				if(cal->getend2() < cal_left2) cal_left2 = cal->getend2();
				if(cal_left2!=left2) continue;


				return cal;
			} 
			else { // ref1 is reference2 of cal
				size_t cal_left1 = cal->getbegin1();
				if(cal->getend1() < cal_left1) cal_left1 = cal->getend1();
				if(cal_left1!=left2) continue;
				if(ref2!=cal->getreference1()) continue;	
				size_t cal_left2 = cal->getbegin2();
				if(cal->getend2() < cal_left2) cal_left2 = cal->getend2();
				if(cal_left2!=left1) continue;
				return cal;
			}
		}*/
		return NULL;
}

std::multimap<size_t, const pw_alignment &> &  overlap_interval_tree::get_als_on_reference(size_t sequence)  {
	return als_on_reference.at(sequence);
}
const std::multimap<size_t, const pw_alignment &> &  overlap_interval_tree::get_als_on_reference_const(size_t sequence)const{
	return als_on_reference.at(sequence);
}
void overlap_interval_tree::get_interval_result_on_interval( size_t& sequence,  size_t& left,  size_t& right,std::vector<const pw_alignment *>& RESULT)const{
//	alind.debug_print();
	std::vector<const pw_alignment *> result;
	alind.search_overlap(sequence,left,right,result);
//	std::cout<< "overlapped als: "<<std::endl;
	for(size_t i =0; i < result.size();i++){
		const pw_alignment* al = result.at(i);
	//	al->print();
		RESULT.push_back(al);
	}
//	return result;
}
void overlap_interval_tree::interval_result_on_point(size_t sequence, size_t & split_point, std::vector<const pw_alignment* > & result)const{
	alind.search_overlap(sequence,split_point,result);
}
std::vector<const pw_alignment *>& overlap_interval_tree::get_interval_result_on_point( size_t& sequence,  size_t& split_point)const{
	std::vector<const pw_alignment *> result;
	alind.search_overlap(sequence,split_point,result);
	return result; // TODO XXX this is wrong
}

// TODO why not used
void overlap_interval_tree::test_all_part()const{
	
	/*std::set < const pw_alignment*, compare_pw_alignment> alignment_parts;
		
		for(size_t i =0 ; i< data.numAlignments(); ++i){
		const pw_alignment & al = data.getAlignment(i);
		size_t al_right1 = al.getend1();
		size_t al_right2 = al.getend2();
		size_t al_left1 = al.getbegin1();
		size_t al_left2 = al.getbegin2();

		if(al_right1 < al_left1){al_right1 = al.getbegin1();
					 al_left1 = al.getend1();
					}
		if(al_right2 < al_left2){al_right2 = al.getbegin2();
					 al_left2= al.getend2();
					}
	//	size_t cal_right1 = 0;
		size_t cal_left1 = al_left1;
		size_t cal_left2 = al_left2;
		do {
			const pw_alignment * cal = get_al_at_left_end(al.getreference1(),al.getreference2(),cal_left1, cal_left2);
			if(cal!=NULL) {
			size_t cal_right1 = cal -> getend1();
			size_t cal_right2 = cal -> getend2();
			cal_left1 = cal -> getbegin1();
			cal_left2 = cal -> getbegin2();
				if(cal_right1 < cal_left1){cal_right1 = cal -> getbegin1();
					                   cal_left1 = cal -> getend1();
				 		}
				if(cal_right2 < cal_left2){cal_right2 = cal -> getbegin2();
			 			 	   cal_left2= cal -> getend2();
						}
			
			cal_left1 = cal_right1 +1;
			cal_left2 = cal_right2 +1;
			alignment_parts.insert(cal);

			if(cal_right1 > al_right1 || cal_right2 > al_right2) {
				std::cout << "Small alignment part is to long" << std::endl;
				exit(1);
			}
			}
			 else { 
			std::cout<< "There is no alignment!"<<std::endl;
			exit(1);
			}


		} while(cal_left1 <= al_right1);
		
	}
		if(alignment_parts.size()!=alignments.size()) {
		std::cout<< " Small parts don't cover the original one!"<<std::endl;
		exit(1);
		}*/
}

void overlap_interval_tree::test_overlap()const{
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it1 = alignments.begin(); it1!=alignments.end(); ++it1){	
		const pw_alignment & alr = *it1;
		const pw_alignment * al1 = &alr;
		size_t l1ref1 = al1->getbegin1();
		size_t r1ref1 = al1->getend1();
		size_t l1ref2 = al1->getbegin2();
		size_t r1ref2 = al1->getend2();
		if(l1ref1>r1ref1){
			l1ref1 = al1->getend1();
			r1ref1 = al1->getbegin1();		
		}
		if(l1ref2>r1ref2){
			l1ref2 = al1->getend2();
			r1ref2 = al1->getbegin2();		
		}
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it2 = alignments.begin(); it2!=alignments.end(); ++it2){
			const pw_alignment & alr2 = *it2;
			const pw_alignment * al2 = &alr2;
			size_t l2ref1 = al2->getbegin1();
			size_t r2ref1 = al2->getend1();
			size_t l2ref2 = al2->getbegin2();
			size_t r2ref2 = al2->getend2();
			if(l2ref1>r2ref1){
				l2ref1 = al2->getend1();
				r2ref1 = al2->getbegin1();		
			}
			if(l2ref2>r2ref2){
				l2ref2 = al2->getend2();
				r2ref2 = al2->getbegin2();		
			}
			if ((l1ref1 < l2ref1 && l2ref1 < r1ref1) || (l1ref1 < l2ref1 && r2ref1 < r2ref1)){
				std::cout<< "There are some overlap!"<<std::endl;
				exit(1);
			}
			else if((l1ref2 < l2ref2 && l2ref2 < r1ref2) || (l1ref2 < l2ref2 && r2ref2 < r2ref2)){
				std::cout<< "There are some overlap!"<<std::endl;
				exit(1);
			}
			else{
				std::cout<< "There is no overlap" << std::endl;
				exit(0);
			}

		}
	}	
}
void overlap_interval_tree::test_no_overlap_between_ccs(const pw_alignment & p, std::set<const pw_alignment*, compare_pointer_pw_alignment> & ccs)const{
	for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.begin();it != ccs.end();it++){
		const pw_alignment * al = *it;
		size_t l1,r1,l2,r2,pl1,pr1,pl2,pr2;
		al->get_lr1(l1,r1);
		al->get_lr2(l2,r2);
		p.get_lr1(pl1,pr1);
		p.get_lr2(pl2,pr2);
		if(p.getreference1() == al->getreference1()){
			if(pl1 <=r1 && pr1 >= l1){
				std::cout << "overlap1! "<<std::endl;
				al->print();
				exit(1);
			}
		}
		if(p.getreference1() == al->getreference2()){
			if(pl1 <=r2 && pr1 >= l2){
				std::cout << "overlap2! "<<std::endl;
				al->print();
				exit(1);
			}
		}
		if(p.getreference2() == al->getreference1()){
			if(pl2 <=r1 && pr2 >= l1){
				std::cout << "overlap3! "<<std::endl;
				al->print();
				exit(1);
			}
		}
		if(p.getreference2() == al->getreference2()){
			if(pl2 <=r2 && pr2 >= l2){
				std::cout << "overlap4! "<<std::endl;
				al->print();
				exit(1);
			}
		}
	}
}
bool overlap_interval_tree::checkAlignments(const pw_alignment & p)const{
	std::set<pw_alignment, compare_pw_alignment>::iterator it=alignments.find(p);
	if(it!=alignments.end()){
		const pw_alignment al = *it;
	//	al.print();
		return true;	
	}
	else return false;	
}
	


size_t overlap_interval_tree::size() const {

//	for(size_t i=0; i<als_on_reference.size(); i++) {
//		std::cout << " ref " << i << " al start/end points: " << als_on_reference.at(i).size() << std::endl;
//	}
	return alignments.size();
}
		
const std::set<pw_alignment, compare_pw_alignment> & overlap_interval_tree::get_all() const {
	return alignments;
}


// true if true partial overlap
bool overlap_interval_tree::check_po(size_t l1, size_t r1, size_t l2, size_t r2) {
	if(r2 >= l1 && l2 <= r1) {
		if(l1!=l2 || r1!=r2) {
			return true;	
		}
	}
	return false;
}

void overlap_interval_tree::test_partial_overlap_set(std::set< const pw_alignment *, compare_pw_alignment> & als) {
	
	std::vector<const pw_alignment *> all;
	for(std::set<const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
		const pw_alignment * al = *it;
		all.push_back(al);
	}
	test_partial_overlap_vec(all);
}

// TODO ???
void overlap_interval_tree::test_partial_overlap_vec(std::vector< const pw_alignment *> & all) {
	for(size_t i=0; i<all.size(); ++i) {
		for(size_t j=i+1; j<all.size(); ++j) {
			const pw_alignment * a = all.at(i);
			const pw_alignment * b = all.at(j);
			size_t al1, ar1, al2, ar2;
			size_t bl1, br1, bl2, br2;
			a->get_lr1(al1, ar1);
			a->get_lr2(al2, ar2);
			b->get_lr1(bl1, br1);
			b->get_lr2(bl2, br2);
			size_t aref1, aref2, bref1, bref2;
			aref1 = a->getreference1();
			aref2 = a->getreference2();
			bref1 = b->getreference1();
			bref2 = b->getreference2();
			if(aref1 == bref1) {
				if(check_po(al1, ar1, bl1, br1)) {
					std::cout << "partial overlap error 1: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}
			
			if(aref1 == bref2) {
				if(check_po(al1, ar1, bl2, br2)) {
					std::cout << "partial overlap error 2: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

			if(aref2 == bref1) {
				if(check_po(al2, ar2, bl1, br1)) {
					std::cout << "partial overlap error 3: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

			if(aref2 == bref2) {
				if(check_po(al2, ar2, bl2, br2)) {
					std::cout << "partial overlap error 4: " << std::endl;
					a->print();
					b->print();
					throw 0;
					exit(1);
				}
			}

		
		}
	}




}
void overlap_interval_tree::test_partial_overlap() const {
	test_al_on_ref();
	std::vector<const pw_alignment *> all;
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = &(*it);
	//	al->print();
		all.push_back(al);
	}

	test_partial_overlap_vec(all);

}
void overlap_interval_tree::test_al_on_ref()const{
	for(size_t i = 0; i < data.numSequences();i++){
		for(std::multimap< size_t, const pw_alignment &>::const_iterator it = als_on_reference.at(i).begin(); it != als_on_reference.at(i).end();it++){
				check_alignment_address(it->second , & it->second);
		}
	}
}
bool overlap_interval_tree::check_alignment_address(const pw_alignment & p, const pw_alignment * p1)const{
	std::set<pw_alignment, compare_pw_alignment>::const_iterator it=alignments.find(p);
	assert(it != alignments.end());
	const pw_alignment * al = &(*it);
//	std::cout << al << " " << p1 << std::endl;
	if( al == p1){
		return true;
	}
	else return false;	
}
#define SPLITPRINT 0
template<typename overlap_type>
splitpoints_interval_tree<overlap_type>::splitpoints_interval_tree(const pw_alignment & p, overlap_type & o, const all_data & d):overl(o), newal(p),data(d), split_points(d.numSequences()) {

}
template<typename overlap_type>
splitpoints_interval_tree<overlap_type>::~splitpoints_interval_tree() {}
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right) {
	std::vector<const pw_alignment*> result;
	overl.get_interval_result_on_interval(sequence,left,right,result);
#if SPLITPRINT
	std::cout << "result size is "<< result.size() << " on seq " << sequence <<std::endl;
#endif
//	if(result.size()>0){
//	std::cout << "result al is: "<<std::endl;
//	const pw_alignment* p = result.at(0);
//	p->print();
//	std::cout << "done! "<<std::endl;
//	}
	for( size_t i =0; i < result.size();i++){
		const pw_alignment * al = result.at(i);
	//	al->print();
		if(al->getreference1()== sequence){
			
			size_t alleft;
			size_t alright;
			al->get_lr1(alleft, alright);
		//	std::cout << alleft << " "<<alright<<std::endl;
			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left,false);
			}
			if(alleft <= right && right < alright) {
				insert_split_point(sequence, right+1,false);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft,false);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1,false);					
			}
		}
	
		if(al->getreference2() == sequence) {
			size_t alleft;
			size_t alright;
			al->get_lr2(alleft, alright);
		//	std::cout << alleft << " "<<alright << " left " << left <<" right " << right <<std::endl;

			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left,false);			
			}
			if(alleft <= right && right < alright) {
		//		std::cout << "alleft <= right && right < alright"<<std::endl;
				insert_split_point(sequence, right+1,false);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft,false);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1,false);
			}
		}
	}
}
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::find_initial_split_points(size_t sequence, size_t left, size_t right) {
	std::vector<const pw_alignment *> result;
	overl.get_interval_result_on_interval(sequence,left,right,result);
	for(size_t i =0; i < result.size();i++){
		const pw_alignment * al = result.at(i);
		if(al->getreference1() == sequence) {
			size_t alleft;
			size_t alright;
			al->get_lr1(alleft, alright);
			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left);
			}
			if(alleft <= right && right < alright) {
				insert_split_point(sequence, right+1);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1);	
			}
		}
		if(al->getreference2() == sequence) {
			size_t alleft;
			size_t alright;
			al->get_lr2(alleft, alright);
			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left);
			}
			if(alleft <= right && right < alright) {
				insert_split_point(sequence, right+1);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1);
			}
		}
	}

}
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::recursive_splits(){




		size_t left1;
		size_t right1;
		size_t left2 ;
		size_t right2;
		newal.get_lr1(left1, right1);
		newal.get_lr2(left2, right2);

#if SPLITPRINT
		std::cout << " Check ref 1 overlaps " << std::endl;
#endif

		find_initial_split_points(newal.getreference1(), left1, right1);

#if SPLITPRINT
		std::cout << " Check ref 2 overlaps " << std::endl;
#endif
	//	std::cout<< "recursive split"<<std::endl;
		find_initial_split_points(newal.getreference2(), left2, right2);
	//	std::cout << "done "<<std::endl;
		if(newal.getreference1()==newal.getreference2()) {
		//	std::cout<< "ref1 == ref2"<<std::endl;
			if(right1>=left2 && left1 < left2){
			//	std::cout<< "here1"<<std::endl;
				insert_split_point(newal.getreference2(),left2);
				insert_split_point(newal.getreference2(),right1+1);
			}
			if(right2>=left1 && left2 < left1){
			//	std::cout<< "here2"<<std::endl;
				insert_split_point(newal.getreference2(),left1);
				insert_split_point(newal.getreference2(),right2+1);
			}
		}

		insert_split_point(newal.getreference1(), left1);
		insert_split_point(newal.getreference1(), right1+1);
		insert_split_point(newal.getreference2(), left2);
		insert_split_point(newal.getreference2(), right2+1);
	//	std::cout<<"the end!"<<std::endl;
	
}
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::nonrecursive_splits(){


		size_t left1;
		size_t right1;
		size_t left2 ;
		size_t right2;
		newal.get_lr1(left1, right1);
		newal.get_lr2(left2, right2);

#if SPLITPRINT
		std::cout << " Check ref 1 overlaps " << std::endl;
#endif

		find_initial_split_points_nonrecursive(newal.getreference1(), left1, right1);

#if SPLITPRINT
		std::cout << " Check ref 2 overlaps " << std::endl;
#endif
	//	std::cout<< "got here!"<<std::endl;
		find_initial_split_points_nonrecursive(newal.getreference2(), left2, right2);
	
		if(newal.getreference1()==newal.getreference2()) {
			if(right1>=left2 && left1 < left2){
				insert_split_point_nonrecursive(newal.getreference2(),left2);
				insert_split_point_nonrecursive(newal.getreference2(),right1+1);
			}
			if(right2>=left1 && left2 < left1){
				insert_split_point_nonrecursive(newal.getreference2(),left1);
				insert_split_point_nonrecursive(newal.getreference2(),right2+1);
			}
		}

		insert_split_point_nonrecursive(newal.getreference1(), left1);
		insert_split_point_nonrecursive(newal.getreference1(), right1+1);
		insert_split_point_nonrecursive(newal.getreference2(), left2);
		insert_split_point_nonrecursive(newal.getreference2(), right2+1);
	


}

template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::insert_split_point_nonrecursive(size_t sequence, size_t position) {
	std::pair<std::set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	
	
	if(insertes.second) {
#if SPLITPRINT
		std::cout << " new split point " << sequence << " at " << position << std::endl;
#endif	
	}
}



template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::insert_split_point(size_t sequence, size_t position, bool recursion) {
	std::pair<std::set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) {
#if SPLITPRINT
		std::cout << " new split point " << sequence << " at " << position << std::endl;
#endif	
		size_t left1;
		size_t right1;
		newal.get_lr1(left1, right1);
		size_t left2;
		size_t right2;
		newal.get_lr2(left2, right2);




		if(newal.getreference1()==sequence) {
		//	std::cout << "on ref 1"<<std::endl;
			if(left1<position && right1>=position) {
				pw_alignment fp;
				pw_alignment sp;
				newal.split(true, position, fp, sp);

				pw_alignment fpe;
				pw_alignment spe;


				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft2= (size_t) -1;
				size_t fperight2 = (size_t) -1;
				size_t speleft2 = (size_t) -1;
				size_t speright2 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
					std::cout << "newal fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
					std::cout << "newal spe " << std::endl;
					spe.print();
					std::cout << std::endl;
#endif
				}


				if(newal.getbegin2() < newal.getend2()) {
						
			//			std::cout << " try ins " << newal.getreference2() << " : " << spleft2 << std::endl;
					if(!sgaps) {
						if(recursion){
							insert_split_point(newal.getreference2(), speleft2);
						}else{
							insert_split_point_nonrecursive(newal.getreference2(), speleft2);
						}
					}
					// additional split point if gaps at split point
					if(!fgaps) {
						if(recursion){
							insert_split_point(newal.getreference2(), fperight2+1);
						}else{
							insert_split_point_nonrecursive(newal.getreference2(), fperight2+1);
						}
					}
				} else {
			//			std::cout << " try ins " << newal.getreference2() << " : " << fpleft2 << std::endl;
					if(!fgaps) {
						if(recursion){
							insert_split_point(newal.getreference2(), fpeleft2);
						}else{
							insert_split_point_nonrecursive(newal.getreference2(), fpeleft2);
						}
					}
					if(!sgaps) {
						if(recursion){
							insert_split_point(newal.getreference2(), speright2+1);
						}else{
							insert_split_point_nonrecursive(newal.getreference2(), speright2+1);

						}
					}

				}
			}
		}
		if(newal.getreference2()==sequence) {
			if(left2<position && right2>=position) {
			//	std::cout<<"HERE!!"<<std::endl;
				pw_alignment fp;
				pw_alignment sp;
				newal.split(false, position, fp, sp);
				pw_alignment fpe;
				pw_alignment spe;
				
				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				size_t fpeleft1 = (size_t) -1;
				size_t fperight1 = (size_t) -1;
				size_t speleft1 = (size_t) -1;
				size_t speright1 = (size_t) -1;
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
					fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
					std::cout << "newal fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
					spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
					std::cout << "newal spe " << std::endl;
					spe.print();
					std::cout << std::endl;
#endif
				}

				
				if(newal.getbegin1() < newal.getend1()) {	
		//			std::cout << " try ins " << newal.getreference1() << " : " << spleft1 << std::endl;
					if(!sgaps) {
						if(recursion){
							insert_split_point(newal.getreference1(), speleft1);
						}else{
							insert_split_point_nonrecursive(newal.getreference1(), speleft1);
						}
					}
					if(!fgaps) {
						if(recursion){
							insert_split_point(newal.getreference1(), fperight1+1);
						}else{
							insert_split_point_nonrecursive(newal.getreference1(), fperight1+1);
						}
					}
				} else {
			//			std::cout << " try ins " << newal.getreference1() << " : " << fpleft1 << std::endl;
					if(!fgaps) {
						if(recursion){
							insert_split_point(newal.getreference1(), fpeleft1);
						}else{
							insert_split_point_nonrecursive(newal.getreference1(), fpeleft1);
						}
					}
					if(!sgaps) {
						if(recursion){
							insert_split_point(newal.getreference1(), speright1+1);
						}else{
							insert_split_point_nonrecursive(newal.getreference1(), speright1+1);
						}
					}
				}
			}
		}

/*		const boost::icl::interval_map<size_t , std::set<const pw_alignment*> > & alonref = overl.get_alignments_interval(sequence);
		boost::icl::interval_map<size_t, std::set<const pw_alignment*> >::const_iterator interval = alonref.find(position);
	if(interval != alonref.end()){
		std::set<const pw_alignment*> als = (*interval).second;
		for(std::set<const pw_alignment*>::iterator it1 = als.begin(); it1 !=als.end(); it1++){
			const pw_alignment *al = *it1;*/

		const std::multimap<size_t, const pw_alignment &> & alonref = overl.get_als_on_reference_const(sequence);
		for(std::multimap<size_t, const pw_alignment & >::const_iterator it = alonref.lower_bound(position); it!=alonref.end(); ++it) {
			const pw_alignment & alr = it-> second;
			const pw_alignment * al = &alr;
				
				
	//			std::cout << " in ins " << sequence << " : " << position << " see: " << std::endl;
	//			al->print();
	//			std::cout << std::endl;


			if(al->getreference1()!=al->getreference2()) {//TODO add the condition that ref1 == ref2 and two refs have overlap with each other
			//	std::cout << "ref1 != ref2" <<std::endl;
				size_t alleft;
				size_t alright;
				al->get_lr_on_reference(sequence, alleft, alright);
				if(alright < position) {
			//			std::cout << " break " << std::endl;
					break;
				}
				if(alleft == position) {
			//			std::cout << " break " << std::endl;
					break;
				} 
				if(alright+1== position) {
			//			std::cout << " break " << std::endl;
					break;
				}
				if(alleft>position){
			//			std::cout << " break " << std::endl;
					break;
				}
			//		al->print(); 
				//	std::cout<<"position"<<position<<std::endl;
				//	std::cout<<"alleft"<<alleft<<std::endl;
				//	std::cout<<"alright"<<alright<<std::endl;
			//	assert(alleft < position && position <= alright);
				pw_alignment fp;
				pw_alignment sp;
				pw_alignment fpe;
				pw_alignment spe;
				/*	for(size_t col = 0; col < al->alignment_length(); col++) {
						char c1;
						char c2;
						al->alignment_col(col, c1, c2);
						std::cout <<col <<"\t"<< c1<<"\t"<<c2<<std::endl;
					}*/
				//	std::cout<<"al length: "<< al->alignment_length() <<std::endl;
				//	if(al->alignment_length()>1){

				al->split_on_reference(sequence, position, fp, sp);
				


				// do splitted parts have only gaps one of their reference sequences
				bool fgaps = onlyGapSample(fp);
				bool sgaps = onlyGapSample(sp);
					
				// find split part alignment ends on reference after removing end gaps
				if(!fgaps) {
					fp.remove_end_gaps(fpe);
#if SPLITPRINT
					std::cout << "al fpe " << std::endl;
					fpe.print();
					std::cout  << std::endl;
#endif
				}
				if(!sgaps) {
					sp.remove_end_gaps(spe);
#if SPLITPRINT
					std::cout << "al spe " << std::endl;
					spe.print();
					std::cout << std::endl;
#endif
				}




					
				if(sp.getreference1()==sequence) {


					
						if(al->getbegin2() < al->getend2()) {
							if(!sgaps) {
								insert_split_point(sp.getreference2(), spe.getbegin2());
							} 
							if(!sgaps) {
								insert_split_point(sp.getreference2(), fpe.getend2()+1);
							}


						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference2(), fpe.getend2());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference2(), spe.getbegin2()+1);
							}
						}
					} else {




						//	std::cout<<"Heya!!!"<<std::endl;
						if(al->getbegin2() < al->getend2()) {
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1());
							}
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1()+1);
							}
						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1()+1);
							}
						}
					}
				} else {
					
					size_t alleft1;
					size_t alright1;
					al->get_lr1(alleft1, alright1);
					size_t alleft2;
					size_t alright2;
					al->get_lr2(alleft2, alright2);
					if(position > alleft1 && position <= alright1) {
		//				std::cout << "inone" << std::endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(true, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(fp);
					bool sgaps = onlyGapSample(sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft2 = (size_t) -1;
					size_t fperight2 = (size_t) -1;
					size_t speleft2 = (size_t) -1;
					size_t speright2 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr2(fpeleft2, fperight2);
#if SPLITPRINT
						std::cout << "al fpe " << std::endl;
						fpe.print();
						std::cout  << std::endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr2(speleft2, speright2);
#if SPLITPRINT
						std::cout << "al spe " << std::endl;
						spe.print();
						std::cout << std::endl;
#endif
					}


						
							
						if(al->getbegin2() < al->getend2()) {
			//				std::cout << " spleft2 " << spleft2 << std::endl;
							if(!sgaps) {
								insert_split_point(sp.getreference2(), speleft2);
							}
							if(!fgaps) {
								insert_split_point(sp.getreference2(), fperight2+1);
							}
						} else {
			//				std::cout << " fpgetend2 " << fp.getend2() << std::endl;
							if(!fgaps) {
								insert_split_point(fp.getreference2(), fpe.getend2());
							}
							if(!sgaps) {
								insert_split_point(fp.getreference2(), spe.getbegin2()+1);
							}

						}
					}
											
					if(position > alleft2 && position <= alright2) {
			//			std::cout << "intwo" << std::endl;
						pw_alignment fp;
						pw_alignment sp;
						al->split(false, position, fp, sp);
					pw_alignment fpe;
					pw_alignment spe;
					
			// do splitted parts have only gaps one of their reference sequences
					bool fgaps = onlyGapSample(fp);
					bool sgaps = onlyGapSample(sp);
					
					// find split part alignment ends on reference after removing end gaps
					size_t fpeleft1 = (size_t) -1;
					size_t fperight1 = (size_t) -1;
					size_t speleft1 = (size_t) -1;
					size_t speright1 = (size_t) -1;
					if(!fgaps) {
						fp.remove_end_gaps(fpe);
						fpe.get_lr1(fpeleft1, fperight1);
#if SPLITPRINT
						std::cout << "al fpe " << std::endl;
						fpe.print();
						std::cout  << std::endl;
#endif
					}
					if(!sgaps) {
						sp.remove_end_gaps(spe);
						spe.get_lr1(speleft1, speright1);
#if SPLITPRINT
						std::cout << "al spe " << std::endl;
						spe.print();
						std::cout << std::endl;
#endif
					}




//							data.alignment_fits_ref(&fp);
//							data.alignment_fits_ref(&sp);
						
						if(al->getbegin1() < al->getend1()) {
							if(!sgaps) {
								insert_split_point(sp.getreference1(), speleft1);
							}
							if(!fgaps){
								insert_split_point(sp.getreference1(), fperight1+1);
							}
						} else {
							if(!fgaps) {
								insert_split_point(sp.getreference1(), fpe.getend1());
							}
							if(!sgaps) {
								insert_split_point(sp.getreference1(), spe.getbegin1()+1);
							}
						}
					//	}
					}
				}
			}
			
		}	

	}

template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::split_all(std::set<pw_alignment, compare_pw_alignment> & remove_alignments, std::vector<pw_alignment> & insert_alignments ){
		for(size_t i = 0; i <data.numSequences();i++){
#if SPLITPRINT
			std::cout << "current number of the split points on sequence " << i << " is "<< split_points.size() << std::endl;
#endif
			for(std::set<size_t>::iterator split = split_points.at(i).begin(); split!= split_points.at(i).end(); ++split){
				size_t splitPoint = *split;
			//	std::cout << "split point is "<< splitPoint << std::endl;
				std::vector<const pw_alignment *> result;
			//	result = overl.get_interval_result_on_point(i,splitPoint);
				overl.interval_result_on_point(i,splitPoint,result);
			//	std::cout << "result size is "<< result.size() << std::endl;
				for(size_t j =0; j < result.size();j++){
					const pw_alignment * p = result.at(j);
					if(p->getreference1() == i) {
						size_t pleft1;
						size_t pright1;
						p->get_lr1(pleft1,pright1);
						if( pleft1<*split && pright1>=*split){
							remove_alignments.insert(*p);
						}
					}
					if(p->getreference2()==i) {
						size_t pleft2;
						size_t pright2;
						p->get_lr2(pleft2,pright2);
						if( pleft2<*split && pright2>=*split){
							remove_alignments.insert(*p);
						}
					}
				}		
			}
		}
	//	std::cout << " in split all "<< remove_alignments.size() << " remove_alignments " << std::endl;
		//Coordinate of the current alignment will be also added to the split points
		std::vector<pw_alignment>  split_pieces;
		splits(newal, split_pieces);
	//	std::cout << "remove als are: " << split_pieces.size() <<std::endl;
		for(std::set<pw_alignment,compare_pw_alignment>::iterator removed = remove_alignments.begin(); removed != remove_alignments.end(); ++removed){
		//	const pw_alignment & p_test = *removed;
		//	p_test.print();
			splits(*removed, split_pieces);			
		}
		std::set<pw_alignment,compare_pw_alignment> inserted_pieces;
		for(size_t i = 0; i<split_pieces.size();i++) {
	//			split_pieces.at(i).print();
#ifndef NDEBUG
				if(!data.alignment_fits_ref(split_pieces.at(i))) {
				//	std::cout<<"fails here!"<<std::endl;
					exit(1);
				}
#endif
			if(!onlyGapSample(split_pieces.at(i))){	
				pw_alignment noendgaps;
				split_pieces.at(i).remove_end_gaps(noendgaps);
#if SPLITPRINT
				std::cout << "ENDGAPS" << std::endl;
				split_pieces.at(i).print();
				std::cout << "TO " << std::endl;
				noendgaps.print();
				std::cout << std::endl;
#endif
				std::set<pw_alignment,compare_pw_alignment>::const_iterator it = inserted_pieces.find(noendgaps);
			//	noendgaps.print();
				if (overl.checkAlignments(noendgaps)) {
				//	std::cout << "check noendgap "<<std::endl;		
				}
				else if ( it != inserted_pieces.end()){
				//	std::cout << "end inserted pieces "<<std::endl;				
				}
				else {
					insert_alignments.push_back(noendgaps);
					inserted_pieces.insert(noendgaps);
				//	std::cout << "insert to the insert_al "<<std::endl;
				}
			}
		}

		if(insert_alignments.size()==0 && remove_alignments.size()==0){//Just added to add independent als to the set of alignments when they dont split and dont add any ins or del. XXX ??
	//		insert_alignments.push_back(newal);
	//		std::cout << "new_al is ignored" <<std::endl;
		}
		
					
	}
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::splits(const pw_alignment & pr,  std::vector<pw_alignment> & insert_alignments){
		pw_alignment p = pr;
#if SPLITPRINT	
		std::cout << "SPL" << std::endl;
		p.print();
		std::cout << std::endl;

		std::cout << " split on " << p.getreference1() << std::endl;
#endif

		pw_alignment p1;
		pw_alignment p2;
		size_t left1;
		size_t right1;	
//		size_t left2;
//		size_t right2;
		p.get_lr1(left1,right1);
//		p->get_lr2(left2,right2);
		for(std::set<size_t>::iterator splitp = split_points.at(p.getreference1()).upper_bound(left1); splitp!= split_points.at(p.getreference1()).end(); splitp++){
			if(right1>=*splitp){
#if SPLITPRINT
				std::cout << "sp " << *splitp << std::endl;
#endif
				p.split(true,*splitp,p1,p2);
				if(p.getbegin1() < p.getend1()) {
					p = p2;

#if SPLITPRINT
					p1.print();
#endif
					if(!onlyGapSample(p1) && p1.alignment_length() >1 &&p1.getbegin1()!=p1.getend1() && p1.getbegin2()!=p1.getend2()) {
						insert_alignments.push_back(p1);
					}
				} else if(p.getbegin1() > p.getend1()) {
					p = p1;

#if SPLITPRINT
					p2.print();
#endif

					if(!onlyGapSample(p2) && p2.alignment_length()>1 && p2.getbegin1()!=p2.getend1() && p2.getbegin2()!=p2.getend2() ){	
						insert_alignments.push_back(p2);
					}
					
				}
			}
			else break;
		}

#if SPLITPRINT
		std::cout << " last part" << std::endl;
		p.print();
		std::cout << std::endl;
#endif

		if(!onlyGapSample(p)&& p.alignment_length()>1 && p.getbegin1()!=p.getend1() && p.getbegin2()!=p.getend2() ){	
			insert_alignments.push_back(p);		
		}
	}
template<typename overlap_type>
bool splitpoints_interval_tree<overlap_type>::onlyGapSample(const pw_alignment & p) const {
		bool gapOnSample1 = true;
		bool gapOnSample2 = true;
		for(size_t i = 0; i< p.alignment_length();i++){
			char p1char = 'X';
			char p2char = 'X';
			p.alignment_col(i, p1char, p2char);	
			if(p1char !='-' ) gapOnSample1 = false;
			if(p2char !='-' ) gapOnSample2 = false;			
			if(! gapOnSample1 && !gapOnSample2) {
				return false;
			}		
		}
		return gapOnSample1 || gapOnSample2;
	}

/*
	std::vector<pw_alignment>  splitpoints::get_insert()const{
		return insert_alignments;
	}
*/




