
#ifndef OVERLAP_CPP
#define OVERLAP_CPP
#include "overlap.hpp"

overlap_interval_tree::overlap_interval_tree(const all_data & d): data(d), alind(d.numSequences()){

	}

overlap_interval_tree::overlap_interval_tree(const overlap_interval_tree & o): data(o.data),alind(o.data.numSequences()){
	assert(alignments.size() == 0); // incomplete copy constructor, does not copy content of object
}

	
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

void overlap_interval_tree::insert_without_partial_overlap(const pw_alignment & p){
//	std::cout << "insert_without_partial_overlap"<<std::endl;
//	alind.insert(&p);
//	std::cout << &p << std::endl;
//	p.print();

// TODO can we avoid double insertions?
//	assert(alignments.find(p) == alignments.end());
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
/*
std::multimap<size_t, const pw_alignment &> &  overlap_interval_tree::get_als_on_reference(size_t sequence)  {
	return als_on_reference.at(sequence);
}
const std::multimap<size_t, const pw_alignment &> &  overlap_interval_tree::get_als_on_reference_const(size_t sequence)const{
	return als_on_reference.at(sequence);
}
*/

void overlap_interval_tree::get_interval_result_on_interval(const size_t& sequence,  const size_t& left,  const size_t& right,std::vector<const pw_alignment *>& result)const{
//	alind.debug_print();
	alind.search_overlap(sequence,left,right,result);
//	std::cout<< "overlapped als: "<<std::endl;
}
void overlap_interval_tree::interval_result_on_point(const size_t & sequence, const size_t & split_point, std::vector<const pw_alignment* > & result)const{
	alind.search_overlap(sequence,split_point,result);
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
	//test_al_on_ref(); TODO other tests?
	std::vector<const pw_alignment *> all;
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = &(*it);
	//	al->print();
		all.push_back(al);
	}

	test_partial_overlap_vec(all);

}
/*
void overlap_interval_tree::test_al_on_ref()const{
	for(size_t i = 0; i < data.numSequences();i++){
		for(std::multimap< size_t, const pw_alignment &>::const_iterator it = als_on_reference.at(i).begin(); it != als_on_reference.at(i).end();it++){
				check_alignment_address(it->second , & it->second);
		}
	}
}
*/
/*
bool overlap_interval_tree::check_alignment_address(const pw_alignment & p, const pw_alignment * p1)const{
	std::set<pw_alignment, compare_pw_alignment>::const_iterator it=alignments.find(p);
	assert(it != alignments.end());
	const pw_alignment * al = &(*it);
//	std::cout << al << " " << p1 << std::endl;
	if( al == p1){
		return true;
	}
	else return false;	
}*/



#define SPLITPRINT 0
template<typename overlap_type>
splitpoints_interval_tree<overlap_type>::splitpoints_interval_tree(const pw_alignment & p, overlap_type & o, const all_data & d):overl(o), newal(p),data(d), split_points(d.numSequences()), multi_mode(false) {

}

template<typename overlap_type>
splitpoints_interval_tree<overlap_type>::splitpoints_interval_tree(const std::vector<const pw_alignment *> & new_als, overlap_type & o, const all_data & d):overl(o), new_als(new_als), data(d), split_points(d.numSequences()), multi_mode(true) {
//	std::cout << " multimode on " << new_als.size() << std::endl;

}

template<typename overlap_type>
splitpoints_interval_tree<overlap_type>::~splitpoints_interval_tree() {}



template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::find_initial_split_points(size_t sequence, size_t left, size_t right, bool recursion) {


	std::vector<const pw_alignment*> result;
	overl.get_interval_result_on_interval(sequence,left,right,result);
#if SPLITPRINT
	std::cout << " Look for initial split points on " << sequence << " from " << left << " to " << right << std::endl;
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

//		std::cout << " overlapper al " << std::endl;
//		al->print();


		if(al->getreference1()== sequence){
			
			size_t alleft;
			size_t alright;
			al->get_lr1(alleft, alright);
		//	std::cout << alleft << " "<<alright<<std::endl;
			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left, recursion);
			}
			if(alleft <= right && right < alright) {
				insert_split_point(sequence, right+1, recursion);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft, recursion);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1, recursion);					
			}
		}
	
		if(al->getreference2() == sequence) {
			size_t alleft;
			size_t alright;
			al->get_lr2(alleft, alright);
		//	std::cout << alleft << " "<<alright << " left " << left <<" right " << right <<std::endl;

			if(alleft < left && left <= alright) {
				insert_split_point(sequence, left, recursion);			
			}
			if(alleft <= right && right < alright) {
		//		std::cout << "alleft <= right && right < alright"<<std::endl;
				insert_split_point(sequence, right+1, recursion);
			}
			if(left < alleft && alleft < right) {
				insert_split_point(sequence, alleft, recursion);
			}
			if(left < alright && alright < right) {
				insert_split_point(sequence, alright+1, recursion);
			}
		}
	}
}


template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::recursive_splits(){
	run_initial(true);
}



template<typename overlap_type> 
void splitpoints_interval_tree<overlap_type>::initial_splits_on_al(const pw_alignment * al, bool recursion) {

	size_t left1;
	size_t right1;
	size_t left2 ;
	size_t right2;
	al->get_lr1(left1, right1);
	al->get_lr2(left2, right2);

// partial overlap caused by borders of an existing alignment within al
#if SPLITPRINT
		std::cout << " Check ref 1 overlaps " <<  std::endl;
#endif

		find_initial_split_points(al->getreference1(), left1, right1, recursion);

#if SPLITPRINT
		std::cout << " Check ref 2 overlaps " << std::endl;
#endif
	//	std::cout<< "got here!"<<std::endl;
		find_initial_split_points(al->getreference2(), left2, right2, recursion);

// partial overlap of al with itself
	if(al->getreference1()==al->getreference2()) {
		if(right1>=left2 && left1 < left2){
			insert_split_point(al->getreference2(),left2, recursion);
			insert_split_point(al->getreference2(),right1+1, recursion);
		}
		if(right2>=left1 && left2 < left1){
			insert_split_point(al->getreference2(),left1, recursion);
			insert_split_point(al->getreference2(),right2+1, recursion);
		}
	}

// partial overlap caused by borders of al within an existing alignment
//	std::cout << " Check borders " << std::endl;
	
	insert_split_point(al->getreference1(), left1, recursion);
	insert_split_point(al->getreference1(), right1+1, recursion);
	insert_split_point(al->getreference2(), left2, recursion);
	insert_split_point(al->getreference2(), right2+1, recursion);

}



template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::run_initial(bool recursion){
	if(!multi_mode) {
		initial_splits_on_al(&newal, recursion);	
	} else {
		for(size_t i=0; i<new_als.size(); ++i) {
			initial_splits_on_al(new_als.at(i), recursion);
		}
	}


/*	std::cout << " SPLITPOINTS INTITIAL , recursion " << recursion << std::endl;
	for(size_t i=0; i<split_points.size(); ++i) {
		if(!split_points.at(i).empty()) {
			std::cout << " ref " << i;
			for(std::set<size_t>::iterator it = split_points.at(i).begin(); it!=split_points.at(i).end(); ++it) {
				std::cout << " " << *it;
			}
			std::cout << std::endl;
		} 
	}
*/
}


template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::nonrecursive_splits(){
	run_initial(false);

}

/*
	A new split point happens on reference al_ref of al, Which other split points are caused by this?

*/
template<typename overlap_type> 
void splitpoints_interval_tree<overlap_type>::split_into_alignment_ref(const pw_alignment * al, size_t al_ref, size_t position, bool recursion) {

//	std::cout << " split into at " << al_ref << " pos " << position << " on al " << std::endl;
//	al->print();

// confusing: splitting into ref1, causes us to split ref2
	size_t l, r; // positions on the reference with the split point
	bool alforward = true; // is al in forward direction on al_other_seq?
	size_t al_ref_seq, al_other_seq;
	if(al_ref) {
		al->get_lr2(l,r);
		if(al->getbegin1() > al->getend1()) alforward = false;
		al_ref_seq = al->getreference2();
		al_other_seq = al->getreference1();
		
	} else {
		al->get_lr1(l,r);
		if(al->getbegin2() > al->getend2()) alforward = false;
		al_ref_seq = al->getreference1();
		al_other_seq = al->getreference2();
	}
		
		//	std::cout << "on ref 1"<<std::endl;
	if(l<position && r>=position) {
		pw_alignment fp;
		pw_alignment sp;
		bool alref_for_split = 1 - (bool) al_ref; // even more confusing: split splits ref 2 if the flag is false
		al->split(alref_for_split, position, fp, sp);
		pw_alignment fpe;
		pw_alignment spe;


		// do splitted parts have only gaps on one of their reference sequences
		bool fgaps = onlyGapSample(fp);
		bool sgaps = onlyGapSample(sp);
					
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


/*
		// print all before starting recursion
		if(alforward) {
						
			if(!sgaps) {
				std::cout << "WILL split " << al_other_seq << " " << speleft << std::endl;
			}
			// additional split point if gaps at split point
			if(!fgaps) {
				std::cout << "WILL split " << al_other_seq << " " << fperight+1 << std::endl;
			}
		} else { // al is backwards on current reference, first part of the split will be after second part
			if(!fgaps) {
				std::cout << "WILL split " << al_other_seq << " " << fpeleft << std::endl;
			}
			if(!sgaps) {
				std::cout << "WILL split " << al_other_seq << " " << speright+1 << std::endl;
			}
		}

*/

		std::set<size_t> split_check;
		if(alforward) {
						
			if(!sgaps) {
				insert_split_point(al_other_seq, speleft, recursion);
				split_check.insert(speleft);
			}
			// additional split point if gaps at split point
			if(!fgaps) {
				insert_split_point(al_other_seq, fperight+1, recursion);
				split_check.insert(fperight+1);
			}
		} else { // al is backwards on current reference, first part of the split will be after second part
			if(!fgaps) {
				insert_split_point(al_other_seq, fpeleft, recursion);
				split_check.insert(fpeleft);
			}
			if(!sgaps) {
				insert_split_point(al_other_seq, speright+1, recursion);
				split_check.insert(fpeleft);
			}
		}

// TODO
//		located_alignment lal(&(overl.data), *al);
	}
}


/*
	new splitpoint put into an alignment, does this imply more splitpoints?

*/
template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::split_into_alignment(const pw_alignment * al, size_t sequence, size_t position, bool recursion) {
	if(al->getreference1()==sequence) {
		split_into_alignment_ref(al, 0, position, recursion);
	}
	if(al->getreference2()==sequence) {
		split_into_alignment_ref(al, 1, position, recursion);
	}
}


template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::insert_split_point(size_t sequence, size_t position, bool recursion) {
	std::pair<std::set<size_t>::iterator, bool> insertes = split_points.at(sequence).insert(position);
	if(insertes.second) { // if we do not already have this split point
#if SPLITPRINT
		std::cout << " new split point " << sequence << " at " << position << std::endl;
#endif	
		if(!recursion) return;

// check if split points creates new split points with the new alignments
//		std::cout << " split into new als " << std::endl;

		if(!multi_mode) {
			split_into_alignment(&newal, sequence, position, recursion);
		} else {
			for(size_t i=0; i<new_als.size(); ++i) {
				split_into_alignment(new_als.at(i), sequence, position, recursion);
			}
		}
// check if split point creates new split points with existing alignments in the overlap class
		std::vector<const pw_alignment *> ol_als;
		overl.interval_result_on_point(sequence, position, ol_als);
		
//		std::cout << " split into " << ol_als.size() << " overlap als " << std::endl;
		for(size_t i=0; i<ol_als.size(); ++i) {
			const pw_alignment * al = ol_als.at(i);
			split_into_alignment(al, sequence, position, recursion);
		}
	}	
}


/*
	TODO at the moment overlap, splitpoints and alignment_index do operation in the order of the total number of sequences
	we can probably get faster and more memory efficient by trying to only do computations on sequences that are currently used


*/

template<typename overlap_type>
void splitpoints_interval_tree<overlap_type>::split_all(std::set<pw_alignment, compare_pw_alignment> & remove_alignments, std::vector<pw_alignment> & insert_alignments ){
		for(size_t i = 0; i <data.numSequences();i++){
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
//		std::cout << " in split all "<< remove_alignments.size() << " remove_alignments " << new_als.size() <<  " mm " << multi_mode<< std::endl;
		std::vector<pw_alignment>  split_pieces;
		if(!multi_mode)
			splits(newal, split_pieces);
		else {

			for(size_t i=0; i<new_als.size(); ++i) {
//				std::cout << " in split all, on new alignment " << i <<std::endl;
				splits(*(new_als.at(i)), split_pieces);
			}

		}
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
					assert(0);
					std::cout<<"fails here!"<<std::endl;
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
		// insert if we did not cut p
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


#endif
