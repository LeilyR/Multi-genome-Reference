#include "alignment_index.hpp"

alignment_index::alignment_index(size_t num_references):trees(num_references) {

}

alignment_index::~alignment_index() {

}

void alignment_index::insert(const pw_alignment * al) {
	size_t l1, r1, l2, r2;
	al->get_lr1(l1, r1);
	trees.at(al->getreference1()).insert(l1, r1, al);
	al->get_lr2(l2, r2);
	trees.at(al->getreference2()).insert(l2, r2, al);
	if(al->getreference1()==al->getreference2()) {
		assert(l1!=l2 || r1!=r2);
	}
}

size_t alignment_index::erase(const pw_alignment * al) {
	size_t l, r;
	al->get_lr1(l, r);
	size_t ref1 = al->getreference1();
	size_t res1 = trees.at(ref1).erase(l, r, al);
	std::cout << "res1 "<< res1 << "tree is "<< std::endl;
	trees.at(ref1).debug_print();

	al->get_lr2(l, r);
	size_t ref2 = al->getreference2();
	size_t res2 = trees.at(ref2).erase(l, r, al);
	std::cout << "res2 "<< res2 << std::endl;
	assert(res1 == res2);
	return res1;
}

void alignment_index::search_overlap(const pw_alignment & al, const size_t & reference, std::vector<const pw_alignment *> & result) {
	size_t fref;
	size_t fleft;
	size_t fright;
	assert(reference<2);
	if(reference == 0) {
		fref = al.getreference1();
		al.get_lr1(fleft, fright);
	} else {
		fref = al.getreference2();
		al.get_lr2(fleft, fright);
	}
	trees.at(fref).overlap_search(fleft, fright, result);	

}
void alignment_index::search_overlap(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result) {
	trees.at(reference).overlap_search(left, right, result);	
}


void alignment_index::super_search_overlap_and_remove(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result, std::multimap<size_t, std::pair<size_t, size_t> > & touched_intervals) {
	trees.at(reference).overlap_search(left, right, result);

	std::cout << "SUPER SEARCH found " << result.size() << " alignments" << std::endl;

	std::multimap<size_t, std::pair<size_t, size_t> > all_intervals; 
	std::set<size_t> used_refs;
	for(size_t i=0; i<result.size(); ++i) {
		const pw_alignment * al = result.at(i);
		al->print(); // TODO 
		std::cout << std::endl;
		size_t l, r;
		al->get_lr1(l, r);
		all_intervals.insert(std::make_pair(al->getreference1(), std::make_pair(l, r) ) );
		used_refs.insert(al->getreference1());
		al->get_lr2(l, r);
		all_intervals.insert(std::make_pair(al->getreference2(), std::make_pair(l, r) ) );
		used_refs.insert(al->getreference2());
	}
	
	for(std::set<size_t>::iterator it = used_refs.begin(); it!=used_refs.end(); ++it) {
		size_t this_ref = *it;
		std::pair< std::multimap<size_t, std::pair<size_t, size_t> >::iterator, std::multimap<size_t, std::pair<size_t, size_t> >::iterator > eqr = all_intervals.equal_range(this_ref);
		tree_type this_tree;
		for(std::multimap<size_t, std::pair<size_t, size_t> >::iterator eit = eqr.first; eit!=eqr.second; ++eit) {
			size_t thisl = eit->second.first;
			size_t thisr = eit->second.second;
		//	std::cout << " r " << this_ref << " l " << thisl << " r " << thisr << std::endl;
			this_tree.insert(thisl, thisr, NULL);
		}
		std::cout << " Ref " << this_ref << " interval tree " << std::endl;
		this_tree.debug_print();

		std::vector<std::pair<size_t, size_t> > this_intervals;
		this_tree.join_intervals(this_intervals);
		std::cout << "touched intervals: "<<std::endl;
		for(size_t i=0; i<this_intervals.size(); ++i) {//XXX Basically this_interval could be already removed from the tree of that ref but we need to find any interval on the tree that has overlap with them.
			std::cout << " ref " << this_ref << " l " << this_intervals.at(i).first << " r " << this_intervals.at(i).second << std::endl;
			touched_intervals.insert(std::make_pair( this_ref, std::make_pair( this_intervals.at(i).first, this_intervals.at(i).second) ) );
		}
	}

	std::cout << " SUPER SEARCH simplified alignments to " << touched_intervals.size() << " intervals on " << used_refs.size() << " references " << std::endl;

	for(size_t i=0; i<result.size(); ++i) {
		const pw_alignment * al = result.at(i);
		std::cout<< "result at "<< i  << " is " << std::endl;
		al->print();
		size_t erase_res = erase(result.at(i));
//		assert(erase_res);
	}


}


void alignment_index::debug_print() const {
	for(size_t i=0; i<trees.size(); ++i) {
		std::cout << " On REF " << i << std::endl;
		trees.at(i).debug_print();
	}
}

#include "intervals.cpp"
