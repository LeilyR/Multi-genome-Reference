#include "alignment_index.hpp"

alignment_index::alignment_index() {}

alignment_index::alignment_index(size_t num_references):trees(num_references), num_threads(1) {

}

/* alignment index can be initialized in parallel.
   we take care that we do not modify the same reference trees in parallel using a lock object for each reference

*/

alignment_index::alignment_index(size_t num_references, size_t num_threads, std::vector<const pw_alignment *> data): trees(num_references), num_threads(num_threads) {
	std::cout << " creating index " << std::endl; 
/*	for(size_t i=0; i<data.size(); ++i) {
		const pw_alignment * al = data.at(i);
		insert(al);
		if(i%1000==0) std::cout << ".";
		if(i%100000==0) std::cout << " " << i << std::endl;
	}
*/


	
	std::vector<std::vector< std::pair< size_t, const pw_alignment *> > > ref_als(num_references);
	for(size_t i=0; i<data.size(); ++i) {
		const pw_alignment * al = data.at(i);
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		ref_als.at(r1).push_back(std::make_pair(0, al));
		ref_als.at(r2).push_back(std::make_pair(1, al));
	}

	std::multimap<size_t, size_t> rsorter;
	for(size_t i=0; i<num_references; i++) {
		rsorter.insert(std::make_pair(ref_als.at(i).size(), i));
	}
	std::vector<size_t> rsorted(rsorter.size());
	size_t num = 0;
	for(std::multimap<size_t, size_t>::iterator it=rsorter.begin(); it!=rsorter.end(); ++it) {
		rsorted.at(num_references - 1 - num) = it->second;
		num++;
	}

#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
	for(size_t k=0; k<num_references; ++k) {
		size_t i = rsorted.at(k);
		clock_t ins_time = clock();
		for(size_t j=0; j<ref_als.at(i).size(); ++j) {
			size_t r = ref_als.at(i).at(j).first;
			const pw_alignment * al = ref_als.at(i).at(j).second; 
			half_insert(al, r);
		}
		ins_time = clock()- ins_time;
		double t = (double)ins_time / (double)CLOCKS_PER_SEC;
		double pt = t / trees.at(i).size();
		double ira = (double) trees.at(i).inside_number() / (double) trees.at(i).size();

/* code for bulk insert method TODO at the moment the bulk insert algorithm is slower	
		std::vector< std::pair< std::pair<size_t, size_t> , const pw_alignment * > > bdata;
		for(size_t j=0; j<ref_als.at(i).size(); ++j) {
			size_t ref = ref_als.at(i).at(j).first;
			const pw_alignment * al = ref_als.at(i).at(j).second; 
			size_t l, r;
			if(ref==0) {
				al->get_lr1(l, r);
			} else {
				al->get_lr2(l, r);
			}
			bdata.push_back(std::make_pair(std::make_pair(l, r), al));
		}
		tree_type btree;
		clock_t bins_time = clock();
		btree.bulk_insert(bdata);
		bins_time = clock() - bins_time;
		double bt = (double)bins_time / (double)CLOCKS_PER_SEC;
		double bpt = bt / btree.size();

		double bira = (double) btree.inside_number() / (double) btree.size();
*/

#pragma omp critical(pr)
{
		std::cout << " ref " << i << " done " << " size " << trees.at(i).size() << " height " << trees.at(i).height() << " levels " << trees.at(i).levels()<< " rebas " << trees.at(i).get_num_reba()<< " moves " << trees.at(i).get_num_moves()<< " list ratio " << trees.at(i).list_ratio()<<  " ins " << ira << " t " << t << " t per al " << pt << std::endl;
//		std::cout << " ref " << i << " bulk " << " size " << btree.size() << " height " << btree.height() << " rebas " << btree.get_num_reba()<< " moves " << btree.get_num_moves()<< " list ratio " << btree.list_ratio() << " ins " << bira<< " t " << bt << " t per al " << bpt << std::endl;
}

	}

/*
	std::vector<omp_lock_t*> reference_locks(num_references);
	std::vector<double> ins_times(num_references, 0);
	for(size_t i=0; i<reference_locks.size(); ++i) {
		omp_lock_t * tl = new omp_lock_t;
		omp_init_lock(tl);
		reference_locks.at(i) = tl;
	}

	int threads_running = 0;

#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
	for(size_t i=0; i<data.size(); ++i) {
		const pw_alignment * al = data.at(i);
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
#pragma omp critical(lock)
{
		omp_set_lock((reference_locks.at(r1)));
}


#pragma omp critical(tr)
{
		threads_running++;
}
		clock_t i_time = clock();
		half_insert(al, 0);
		i_time = clock() - i_time;
		ins_times.at(r1)+=((double)i_time/CLOCKS_PER_SEC);

#pragma omp critical(tr)
{
		threads_running--;
}
#pragma omp critical(unlock)
{
		omp_unset_lock(reference_locks.at(r1));
}

#pragma omp critical(lock)
{
		omp_set_lock((reference_locks.at(r2)));
}
#pragma omp critical(tr)
{
		threads_running++;
}
		clock_t i_time2 = clock();
		half_insert(al, 1);
		i_time2 = clock() - i_time2;
		ins_times.at(r2)+=((double)i_time2/CLOCKS_PER_SEC);
#pragma omp critical(tr)
{
		threads_running--;
}

#pragma omp critical(unlock)
{
		omp_unset_lock(reference_locks.at(r2));
}

#pragma omp critical(print)
{
		if(i%10==0) std::cout << "." << std::flush;
		if(i%1000==0 || i+1 == data.size()) { 
			std::cout << " " << i << " r1 " << r1 << " r2 " << r2 << " threads running "<< threads_running<< std::endl <<  std::flush;

			for(size_t j=0; j<trees.size(); ++j) {
				if(trees.at(j).size()>0) {
					double it = ins_times.at(j);
					double wit = it / (double) trees.at(j).size();
					std::cout << " tree "<< j << " size " << trees.at(j).size() << " height " << trees.at(j).height() << " list ratio " << trees.at(j).list_ratio() << " tot time " << it << " time per ins " << wit << std::endl;

				}
			}
		}

}
	}
	for(size_t i=0; i<reference_locks.size(); ++i) {
		omp_destroy_lock((reference_locks.at(i)));
		delete reference_locks.at(i);
	}
*/
	

}
alignment_index::~alignment_index() {

}

void alignment_index::insert(const pw_alignment * al) {
	size_t l1, r1, l2, r2;
	al->get_lr1(l1, r1);
	trees.at(al->getreference1()).insert(l1, r1, al);
	al->get_lr2(l2, r2);
	trees.at(al->getreference2()).insert(l2, r2, al);
}
void alignment_index::half_insert(const pw_alignment * al, size_t reference) {
	assert(reference<2);
	size_t l1, r1, l2, r2;
	if(reference==0) {
	al->get_lr1(l1, r1);
	trees.at(al->getreference1()).insert(l1, r1, al);
	} else {
	al->get_lr2(l2, r2);
	trees.at(al->getreference2()).insert(l2, r2, al);
	}
}

size_t alignment_index::erase(const pw_alignment * al) {
	size_t l, r;
	al->get_lr1(l, r);
	size_t ref1 = al->getreference1();
	size_t res1 = trees.at(ref1).erase(l, r, al);
//	std::cout << "res1 "<< res1 << "tree is "<< std::endl;
//	trees.at(ref1).debug_print();

	al->get_lr2(l, r);
	size_t ref2 = al->getreference2();
	size_t res2 = trees.at(ref2).erase(l, r, al);
	assert(res1 == res2);
	return res1;
}

void alignment_index::search_overlap(const pw_alignment & al, const size_t & reference, std::vector<const pw_alignment *> & result) const {
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
void alignment_index::search_overlap(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result) const {
	trees.at(reference).overlap_search(left, right, result);	
}


void alignment_index::search_overlap_fraction(const size_t & reference, const size_t & left, const size_t & right, const double & fraction, std::vector<const pw_alignment *> & result) const {
	trees.at(reference).overlap_fraction_search(left, right, fraction, result);	
	
}

void alignment_index::super_search_overlap_and_remove(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result, std::multimap<size_t, std::pair<size_t, size_t> > & touched_intervals) {

	// find result on a single tree
	trees.at(reference).overlap_search_erase(left, right, result);
	std::cout << " SSOR search on ref " << reference << " tree has " << trees.at(reference).size() << " intvs, found " << result.size() << std::endl << std::flush;


	
	// find all intervals touched on other references
	std::multimap<size_t, std::pair<size_t, size_t> > all_intervals; 
	std::multimap<size_t, std::pair< std::pair<size_t, size_t>, const pw_alignment*> > erase_data;
	std::set<size_t> used_refs;


	for(size_t i=0; i<result.size(); ++i) {
		const pw_alignment * al = result.at(i);
		size_t l1, r1;
		al->get_lr1(l1, r1);
		all_intervals.insert(std::make_pair(al->getreference1(), std::make_pair(l1, r1) ) );
		used_refs.insert(al->getreference1());
		size_t l2, r2;
		al->get_lr2(l2, r2);
		all_intervals.insert(std::make_pair(al->getreference2(), std::make_pair(l2, r2) ) );
		used_refs.insert(al->getreference2());
		size_t ref1 = al->getreference1();
		size_t ref2 = al->getreference2();

		if(ref1 != reference) {
//			trees.at(ref1).erase(l1, r1, al);
			erase_data.insert(std::make_pair(ref1, std::make_pair(std::make_pair(l1, r1), al) ));
		} else if(r1 >= left && l1 <= right) {
			// this alignment was already erased be overlap_search_remove
		} else {
			erase_data.insert(std::make_pair(ref1, std::make_pair(std::make_pair(l1, r1), al) ));
//			trees.at(ref1).erase(l1, r1, al);
		}
		if(ref2 != reference) {
//			trees.at(ref2).erase(l2, r2, al);
			erase_data.insert(std::make_pair(ref2, std::make_pair(std::make_pair(l2, r2), al) ));
		} else if(r2 >= left && l2 <= right) {
			// this alignment was already erased be overlap_search_remove
		} else {
//			trees.at(ref2).erase(l2, r2, al);
			erase_data.insert(std::make_pair(ref2, std::make_pair(std::make_pair(l2, r2), al) ));
		}
	}

	// erase all alignments on other references
	for(std::set<size_t>::iterator it = used_refs.begin(); it!=used_refs.end(); ++it) {
		std::vector<std::pair< std::pair<size_t, size_t>, const pw_alignment*> > erase_here;
		std::pair< std::multimap<size_t, std::pair< std::pair<size_t, size_t>, const pw_alignment*> >::iterator,
			std::multimap<size_t, std::pair< std::pair<size_t, size_t>, const pw_alignment*> >::iterator > eqr = erase_data.equal_range(*it);
		for(std::multimap<size_t, std::pair< std::pair<size_t, size_t>, const pw_alignment*> >::iterator eit = eqr.first; eit!= eqr.second; ++eit) {
			erase_here.push_back(eit->second);
		}
		if(!erase_here.empty()) {
			std::cout << " on ref " << *it << " erasing " << erase_here.size() << " nodes " << std::endl << std::flush;  // TODO
		}
		trees.at(*it).erase_vector(erase_here);
	}



	// simplify results on other references
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
//		std::cout << " Ref " << this_ref << " interval tree " << std::endl;
//		this_tree.debug_print();

		std::vector<std::pair<size_t, size_t> > this_intervals;
		this_tree.join_intervals(this_intervals);
//		std::cout << "touched intervals: "<<std::endl;
		for(size_t i=0; i<this_intervals.size(); ++i) {
//			std::cout << " ref " << this_ref << " l " << this_intervals.at(i).first << " r " << this_intervals.at(i).second << std::endl;
			touched_intervals.insert(std::make_pair( this_ref, std::make_pair( this_intervals.at(i).first, this_intervals.at(i).second) ) );
		}
		if(!this_intervals.empty()) {
			std::cout << " ref " << *it << " was touched by  " << this_tree.size() << " results alignments. We simplified this to " << this_intervals.size() << " intervals " << std::endl << std::flush;
		}
	}
}


void alignment_index::super_search_overlap_no_remove(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result, std::multimap<size_t, std::pair<size_t, size_t> > & touched_intervals) const {



}


void alignment_index::debug_print() const {
//	for(size_t i=0; i<trees.size(); ++i) {
			//	trees.at(i).debug_print();
//	}
}

#include "intervals.cpp"
