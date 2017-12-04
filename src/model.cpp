#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP
#define PRINT 0



/*
	TODO we can make the slow part even faster:
	As the repetitive part is already done, new insertions will probably not cut the repetitve part. Keep previous content fixed in the overlap class. A function to insert new pieces already cut in a way that they fit to the repetitve pieces


*/
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::compute(overlap_type & o) {

//	compute_simple_lazy_splits(o);
//	compute_simple(o);
	compute_vcover_clarkson(o);

}
/*
input:
	all_ins/all_rem contains modifications that were previously done to ovrlp
	this_ins/rem contains new modifications that were just done to overlp

output:
	consolidated modifications in all_ins/rem
	

*/
template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::insert_alignment_sets(std::set<pw_alignment, compare_pw_alignment> & all_ins, std::set<pw_alignment, compare_pw_alignment> & all_rem, std::vector<pw_alignment> & this_ins, std::vector<pw_alignment> & this_rem) {
	
	for(size_t i=0; i<this_ins.size(); i++) {
		const pw_alignment & al = this_ins.at(i);
		std::set<pw_alignment>::iterator findr = all_rem.find(al);
		if(findr!=all_rem.end()) {
// al was previously removed. Now it was inserted again
			all_rem.erase(al);
		} else {
			all_ins.insert(al);
		}
	}

	for(size_t i=0; i<this_rem.size(); i++) {
		const pw_alignment & al = this_rem.at(i);
		std::set<pw_alignment>::iterator findr = all_ins.find(al);
		if(findr!=all_ins.end()) {
// al was previously inserted. Now it was removed again
			all_ins.erase(al);
		} else {
			all_rem.insert(al);
		}

	}
}




/*
	Insert remove dynamics:
	we never want to insert an alignment that was just removed (in same function level)
	alignments that were just inserted, might be removed by the next call. In that case we can permanentely delete

*/
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::insert_alignment_sets(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_ins, std::set<pw_alignment, compare_pw_alignment> & all_rem, std::vector<pw_alignment> & this_ins, std::vector<pw_alignment> & this_rem) {
	
	for(size_t i=0; i<this_ins.size(); i++) {
		all_ins.insert(this_ins.at(i));
	}

	for(size_t i=0; i<this_rem.size(); i++) {
		std::set<pw_alignment>::iterator findr = all_ins.find(this_rem.at(i));
		if(findr!=all_ins.end()) {
			const pw_alignment & al = this_rem.at(i);
			all_ins.erase(al);
			all_rem.erase(al);
//	std::cout << "is gonna be removed3: "<<&al <<std::endl;
	//	al.print();

			ovrlp.remove_alignment(al); 
		} else {
			all_rem.insert(this_rem.at(i));
		}

	}
}

template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::local_undo(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_inserted, std::set<pw_alignment , compare_pw_alignment> & all_removed) {
	for(std::set< pw_alignment , compare_pw_alignment >::iterator it = all_inserted.begin(); it!=all_inserted.end(); it++) {
		const pw_alignment & al = *it;
	//	std::cout << "is gonna be removed2: "<<&al <<std::endl;
	//	al.print();
		ovrlp.remove_alignment(al);
	}

	for(std::set< pw_alignment, compare_pw_alignment >::iterator it = all_removed.begin(); it!= all_removed.end(); it++) {
		const pw_alignment & ral = *it;
		ovrlp.insert_without_partial_overlap(ral);

// slow test
//		ovrlp.test_partial_overlap();


	}

}


template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::local_undo(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_inserted, std::set<pw_alignment , compare_pw_alignment> & all_removed) {
	for(std::set< pw_alignment , compare_pw_alignment >::iterator it = all_inserted.begin(); it!=all_inserted.end(); it++) {
		const pw_alignment & al = *it;
	//	std::cout << "is gonna be removed2: "<<&al <<std::endl;
	//	al.print();
		ovrlp.remove_alignment(al);
	}

	for(std::set< pw_alignment, compare_pw_alignment >::iterator it = all_removed.begin(); it!= all_removed.end(); it++) {
		const pw_alignment & ral = *it;
		ovrlp.insert_without_partial_overlap(ral);
	}

}

template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::all_push_back(std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, std::set<pw_alignment , compare_pw_alignment> & all_inserted, std::set< pw_alignment, compare_pw_alignment> & all_removed ) {
	for(std::set<pw_alignment>::iterator it = all_inserted.begin(); it!=all_inserted.end(); it++) {
		const pw_alignment & p = *it;
		inserted_alignments.push_back(p);
	}

	for(std::set<pw_alignment>::iterator it = all_removed.begin(); it!= all_removed.end(); it++) {
		const pw_alignment & p = *it;
		removed_alignments.push_back(p);
	}

}

template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::all_push_back(std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, std::set<pw_alignment , compare_pw_alignment> & all_inserted, std::set< pw_alignment, compare_pw_alignment> & all_removed ) {
	for(std::set<pw_alignment>::iterator it = all_inserted.begin(); it!=all_inserted.end(); it++) {
		const pw_alignment & p = *it;
		inserted_alignments.push_back(p);
	}

	for(std::set<pw_alignment>::iterator it = all_removed.begin(); it!= all_removed.end(); it++) {
		const pw_alignment & p = *it;
		removed_alignments.push_back(p);
	}

}


template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::lazy_split_full_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment & alin, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {

//	std::cout << " full insert step on level " << level  << std::endl;

	size_t orig_ovrlp_size = ovrlp.get_all().size();

	std::set<pw_alignment, compare_pw_alignment> remove_als; 
	std::vector<pw_alignment> insert_als;
	splitpoints_interval_tree<overlap_type> spl2(alin, ovrlp, data);
	spl2.recursive_splits();
//	std::cout<< "inja!"<<std::endl;

	spl2.split_all(remove_als, insert_als);
//	std::cout<< "inja1!"<<std::endl;

	size_t inserted_counter = 0;
	if(remove_als.empty()) {
// no remove als, we insert all pieces with positive gain
		local_gain = 0;
		for(size_t i=0; i<insert_als.size(); i++) {
			double g1;
			double g2;
			model.gain_function(insert_als.at(i), g1, g2);
			double avg = (g1 + g2) / 2.0;
		//	std::cout << "avg "<< avg <<std::endl;
			if(avg > 0) {
				ovrlp.insert_without_partial_overlap(insert_als.at(i));

//  
//				ovrlp.test_partial_overlap();

				inserted_alignments.push_back(insert_als.at(i));
				local_gain+=avg;
				inserted_counter++;
			}
		}		
	//	std::cout << "full insert: " << level << " on al length " << alin.alignment_length() << " finally inserted pieces: " << inserted_counter << " real local gain " << local_gain << std::endl;		
	// insert pieces for recursion:
	} else {
		std::set< pw_alignment , compare_pw_alignment> all_removed;
		std::set<pw_alignment, compare_pw_alignment> all_inserted;
		double lost_gain = 0;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); it++) {
			const pw_alignment & ral = *it;
			double gain1, gain2;
			model.gain_function(ral, gain1, gain2);
			double rav_gain = (gain1 + gain2)/2;
			lost_gain += rav_gain;
		//	std::cout << "is gonna be removed1: " << &ral<<std::endl;
		//	ral.print();

			ovrlp.remove_alignment(ral);
			all_removed.insert(ral);
		//			std::cout << "REM2 " << std::endl;
		//	       		ral->print();
		//			std::cout << std::endl;	
		}

		double max_parts_gain = 0;
		std::multimap<double, const pw_alignment *> ordered_parts;
		for(size_t i=0; i<insert_als.size(); i++) {
			double gain1, gain2;
			model.gain_function(insert_als.at(i), gain1, gain2);
			double iav_gain = (gain1 + gain2)/2;
			if(iav_gain > 0) {
				ordered_parts.insert(std::make_pair(iav_gain, & insert_als.at(i)));
				max_parts_gain+=iav_gain;
			}

		}
	//	std::cout << "2nd split: " << level << " on al length " << alin.alignment_length() <<  " new remove " << remove_als.size() << " lost gain " << lost_gain <<
	//	       	" new insert " << ordered_parts.size() << " gain "<<  max_parts_gain << std::endl;

		if(max_parts_gain > lost_gain) {
// we could get positive gain from recursive calls
			double sum_of_gain = 0;
			size_t alnum = 0;


			for(std::multimap<double, const pw_alignment *>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); ++it) {
				double thisgain;
				std::vector<pw_alignment> this_inserted;
				std::vector<pw_alignment> this_removed;
				const pw_alignment & thisal = *it->second;
				std::vector<pw_alignment> alv;
				alv.push_back(thisal);

	//		std::cout <<"full insert function level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << std::endl; 
			// this function only changes ovrlp if this was locally good. We may need to undo the changes if they were globally bad
			// we go back to normal non-recursive split function to avoid uneccessary splitting
				lazy_split_insert_step(ovrlp, level + 1, rec_calls,  alv, this_inserted, this_removed, thisgain);
				insert_alignment_sets(all_inserted, all_removed, this_inserted, this_removed);
				sum_of_gain += thisgain;

				alnum++;
			}
			if(sum_of_gain > lost_gain) {
// recursive calls did result in positive gain
// take this step
				all_push_back(inserted_alignments, removed_alignments, all_inserted, all_removed);
				local_gain = sum_of_gain - lost_gain;

				assert(ovrlp.get_all().size() + removed_alignments.size() == orig_ovrlp_size + inserted_alignments.size());

				return;
			}

		} 

// we did not get positive gain, local undo
		local_gain = 0;
		local_undo(ovrlp, all_inserted, all_removed);
	}
	assert(ovrlp.get_all().size() + removed_alignments.size() == orig_ovrlp_size + inserted_alignments.size());
}

template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::lazy_split_full_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment & alin, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {

//	std::cout << " full insert step on level " << level  << std::endl;

	std::set<pw_alignment, compare_pw_alignment> remove_als; 
	std::vector<pw_alignment> insert_als;
	splitpoints_interval_tree<overlap_type> spl2(alin, ovrlp, data);
	spl2.recursive_splits();
//	std::cout<< "inja!"<<std::endl;

	spl2.split_all(remove_als, insert_als);
//	std::cout<< "inja1!"<<std::endl;

	size_t inserted_counter = 0;
	if(remove_als.empty()) {
		local_gain = 0;
		for(size_t i=0; i<insert_als.size(); i++) {
			double g1;
			double g2;
			common_model.gain_function(insert_als.at(i), g1, g2);
			double avg = (g1 + g2) / 2.0 - base_cost;
		//	std::cout << "avg "<< avg <<std::endl;
			if(avg > 0) {
				ovrlp.insert_without_partial_overlap(insert_als.at(i));
				inserted_alignments.push_back(insert_als.at(i));
				local_gain+=avg;
				inserted_counter++;
			}
		}		
	//	std::cout << "full insert: " << level << " on al length " << alin.alignment_length() << " finally inserted pieces: " << inserted_counter << " real local gain " << local_gain << std::endl;		
	// insert pieces for recursion:
	} else {
		std::set< pw_alignment , compare_pw_alignment> all_removed;
		double lost_gain = 0;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); it++) {
			const pw_alignment & ral = *it;
			double gain1, gain2;
			common_model.gain_function(ral, gain1, gain2);
			double rav_gain = (gain1 + gain2)/2 - base_cost;
			lost_gain += rav_gain;
		//	std::cout << "is gonna be removed1: " << &ral<<std::endl;
		//	ral.print();

			ovrlp.remove_alignment(ral);
			all_removed.insert(ral);
		//			std::cout << "REM2 " << std::endl;
		//	       		ral->print();
		//			std::cout << std::endl;	
		}

		double max_parts_gain = 0;
		std::multimap<double, const pw_alignment *> ordered_parts;
		for(size_t i=0; i<insert_als.size(); i++) {
			double gain1, gain2;
			common_model.gain_function(insert_als.at(i), gain1, gain2);
			double iav_gain = (gain1 + gain2)/2 - base_cost;
			if(iav_gain > 0) {
				ordered_parts.insert(std::make_pair(iav_gain, & insert_als.at(i)));
				max_parts_gain+=iav_gain;
			}

		}
	//	std::cout << "2nd split: " << level << " on al length " << alin.alignment_length() <<  " new remove " << remove_als.size() << " lost gain " << lost_gain <<
	//	       	" new insert " << ordered_parts.size() << " gain "<<  max_parts_gain << std::endl;

		double sum_of_gain = 0;
		size_t alnum = 0;

		std::set<pw_alignment, compare_pw_alignment> all_inserted;

		for(std::multimap<double, const pw_alignment *>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); ++it) {
			double thisgain;
			std::vector<pw_alignment> this_inserted;
			std::vector<pw_alignment> this_removed;
			const pw_alignment & thisal = *it->second;

	//		std::cout <<"full insert function level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << std::endl; 
			// this function only changes ovrlp if this was locally good. We may need to undo the changes if they were globally bad
			lazy_split_insert_step(ovrlp, level + 1, rec_calls,  thisal, this_inserted, this_removed, thisgain);
			insert_alignment_sets(ovrlp, all_inserted, all_removed, this_inserted, this_removed);
			sum_of_gain += thisgain;

			alnum++;
		}


		if(lost_gain > sum_of_gain) {
		// Undo whole function call if it did not pay off
			local_gain = 0;
			local_undo(ovrlp, all_inserted, all_removed);
		
		} else {
		// take this step
			all_push_back(inserted_alignments, removed_alignments, all_inserted, all_removed);
			local_gain = sum_of_gain - lost_gain;
		}
	}
}


template<typename T, typename overlap_type>
bool initial_alignments_from_groups<T,overlap_type>::try_rnodes_unselect(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> insert_gain, const std::vector<double> & remove_gain, const std::set<size_t> & rem_to_unselect, vector<bool> & ins_selected, vector<bool> & rem_selected, double & new_gain) const {
	
	std::set<size_t> ins_unsel; // if we unselect rem_to_unselect, we have to unselect ins_unsel
	for(std::set<size_t>::const_iterator it = rem_to_unselect.begin(); it!=rem_to_unselect.end(); it++) {
		size_t rem_uns = *it;
		assert(rem_selected.at(rem_uns));
		new_gain += remove_gain.at(rem_uns);
		const std::set<size_t> & ru_neigh = rem_to_ins.at(rem_uns);
		for(std::set<size_t>::const_iterator itt = ru_neigh.begin(); itt!=ru_neigh.end(); ++itt) {
			if(ins_selected.at(*itt)) {
				new_gain -= insert_gain.at(*itt);
//				std::cout << " rem " << rem_uns << " has neigh " << *itt << std::endl;
				ins_unsel.insert(*itt);
			}
		}
	}
	std::set<size_t> ins_unsel_n; // remove nodes neighboring ins_unsel, we may be able to also unselect them
	for(std::set<size_t>::const_iterator it = ins_unsel.begin(); it!=ins_unsel.end(); ++it) {
		size_t iu = *it;
		const std::set<size_t> & iun = ins_to_rem.at(iu);
		for(std::set<size_t>::const_iterator itt = iun.begin(); itt!=iun.end(); ++itt) {
			size_t iuni = *itt;
			if(rem_selected.at(iuni)) {
				std::set<size_t>::const_iterator find_iuni = rem_to_unselect.find(iuni);
				if(find_iuni == rem_to_unselect.end()) {
					ins_unsel_n.insert(iuni);
//					std::cout << " ins to unsel " << iu << " has neigh " << iuni << std::endl;
				}
			}
		}

	}

	std::set<size_t> ins_unsel_n_unsel; // subsect of ins_unsel_n which can actually be unselected
	for(std::set<size_t>::const_iterator it = ins_unsel_n.begin(); it!=ins_unsel_n.end(); ++it) {
		size_t iu = *it;
		bool canunsel = true;
		const std::set<size_t> & iun = rem_to_ins.at(iu);
		for(std::set<size_t>::const_iterator itt = iun.begin(); itt!=iun.end(); ++itt) {
			size_t p = *itt;
			if(ins_selected.at(p)) {
				std::set<size_t>::const_iterator findp = ins_unsel.find(p);
				if(findp == ins_unsel.end()) {
					canunsel = false; // we have to keep remove alignment iu, because it is caused by insert_alignment p
					break;
				}

			}

		}
		if(canunsel) {
			ins_unsel_n_unsel.insert(iu);
			new_gain += remove_gain.at(iu);
		}
	}
 
	if(new_gain > 0) {
		for(std::set<size_t>::const_iterator it = rem_to_unselect.begin(); it!=rem_to_unselect.end(); it++) {
			rem_selected.at(*it) = false;
		}
		for(std::set<size_t>::const_iterator it = ins_unsel.begin(); it!=ins_unsel.end(); ++it) {
			ins_selected.at(*it) = false;
		}
		for(std::set<size_t>::const_iterator it = ins_unsel_n_unsel.begin(); it!=ins_unsel_n_unsel.end(); ++it) {
			rem_selected.at(*it) = false;
		}
		
		return true;
	} else {
		new_gain = 0;
		return false;
	}

}


/*
 Compute current upper bounds of gains 
	a) sum of possible gain of all remaining (free & !taken) remove alignments, 
	b) sum of gain of all insert alignments which are not yet taken

*/
template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::ub_gain(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, 
std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & ub) const {
	double uba = 0;

	for(size_t i=0; i<rem_to_ins.size(); ++i) {
	// if i is not taken alreay and free to take in that subtree
		if(rem_free.at(i) && !rem_taken.at(i)) {
			double rgain = remove_gain.at(i);
			double igain_possible = 0;
			const std::set<size_t> & ins_of_i = rem_to_ins.at(i); // taking rem i, could make it possible to take all these ins alignments
			for(std::set<size_t>::const_iterator it = ins_of_i.begin(); it!=ins_of_i.end(); ++it) {
				size_t j=*it;
				if(!ins_taken.at(j)) {
					igain_possible += insert_gain.at(j);		
				}
			}
			double pgain = igain_possible - rgain;
//			std::cout << "in ub r al " << i << " possible gain " << pgain << std::endl;
			if(pgain > 0) {
				uba += pgain;
			}
		}
	}

	double ubb = 0;
	for(size_t i=0; i<insert_gain.size(); ++i) {
		if(!ins_taken.at(i)) {
			const std::set<size_t>  irem = ins_to_rem.at(i);
			double igain = insert_gain.at(i);
			bool cantakei = true; 
			for(std::set<size_t>::iterator it = irem.begin(); it!=irem.end(); ++it) {
				size_t j = *it;
				if(!rem_taken.at(j)) {
					if(!rem_free.at(j)) {
						cantakei = false;
					} 
				}
			}
//			std::cout << " for ins "<< i << " possible gain is " << igain << std::endl;
			if(cantakei) {
				ubb+=igain;
			}
		}
	}

//	std::cout << " uba " << uba << " ubb " << ubb << std::endl;
	ub = uba;
	if(ubb < ub) ub = ubb;

}






/*
	which insert alignments will directly pay of (insert gain higher than all removed gain caused by it)

	as long as the data changes by this procedure we go on

*/
template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::easy_insert(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, 
std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & extra_gain) const {
	bool go_on = true;
	while(go_on) {
		go_on = false;

		for(size_t i=0; i<ins_to_rem.size(); ++i) {
			if(!ins_taken.at(i)) {
				const std::set<size_t>  irem = ins_to_rem.at(i);
				double igain = insert_gain.at(i);
				for(std::set<size_t>::iterator it = irem.begin(); it!=irem.end(); ++it) {
					size_t j = *it;
					if(!rem_taken.at(j)) {
						if(rem_free.at(j)) {
							igain -= remove_gain.at(j);
						} else {
							igain = -std::numeric_limits<double>::max();
							break;
						}
					}
				}
//				std::cout << " for ins " << i << " new direct gain is " << igain << std::endl;
				if(igain > 0) {
					ins_taken.at(i) = true;
					for(std::set<size_t>::iterator it = irem.begin(); it!=irem.end(); ++it) {
						size_t j = *it;
						rem_taken.at(j) = true;
						rem_free.at(j) = false;
					}
					go_on = true;
//					std::cout << " TAKE AND GO ON " << std::endl;
					extra_gain += igain;
				}
			}
		}
	}
//	std::cout << " easy insert extra gain " << extra_gain << std::endl;

}



/*
	chooses free rem alignment with best possible gain
	then selects all other alignments that increase gain


	in: old real_gain, out: new real_gain, possible_gain is in addition to new real_gain

	we change free vector so that we will not try the same things twice



1) which remove alignment has the highest possible gain (chosen)

2) insert chosen into ins_taken and rem_taken, remove it from rem_free

3) compute new gain

*/
template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::take_highest_possible(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, 
std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & gain 
) const {
// search for highest possible gain
	size_t chosen = rem_to_ins.size();
	double best_possible = 0;

// compute possible gain for each remove al, take best one
	for(size_t i=0; i<rem_to_ins.size(); ++i) {
		// if i is not taken alreay and free to take in that subtree
		if(rem_free.at(i) && !rem_taken.at(i)) {
			double rgain = remove_gain.at(i);
			double igain_possible = 0;
			const std::set<size_t> & ins_of_i = rem_to_ins.at(i); // taking rem i, could make it possible to take all these ins alignments
			for(std::set<size_t>::const_iterator it = ins_of_i.begin(); it!=ins_of_i.end(); ++it) {
				size_t j=*it;
				if(!ins_taken.at(j)) {
					igain_possible += insert_gain.at(j);		
				}
			}
			double pgain = igain_possible - rgain;
//			std::cout << " r al " << i << " possible gain " << pgain << std::endl;
			if(pgain > best_possible) {
				best_possible = pgain;
				chosen = i;
			}
		}
	}
// choose the best
	if(chosen < rem_to_ins.size()) {
		// take chosen remove alignment
		rem_free.at(chosen) = false;
		rem_taken.at(chosen) = true;
//		std::cout << " take " << chosen << std::endl;
		gain -= remove_gain.at(chosen);
	}


}


template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::bb_step(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, const std::vector<bool> & rem_free, const std::vector<bool> & ins_taken, const std::vector<bool> & rem_taken, const double & real_gain, double & best_gain, std::vector<bool> & best_ins, std::vector<bool> & best_rem) const {


// storage for both recursive calls
	std::vector<bool> take_ins = ins_taken;
	std::vector<bool> notake_ins = ins_taken;
	std::vector<bool> take_rem = rem_taken;
	std::vector<bool> notake_rem = rem_taken;
	std::vector<bool> both_remfree = rem_free;


	double take_gain = real_gain;
	double notake_gain = real_gain;
	double take_extra_gain = 0;
	double notake_extra_gain = 0;
	double take_ub;
	double notake_ub;

// select best remove alignment (highest possible gain), this modifies this_rem_free, to never make the same choice in this recursion subtree
	take_highest_possible(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, take_ins, take_rem, take_gain);

// do easy insertions in both
	easy_insert(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, take_ins, take_rem, take_extra_gain);
//	std::cout << " easy ins in take: " << take_extra_gain << std::endl;
	easy_insert(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, notake_ins, notake_rem, notake_extra_gain);
//	std::cout << " easy ins in notake: " << notake_extra_gain << std::endl;
	take_gain += take_extra_gain;
	notake_gain += notake_extra_gain;

// compute upper bounds 
	ub_gain(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, take_ins, take_rem, take_ub);
	ub_gain(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, notake_ins, notake_rem, notake_ub);


// found an optimum
	if(take_gain > best_gain) {
		best_gain = take_gain;
		best_ins = take_ins;
		best_rem = take_rem;
	}
	if(notake_gain > best_gain) {
		best_gain = notake_gain;
		best_ins = notake_ins;
		best_rem = notake_rem;
	}	



// current upper bound more than lower bound (best previous)
	if(take_gain + take_ub > best_gain) {
		
// take current remove alignment with highest possible gain
		bb_step(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, take_ins, take_rem, take_gain, best_gain, best_ins, best_rem);

	}

	if(notake_gain + notake_ub > best_gain) {
// never take current remove alignment with highest possible gain
		bb_step(ins_to_rem, rem_to_ins, insert_gain, remove_gain, both_remfree, notake_ins, notake_rem, notake_gain, best_gain, best_ins, best_rem);
	}

}




/*
	TODO 
this algorithm is not that good. I think it only gets called on small problems. If not it could be too small. 
can we prove that the problem is NP-complete?
Or is there an easy algorithm to solve it and I have no idea of it.
To me it seems  like something that should be known from operations research


*/



template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::select_from_groups(const std::vector<std::vector<pw_alignment> > & insert_als, const std::vector<std::vector<pw_alignment> > & rem_als_per_ins, std::vector<pw_alignment> & result_ins, std::vector<pw_alignment> & result_rem, double & total_gain) const {
// sort all rem alignments:

	clock_t st = clock();

	clock_t sort_rem_time = clock();
	std::map<pw_alignment, size_t, compare_pw_alignment> sorted_rem; // each remove alignment to a unique index in all_rem_als
	std::vector<const pw_alignment *> all_rem_als;
	for(size_t i=0; i<rem_als_per_ins.size(); ++i) {
		for(size_t j=0; j<rem_als_per_ins.at(i).size(); ++j) {
			const pw_alignment & al = rem_als_per_ins.at(i).at(j);
			std::map<pw_alignment, size_t, compare_pw_alignment>::iterator findal = sorted_rem.find(al);
			if(findal==sorted_rem.end()) {
				sorted_rem.insert(std::make_pair(al, sorted_rem.size()));
				all_rem_als.push_back(&al);
			}
		}
	}
	sort_rem_time = clock() - sort_rem_time;
	std::cout << " TIME select from groups, sort remove alignments " << (double)sort_rem_time/CLOCKS_PER_SEC << std::endl;

	clock_t sr_gain = clock();

// compute all gains
	std::vector<double> insert_gain(insert_als.size());
	for(size_t i=0; i<insert_als.size(); ++i) {
		double igain = 0;
		for(size_t j=0; j<insert_als.at(i).size(); ++j) {
			double g1, g2;
			model.gain_function(insert_als.at(i).at(j), g1, g2);
			double gain = (g1+g2)/2.0;
			if(gain < 0) gain = 0;
			igain += gain;
		}
		insert_gain.at(i) = igain;
//		std::cout << " i " << i << " g " << igain << std::endl;
	}
	std::vector<double> rem_gain(all_rem_als.size());
	for(size_t i=0; i<rem_gain.size(); ++i) {
		double g1, g2;
		model.gain_function(*(all_rem_als.at(i)), g1, g2);
		double gain = (g1+g2)/2.0;
		if(gain <0) gain = 0;
		rem_gain.at(i) = gain;
//		std::cout << " r " << i<< " g " << gain << std::endl;
	}

	sr_gain = clock() - sr_gain;
	std::cout << " TIME select from groups compute gains " << (double)sr_gain/CLOCKS_PER_SEC << std::endl;

	clock_t graph_time = clock();

// bipartite dependency graph
	std::vector<std::set<size_t> > ins_to_rem(insert_als.size());
	std::vector<std::set<size_t> > rem_to_ins(all_rem_als.size());
	size_t nume = 0;
	for(size_t i=0; i<rem_als_per_ins.size(); ++i) {
		for(size_t j=0; j<rem_als_per_ins.at(i).size(); ++j) {
			const pw_alignment & al = rem_als_per_ins.at(i).at(j);
			std::map<pw_alignment, size_t, compare_pw_alignment>::iterator findal = sorted_rem.find(al);
			assert(findal!=sorted_rem.end());
			size_t al_index = findal->second;
			ins_to_rem.at(i).insert(al_index);
			rem_to_ins.at(al_index).insert(i);
			nume++;

//			std::cout << " ins " << i << " rem " << al_index << std::endl;
		}
	}
//	std::cout << " we have " << nume << " edges " << std::endl;

	graph_time = clock() - graph_time;
	std::cout << " TIME select groups make graph " << (double)graph_time /CLOCKS_PER_SEC << std::endl;
// select nothing initially
	vector<bool> ins_selected(insert_als.size(), 0);
	vector<bool> rem_selected(all_rem_als.size(), 0);
// backtracker is allowed to use all remove alignments
	vector<bool> rem_free(all_rem_als.size(), 1);
	double gain = 0;
	double best_gain = 0;
	std::vector<bool> best_ins(insert_als.size(), 0);
	std::vector<bool> best_rem(all_rem_als.size(), 0);

	clock_t rec_time = clock();
	easy_insert(ins_to_rem, rem_to_ins, insert_gain, rem_gain, rem_free, ins_selected, rem_selected, gain);
//	std::cout << " gain after initial easy insert " << gain << std::endl;
// run backtracking algorithm to find best selection
	bb_step(ins_to_rem, rem_to_ins, insert_gain, rem_gain, rem_free, ins_selected, rem_selected, gain, best_gain, best_ins, best_rem);
	
	rec_time = clock() - rec_time;

	std::cout << " TIME select groups recursions " << (double) rec_time/CLOCKS_PER_SEC << std::endl;

// check result
	clock_t check_time = clock();
	gain = 0;
	for(size_t i=0; i<ins_selected.size(); ++i) {
		if(best_ins.at(i)) {
			gain+=insert_gain.at(i);
			for(size_t j=0; j<insert_als.at(i).size(); ++j) {
				result_ins.push_back(insert_als.at(i).at(j));
			}
		}
	}
	for(size_t i=0; i<rem_selected.size(); ++i) {
		if(best_rem.at(i)) {
			gain-=rem_gain.at(i);
			result_rem.push_back(*(all_rem_als.at(i)));
		}
	}
	check_time = clock() - check_time;
	std::cout << " TIME select groups check " << (double) check_time / CLOCKS_PER_SEC << std::endl;

//	std::cout << "XXX i taken " << result_ins.size() << " r taken " << result_rem.size() << " gain check " << gain << " " << best_gain << std::endl;

/*

	std::map<pw_alignment, size_t, compare_pw_alignment> sorted_rem; // each remove alignment to a unique index in all_rem_als
	std::vector<const pw_alignment *> all_rem_als;
	for(size_t i=0; i<rem_als_per_ins.size(); ++i) {
		for(size_t j=0; j<rem_als_per_ins.at(i).size(); ++j) {
			const pw_alignment & al = rem_als_per_ins.at(i).at(j);
			std::map<pw_alignment, size_t, compare_pw_alignment>::iterator findal = sorted_rem.find(al);
			if(findal==sorted_rem.end()) {
				sorted_rem.insert(std::make_pair(al, sorted_rem.size()));
				all_rem_als.push_back(&al);
			}
		}
	}
	
// compute all gains
	std::vector<double> insert_gain(insert_als.size());
	for(size_t i=0; i<insert_als.size(); ++i) {
		double igain = 0;
		for(size_t j=0; j<insert_als.at(i).size(); ++j) {
			double g1, g2;
			model.gain_function(insert_als.at(i).at(j), g1, g2);
			double gain = (g1+g2)/2.0;
			igain += gain;
		}
		insert_gain.at(i) = igain;
	}
	std::vector<double> rem_gain(all_rem_als.size());
	for(size_t i=0; i<rem_gain.size(); ++i) {
		double g1, g2;
		model.gain_function(*(all_rem_als.at(i)), g1, g2);
		double gain = (g1+g2)/2.0;
		rem_gain.at(i) = gain;
	}

// bipartite dependency graph
	std::vector<std::set<size_t> > ins_to_rem(insert_als.size());
	std::vector<std::set<size_t> > rem_to_ins(all_rem_als.size());
	for(size_t i=0; i<rem_als_per_ins.size(); ++i) {
		for(size_t j=0; j<rem_als_per_ins.at(i).size(); ++j) {
			const pw_alignment & al = rem_als_per_ins.at(i).at(j);
			std::map<pw_alignment, size_t, compare_pw_alignment>::iterator findal = sorted_rem.find(al);
			assert(findal!=sorted_rem.end());
			size_t al_index = findal->second;
			ins_to_rem.at(i).insert(al_index);
			rem_to_ins.at(al_index).insert(i);

//			std::cout << " ins " << i << " rem " << al_index << std::endl;
		}
	}


	vector<bool> ins_selected(insert_als.size(), 1);
	vector<bool> rem_selected(all_rem_als.size(), 1);

// At first, we select all insert alignments
	total_gain =0;
	for(size_t i=0; i<insert_als.size(); ++i) {
		total_gain += insert_gain.at(i);

//		std::cout << " i gain " << insert_gain.at(i) << std::endl;
	}
// This causes that we need to remove all remove alignments
	for(size_t i=0; i<all_rem_als.size(); ++i) {
		total_gain -= rem_gain.at(i);

//		std::cout << " r gain " << rem_gain.at(i) << std::endl;
	}

//	std::cout << " in select from groups " << total_gain << " total gain from " << insert_als.size() << " ins groups and " << all_rem_als.size() << " rem alignments " << std::endl;
// can we now increase total gain by not selecting some individual or pairs of remove alignments?
// this is an approximation, the correct optimal solution would be to iterate over all subsets which is too slow


	bool changed = true;
	while(changed) { // if we managed to remove something, we will try all options again
		changed = false;
		// we try to remove a single remove alignment (and all insert alignments that can be combined with it to get positive gain)
		for(size_t i=0; i<all_rem_als.size(); ++i) {
			if(rem_selected.at(i)) {
//				std::cout << " try unselect " << i << std::endl;
				std::set<size_t> torem;
				torem.insert(i);
				double thisgain;
				bool res = try_rnodes_unselect(ins_to_rem, rem_to_ins, insert_gain, rem_gain, torem, ins_selected, rem_selected, thisgain);
				if(res) {
//					std::cout << " res " << thisgain << std::endl;
					total_gain += thisgain;
					changed = true;
				}
			}
		}


		// we try to use a pair of remove alignments (and all insert alignments that can be combined with them to get positve gain)
		for(size_t i=0; i<all_rem_als.size(); ++i) {
			for(size_t j=i+1; j<all_rem_als.size(); ++j) {
				if(rem_selected.at(i) && rem_selected.at(j)) {
//				std::cout << " try unselect " << i<< " "<< j<< std::endl;
					std::set<size_t> torem;
					torem.insert(i);
					torem.insert(j);
					double thisgain;
					bool res = try_rnodes_unselect(ins_to_rem, rem_to_ins, insert_gain, rem_gain, torem, ins_selected, rem_selected, thisgain);
					if(res) {
//						std::cout << "pair res " << thisgain << std::endl;
						total_gain += thisgain;
						changed = true;
					}
				}				
			}
		}

		// triplet as well (no more than 3 for speed)
		for(size_t i=0; i<all_rem_als.size(); ++i) {
			for(size_t j=i+1; j<all_rem_als.size(); ++j) {
				for(size_t k=j+1; k<all_rem_als.size(); ++k) {
					if((rem_selected.at(i) && rem_selected.at(j)) && rem_selected.at(k)) {
//						std::cout << " try unselect " << i <<" "<< j << " " << k<< std::endl;
						std::set<size_t> torem;
						torem.insert(i);
						torem.insert(j);
						torem.insert(k);
						double thisgain;
						bool res = try_rnodes_unselect(ins_to_rem, rem_to_ins, insert_gain, rem_gain, torem, ins_selected, rem_selected, thisgain);
						if(res) {
//							std::cout << "triple res " << thisgain << std::endl;
							total_gain += thisgain;
							changed = true;
						}
					}
				}
			}
		}


	}


	double check_gain = 0;
	for(size_t i=0; i<rem_selected.size(); ++i) {
		if(rem_selected.at(i)) {
			result_rem.push_back(*(all_rem_als.at(i)));
			check_gain -= rem_gain.at(i);
		}
	}
	size_t igroup_sel = 0;
	for(size_t i=0; i<insert_als.size(); ++i) {
		if(ins_selected.at(i)) {
			igroup_sel++;
			for(size_t j=0; j<insert_als.at(i).size(); ++j) {
				result_ins.push_back(insert_als.at(i).at(j));
			}
			check_gain += insert_gain.at(i);
		}
	}

//	std::cout << " select from groups result " << total_gain << " " << check_gain<< " total gain from " << igroup_sel << " ins groups and " << result_rem.size() << " rem alignments " << std::endl;
*/
	st = clock() - st;
	select_groups_time += (double) st / CLOCKS_PER_SEC;
}







template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::lazy_split_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const std::vector<pw_alignment>  & als, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {
	rec_calls++; // count number of calls
#if PRINT
//	std::cout<< "al in lazy split recursion level " << level << " we have " << als.size() << " alignments"<< std::endl;
	for(size_t i=0; i<als.size(); ++i) {
		als.at(i).print();
	}
#endif
	clock_t group_gain_time = clock();

	size_t orig_ovrlp_size = ovrlp.get_all().size();
// slow test
//	ovrlp.test_partial_overlap();


	std::vector<pw_alignment> insert_als;// insert and remove for next recursion
	std::vector<pw_alignment> remove_als;

// at first, we try to insert all at once
	double group_gain_in = 0;
	std::vector< const pw_alignment *> alsp(als.size());
	for(size_t i=0; i<als.size(); ++i) {
		double g1, g2;
		model.gain_function(als.at(i), g1, g2);
		group_gain_in += (g1+g2)/2.0;
		alsp.at(i) = &(als.at(i));
	}


	group_gain_time = clock() - group_gain_time;
	std::cout << " TIME group gain " << (double)group_gain_time/CLOCKS_PER_SEC<<std::endl;

	clock_t split_time = clock();

	splitpoints_interval_tree<overlap_type> spl1(alsp, ovrlp, data);
	spl1.nonrecursive_splits();
	std::set<pw_alignment, compare_pw_alignment> als_remove_als_set; 
	std::cout << "spl1"<<std::endl;
	spl1.split_all(als_remove_als_set, insert_als);
	
	split_time = clock() - split_time;
	std::cout << " TIME split " << (double)split_time/CLOCKS_PER_SEC << std::endl;
	
	clock_t ins_rem_time = clock();

	double als_rem_gain = 0;
	double als_ins_gain = 0;
	for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = als_remove_als_set.begin(); it!=als_remove_als_set.end(); ++it) {
		remove_als.push_back(*it);
		const pw_alignment & ral = *it;
		double g1, g2;
		model.gain_function(ral, g1, g2);
		als_rem_gain += (g1+g2)/2.0;	
	}
	for(size_t i=0; i<insert_als.size(); ++i) {
		const pw_alignment & ial = insert_als.at(i);
		double g1, g2;
		model.gain_function(ial, g1, g2);
		als_ins_gain += (g1+g2)/2.0;
	}
	double group_gain = 0;
// how much gain can we keep, when we take the entire group
	double als_total_gain = als_ins_gain - als_rem_gain;

	ins_rem_time = clock() - ins_rem_time;
	std::cout << " TIME ins rem gains " << (double) ins_rem_time/CLOCKS_PER_SEC << std::endl;

//	std::cout <<" l " << level << " groups size " << als.size()<<   " gain in " << group_gain_in << " split groups gain " << als_total_gain << std::endl;
// if we keep less than 90%, we try to improve it using the slow select_from_groups function which takes only some of the alignments of the group
	if(group_gain_in * 0.9 > als_total_gain) {
		std::cout << " REDO " << std::endl;
		
		clock_t groups_time = clock();

		// remove results from group split step
		insert_als.clear();
		remove_als.clear();

		// inserting i from als into overlap would lead to remove remove_groups.at(i) from overlap and insert insert_groups.at(i) into overlap
		std::vector<std::vector<pw_alignment> > insert_groups;
		std::vector<std::vector<pw_alignment> > remove_groups;
		for(size_t i=0; i<als.size(); ++i) {

			const pw_alignment & al = als.at(i);
	
			splitpoints_interval_tree<overlap_type> spl(al, ovrlp, data);
			// nonrecursive splits only eliminate direct partial overlap, we may need to go on in further steps of lazy_split_insert
			spl.nonrecursive_splits();
			// sets of alignments that need to be removed and inserted if we want the current alignment 
			std::set<pw_alignment, compare_pw_alignment> al_remove_als_set; 
			std::vector<pw_alignment> al_insert_als;
			std::vector<pw_alignment> al_remove_als;
			spl.split_all(al_remove_als_set, al_insert_als);
			for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = al_remove_als_set.begin(); it!=al_remove_als_set.end(); ++it) {
				al_remove_als.push_back(*it);
			}
			insert_groups.push_back(al_insert_als);
			remove_groups.push_back(al_remove_als);
		}
	
		groups_time = clock() - groups_time;
		std::cout << " TIME make groups " << (double) groups_time / CLOCKS_PER_SEC << std::endl;

		clock_t sgroups_time = clock();

		// select a subset of als such that inserting it into o, yields the most gain
		select_from_groups(insert_groups, remove_groups, insert_als, remove_als, group_gain);

		sgroups_time = clock() - sgroups_time;
		std::cout << " TIME select groups " << (double) sgroups_time / CLOCKS_PER_SEC << std::endl;
//		std::cout << " select groups gain " << group_gain << std::endl;
	} else {
		group_gain = als_total_gain;
	}
//	std::cout << " GROUP gain " << group_gain <<  std::endl;


//	std::cout << " ins "<< insert_als.size() << " rem " << remove_als.size() << std::endl;

	local_gain = 0;

// we continue if information gain is possible from the current group
	if(group_gain > 0) {

	
#if PRINT
		std::cout <<"initial: level " << level << " positive gain: " << group_gain << " split res rem " << remove_als.size() << " ins " << insert_als.size() << std::endl;
#endif
		if(remove_als.empty() ) {
#if PRINT
			std::cout << " no al to delete "<<std::endl;
#endif

			clock_t oparts_time = clock();
//	Nothing to remove: Order parts by gain and insert into overlap using fully recursive alignment splits
			std::multimap<double, const pw_alignment *> ordered_parts;
			double max_parts_gain = 0;
			for(size_t i=0; i<insert_als.size(); i++) {
				const pw_alignment * al = & insert_als.at(i);
				double g1, g2;
				model.gain_function(insert_als.at(i), g1, g2);
				double iav_gain = (g1 + g2)/2;
				if(iav_gain > 0) {
					max_parts_gain += iav_gain;
					ordered_parts.insert(std::make_pair(iav_gain, al));
				}
//				std::cout << "INS  gain " <<iav_gain << std::endl;
//				insert_als.at(i).print();
//				std::cout << std::endl;
			}
			double children_gain = 0;
			std::set<pw_alignment, compare_pw_alignment> l_inserted;
			std::set<pw_alignment, compare_pw_alignment> l_removed;

			oparts_time = clock() - oparts_time;
			std::cout << " TIME no remove order parts " << (double) oparts_time/CLOCKS_PER_SEC << std::endl;

			for(std::multimap<double, const pw_alignment *>::reverse_iterator rit = ordered_parts.rbegin(); rit!=ordered_parts.rend(); ++rit) {
				pw_alignment nal = *(rit->second);
				double lgain = 0;
				std::vector<pw_alignment> local_inserted; // these alignments were actually inserted/removed by child calls
				std::vector<pw_alignment> local_removed;
				// full split step on single alignments with recursive split points that guarantee to eliminate all partial overlap
				lazy_split_full_insert_step(ovrlp, level, rec_calls, nal, local_inserted, local_removed, lgain);
				children_gain += lgain;
				insert_alignment_sets(l_inserted, l_removed, local_inserted, local_removed);
			}
			all_push_back(inserted_alignments, removed_alignments, l_inserted, l_removed);
			local_gain = children_gain;
			assert(local_gain>=0); // remove_als was empty, therefore no loss of gain, no local_undo
			assert(ovrlp.get_all().size() + removed_alignments.size() == orig_ovrlp_size + inserted_alignments.size());
			return;
		}


// remove all

		clock_t rtime = clock();
		std::set<pw_alignment, compare_pw_alignment> all_inserted;
		std::set<pw_alignment, compare_pw_alignment> all_removed;
		size_t remove_count = 0;
		double removed_gain = 0;
		for(size_t i=0; i<remove_als.size(); ++i) {
			const pw_alignment & pwa = remove_als.at(i);
			all_removed.insert(pwa); 
			double g1, g2;
			model.gain_function(pwa, g1, g2);
			double gain = (g1+g2)/2.0;
			removed_gain+=gain;
	//			std::cout <<"level " << level << "REMOVE NOW: " << std::endl;
	//			pwa->print();
	//			std::cout << " ovrlp size " << ovrlp.size() << std::endl;
		//		std::cout << "is gonna be removed5: " << &pwa<<std::endl;
		//		pwa.print();

			ovrlp.remove_alignment(pwa);

		
			remove_count++;	
		}

		rtime = clock() - rtime;
		std::cout << " TIME remove " << (double) rtime /CLOCKS_PER_SEC << std::endl;
// 
//	ovrlp.test_partial_overlap();

// recursive insert
		double next_group_gain = 0;
		std::vector<pw_alignment> this_inserted;
		std::vector<pw_alignment> this_removed;
		lazy_split_insert_step(ovrlp, level + 1, rec_calls,  insert_als, this_inserted, this_removed, next_group_gain);
		insert_alignment_sets(all_inserted, all_removed, this_inserted, this_removed);

		if(next_group_gain > removed_gain) {
			all_push_back(inserted_alignments, removed_alignments, all_inserted, all_removed);

			local_gain = next_group_gain - removed_gain;

		} else {
//			std::cout << " gr gain was " << next_group_gain << " removed gain " << removed_gain << std::endl;
//			std::cout << " UNDO ins " << all_inserted.size() << " rem " << all_removed.size() << std::endl;
			local_undo(ovrlp, all_inserted, all_removed);
			local_gain = 0;
		}
	}

// TODO exclude double insert into overlp events to make this assertion work
//	assert(ovrlp.get_all().size() + removed_alignments.size() == orig_ovrlp_size + inserted_alignments.size());

}


/*

insert al recursively
each insert operation causes the removal of some alignments from overlap
then we split al and the removed alignments (without induced splits)
the gain of all new pieces with positive gain - the gain of all removed pieces is an upper bound of 
the information gain of inserting al. If this upper bound is positive we continue and try to recursively insert all pieces
Afterwards, we have to compute the actual gain. If it is negative we undo the local recursive insert operation

*/
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::lazy_split_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment  & al, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {
// start with computing gain of al
	rec_calls++; // count number of calls
	double gain1;
	double gain2;
	common_model.gain_function(al, gain1, gain2);
	double av_al_gain = (gain1 + gain2) / 2 - base_cost;
	// we continue if information gain is possible from the current alignment
	local_gain = 0;
#if PRINT
	std::cout<< "al in lazy split recursion level " << level << std::endl;
	al.print();
#endif
//	std::cout << "av_al_gain "<< av_al_gain <<std::endl;
	if(av_al_gain > 0) {
		splitpoints_interval_tree<overlap_type> spl(al, ovrlp, data);
		spl.nonrecursive_splits();
		// sets of alignments that need to be removed and inserted if we want the current alignment 
		std::set<pw_alignment, compare_pw_alignment> remove_als; // remove alignments are pointers to objects contained in the overlap structure
		std::vector<pw_alignment> insert_als;
	//	std::cout << "split all is called 1! "<<std::endl;

		spl.split_all(remove_als, insert_als);
#if PRINT
		std::cout <<"initial: level " << level << " positive gain: " << av_al_gain << " split res rem " << remove_als.size() << " ins " << insert_als.size() << std::endl;
#endif
		// we evaluate an upper bound of information gain for the current alignment
		double lost_information_gain = 0;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			const pw_alignment & ral = *it;
			common_model.gain_function(ral, gain1, gain2);
			double rav_gain = (gain1 + gain2)/2 - base_cost;
			lost_information_gain += rav_gain;

//			std::cout << "REM " << std::endl;
//		       	ral->print();
//			std::cout << std::endl;	
		}
		std::multimap<double, const pw_alignment *> ordered_parts;
		double max_parts_gain = 0;
		for(size_t i=0; i<insert_als.size(); i++) {
			const pw_alignment * al = & insert_als.at(i);
			common_model.gain_function(insert_als.at(i), gain1, gain2);
			double iav_gain = (gain1 + gain2)/2 - base_cost;
			if(iav_gain > 0) {
				max_parts_gain += iav_gain;
				ordered_parts.insert(std::make_pair(iav_gain, al));
			}
		//	std::cout << "INS  gain " <<iav_gain << std::endl;
		//	insert_als.at(i).print();
		//	std::cout << std::endl;
		}

//		std::cout << "splits: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " causes to remove " << remove_als.size() << " alignments with " << lost_information_gain << 
//			" gain. We could insert " << ordered_parts.size() << " small pieces. Upper bound of gain is " << max_parts_gain << std::endl;

		// here: we actually decide to insert a small part (in the next function we check more carefully for indirectly induced splits)
		if(remove_als.empty() && (ordered_parts.size() == 1)) {
#if PRINT
			std::cout << "Only 1 al to insert and no al to delete "<<std::endl;
#endif
		/*	for(std::multimap<double , const pw_alignment &>::iterator it = ordered_parts.begin(); it!= ordered_parts.end();it++){
				std::pair<std::multimap<double, const pw_alignment&>::iterator, std::multimap<double, const pw_alignment&>::iterator  > r2 = ordered_parts.equal_range(it->first);
				for(std::multimap<double, const pw_alignment &>::iterator it1 = r2.first; it1!=r2.second; ++it1) {
					const pw_alignment & al = it1->second;
					al.print();
					std::cout << &al <<std::endl;
				}
			}*/
			const pw_alignment & nal = *ordered_parts.begin()->second;
		//	nal.print();
			// full split step, we can just copy local gain
			lazy_split_full_insert_step(ovrlp, level, rec_calls, nal, inserted_alignments, removed_alignments, local_gain);
			return;
		}





		// if it is possible to increase information gain by inserting this alignment, we 
		// recursively insert all small pieces
		// we start with the high gain pieces, because they have more potential to loose information gain, which can be used to abort the current insert
		double upper_bound_of_gain = max_parts_gain - lost_information_gain;
		
		if(upper_bound_of_gain > 0) {
			double all_lost_gain = 0; // difference between gain realized by inserting a piece versus the gain of that piece
			double total_recursive_gain = 0; // total gain of all subtrees processed so far

			// all_ins/rem accumulated over all recursive steps below current. If we take current step, those are tranferred to inserted_alignments/removed_alignments
			std::set<pw_alignment, compare_pw_alignment> all_inserted;
			std::set<pw_alignment, compare_pw_alignment> all_removed;
//			std::cout << "size of all inserted "<<all_inserted.size()<<std::endl;
			// remove before recursion
			size_t remove_count = 0;
			for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); it++) {
				const pw_alignment & pwa = *it;
				all_removed.insert(pwa); 			
	//			std::cout <<"level " << level << "REMOVE NOW: " << std::endl;
	//			pwa->print();
	//			std::cout << " ovrlp size " << ovrlp.size() << std::endl;
		//		std::cout << "is gonna be removed5: " << &pwa<<std::endl;
		//		pwa.print();

				ovrlp.remove_alignment(pwa);


				remove_count++;
			}

	//		std::cout << " before recursion, we have removed " << remove_count << " alignments with gain " << lost_information_gain <<  " upper bound of gain is " << upper_bound_of_gain << std::endl;
			size_t alnum = 0;
			for(std::map<double, const pw_alignment *>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); it++) {
				double ub_of_thisgain = it->first;
				double thisgain;
				std::vector<pw_alignment> this_inserted;
				std::vector<pw_alignment> this_removed;
				const pw_alignment & thisal = *it->second;
//				std::cout <<"this al "<<std::endl;
//				thisal.print();
			//	std::cout << std::endl;

	//			std::cout <<"start recursion: level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << std::endl; 
				// this function only changes ovrlp if this was locally good. We may need to undo the changes if they were globally bad
				lazy_split_insert_step(ovrlp, level + 1, rec_calls,  thisal, this_inserted, this_removed, thisgain);
				double lost_gain = ub_of_thisgain - thisgain; // we have stepped down that much from the previous upper bound of gain
				all_lost_gain += lost_gain;
				total_recursive_gain += thisgain;

		//		for(size_t i =0; i < this_inserted.size();i++){
				//	this_inserted.at(i).print();
		//		}
		//		std::cout <<"end recursion: level "<< level << " end rec "<< alnum << " on al length " << thisal.alignment_length() << " subtree local gain " << thisgain << " ub was " << ub_of_thisgain << " lost gain is " << lost_gain << " total lost gain is " << all_lost_gain << " total gain of this level until now: " << total_recursive_gain <<  " ub is " << upper_bound_of_gain << std::endl; 
			
				insert_alignment_sets(ovrlp, all_inserted, all_removed, this_inserted, this_removed);

				// if we have lost more than we may gain, we undo the entire recursive insertion subtree from here
				if(all_lost_gain > upper_bound_of_gain) {
			//		std::cout << "local undo "<<std::endl;
					local_undo(ovrlp, all_inserted, all_removed);


		//			std::cout << "undo subtree: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " ABORT because of low upper bound of gain" << std::endl;
					// gain of this subtree is 0 because we do not insert
					local_gain = 0;


					return;
				} 
				alnum++;
			} // for ordered_parts

			// we take all alignments underneath of current recursion level if they have more gain than all the one that we had to remove for the current split operation
			if(total_recursive_gain > lost_information_gain) {
				local_gain = total_recursive_gain - lost_information_gain;

				all_push_back(inserted_alignments, removed_alignments, all_inserted, all_removed);
		//		std::cout << "all push back"<<std::endl;

//	std::cout << "take subtree: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " causes to remove " << all_removed.size() << " alignments with " << lost_information_gain << 
//			" gain. We want to insert " << all_inserted.size() << " small pieces. TOTAL gain is " << local_gain << std::endl;

				return;
			
			} else {
				local_gain = 0;
				local_undo(ovrlp, all_inserted, all_removed);
			}
		} // upper bound gain pos

	


	} // pos gain

}

/* TODO
	maybe we can increase the amount of gain taken by 
	integrating all the groups that were not taken into an empty overlap structure and then inserting its groups in the real o
*/

template<typename T, typename overlap_type>
void initial_alignments_from_groups<T,overlap_type>::compute(overlap_type & o, std::string ihead,  std::string & info){
// TODO why does compute take so long on leftover alignments?
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;

	double longtime = 0;
	size_t longtime_at = 0;
	double long_groupstime = 0;

	size_t long_groupsize = 0;
	size_t long_groupins = 0;
	size_t long_grouprem = 0;
	double long_groupgain = 0;

	std::cout << " IAS COMPUTE START, we have " << sorted_original_als.size() << " groups " << std::endl;
	for (size_t i=0; i<sorted_original_als.size(); ++i) {

		clock_t ins_time = clock();
	
		std::vector<pw_alignment> als;
		for(size_t j=0; j<sorted_original_als.at(i).size(); ++j) {
			const pw_alignment & al = *(sorted_original_als.at(i).at(j));
			als.push_back(al);
		}
		double gain_of_als = 0;


		



#if PRINT
		std::cout << "on alignment group " << i << " with " << als.size() << " alignments"<< std::endl;
		cout << " start recursive lazy insert tree on " << endl;
#endif
		for(size_t j=0; j<als.size(); ++j) {
			double g1, g2;
			model.gain_function(als.at(j), g1, g2);
			double gain = (g1+g2)/2.0;
			gain_of_als+=gain;	
#if PRINT
			(als.at(j)).print();
			std::cout << std::endl;
#endif

		}
#if PRINT
		std::cout << "group gain " << gain_of_als << std::endl;
#endif

/*		std::cout << " Current overlap content " << std::endl;
		const std::set<pw_alignment, compare_pw_alignment> oar = o.get_all();
		for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = oar.begin(); it!=oar.end(); ++it) {
			it->print();
		}
*/

		clock_t ltime = clock();

		std::vector<pw_alignment> ins_als;
		std::vector<pw_alignment> rem_als;
		size_t count_rec_calls = 0;
		double gain_taken = 0;
		select_groups_time = 0;
		lazy_split_insert_step(o,0, count_rec_calls, als, ins_als, rem_als, gain_taken);//Splits alignments that have overlap with 'al'
		
		ltime = clock() - ltime;
		double lts = (double)ltime / CLOCKS_PER_SEC;

		if(lts > 10) {
			longtime = lts;
			longtime_at = i;
			long_groupstime = select_groups_time;
			long_groupsize = als.size();
			long_groupins = ins_als.size();
			long_grouprem = rem_als.size();
			long_groupgain = gain_of_als;
		}

		std::stringstream str;
		str <<ihead << ", ias compute " << i << "of "<< sorted_original_als.size() << " we had " << als.size() << " als in " << ins_als.size() << " out " << rem_als.size() <<  " gain " << gain_taken << " of " << gain_of_als << " time " << lts << " last long time " << longtime << " at " << longtime_at<< " in select gr " << long_groupstime << " size " << long_groupsize << " ins " << long_groupins << " rem " << long_grouprem << " allgain " << long_groupgain<<  " in overlap " << o.get_all().size() ;
#pragma omp critical(scheduling)
{
		info = str.str();
}


		double rgain = 0;
		double igain = 0;
		double gain1, gain2;
		for(size_t j=0; j<ins_als.size(); ++j) {
			model.gain_function(ins_als.at(j), gain1, gain2);
			double vgain = (gain1+gain2)/2;
			igain+=vgain;
		}
		for(size_t j=0; j<rem_als.size(); ++j) {
			model.gain_function(rem_als.at(j), gain1, gain2);
			double vgain = (gain1+gain2)/2;
			rgain+=vgain;
		}
		

		double checkgain = igain - rgain;


		if(!ins_als.empty()) {
			used++;
			total_gain+=gain_of_als;
			pcs_ins+=ins_als.size();
			pcs_rem+=rem_als.size();
			std::cout << " alignment taken "<< std::endl;
		} else {
			not_used++;
			std::cout << " alignment not taken " << std::endl;
		}
		
		std::cout << " on alignments " << i << " number " << als.size() << " gain " << gain_taken <<" of " << gain_of_als << " recursive lazy split insert results:" << std::endl;
		std::cout << " al " << i << " gain in " << igain << " gain out " << rgain << " check gain (ins - rem) " << checkgain << std::endl;
		std::cout << " out: " << rem_als.size() << " gain " << rgain << " in: " << ins_als.size() << " gain " << igain << std::endl;
		std::cout << " total gain until here: " << total_gain << " needed " << count_rec_calls<< " recursive lazy split insert calls " << std::endl;


		ins_time = clock() - ins_time;
		std::cout << " INSERT "<< i << " time " << (double)ins_time/CLOCKS_PER_SEC<< " lazy time " << (double)ltime/CLOCKS_PER_SEC << std::endl;
//		std::cout << " CHECK PO on " << o.get_all().size() << std::endl;
//		o.test_partial_overlap();
		


	}
	result_gain = total_gain;
	used_alignments = used;
	std::cout << "Input size " << sorted_original_als.size() << " max gain " << max_gain << std::endl;
	std::cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << std::endl;
	std::cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << std::endl;




}


template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::compute_simple_lazy_splits(overlap_type & o){
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;

	for (size_t i=0; i<sorted_original_als.size(); ++i) {//XXX Can we optimize it a bit by not walking on all these als? Since we may already cut them.
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;
#if PRINT
		std::cout << "on alignment " << i << std::endl;
		cout << " start recursive lazy insert tree on " << endl;
	//	al->print();
	//	cout << endl;
#endif

		std::vector<pw_alignment> ins_als;
		std::vector<pw_alignment> rem_als;
		size_t count_rec_calls = 0;
		lazy_split_insert_step(o,0, count_rec_calls, *al, ins_als, rem_als, gain_of_al);//Splits alignments that have overlap with 'al'
		
	//	double gain1;
	//	double gain2;
	//	common_model.gain_function(*al, gain1, gain2);
	//	double vgain = (gain1 + gain2) / 2 - base_cost;


/*//XXX Commented double rgain = 0;
		double igain = 0;
		for(size_t j=0; j<ins_als.size(); ++j) {
			common_model.gain_function(ins_als.at(j), gain1, gain2);
			vgain = (gain1+gain2)/2 - base_cost;
			igain+=vgain;
		}
		for(size_t j=0; j<rem_als.size(); ++j) {
			common_model.gain_function(rem_als.at(j), gain1, gain2);
			vgain = (gain1+gain2)/2 - base_cost;
			rgain+=vgain;
		}
		

		double checkgain = igain - rgain;*/

		if(!ins_als.empty()) {
			used++;
			total_gain+=gain_of_al;
			pcs_ins+=ins_als.size();
			pcs_rem+=rem_als.size();
		//	std::cout << " alignment taken "<< std::endl;
		} else {
			not_used++;
	//		std::cout << " alignment not taken " << std::endl;
		}
		
//		std::cout << " on alignment " << i << " length " << al.alignment_length() << " gain " << vgain << " recursive lazy split insert results:" << std::endl;
//		std::cout << " al " << i << " gain in " << vgain << " gain out " << gain_of_al << " check gain (ins - rem) " << checkgain << std::endl;
//		std::cout << " out: " << rem_als.size() << " gain " << rgain << " in: " << ins_als.size() << " gain " << igain << std::endl;
//		std::cout << " total gain until here: " << total_gain << " needed " << count_rec_calls<< " recursive lazy split insert calls " << std::endl;
		
		
	}
	result_gain = total_gain;
	used_alignments = used;
//	std::cout << "Input size " << sorted_original_als.size() << " max gain " << max_gain << std::endl;
//	std::cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << std::endl;
//	std::cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << std::endl;



}



/**
	greedy test function 
**/
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::compute_simple(overlap_type & o) {
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;
//	std::cout<< "sorted alignment size" << sorted_original_als.size()<<std::endl;	
	for(size_t i=0; i<sorted_original_als.size(); i++) {
		const pw_alignment & al = sorted_original_als.at(i);
		double gain_of_al = 0;

		// TODO remove
		double gain1, gain2;
		common_model.gain_function(al, gain1, gain2);
		gain1-=base_cost;
//		std::cout << std::endl<<"at alignment " << i << " length " << al->alignment_length() << " al base gain " << gain1 << std::endl;
		//al->print();
//		std::cout << std::endl;


		std::set<pw_alignment, compare_pw_alignment> remove_als;
		std::vector<pw_alignment> insert_als;

//		// TODO
//		std::cout << " splitpoints on: " << al << " at " << i <<  " size " <<sorted_original_als.size()<< std::endl; 
//		al->print();
//		std::cout << std::endl;	
//		o.print_all_alignment();


		splitpoints_interval_tree<overlap_type> spl(al, o, data);
		spl.recursive_splits();
	//	std::cout << "split all is called! "<<std::endl;
		spl.split_all(remove_als, insert_als);
		std::vector<double> insert_gains(insert_als.size(), 0);

		for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			double g1;
			double g2;
			common_model.gain_function(*it, g1, g2);
		//	if(g2<g1) g1 = g2;
			g1-=base_cost;
	//		std::cout << "r " << (*it)->alignment_length() << " g " << g1 << std::endl;
			gain_of_al -= g1;

		//	(*it)->print();
		//	std::cout << std::endl;

		}	
		for(size_t j=0; j<insert_als.size(); ++j) {
			
			double g1;
			double g2;
			common_model.gain_function(insert_als.at(j), g1, g2);
		//	if(g1<g2) g1 = g2;
			g1-=base_cost;
			insert_gains.at(j) = g1;
			if(g1>0) {
		//		std::cout << "i " << (insert_als.at(j)).alignment_length() << " g " << g1 << std::endl;
				gain_of_al += g1;

			//	insert_als.at(j).print();
			//	std::cout << std::endl;
			}

		}
//		std::cout << " al " << i << " rem " << remove_als.size() << " insert " << insert_als.size() << " gain " << gain_of_al << std::endl;	
		if(gain_of_al>=0) {
			used++;
			for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			//	std::cout << "is gonna be removed4: "  << &*it<<std::endl;
				pw_alignment &al=*it;
			//	al.print();
				o.remove_alignment(*it);
				pcs_rem++;
			}	
			for(size_t j=0; j<insert_als.size(); j++) {
				if(insert_gains.at(j)>0) {
				o.insert_without_partial_overlap(insert_als.at(j));
				pcs_ins++;
				}
			}
// slow check
//			o.test_partial_overlap();
			total_gain+=gain_of_al;

		} else {
			not_used++;
		}

	
	}
	result_gain = total_gain;
//	std::cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << std::endl;
//	std::cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << std::endl;

}



template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::search_area(size_t from, size_t to, const std::multimap<size_t, size_t> & data, std::set<size_t> & res) {
	std::multimap<size_t, size_t>::const_iterator search_begin = data.lower_bound(from);
	std::multimap<size_t, size_t>::const_iterator search_end = data.upper_bound(to);

//	cout << " search fr " << from << " to " << to << endl;
	for(std::multimap<size_t, size_t>::const_iterator it = search_begin; it!= search_end; ++it) {
		res.insert(it->second);
//		cout << " found " << it->first << " al " << it->second << endl;
	}
}


template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::overlap_graph(std::vector<std::set<size_t> > & edges) {

	std::vector<std::multimap<size_t, size_t> > all_left_pos(data.numSequences()); // reference -> alignment left pos, alignment index

	// all alignment left pos into a vector of multimaps
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		size_t left, right;
		al->get_lr1(left, right);
		all_left_pos.at(r1).insert(std::make_pair(left, i));
		al->get_lr2(left, right);
		all_left_pos.at(r2).insert(std::make_pair(left, i));
	}

	edges = std::vector<std::set<size_t> >(sorted_original_als.size()); // has to be completed with other direction

	// for each alignment: if there is overlap, we will find a left pos of another alignment between it left and right (or the other way around)
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		size_t left1, right1, left2, right2;
	//	cout << " on al " << i << endl;

		al->get_lr1(left1, right1);
		search_area(left1, right1, all_left_pos.at(r1), edges.at(i));
		al->get_lr2(left2, right2);
		search_area(left2, right2, all_left_pos.at(r2), edges.at(i));
	}

	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		// self edge remove:
		edges.at(i).erase(i);
		// insert reversed edges
		for(std::set<size_t>::iterator it = edges.at(i).begin(); it!=edges.at(i).end(); ++it) {
			size_t j = *it;
			edges.at(j).insert(i);
		}
	}

	size_t edge_count = 0;
	for(size_t i=0; i<edges.size(); ++i) {
		edge_count += edges.at(i).size();
	}

//	cout << " on a graph with " << edges.size() << " nodes and " << edge_count/2 << " edges" << endl;


}

template<typename T, typename overlap_type>

void initial_alignment_set<T,overlap_type>::remove_val(size_t index, std::map<size_t, double> & index_to_weight,  std::multimap<double, size_t> & weight_to_index) {
	std::map<size_t, double>::iterator findi = index_to_weight.find(index);
	assert(findi!=index_to_weight.end());
	double weight = findi->second;
	index_to_weight.erase(findi);
	std::pair<std::multimap<double, size_t>::iterator, std::multimap<double, size_t>::iterator > eqr = weight_to_index.equal_range(weight);
	for(std::multimap<double, size_t>::iterator it = eqr.first; it!=eqr.second; ++it) {
		size_t thisindex = it->second;
		if(thisindex == index) {
			weight_to_index.erase(it);
			return;
		}
	}
	assert(false);
}

/* 
	node: index was added to vertex cover.

	here we remove it from maps and update the weights

*/
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::update_remove(size_t index, std::vector<std::set<size_t> > & edges, std::map<size_t, double> & index_to_weight, std::multimap<double, size_t> & weight_to_index) {
	assert(index_to_weight.size() == weight_to_index.size());
	size_t oldsize = index_to_weight.size();
	std::map<size_t, double>::iterator findi = index_to_weight.find(index);
	assert(findi!=index_to_weight.end());
	double iweight = findi->second;

	remove_val(index, index_to_weight, weight_to_index);
	assert(index_to_weight.size() == oldsize -1 );
	assert(weight_to_index.size() == oldsize -1 );

	std::vector<size_t> iedges(edges.at(index).size());
	std::vector<double> weights(edges.at(index).size());
	std::vector<size_t> degrees(edges.at(index).size());
	size_t at = 0;
	for(std::set<size_t>::iterator it = edges.at(index).begin(); it!=edges.at(index).end(); ++it) {
		size_t j = *it;
		std::map<size_t, double>::iterator findj = index_to_weight.find(j);
		if(findj!=index_to_weight.end()) {
			iedges.at(at) = findj->first;
			weights.at(at) = findj->second;
			degrees.at(at) = edges.at(j).size();
		//	cout << " find neighbor " << findj->first << " cweight " << findj->second << " degree " << edges.at(j).size() << endl;
			at++;
		} else {
			assert(0);
		}
	}
//	cout << " remove node had weight " << iweight << " deg " << edges.at(index).size() << endl;
	for(size_t j=0; j<iedges.size(); ++j) {
		size_t oldjdegree = degrees.at(j);
		double oldjweight = weights.at(j)*oldjdegree;
//		cout << " update jnode " << iedges.at(j) << " oldweight " << oldjweight << endl;
		double newjweight = oldjweight - iweight;
//		cout << " newjweight " << newjweight;
		if(newjweight < 0) newjweight = 0;
		if(oldjdegree > 1) {
			newjweight /= (oldjdegree - 1);
		} else {
			newjweight = 0;
		}
//		cout << " newjweight " << newjweight << endl;
		remove_val(iedges.at(j), index_to_weight, weight_to_index);
		index_to_weight.insert(std::make_pair(iedges.at(j), newjweight));
		weight_to_index.insert(std::make_pair(newjweight, iedges.at(j)));
		// remove edge to j:
		edges.at(index).erase(iedges.at(j));
		edges.at(iedges.at(j)).erase(index);
	}

}

template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::compute_vcover_clarkson(overlap_type & o) {
	std::vector<std::set<size_t> > edges;
	overlap_graph(edges);


	// initialize the clarkson weights (wc(v)/dc(v) in terminology of the original paper)
	std::map<size_t, double> index_to_weight;
	std::multimap<double, size_t> weight_to_index;
	vector<double> orig_weights(sorted_original_als.size());
	double total_in_weight = 0;
# pragma omp parallel for num_threads(num_threads)
	for(size_t i=0; i<edges.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain1, gain2;
		common_model.gain_function(*al, gain1, gain2);
		double vgain = (gain1+gain2)/2 - base_cost;
	//	cout << " al " << i << " gain " << vgain <<  " at " << al << endl;
		assert(vgain >= 0.0);
#pragma omp critical(gain)
{
		total_in_weight += vgain;
		orig_weights.at(i) = vgain;
		double weight = vgain / edges.at(i).size();
		index_to_weight.insert(std::make_pair(i, weight));
		weight_to_index.insert(std::make_pair(weight, i));
}	
	}


	std::set<size_t> removed; // nodes removed = vertex cover 
	double orig_weight_removed = 0;

	while(!weight_to_index.empty()) {
		std::multimap<double, size_t>::iterator minit = weight_to_index.begin();
		size_t minind = minit->second;
		double clarkson_weight = minit->first;
//		if(clarkson_weight <= 0) {
//			break;
//		}

		if(edges.at(minind).size()>0) {
			removed.insert(minind);
			orig_weight_removed += orig_weights.at(minind);
		}
//		cout << " remove: " << minind << " now removed: " << removed.size() << endl;
//		cout << "orig weight " << orig_weights.at(minind) << " clarkson weight " << minit->first << endl;
		update_remove(minind, edges, index_to_weight, weight_to_index);
//		cout << " size " << weight_to_index.size() << endl;


		size_t numedges = 0;
		for(size_t i=0; i<edges.size(); ++i) {
			numedges+=edges.at(i).size();
		}
//		cout << " num edges now " << numedges << endl;

		if(numedges==0) {
			break;
		}

	
	}

	
	std::vector<const pw_alignment*> backup_soa(sorted_original_als);

	// add independent set
	size_t at = 0;
	double weight_in_independent_set = 0;
	std::set<const pw_alignment*> independent_als;
	for(size_t i=0; i<sorted_original_als.size(); i++) {
		std::set<size_t>::iterator findi = removed.find(i);
		if(findi == removed.end()) {
			sorted_original_als.at(at) = backup_soa.at(i);
			independent_als.insert(backup_soa.at(i));
			weight_in_independent_set += orig_weights.at(i);
			at++;
		//	cout << " keep " << i << " weight " << orig_weights.at(i) << " degree " << edges.at(i).size() << endl;
		}
	}
	size_t indep_set_size = at;
	double remainder_weight = 0;
	find_als_weight(independent_als, removed, backup_soa,at);
	// removed alignments after the others
/*	for(std::set<size_t>::iterator it = removed.begin(); it!=removed.end(); ++it) {
		//XXX Instead we added the new oredering based on the gain value of each alignment in removed set divided by number of independent als which have overlap with that al. Alignments with the bigger weight go first
	//	sorted_original_als.at(at) = backup_soa.at(*it);
		remainder_weight += orig_weights.at(*it);
	//	at++;
	}*/
	assert(at == sorted_original_als.size());

//	cout << " on " << backup_soa.size() << " alignments with total gain of " << total_in_weight << endl;
//	cout << " on " << backup_soa.size() << " alignments with total gain of " << max_gain << endl;
//	cout << " found an INDEPENDENT SET of " << indep_set_size << " alignments with total gain " << weight_in_independent_set <<endl;
//	cout << " remainder contains " << removed.size() << " alignments with total gain of " << remainder_weight << endl;

/*
	cout << endl << "New ordering: " << endl;
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		cout << " " <<i << " length " << al->alignment_length() << " ";
		if(i<indep_set_size) {
			cout << "+";
		} else {
			cout << "-";
		}
		cout << endl;
	
	}
*/	
//	std::cout << "compute lazy split : "<<std::endl;
	compute_simple_lazy_splits(o);



}

template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::find_als_weight(std::set<const pw_alignment*>& independent_set, std::set<size_t>& removed_set, std::vector<const pw_alignment*>& backup, size_t & at){
//	std::cout << "find weights of als "<<std::endl;
	alignment_index alind(data.numSequences());
	for(std::set<const pw_alignment*>::iterator it = independent_set.begin();it != independent_set.end();it++){
		const pw_alignment * p = *it;
		alind.insert(p);
	}
	std::multimap<double, size_t> sorter; // first is the weight, second is al index
	for(std::set<size_t>::iterator it = removed_set.begin();it!=removed_set.end();it++){
		const pw_alignment * al = backup.at(*it);
		std::vector<const pw_alignment *> result1;
		std::vector<const pw_alignment *> result2;
		alind.search_overlap(*al, 0, result1);
		alind.search_overlap(*al, 1, result2);
	//	assert(result1.size()!=0 || result2.size()!=0);
		double g1 ,g2;
		common_model.gain_function(*al,g1,g2);
		double av_gain = (g1+g2)/2;
		double weight = 0;
		if(result1.size()+result2.size()>0){
			weight = (av_gain/(result1.size()+result2.size()));
		}else{
			weight = (av_gain/1);
		}
		sorter.insert(std::make_pair(weight, *it));
	}
	for(std::multimap<double, size_t>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); rit++) {
		sorted_original_als.at(at) = backup.at(rit->second);
		at++;
	}
}
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::finding_repeats(std::set< const pw_alignment* , compare_pointer_pw_alignment> & cc,std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> >& overlapped, std::set< const pw_alignment* , compare_pointer_pw_alignment> & less_redundent){//XXX It can be written in a more efficient way, maybe i need to somehow integrate it with the old code
		alignment_index alind1(data.numSequences());
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc.begin(); it!=cc.end();it++){
			const pw_alignment * p = *it;
			alind1.insert(p);
		}
		std::set<const pw_alignment* , compare_pointer_pw_alignment> seen;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc.begin(); it!=cc.end();it++){	
			std::set<const pw_alignment*, compare_pointer_pw_alignment> highly_overlapped;
			const pw_alignment * p = *it;
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it1 = seen.find(p);
			if(it1 == seen.end()){
				highly_overlapped.insert(p);
				seen.insert(p);
				size_t left1,left2,right1,right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
				size_t length = p->alignment_length();
				std::vector<const pw_alignment*> results;
				alind1.search_overlap(*p, 0,results);
			//	find_overlap_rate(left1,right1,results,)
				for(size_t i =0; i < results.size();i++){
					const pw_alignment * al = results.at(i);
						size_t l1,l2,r1,r2;
						al->get_lr1(l1,r1);
						al->get_lr2(l2,r2);
						if(l1<right1 && right1< r1){
							double overlap = (right1-l1)/(double(length+al->alignment_length()));
							if(overlap>0.95){
								highly_overlapped.insert(al);//TODO maybe better check for 95% of this al here
								seen.insert(al);
							}
						}
						else if(left1<r1 && right1 > r1){
							double overlap = (left1-r1)/(double)(length+al->alignment_length());
							if((overlap*100)>0.95){
								highly_overlapped.insert(al);
								seen.insert(al);
							}


						}
						else if(l2<right1 && right1< r2){
							double overlap = (right1-l2)/(double)(length+al->alignment_length());
							if((overlap*100)>0.95){
								highly_overlapped.insert(al);
								seen.insert(al);
							}

						}
						else if(left1<r2 && right1 > r2){
							double overlap = (left1-r2)/(double)(length+al->alignment_length());
							if((overlap*100)>0.95){
								highly_overlapped.insert(al);
								seen.insert(al);
							}

						}	
				}
				if(highly_overlapped.size() != 1){
					overlapped.push_back(highly_overlapped);
				}
			}
		}
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc.begin(); it!=cc.end();it++){	
			const pw_alignment * p = *it;
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it1 = seen.find(p);
			if(it1 == seen.end()){
				less_redundent.insert(p);
			}
		}
	}

template<typename overlap_type>
void compute_cc_with_interval_tree<overlap_type>::add_the_interval(const pw_alignment* p){
	size_t ref1 = p->getreference1();
	size_t ref2 = p->getreference2();
	size_t l1,l2,r1,r2;
	p->get_lr1(l1,r1);
	p->get_lr2(l2,r2);

	intervals.at(ref1).push_back(Interval<const pw_alignment*>(l1,r1,p));
	intervals.at(ref2).push_back(Interval<const pw_alignment*>(l2,r2,p));
/*	AlignmentInterval itv1(p,0);
	std::cout<< &itv1<<std::endl;
	itv1.Print();
	IntervalTree & tree1 = trees.at(ref1);
	tree1.Insert(&itv1);
	AlignmentInterval itv2(p,1);
	itv2.Print();
	IntervalTree& tree2 = trees.at(ref2);
	tree2.Insert(&itv2);*/

}
template<typename overlap_type>
void compute_cc_with_interval_tree<overlap_type>::compute(std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> > & ccs){
//	IntervalTree::findOverlapping
	std::set <const pw_alignment*, compare_pointer_pw_alignment> seen;
	std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
	for(size_t i = 0; i < alignments.size();i++) {
	//	std::cout << "i "<< i << std::endl;
		const pw_alignment * al = alignments.at(i);
	//	al->print();
		std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) {
			std::set< const pw_alignment*, compare_pointer_pw_alignment> cc;
			get_cc(*al, cc, seen);
			sorter.insert(std::make_pair(cc.size(), cc));
			std::cout << "cc size is "<< cc.size()<<std::endl;
		}
	}
	for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
		ccs.push_back(it->second);
	}	std::cout << "ccs size is "<< ccs.size()<<std::endl;
}
template<typename overlap_type>
void compute_cc_with_interval_tree<overlap_type>::get_cc(const pw_alignment & al, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment*, compare_pointer_pw_alignment> & seen) {
	std::vector<size_t> left(2);
	std::vector<size_t> right(2);
	al.get_lr1(left.at(0), right.at(0));
	al.get_lr2(left.at(1), right.at(1));
	std::vector<size_t>reference(2);
	reference.at(0) = al.getreference1();
	reference.at(1) = al.getreference2();
//#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < 2;i++){
	//	std::cout << "on reference "<< i << std::endl;
	//	std::cout << reference.at(i) << " " << left.at(i)<< " "<< right.at(i)<<std::endl;
		cc_step(reference.at(i), left.at(i), right.at(i), cc, seen);	
	}	

}

template<typename overlap_type>
void compute_cc_with_interval_tree<overlap_type>::cc_step(size_t ref, size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment* , compare_pointer_pw_alignment>  & seen ) {
	std::set <const pw_alignment * , compare_pointer_pw_alignment> seen1;
	std::vector<Interval<const pw_alignment*> > ov;
	trees.at(ref).findOverlapping(left, right, ov);
//	std::cout << "ref " << ref << " left "<< left << " right "<<right <<std::endl;
	for(size_t i =0;i < ov.size();i++){
		const pw_alignment* al = ov.at(i).GetId();
	//	ov.at(i).PrintInterval();
		std::set <const pw_alignment* , compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()){ // if current al not contained in any connected component
			size_t aleft1, aright1,aleft2, aright2;	
			al->get_lr1(aleft1, aright1);//current alignment
			al->get_lr2(aleft2, aright2);
			seen.insert(al);
			seen1.insert(al);
			cc.insert(al);
		}else { }
	}
	//TODO remove seen1 objects
	std::vector<const pw_alignment* > seen2;
	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		const pw_alignment * al = *it;
		seen2.push_back(al);
	//	Interval<const pw_alignment> itv1(l1,r1,ref);
	}
//	std::cout << "seen2 size "<< seen2.size()<<std::endl;
//#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i< seen2.size(); i++){
		const pw_alignment * al;
		al = seen2.at(i);
		size_t aleft, aright;
		al->get_lr2(aleft, aright);
		cc_step(al->getreference2(), aleft, aright, cc, seen);
		al->get_lr1(aleft, aright);
		cc_step(al->getreference1(), aleft, aright, cc, seen);
	} 
}

template<typename overlap_type>
void compute_cc_avl<overlap_type>::compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs) {
//	alind.debug_print(); 
	std::cout << "compute CC on " << alignments.size() << std::endl;
	std::set <const pw_alignment*, compare_pointer_pw_alignment> seen;
	std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
	for(std::set<const pw_alignment *, compare_pointer_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
	//	std::cout << "compute_cc" <<std::endl;
		const pw_alignment * al = *it;
	//	al->print();
		std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
//		std::cout << " seen " << seen.size() << std::endl;
		if(seenal == seen.end()) {
		//	std::cout << " getcc" << std::endl;
			std::set< const pw_alignment*, compare_pointer_pw_alignment> cc;
			get_cc(*al, cc, seen);
		//	std::cout << "FOUND CC size " << cc.size() << std::endl;
		//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::const_iterator it = cc.begin(); it != cc.end();it++){
		//		const pw_alignment * test = *it;
		//		test->print();
		//	}
			sorter.insert(std::make_pair(cc.size(), cc));
		}	
	}
	for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
		ccs.push_back(it->second);
	}
}
template<typename overlap_type>
void compute_cc_avl<overlap_type>::get_cc(const pw_alignment & al , std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment*, compare_pointer_pw_alignment> & seen) {
	std::vector<size_t> left(2);
	std::vector<size_t> right(2);
	al.get_lr1(left.at(0), right.at(0));
	al.get_lr2(left.at(1), right.at(1));
	std::vector<size_t>reference(2);
	reference.at(0) = al.getreference1();
	reference.at(1) = al.getreference2();
#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < 2;i++){
	//	std::cout << "on reference "<< i << std::endl;
	//	std::cout << reference.at(i) << " " << left.at(i)<< " "<< right.at(i)<<std::endl;
		cc_step(reference.at(i), left.at(i), right.at(i), cc, seen);	
	}	
}

template<typename overlap_type>
void compute_cc_avl<overlap_type>::cc_step(size_t current , size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment* , compare_pointer_pw_alignment>  & seen ) {
//	std::cout<< "ref " << current << " l " << left << " r "<< right<<std::endl;
	vector<const pw_alignment *> seen1;
	std::multimap<size_t, std::pair<size_t, size_t> > touched_intervals;
	std::vector<const pw_alignment* > seen2;
#pragma omp critical(seen)
{
//	alind->super_search_overlap_and_remove(current, left, right, seen1, touched_intervals);
	alind->search_overlap(current, left, right, seen1);
//	std::cout << "seen1 "<<seen1.size() << " " << touched_intervals.size()<<std::endl;
	for(size_t i=0; i<seen1.size(); ++i) {
		const pw_alignment * seen_p = seen1.at(i);
		std::set <const pw_alignment* , compare_pointer_pw_alignment>::iterator find_seen = seen.find(seen_p);
		if(find_seen == seen.end()){ 
		//	seen.insert(seen1.at(i)); // TODO for now we keep seen, but it should not be necessary as all alignments are deleted from the index
			cc.insert(seen1.at(i));
			seen2.push_back(seen1.at(i));
	                //size_t l1, r1, l2, r2;
                	//seen_p->get_lr2(l2, r2);
                	//cc_step(seen_p->getreference2(), l2, r2, cc, seen);
                	//seen_p->get_lr1(l1, r1);
                	//cc_step(seen_p->getreference1(), l1, r1, cc, seen);
			seen.insert(seen1.at(i)); // TODO for now we keep seen, but it should not be necessary as all alignments are deleted from the index
		}
	}
}
       // std::cout << "seen2 size "<< seen2.size()<<std::endl;
        for(size_t i =0; i< seen2.size(); i++){
                const pw_alignment * al;
                al = seen2.at(i);
                size_t aleft, aright;
                al->get_lr2(aleft, aright);
                cc_step(al->getreference2(), aleft, aright, cc, seen);
                al->get_lr1(aleft, aright);
                cc_step(al->getreference1(), aleft, aright, cc, seen);
        }
//	for(std::multimap<size_t, std::pair<size_t, size_t> >::iterator it= touched_intervals.begin(); it!=touched_intervals.end(); ++it) {
//		size_t ref = it->first;
//		size_t l = it->second.first;
//		size_t r = it->second.second;
//		std::cout<< "ref " << ref << " l " << l << " r "<< r <<std::endl;
//	}

	/*for(std::multimap<size_t, std::pair<size_t, size_t> >::iterator it= touched_intervals.begin(); it!=touched_intervals.end(); ++it) {
		size_t ref = it->first;
		size_t l = it->second.first;
		size_t r = it->second.second;
	//	std::cout<< " touched ref " << ref << " l " << l << " r "<< r <<std::endl;

		cc_step(ref, l, r, cc, seen);
	}*/
}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::add_on_intmap(const pw_alignment * al){
	size_t ref1 = al->getreference1();
	size_t ref2 = al->getreference2();
	std::set<size_t> id;
	id.insert(alignments.size()-1);
//	std::cout<< "entering id "<< alignments.size()-1 <<std::endl;
	size_t l1,l2,r1,r2;
	al->get_lr1(l1,r1);
	al->get_lr2(l2,r2);
//	std::cout << l1 << " " << r1 << " "<< " " << l2 << " "<< r2<<std::endl;
	boost::icl::discrete_interval<size_t> bounds_int1 = boost::icl::construct<boost::icl::discrete_interval<size_t> >(l1,r1+1); //XXX Since right is excluded in boost intervals, i added +1.

	boost::icl::discrete_interval<size_t> bounds_int2 = boost::icl::construct<boost::icl::discrete_interval<size_t> >(l2,r2+1);

	als_on_reference.at(ref1).add(make_pair(bounds_int1,id));
//	all_intervals.at(ref1).insert(bounds_int1);
//	all_intervals.at(ref2).insert(bounds_int2);
	als_on_reference.at(ref2).add(make_pair(bounds_int2,id));
}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::remove_from_intmap(size_t & id, size_t & ref1, size_t & ref2, size_t& left1, size_t & right1, size_t & left2, size_t & right2){
	std::cout << "id " << id << " ref1 "<< ref1 << " left 1 "<< left1 << " right1 "<<right1 <<std::endl;
	for(boost::icl::interval_map<size_t, std::set<size_t> >::iterator it = als_on_reference.at(ref1).begin(); it != als_on_reference.at(ref1).end(); ++it){//TODO change it to the equal_range
//	for(boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator it = als_on_reference.at(ref1).find(left1); it != als_on_reference.at(ref1).end(); ++it){//TODO change it to the equal_range
		//Check all the alignments in this interval
		boost::icl::discrete_interval<size_t> itv  = (*it).first;
if(lower(itv)>=left1 && lower(itv) <= right1){
		std::set<size_t> overlaps_al = (*it).second;
	//	std::cout << "in interval " << itv << std::endl;
	//	std::cout<< "size of set "<< overlaps_al.size() << std::endl;
	//	if(boost::icl::lower(itv) > right1){
	//		std::cout << "break! " << lower(itv) << " > " << right1 <<std::endl;
	//		break;
	//	}
		std::set<size_t>::iterator setit = (*it).second.find(id);
		if(setit != (*it).second.end()){
			(*it).second.erase(setit);
		}
	//	TODO can not change the const set!
	//	(*it).second = overlaps_al;
}
if(boost::icl::lower(itv) > right1){
			std::cout << "break! " << lower(itv) << " > " << right1 <<std::endl;
			break;
}

	}
	std::cout << "id " << id << " ref2 " << ref2 << " left 2 "<< left2 << " right2 "<<right2 <<std::endl;
	for(boost::icl::interval_map<size_t, std::set<size_t> >::iterator it = als_on_reference.at(ref2).begin(); it != als_on_reference.at(ref2).end(); ++it){//TODO change it to the equal_range
//	for(boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator it = als_on_reference.at(ref2).find(left2); it != als_on_reference.at(ref2).end(); ++it){//TODO change it to the equal_range
		//Check all the alignments in this interval
		boost::icl::discrete_interval<size_t> itv  = (*it).first;
if(lower(itv)>=left2 && lower(itv) <= right2){

		std::set<size_t> overlaps_al = (*it).second;
	//	std::cout << "in interval " << itv << std::endl;
	//	std::cout<< "size of set "<< overlaps_al.size() << std::endl;
//		if(boost::icl::lower(itv) > right2){
//			std::cout << "break! " << lower(itv) << " > " << right2 <<std::endl;
//			break;
//		}
		std::set<size_t>::iterator setit = (*it).second.find(id);
		if(setit != (*it).second.end()){
			(*it).second.erase(setit);
}
	//	TODO can not change the const set!
	//	(*it).second = overlaps_al;
}
if(boost::icl::lower(itv) > right2){
			std::cout << "break! " << lower(itv) << " > " << right2 <<std::endl;
			break;
}

	}



}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::compute_test(std::vector<std::set<size_t> >& ccs){
	std::set<size_t> seen;
	for(std::map<size_t , std::pair<size_t,size_t> >::iterator it = id_and_bounds.begin(); it!= id_and_bounds.end();it++){
		size_t current = it->first;
		std::set<size_t>::iterator setit = seen.find(current);
		if(setit == seen.end()){
			std::set<size_t> cc;
			cc_step_current(current,it->second.first,it->second.second,cc,seen);
			ccs.push_back(cc);
		}
	}
	boost::icl::discrete_interval<size_t> itv = boost::icl::construct<boost::icl::discrete_interval<size_t> >(0,3);
//	boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator it =als_on_ref.find(itv);
	als_on_ref-= itv;
	std::cout << "finally "<< std::endl;
	for(boost::icl::interval_map<size_t, std::set<size_t> >::iterator it =als_on_ref.begin();it!=als_on_ref.end();it++){
		boost::icl::discrete_interval<size_t> itv  = (*it).first;
		std::cout << "from "<< first(itv)   << " to " << last(itv) <<std::endl;
		std::set<size_t> overlaps_al = (*it).second;
		std::cout << "in interval " << itv << std::endl;
		for(std::set<size_t>::iterator it1 = overlaps_al.begin(); it1 != overlaps_al.end(); it1++){
			std::cout << *it1 <<std::endl;
		}
	}

}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::cc_step_current(size_t & current, size_t & left, size_t & right, std::set<size_t>& cc, std::set<size_t> & seen){

	boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator searchbegin;
	boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator searchend;
	std::cout<<"from "<< left << " to "<< right <<std::endl;
	searchbegin = als_on_ref.find(left);
	searchend = als_on_ref.find(right);
	std::set<size_t> seen1;
	for(boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator it = searchbegin; it != searchend; ++it){
		//Check all the alignments in this interval
		boost::icl::discrete_interval<size_t> itv  = (*it).first;
		std::set<size_t> overlaps_al = (*it).second;
		std::cout << "in interval " << itv << std::endl;
		std::cout<< "size of set "<< overlaps_al.size() << std::endl;
		for(std::set<size_t>::iterator setit = overlaps_al.begin(); setit != overlaps_al.end();setit++){
			std::set<size_t>::iterator it_set = seen.find(*setit);
			if(it_set==seen.end()){
				seen.insert(*setit);
				seen1.insert(*setit);
				cc.insert(*setit);
			}
		}
	}
	std::vector<size_t> seen2;
	for(std::set<size_t>::iterator it = seen1.begin(); it != seen1.end(); it++){
		seen2.push_back(*it);
	}
	for(size_t i = 0; i < seen2.size();i++){
		std::map<size_t, std::pair<size_t,size_t> >::iterator boundries= id_and_bounds.find(seen2.at(i));
		cc_step_current(seen2.at(i), boundries->second.first,boundries->second.second,cc,seen);
	}


}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs) {//Filling in the ccs
	std::set <const pw_alignment*, compare_pointer_pw_alignment> seen;
	std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
	for(size_t i = 0; i < alignments.size();i++) {
	//	std::cout << "i "<< i << std::endl;
		const pw_alignment * al = alignments.at(i);
	//	al->print();
		std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) {
			std::set< const pw_alignment*, compare_pointer_pw_alignment> cc;
			get_cc(*al, cc, seen);
			sorter.insert(std::make_pair(cc.size(), cc));
			std::cout << "cc size is "<< cc.size()<<std::endl;
		}
	}
	for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
		ccs.push_back(it->second);
	}
	std::cout << "ccs size is "<< ccs.size()<<std::endl;

}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::get_cc(const pw_alignment & al, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment*, compare_pointer_pw_alignment> & seen) {
	std::vector<size_t> left(2);
	std::vector<size_t> right(2);
	al.get_lr1(left.at(0), right.at(0));
	al.get_lr2(left.at(1), right.at(1));
	std::vector<size_t>reference(2);
	reference.at(0) = al.getreference1();
	reference.at(1) = al.getreference2();
//#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i < 2;i++){
	//	std::cout << "on reference "<< i << std::endl;
	//	std::cout << reference.at(i) << " " << left.at(i)<< " "<< right.at(i)<<std::endl;
		cc_step(reference.at(i), left.at(i), right.at(i), cc, seen);	
	}	
}
template<typename overlap_type>
void compute_cc_with_icl<overlap_type>::cc_step(size_t ref, size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment* , compare_pointer_pw_alignment>  & seen ) {
	std::cout << "seen size "<<seen.size()<<std::endl;
	// search bounds (where could other alignments with overlap with the current one start or end)
	boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator searchbegin;
	searchbegin = als_on_reference.at(ref).find(left);
	std::cout<< "left " << left << " right " << right <<std::endl;
//	if(searchbegin == als_on_reference.at(ref).end()){
//		std::cout<< "reached end of the map! "<<std::endl;
//		std::pair<boost::icl::interval_set<size_t>::const_iterator, boost::icl::interval_set<size_t>::const_iterator> itRes = all_intervals.at(ref).equal_range(boost::icl::discrete_interval<size_t>::closed(left, right));
//		for(boost::icl::interval_set<size_t>::const_iterator it = itRes.first; it != itRes.second; ++it){
//			assert(it != all_intervals.at(ref).end());
//			boost::icl::discrete_interval<size_t> itv  = *it;
//			std::cout << "iterval is "<< itv << std::endl;
//			searchbegin = als_on_reference.at(ref).find(lower(itv));
//			std::cout << "lower " << lower(itv)<<std::endl;
//			if(searchbegin != als_on_reference.at(ref).end() && lower(itv) <=right && lower(itv)>=left){
//				break;
//			}
//			else{
//				searchbegin = als_on_reference.at(ref).end();
//				std::cout<< "skip the next for loop"<<std::endl;
//			}
//		}
//	}
	std::set <const pw_alignment * , compare_pointer_pw_alignment> seen1;


	// search for overlap first, then do all recursive calls after overlapping alignments have been put to seen 
	// this reduces the maximal recursion level
//	size_t numseen = 0;
//	for(boost::icl::interval_map<size_t, std::set<size_t> >::iterator it =  als_on_reference.at(ref).begin(); it != als_on_reference.at(ref).end(); ++it){
//		boost::icl::discrete_interval<size_t> itv  = (*it).first;
//		if(lower(itv) >= left && lower(itv) <= right){
//			std::cout << "in interval_test " << itv << std::endl;
//		}
//	}
	for(boost::icl::interval_map<size_t, std::set<size_t> >::const_iterator it = searchbegin; it != als_on_reference.at(ref).end(); ++it){//TODO change it to the equal_range
		//Check all the alignments in this interval
		boost::icl::discrete_interval<size_t> itv  = (*it).first;
		std::set<size_t> overlaps_al = (*it).second;
		std::cout << "in interval " << itv << std::endl;
	//	std::cout<< "lower is "<< boost::icl::lower(itv)<<std::endl;
		std::cout<< "size of set "<< overlaps_al.size() << std::endl;
		if(boost::icl::lower(itv) > right){
			std::cout << "break! " << lower(itv) << " > " << right <<std::endl;
			break;
		}
		std::set<size_t> intermediate;
		for(std::set<size_t>::iterator setit = overlaps_al.begin(); setit != overlaps_al.end();setit++){
			const pw_alignment* al = alignments.at(*setit);
		//	al->print();
//#pragma omp critical(seen)
//{
			std::set <const pw_alignment* , compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
			if(seenal == seen.end()){ // if current al not contained in any connected component
				size_t aleft1, aright1,aleft2, aright2;	
				al->get_lr1(aleft1, aright1);//current alignment
				al->get_lr2(aleft2, aright2);
				seen.insert(al);
				seen1.insert(al);
				cc.insert(al);
				intermediate.insert(*setit);
//(*it).second.erase(*setit);
			//	std::cout << "left "<< left << " right "<<right <<std::endl;
				assert((aright1 >= left && aleft1 <= right)||(aright2 >= left && aleft2 <= right));
			}
			 else {
			//	std::cout << "numseen "<<numseen << std::endl;
			//	numseen++;
			}
//}
		}
		std::cout<< "size of intermediate " << intermediate.size()<<std::endl;
//		for(std::set<size_t>::iterator it = intermediate.begin(); it!=intermediate.end();it++){
//			std::cout<< *it <<std::endl;
//			std::set<size_t>::iterator it1 = overlaps_al.find(*it);
//			const pw_alignment* al = alignments.at(*it1);
//			size_t aleft1, aright1,aleft2, aright2,ref1,ref2;	
//			al->get_lr1(aleft1, aright1);
//			al->get_lr2(aleft2, aright2);
//			ref1 = al->getreference1();
//			ref2 = al->getreference2();
//			size_t id = *it;
//			remove_from_intmap(id , ref1, ref2, aleft1, aright1, aleft2, aright2);
//		//	overlaps_al.erase(it1);//I should be able to remove this!
//		}
		std::cout<< "After removing "<<std::endl;
	}
//	boost::icl::discrete_interval<size_t> bounds = boost::icl::construct<boost::icl::discrete_interval<size_t> >(left,right);
//	als_on_reference.at(ref) -= bounds;
//		RIGHT = right;
//		REFERENCE = ref;
	std::vector<const pw_alignment* > seen2;
//#pragma omp parallel for num_threads(num_threads)
	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		const pw_alignment * al = *it;
		seen2.push_back(al);
	}
	std::cout << "seen2 size "<< seen2.size()<<std::endl;
//	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
//		const pw_alignment * al = *it;
	//	std::cout << "seen1 "<<std::endl;
//		al.print();
//#pragma omp parallel for num_threads(num_threads)
	for(size_t i =0; i< seen2.size(); i++){
		const pw_alignment * al;
		al = seen2.at(i);
		size_t aleft, aright;
		al->get_lr2(aleft, aright);
		cc_step(al->getreference2(), aleft, aright, cc, seen);
		al->get_lr1(aleft, aright);
		cc_step(al->getreference1(), aleft, aright, cc, seen);
	} 
/*	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		const pw_alignment & al = *it;
		std::cout << "seen2 "<<std::endl;
		al.print();
		size_t aleft, aright;	
		al.get_lr1(aleft, aright);
		cc_step(al.getreference1(), aleft, aright, cc, seen);
	} */	
//	boost::icl::discrete_interval<size_t> bounds = boost::icl::construct<boost::icl::discrete_interval<size_t> >(left,right);
//	als_on_reference.at(ref) -= bounds;
	
}
template<typename overlap_type>
compute_cc<overlap_type>::compute_cc(const all_data & dat): als_on_reference(dat.numSequences()), last_pos(dat.numSequences(), 0) {
	max_al_ref_length = 0;
	for(size_t i=0; i<dat.numAlignments(); ++i) {
		const pw_alignment & a = dat.getAlignment(i);
		alignments.insert(&a);
		add_on_mmaps(&a);
	}

}
template<typename overlap_type>
void compute_cc<overlap_type>::add_on_mmaps(const pw_alignment * pwa) {
	size_t ref1 = pwa->getreference1();
	size_t ref2 = pwa->getreference2();
/*	
	std::cout << "INS " <<std::endl;
	pwa->print();
	std::cout << std::endl;	
*/
	als_on_reference.at(ref1).insert(std::make_pair(pwa->getbegin1(), pwa));
	als_on_reference.at(ref1).insert(std::make_pair(pwa->getend1(), pwa));
	als_on_reference.at(ref2).insert(std::make_pair(pwa->getbegin2(), pwa));
	als_on_reference.at(ref2).insert(std::make_pair(pwa->getend2(), pwa));

	size_t left, right;
	pwa->get_lr1(left, right);
	size_t length = right - left + 1;
	if(length > max_al_ref_length) max_al_ref_length = length;
	if(right > last_pos.at(ref1)) last_pos.at(ref1) = right;
	pwa->get_lr2(left, right);
	length = right - left + 1;
	if(length > max_al_ref_length) max_al_ref_length = length;
	if(right > last_pos.at(ref2)) last_pos.at(ref2) = right;

}
template<typename overlap_type>
void compute_cc<overlap_type>::remove_on_mmaps(const pw_alignment * al) {
	size_t ref1 = al->getreference1();
	size_t left, right;
	al->get_lr1(left, right);

	std::pair<std::multimap<size_t, const pw_alignment*>::iterator, std::multimap<size_t, const pw_alignment*>::iterator > l1 = 
		als_on_reference.at(ref1).equal_range(left);
	for(std::multimap<size_t, const pw_alignment*>::iterator it = l1.first; it!=l1.second; ++it) {
		if(al == it->second) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}
	std::pair<std::multimap<size_t, const pw_alignment*>::iterator, std::multimap<size_t, const pw_alignment*>::iterator  > r1 = 
		als_on_reference.at(ref1).equal_range(left);//Changed them to right!
	for(std::multimap<size_t, const pw_alignment*>::iterator it = r1.first; it!=r1.second; ++it) {
		if(al == it->second) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}

	size_t ref2 = al->getreference2();
	al->get_lr2(left, right);
	std::pair<std::multimap<size_t, const pw_alignment*>::iterator, std::multimap<size_t, const pw_alignment*>::iterator  > l2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(std::multimap<size_t, const pw_alignment*>::iterator it = l2.first; it!=l2.second; ++it) {
		if(al==it->second) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}
	std::pair<std::multimap<size_t, const pw_alignment*>::iterator, std::multimap<size_t, const pw_alignment*>::iterator  > r2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(std::multimap<size_t, const pw_alignment*>::iterator it = r2.first; it!=r2.second; ++it) {
		if(al==it->second) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}

}
template<typename overlap_type>
void compute_cc<overlap_type>::compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs) {//Filling in the ccs
	std::set <const pw_alignment*, compare_pointer_pw_alignment> seen;
	std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
	for(std::set<const pw_alignment *, compare_pointer_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
	//	std::cout << "compute_cc" <<std::endl;
		const pw_alignment * al = *it;
	//	al->print();
		std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
		std::cout << " seen " << seen.size() << std::endl;
		if(seenal == seen.end()) {
		//	std::cout << " getcc" << std::endl;
			std::set< const pw_alignment*, compare_pointer_pw_alignment> cc;
			get_cc(*al, cc, seen);
			std::cout << "FOUND CC size " << cc.size() << std::endl;
			sorter.insert(std::make_pair(cc.size(), cc));
		}	
	}
	for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
		ccs.push_back(it->second);
	}



}
template<typename overlap_type>
void compute_cc<overlap_type>::get_cc(const pw_alignment & al, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment*, compare_pointer_pw_alignment> & seen) {
/*	size_t left, right; //XXX The commeted part is the original code. I change it in a way to do them in parallel
	al.get_lr1(left, right);
//	std::cout << "left1 "<< left << "right1 " << right <<std::endl;
	cc_step(al.getreference1(), left, right, cc, seen);
	al.get_lr2(left, right);
//	std::cout << "left2 "<< left << "right2 " << right <<std::endl;
	cc_step(al.getreference2(), left, right, cc, seen);*/
//XXX New code:
	std::vector<size_t> left(2);
	std::vector<size_t> right(2);
	al.get_lr1(left.at(0), right.at(0));
	al.get_lr1(left.at(1), right.at(1));
	std::vector<size_t>reference(2);
	reference.at(0) = al.getreference1();
	reference.at(1) = al.getreference2();
#pragma omp parallel for num_threads(numThreads)
	for(size_t i =0; i < 2;i++){
		std::cout<< "i" << i << std::endl;
		cc_step(reference.at(i), left.at(i), right.at(i), cc, seen);	
	}
	
}




// TODO further improvements to this function are possible if we store intervals on the references in which all alignments were already processed
template<typename overlap_type>
void compute_cc<overlap_type>::cc_step(size_t ref, size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment* , compare_pointer_pw_alignment>  & seen ) {
//	boost::icl::interval_map<size_t, size_t> intervals;

//	std::cout << "beginning of cc step "<<std::endl;
	// search bounds (where could other alignments which overlap with the current one start or end)
	// leftbound: all alignments starting at leftbound or ealier either have the end in the search interval or no overlap with the search interval
//	std::cout << " cc step " << ref << " from " << left << " to " << right << " seen is " << seen.size() << " we are on " << alignments.size() << " alignments" <<  std::endl; 
	std::multimap<size_t, const pw_alignment*>::iterator searchbegin;
	std::multimap<size_t, const pw_alignment*>::iterator searchend;
	if(right > max_al_ref_length) {
		size_t leftbound = right - max_al_ref_length;
		searchbegin = als_on_reference.at(ref).upper_bound(leftbound);
	} else {
		searchbegin = als_on_reference.at(ref).begin();
	}

	// rightbound: all alignments ending at rightbound or later either have their start in the search interval or no overlap with the search interval
	if(left + max_al_ref_length < last_pos.at(ref)) {
		size_t rightbound = left + max_al_ref_length;
		searchend = als_on_reference.at(ref).lower_bound(rightbound);
		searchend++; // make end exclusive

	} else {
		searchend = als_on_reference.at(ref).end();
	}

	std::set <const pw_alignment * , compare_pointer_pw_alignment> seen1;
//	std::set <pw_alignment, compare_pw_alignment> seen2;


	// search for overlap first, then do all recursive calls after overlapping alignments have been put to seen 
	// this reduces the maximal recursion level
	size_t numseen = 0;
	for(std::multimap<size_t, const pw_alignment* >::iterator it = searchbegin; it!=searchend; ++it) {
		const pw_alignment * al = it->second;
	//	std::cout<<"search"<<std::endl;
	//	al->print();
#pragma omp critical(seen)
{
		std::set <const pw_alignment* , compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) { // if current al not contained in any connected component
		//	std::cout << " not seen" << std::endl;
			size_t aleft, aright;	
	//		size_t leftmost_point_of_al_on_ref = numeric_limits<size_t>::max(); 
			if(al->getreference1()==ref) {
				al->get_lr1(aleft, aright);//current alignment
		//		std::cout << "alleft "<< aleft << "alright "<<aright <<std::endl;
				if(aright >= left && aleft <= right) {
				//	std::cout << "on ref 1 "<<std::endl;
					seen.insert(al);
					seen1.insert(al);
					cc.insert(al);
			//		std::cout << "ovlr " << cc.size() << " "  << seen.size() <<  " ref "<< ref << " : " << left << " " << right << " ovrlaps " << std::endl;
				//	al.print();
				}
		//		if(aleft < leftmost_point_of_al_on_ref ) {
		//			 leftmost_point_of_al_on_ref = aleft;
		//		}
			}
			if(al->getreference2()==ref) {
				al->get_lr2(aleft, aright);
		//		std::cout << "alleft "<< aleft << "alright "<<aright <<std::endl;
				if(aright >= left && aleft <= right) {
				//	std::cout<<"on ref 2 "<<std::endl;
			//		al->print();
					seen.insert(al);
					seen1.insert(al);
					std::set <const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it1 = cc.find(al);
				//	if(it1 != cc.end()){
				//		std::cout << "already there! "<<std::endl;
				//	}
					cc.insert(al);
				//	std::cout << "ovlr " << cc.size() << " "  << seen.size() << " ref "<< ref << " : " << left << " " << right << " ovrlaps " << std::endl;
				//	al.print();
				//	al->get_lr1(aleft, aright);
				//	cc_step(al->getreference1(), aleft, aright, cc, seen);
				}
		//		if(aleft < leftmost_point_of_al_on_ref ) {
		//			 leftmost_point_of_al_on_ref = aleft;
		//		}

			}
		//	if(leftmost_point_of_al_on_ref > rightbound) {
		//		break;
		//	}
		} else {
			numseen++;
			std::cout << numseen << std::endl;
		}
}
	}
//	std::cout << " found overlap with " << seen1.size() << " alignments, already seen before: " << numseen << std::endl;


	// now remove all seen alignments to be faster 
//	std::cout << "remove als"<<std::endl;
//	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
//		const pw_alignment * al = *it;
	//	al->print();
//		remove_on_mmaps(al);
//	}
/*	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		remove_on_mmaps(*it);
	}*/

	size_t debugsum = 0;
	for(size_t i=0; i<als_on_reference.size(); ++i) {
		debugsum+=als_on_reference.at(i).size();
	}
//	std::cout << " mmaps length " <<debugsum << std::endl;


	std::vector<const pw_alignment* > seen2;
//#pragma omp parallel for //TODO
	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		const pw_alignment * al = *it;
		seen2.push_back(al);
	}
	std::cout << "seen size " << seen2.size() <<std::endl;
//	for(std::set<const pw_alignment * , compare_pointer_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
//		const pw_alignment * al = *it;
	//	std::cout << "seen1 "<<std::endl;
//		al.print();
//#pragma omp parallel for num_threads(numThreads)
	for(size_t i =0; i< seen2.size(); i++){
		const pw_alignment * al;
		al = seen2.at(i);
		al->print();
		size_t aleft, aright;
	//	if(i%2==0) {	
		std::cout<< "its second ref "<<std::endl;
		al->get_lr2(aleft, aright);
		cc_step(al->getreference2(), aleft, aright, cc, seen);
	//	} else {
		std::cout<< "its first ref "<<std::endl;
		al->get_lr1(aleft, aright);
		cc_step(al->getreference1(), aleft, aright, cc, seen);
	//	}
	//	size_t tNum = omp_get_thread_num();
	//	std::cout << "num of threads "	 << tNum <<std::endl;
	} 
/*	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		const pw_alignment & al = *it;
		std::cout << "seen2 "<<std::endl;
		al.print();
		size_t aleft, aright;	
		al.get_lr1(aleft, aright);
		cc_step(al.getreference1(), aleft, aright, cc, seen);
	} */	
	
}
/*
template<typename tmodel>
clustering<tmodel>::clustering(overlap & o, all_data & d,tmodel & m):overl(o),data(d),model(m),als_on_ref(data.numSequences()),gain(data.numSequences(),std::vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),ava(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),res(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))){
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * a = &(data.getAlignment(i));
			alignments.insert(a);
		}

	}
template<typename tmodel>
clustering<tmodel>::~clustering(){}

template<typename tmodel>
void clustering<tmodel>::als_on_reference(const pw_alignment * p) {
	size_t ref1 = p->getreference1();
	size_t ref2 = p->getreference2();

	als_on_ref.at(ref1).insert(make_pair(p->getbegin1(), p));
	als_on_ref.at(ref1).insert(make_pair(p->getend1(), p));
	als_on_ref.at(ref2).insert(make_pair(p->getbegin2(), p));
	als_on_ref.at(ref2).insert(make_pair(p->getend2(), p));*

	}

//making a matrix that represents gain values:(First i did everything considering references on each alignemnt then i thought i may need to specify their position on the reference 'cus otherwise i couldnt clearly present the specific piece of alignemnt on the reference sequence. But then adding positions made things complicated and it is not working properly.I mean i had no idea how i can initialize the gain matrix and the way it works now it is so slow.)
template<typename tmodel>
void clustering<tmodel>::calculate_similarity(){
		double g1;
		double g2;
		for (size_t i = 0; i< data.numSequences();i++){
			als_on_ref.at(i) = overl.get_als_on_reference(i);
			for(std::multimap<size_t, pw_alignment*>::iterator it =als_on_ref.at(i).begin(); it!=als_on_ref.at(i).end(); ++it){
				const pw_alignment * p = it-> second;
				//p->print();
				size_t position = it->first;
				std::cout<<"position: "<<position<<std::endl;
				model.gain_function(*p,g1,g2);
				gain.at(p->getreference1()).at(p->getreference2()).at(position) += g1;//position has been used to specify the part of reference contains P's sample
				gain.at(p->getreference2()).at(p->getreference1()).at(position) += g2;
				std::cout<<"gain2: "<< gain.at(p->getreference2()).at(p->getreference1()).at(position)<<std::endl;
			}
	
		}
	}

template<typename tmodel>
	void clustering<tmodel>::update_values(){
	//	size_t iteration = 100;
		//for(size_t l=0; l<iteration;l++){
		for(size_t i= 0; i<data.numSequences();i++){
			double max = 0;
			double max1 = 0;	
			for(size_t j = 0; j< data.numSequences();j++){
				for(size_t m= 0; m <gain[i][j].size();m++ ){
					for(size_t n = 0; n < gain[i][j].size();n++ ){
						if(n!=m){
							double sum = ava.at(i).at(j).at(n) + gain.at(i).at(j).at(n);
							if(sum>max){
								max = sum;
							}	
						}else continue;
					}
				res.at(i).at(j).at(m) = gain.at(i).at(j).at(m)-max;
				for(size_t k= 0; k <gain[i][j].size();k++ ){
					if(k!=m ){
						if(res.at(i).at(j).at(k)> 0){
							max1 += res.at(i).at(j).at(k);
						}else{	
							max1 += 0;
						}
					}else continue;
				}
				if(i==j){
					ava.at(i).at(j).at(m) = max1;
				}else{
					double value = res.at(j).at(j).at(m) + max1;
					if(value < 0){
						ava.at(i).at(j).at(m) = value;
					}else {
						ava.at(i).at(j).at(m) = 0;
					}
				}
			}
			for(size_t o=0;o<ava[i][j].size();o++){
			std::cout<<"ava"<<	ava.at(i).at(j).at(o);
			}
		}
		}
//	}
}

template<typename tmodel>
	void clustering<tmodel>::update_clusters(size_t acc){
		size_t iteration = 100;
		double damp_value = 0.6;
		std::vector<size_t> examplar(data.numAcc(),0); //examplars of the class of acc, acc will be the center.
		for(size_t i = 0; i<data.numAcc(); i++){
			double value = res.at(acc).at(i) + ava.at(acc).at(i);
			for(size_t j = 0; j< data.numAcc(); j++){
				double max = res.at(j).at(i) + ava.at(j).at(i);
				if (max > value){
					break;
				}else {
				
				examplar.push_back(i);
				
				}

			}
		}
		
	}

template<typename tmodel>
	void clustering<tmodel>::update_clusters(){
		std::vector<vector<vector<double> > >examplar;
		std::vector<size_t> center;
		for(size_t i = 0; i<data.numSequences(); i++){
			for(size_t n = 0 ; n < gain[i][i].size();n++){
				examplar.at(i).at(i).at(n) = res.at(i).at(i).at(n) + ava.at(i).at(i).at(n);
				if (examplar.at(i).at(i).at(n)>0){
					center.push_back(i);
				}
			}
		//	std::cout<<"center at: "<< i << std::endl;

		}
		std::vector<size_t> idx(data.numSequences(),0);
		for(size_t i = 0; i<data.numSequences(); i++){
			double max = -HUGE_VAL;
			for(size_t j = 0; j<center.size(); j++){
				size_t c = center.at(j);
			for(size_t n = 0 ; n < gain[i][c].size();n++){
				if(gain.at(i).at(c).at(n)>max){
					max = gain.at(i).at(c).at(n);
					idx.at(i) = c;
				}
			}
}
		}
		for (size_t k = 0;k < data.numSequences();k++ ){
	//		std::cout << "center of "<< k << " is "<<idx.at(k)<<std::endl;
		}

	}

*/



template<typename tmodel, typename overlap_type>
void affpro_clusters<tmodel,overlap_type>::add_alignment(const pw_alignment & al) {
	// Get identifiers for both parts of the pairwise alignment
	std::stringstream sstr1;
//	std::cout<<"data1 ad in add_al: "<< & dat << std::endl;	
	size_t left1, right1;
	al.get_lr1(left1, right1);
//	if(al.getbegin1()<al.getend1()){
//		sstr1 << 0 <<":"<< al.getreference1()<<":"<<left1;
//	}else{
//		sstr1 << 1 <<":"<< al.getreference1()<<":"<<left1;		
//	}
	std::stringstream sstr2;
	size_t left2, right2;
	al.get_lr2(left2, right2);
//	if(al.getbegin2()<al.getend2()){
//		sstr2 << 0 <<":"<< al.getreference2()<<":"<<left2;
//	}else{
//		sstr2 << 1 <<":"<< al.getreference2()<<":"<<left2;		
//	}
	sstr1<< al.getreference1()<<":"<<left1;
	sstr2<< al.getreference2()<<":"<<left2;	
	std::string ref1name = sstr1.str();
	std::string ref2name = sstr2.str();

	// Get data indices for both alignment parts/ sequence pieces
	size_t ref1idx = 0;
	size_t ref2idx = 0;
	std::map<std::string, size_t>::iterator find1 = sequence_pieces.find(ref1name);
	std::map<std::string, size_t>::iterator find2 = sequence_pieces.find(ref2name);

	if(find1 == sequence_pieces.end()) {
		ref1idx = sequence_pieces.size();
		sequence_pieces.insert(std::make_pair(ref1name, ref1idx));
		sequence_names.push_back(ref1name);
		sequence_lengths.push_back(right1 - left1 + 1);
	} else {
		ref1idx = find1->second; 
	}
	if(find2 == sequence_pieces.end()) {
		ref2idx = sequence_pieces.size();
		sequence_pieces.insert(std::make_pair(ref2name, ref2idx));
		sequence_names.push_back(ref2name);
		sequence_lengths.push_back(right2 - left2 + 1);
	} else {
		ref2idx = find2->second;
	}
	
//	std::cout<<"data2 ad in add_al: "<< & dat << std::endl;	
//	std::cout << "ref1 id" << ref1idx <<" ref2 id "<< ref2idx<<std::endl;
	// enlarge similarity matrix
	size_t max = ref1idx;
	if(ref2idx > max) {
		max = ref2idx;
	}
	max++;
	if(max > simmatrix.size()) {
		for(size_t i=0; i<simmatrix.size(); ++i) {
			simmatrix.at(i).resize(max, -HUGE_VAL);
		}
		simmatrix.resize(max, std::vector<double>(max, -HUGE_VAL));
	}
	double c1;
	double c2;
	double m1; 
	double m2;
//	al->print();
//	std::cout<<"data3 ad in add_al: "<< & dat << std::endl;	
//	dat.numAcc();
//	std::cout << " dat adress " << & dat<< std::endl;
	model.cost_function(al, c1, c2, m1, m2);
	// preferences
	simmatrix.at(ref1idx).at(ref1idx) = -c1 - base_cost;
	simmatrix.at(ref2idx).at(ref2idx) = -c2 - base_cost;
	// similarity (i,j): i has exemplar j, this means modify j to i?
	simmatrix.at(ref2idx).at(ref1idx) = -m2;
	simmatrix.at(ref1idx).at(ref2idx) = -m1;
	double random_number;
	for(size_t i=0; i<simmatrix.size(); i++){//This for loops was added because we observed that AFP has issue in choosing the best center when there are several centers that return the same value.
//TODO replace it with mt 19937 random number generator!
		random_number = std::rand() % 10000;
		random_number = 1 / random_number;
		simmatrix.at(i).at(i) -= random_number;	
	}

}

	template<typename tmodel,typename overlap_type>
	size_t affpro_clusters<tmodel,overlap_type>::get_sequence_length(size_t ref_idx)const{	
		return 	sequence_lengths.at(ref_idx);
	}
	template<typename tmodel,typename overlap_type>
	void affpro_clusters<tmodel,overlap_type>::make_an_alignment(std::string & center , std::string & member , pw_alignment & p){
	//	std::map<std::string, std::vector<pw_alignment> >::iterator it = local_al_in_a_cluster.find(center);
	//	if(it == local_al_in_a_cluster){
	//		local_al_in_a_cluster.insert(make_pair(center, std::vector<pw_alignment>()));
	//		it = local_al_in_a_cluster.find(center);
	//	}
		std::vector<std::string> center_parts;
		strsep(center, ":" , center_parts);
	//	unsigned int dir = atoi(center_parts.at(0).c_str());
		unsigned int ref = atoi(center_parts.at(0).c_str());
		unsigned int left = atoi(center_parts.at(1).c_str());
		std::vector<std::string> member_parts;
		strsep(member, ":" , member_parts);
	//	unsigned int mem_dir = atoi(member_parts.at(0).c_str());
		unsigned int mem_ref = atoi(member_parts.at(0).c_str());
		unsigned int mem_left = atoi(member_parts.at(1).c_str());
		// look at all input alignments and determine which cluster it belongs to
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it1 = all_als.begin(); it1!=all_als.end(); ++it1){
			const pw_alignment * al = *it1;
		/*	unsigned int al_dir1;
			if(al->getbegin1() < al->getend1()){
				al_dir1 =0;
			}else{
				al_dir1 = 1;
			}
			unsigned int al_dir2;
			if(al->getbegin2() < al->getend2()){
				al_dir2 =0;
			}else{
				al_dir2 = 1;
			}*/
			size_t ref1 = al->getreference1();
			size_t ref2 = al->getreference2();
			size_t left1;
			size_t left2;
			size_t right1;
			size_t right2;
			al->get_lr1(left1,right1);
			al->get_lr2(left2,right2);
			if(ref1 == ref && left1 == left && ref2 == mem_ref && left2 == mem_left){
				p = *al;
			//	p.print();
			//	it->second.push_back(*al);
				double gain = model.get_the_gain(*al, center);
			//	assert(gain > 0);
				break;
			}
			if(ref2 == ref && left2 == left && ref1 == mem_ref && left1 == mem_left ){
				p = *al;
			//	it->second.push_back(*al);
				double gain = model.get_the_gain(*al, center);
			//	p.print();
			//	assert(gain > 0);
				break;
			}
		}
	}
	template<typename tmodel,typename overlap_type>
	void affpro_clusters<tmodel,overlap_type>::write_clusters(std::map<std::string, std::vector<std::string> > & cluster_result, std::map<std::string, std::vector<pw_alignment> > & local_al_in_a_cluster){
	//Go through alignments and members_of_clusters
		for(std::map<std::vector<std::string>,pw_alignment>::iterator it = alignments.begin(); it != alignments.end(); it++){
		//	std::cout << "cent is " << it->first.at(0) << std::endl;

			std::map<std::string, std::vector<pw_alignment> >::iterator al = local_al_in_a_cluster.find(it->first.at(0));
			if(al == local_al_in_a_cluster.end()){
				local_al_in_a_cluster.insert(std::make_pair(it->first.at(0), std::vector<pw_alignment>()));
				al= local_al_in_a_cluster.find(it->first.at(0));
			}
			al->second.push_back(it->second);
		}
		for(std::map<std::string, std::string>::iterator it= members_of_clusters.begin(); it != members_of_clusters.end(); it++){
		//	std::cout << "center is " << it->second << std::endl;
			std::map<std::string, std::vector<std::string> >::iterator cl = cluster_result.find(it->second);
			assert(cl != cluster_result.end());
		//	if(cl == cluster_result.end()){
		//		cluster_result.insert(std::make_pair(it->second,std::vector<std::string>()));
		//		cl=cluster_result.find(it->second);
		//	}
			cl->second.push_back(it->first);
		}
	}
	template<typename tmodel,typename overlap_type>
	void affpro_clusters<tmodel,overlap_type>::check_redundancy(){
		//check if a single piece of sequence with two different direction is assigned to two different clusters as members.
		std::vector<std::string> members;
		std::set<std::string> centers;
		for(std::map<std::string, std::string>::iterator it = members_of_clusters.begin() ; it!= members_of_clusters.end(); it++){
			members.push_back(it->first);
			centers.insert(it->second);
		}
		for(size_t i = 0; i < members.size()-1; i++){
			std::string current = members.at(i);
			for(size_t j = i+ 1; j < members.size();j++){
				if(members.at(j)==current){
					std::cout << "redundancy! " << current << " "<<members.at(j)<<std::endl;
					exit(1);
				}
				std::string rev = make_reverse(current);
				if(members.at(j)== rev){
					std::cout << "redundancy! " << current << " " << rev << std::endl;
					exit(1);
				}
				
			}
		}
		//A dirty solution to solve the problem of having the same piece both as a member and as a center
		for(std::set<std::string>::iterator it = centers.begin() ; it != centers.end() ;it++){
			std::string current = *it;
			for(size_t i =0; i < members.size();i++){
				if(members.at(i)==current){
					std::cout << "Warning!! A piece of sequence classified twice!" <<std::endl;
					exit(1);
				}
				std::string rev = make_reverse(current);
				if(members.at(i)==rev){//A piece of sequence happened both as a member and a center in two different direction
					//The one that is a member is removed!
					std::map<std::string, std::string>::iterator it1 = members_of_clusters.find(members.at(i));
					assert(it1 != members_of_clusters.end());
					std::vector<std::string> temp;
					temp.push_back(it1->second);
					temp.push_back(it1->first);
					std::map<std::vector<std::string>, pw_alignment>::iterator al = alignments.find(temp);
					assert(al != alignments.end());
					alignments.erase(al);
					members_of_clusters.erase(it1);
				}
			}
			
		}
		
		

	}
//This is the old suffix tree class, I am not using it anymore
/*
	suffix_tree::suffix_tree(all_data & d, finding_centers & c):data(d),centers(c), successive_centers(data.numSequences()), suffixes(data.numSequences()){
		size_t numberOfPowers = 32;
		powerOfTwo = std::vector<size_t>(numberOfPowers, 1);
		for(size_t i=1; i< numberOfPowers; ++i) {
			powerOfTwo.at(i) = powerOfTwo.at(i-1)*2;	
		}

	}
	suffix_tree::~suffix_tree(){
	
	}
	void suffix_tree::create_suffix(size_t seq_id){
//		std::cout<< "seq id " << seq_id << std::endl;
		std::vector<std::vector<std::vector<size_t> > > AllSuffixes;
		for(std::map<size_t , std::vector<size_t> >::iterator it=successive_centers.at(seq_id).begin(); it != successive_centers.at(seq_id).end(); it++){
			std::vector<std::vector<size_t> > ItsSuffix;
			for(size_t j =0; j < it->second.size(); j++){
				std::vector<size_t> suffix;
				for(size_t k =j; k < it->second.size();k++){
					suffix.push_back(it->second.at(k));
				}
				ItsSuffix.push_back(suffix);
			}
			AllSuffixes.push_back(ItsSuffix);
		}
		size_t counter = 0;
		for(std::map<size_t , std::vector<size_t> >::iterator it=successive_centers.at(seq_id).begin(); it != successive_centers.at(seq_id).end(); it++){
			size_t last_center = it->second.at(it->second.size()-1);
			for(size_t j =0; j < it->second.size()-1; j++){
				if(it->second.at(j)==last_center){
					for(size_t m =0; m < AllSuffixes.at(counter).size();m++){
						std::vector<size_t> new_suffix = AllSuffixes.at(counter).at(m);
						new_suffix.push_back(powerOfTwo.at(31));
						AllSuffixes.at(counter).at(m)= new_suffix;
					//	std::cout<< "new suffix "<<std::endl;
					//	for(size_t n =0; n < new_suffix.size();n++){
					//		std::cout<< new_suffix.at(n)<< " ";
					//	}
					//	std::cout << " " << std::endl;
					}
					std::vector<size_t> new_suffix;
					new_suffix.push_back(powerOfTwo.at(31));
					AllSuffixes.at(counter).push_back(new_suffix);
					break;
				}
			}			
		//	std::cout << "counter "<<counter <<std::endl;
			counter = counter + 1;
		}
		suffixes.at(seq_id)=(AllSuffixes);
//		if(successive_centers.at(seq_id).size() > 1){//Adding an extar char at the end of each suffix if the last center happens more than once!
//			size_t last_center = successive_centers.at(seq_id).at(successive_centers.at(seq_id).size()-1);
//			for(size_t i =0; i < successive_centers.at(seq_id).size()-1; i ++){
//				size_t current_center = successive_centers.at(seq_id).at(i);
//			//	std::cout<< " current " << current_center << " last center "<< last_center<<std::endl;
//				if(current_center == last_center){
//					for(size_t j =0;j < suffixes.at(seq_id).size();j++){
//						std::vector<size_t> new_suffix = suffixes.at(seq_id).at(j);
//						new_suffix.push_back(powerOfTwo.at(31));
//						suffixes.at(seq_id).at(j) = new_suffix;
//					}
//					std::vector<size_t> new_suffix;
//					new_suffix.push_back(powerOfTwo.at(31));
//					suffixes.at(seq_id).push_back(new_suffix);
//					break;
//				}
//			}
//		}
//		if(successive_centers.size()==1){
//			for(size_t j =0;j < suffixes.at(seq_id).size();j++){
//				std::vector<size_t> new_suffix = suffixes.at(seq_id).at(j);
//						new_suffix.push_back(powerOfTwo.at(31));
//						suffixes.at(seq_id).at(j) = new_suffix;
//			}
//			std::vector<size_t> new_suffix;
//			new_suffix.push_back(powerOfTwo.at(31));
//			suffixes.at(seq_id).push_back(new_suffix);
//		}
//		std::cout<< "seq id "<< seq_id <<std::endl;
//		if(suffixes.at(seq_id).size() > 0){
//			std::cout << " suffixes are " << std::endl;
//			for(size_t i =0; i < suffixes.at(seq_id).size();i++){
//				std::vector<std::size_t> suf = suffixes.at(seq_id).at(i);
//				for( size_t j =0; j < suf.size() ; j ++){
//					std::cout << suf.at(j)<< " ";
//				}
//					std::cout<< " " << std::endl;
//			}
//		}
//		std::cout<< "all the suffixes: "<<std::endl;
//		for(size_t i =0; i < suffixes.at(seq_id).size(); i++){
//			for(size_t j =0; j < suffixes.at(seq_id).at(i).size(); j++){
//				for(size_t  k =0; k < suffixes.at(seq_id).at(i).at(j).size();k++){
//					std::cout << suffixes.at(seq_id).at(i).at(j).at(k) << " ";
//				}
//				std::cout<< " " <<std::endl;
//			}
//		}
	}
	void suffix_tree::update_successive_centers(std::vector<size_t> & highest_path,size_t & index ,size_t & seq){
		for(std::map<size_t, std::vector<size_t> >::iterator it = successive_centers.at(seq).begin();it != successive_centers.at(seq).end();it++){
			std::vector<size_t> centers = it->second;
			for(size_t k =0; k < centers.size();k++){
				if(centers.at(k) == highest_path.at(0)&& ((centers.size()-k)>= highest_path.size())){
					size_t first_common_index;
					std::vector<size_t> common;
					for(size_t i = 0; i < highest_path.size();i++){
						first_common_index = k;
						if(centers.at(i+k)== highest_path.at(i)){
							common.push_back( highest_path.at(i));
						}else break;
					}
					if(common.size() == highest_path.size()){
						std::vector<size_t> new_center;
						for(size_t m = 0; m < first_common_index ;m++){
							new_center.push_back(centers.at(m));
						}
						new_center.push_back(index);
						for(size_t m = first_common_index + common.size(); m < centers.size();m++){
							new_center.push_back(centers.at(m));
						}
						centers=new_center;
					}
				}
			}
			it->second = centers;
		//	for(size_t i = 0; i < it->second.size(); i++){
		//		std::cout << it->second.at(i) << std::endl;
		//	}
		}
	}
	void suffix_tree::find_parent(size_t & node_index, size_t & parent){
		for(std::multimap<size_t,size_t>::iterator par =nodes_relation.begin();par!=nodes_relation.end();par++){
			size_t common_par = par->first;
			pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(common_par);
			for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				if(it1->second == node_index){
					parent = common_par;	
				//	std::cout<< "parent " << parent <<std::endl;
					break;
				}
			}
		}
	}
	void suffix_tree::delete_relation(size_t & parent, size_t & child_node){
		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(parent);
			for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				if(it1->second == child_node){
					nodes_relation.erase(it1);
				//	std::cout<<"child_node:"<< child_node<< " " << it1->second << "parent " << parent <<std::endl;
					break;
				}
			}
			//std::cout<<"nodes relation: "<<std::endl;
			//for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
			//	std::cout << it->first <<" "<< it->second << std::endl;
			//}
	}
	void suffix_tree::get_first_parent(std::vector<size_t> & current , std::vector<size_t> & first_parent){
		for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::vector<size_t> parent = it->first;
			if(current.at(0) == parent.at(0)){
				first_parent = it->first;
		//		std::cout<< "parent index: "<< it->second<<std::endl;
				break;
			}
		}
	//	std::cout<< "first parent in read function:" <<std::endl;
	///	for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
	//		std::vector<size_t> f_parent = it->first;
		//	std::cout << "node index is " << it->second << " ";
		//	for( size_t j =0; j < f_parent.size() ; j ++){
		//		std::cout << f_parent.at(j)<< " ";
		//	}
		//	std::cout<< " " << std::endl;
//		}
//		std::cout<< "all the nodes:"<<std::endl;
//		for(size_t i = 0; i < nodes.size(); i++){
//			std::string node = nodes.at(i);
//			std::cout << "node at " << i << " is ";
//			for(size_t j =0; j < node.size();j++){
//				std::cout << int(node.at(j)) << " " ;
//			}
//			std::cout << " " <<std::endl;
//		}
	}
	void suffix_tree::check_first_nodes_after_root(std::vector<size_t> & current, size_t & node_index){
		std::map<size_t,size_t>::iterator it = FirstCenterOnFirstNodes.find(current.at(0));
		if(it!=FirstCenterOnFirstNodes.end()){
			node_index = it->second;
		}
	}
	void suffix_tree::find_child_nodes(size_t & index, vector<size_t>& childs){
		childs.clear();
		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			childs.push_back(it1->second);
		//	std::cout << "index "<< index << "children: " << it1->second << std::endl;
		}
	}
	void suffix_tree::first_parent_index(size_t & current_index, size_t & first_parent){
			size_t index = 0;
			for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
				if(it->second <= current_index && index <= it->second){
					index = it->second;
					first_parent = it->second;
				}
			}		
	}
	void suffix_tree::create_tree(std::vector<size_t> & center_with_highest_gain, size_t & highest_index ){//make a tree dependent to high_gain. For the first time high_gain == 0 and then it is replaced by new one from merging. 
		nodes.clear();
		nodes_relation.clear();
		firstParent.clear();
		for(size_t seq_id =0; seq_id < data.numSequences(); seq_id++){
			if(center_with_highest_gain.size()==0){
				std::cout << "highest is zero! "<<std::endl;
				successive_centers.at(seq_id) = centers.get_center(seq_id);
				create_suffix(seq_id);
			}else{
			//	std::cout << "when highest is " << highest_index << std::endl;
				update_successive_centers(center_with_highest_gain,highest_index,seq_id);
				create_suffix(seq_id);
			}
			if(suffixes.at(seq_id).size() != 0){
				for(size_t i = 0; i < suffixes.at(seq_id).size(); i++){
				for(size_t q = 0; q < suffixes.at(seq_id).at(i).size();q++){
				//	std::cout << "suffixes at " << i  << " at " << q <<std::endl;
					std::vector<size_t> first_parent;
					get_first_parent(suffixes.at(seq_id).at(i).at(q),first_parent);
				//	size_t FirstParentNodeIndex = nodes.size();
				//	check_first_nodes_after_root(suffixes.at(seq_id).at(i).at(q),FirstParentNodeIndex);
				//	std::cout<< "f_parent "<<std::endl;
				//	for(size_t j =0; j < first_parent.size(); j++){
				//		std::cout << first_parent.at(j)<< " ";
				//	}
				//	std::cout << "" << std::endl;
					if(first_parent.size()==0){//if there is no first parent starts with the current suffix
				//	if(FirstParentNodeIndex==nodes.size()){
					//	std::cout << "if first parent is empty"<<std::endl;
						nodes.push_back(suffixes.at(seq_id).at(i).at(q));
						firstParent.insert(make_pair(suffixes.at(seq_id).at(i).at(q),nodes.size()-1));	
					//	FirstCenterOnFirstNodes.insert(make_pair(suffixes.at(seq_id).at(i).at(q).at(0),nodes.size()-1));//TODO should be updated!
					//	IndexOfFirstNodesAfterRoot.insert(nodes.size()-1);//TODO should be updated
					}else{
					//	std::vector<size_t> first_parent =  nodes.at(FirstParentNodeIndex);
					//	std::cout << "else "<<std::endl;// there is already a first parent starts with the current suffix
						std::vector<size_t> common_part;
						size_t shorter_length = suffixes.at(seq_id).at(i).at(q).size();
						if(first_parent.size()<= shorter_length){
							shorter_length = first_parent.size();
						}
						for(size_t j = 0; j < shorter_length ;  j++){
							if(first_parent.at(j)==suffixes.at(seq_id).at(i).at(q).at(j)){
								common_part.push_back(first_parent.at(j));
							}else break;
						}
					//	std::cout << "common part " << std::endl; 
					//	for(size_t j =0; j < common_part.size();j++){
					//		std::cout<< common_part.at(j)<< " ";
					//	}
					//	std::cout << " " << std::endl;
						std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(first_parent);
						assert(it!=firstParent.end());
						size_t node_index = it->second;//It is updated later on in a way that always is equal to the parent node index
						std::vector<size_t> current_parent = first_parent;//it will be updated
						std::vector<size_t> current_string = suffixes.at(seq_id).at(i).at(q);//it will be updated
						bool making_tree = true;
						while(making_tree == true){						
							//The first case:
							if(shorter_length == current_parent.size() && common_part.size() == shorter_length){//Current parent is shorter than the current suffix, and the current suffix contains all of it
								assert(nodes.at(node_index) == common_part);
							//	std::cout << "first case! "<<std::endl;
								std::vector<size_t> other;
								for(size_t j = common_part.size(); j < current_string.size(); j++){
									other.push_back(current_string.at(j));
								}
							//	std::cout << "other " <<std::endl;
							//	for(size_t j =0 ; j < other.size(); j++){
							//		std::cout << other.at(j) << " " ;
							//	}
							//	std::cout << " "<<std::endl;
								std::vector<size_t> childs;
								find_child_nodes(node_index,childs);
							//	std::cout << "child size "<<childs.size()<<std::endl;
								if(childs.size() == 0){
									std::cout<< "node index "<< node_index << std::endl;
							//		std::cout << "if it has no child node"<<std::endl;
									if(node_index == nodes.size()-1){//If it is the last node on the tree
										if(other.size() != 0){
											nodes.push_back(other);
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.push_back(new_suffix);
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
										}
										else{//when the entire current_string is on the parent node
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											if(current_parent.at(current_parent.size()-1) != powerOfTwo.at(31)){
												nodes.push_back(new_suffix);
												nodes_relation.insert(make_pair(node_index,node_index+1));
										//		std::cout << "# is pushed back!" <<std::endl;
											}
											else{ 
												std::cout<< "we are here ! 35" << std::endl;
												making_tree = false;
												break;
											}
										}
									}
									else if(node_index == nodes.size()-2){
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										std::vector<size_t> new_suffix;
										new_suffix.push_back(powerOfTwo.at(31));
										if(other.size() != 0){
											nodes.at(node_index +1) =other;
											nodes.push_back(new_suffix);
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
											std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(last_node); 
											std::map<size_t,size_t>::iterator it1 = FirstCenterOnFirstNodes.find(last_node.at(0));
											if(it != firstParent.end() && it->second == node_index+1){
												it->second = node_index+3;
												it1->second = node_index+3;
											//	std::set<size_t>::iterator it2 = IndexOfFirstNodesAfterRoot.find(node_index+1);
											//	assert(it2 != IndexOfFirstNodesAfterRoot.end());
											//	IndexOfFirstNodesAfterRoot.erase(it2);
											//	IndexOfFirstNodesAfterRoot.insert(node_index+3);
											}
										}else{
											if(current_parent.at(current_parent.size()-1) != powerOfTwo.at(31)){
												nodes.at(node_index +1) = new_suffix;
												nodes_relation.insert(make_pair(node_index,node_index+1));
												std::map<std::vector<size_t>, size_t>::iterator it_1 = firstParent.find(last_node);
												std::map<size_t,size_t>::iterator it1 = FirstCenterOnFirstNodes.find(last_node.at(0)); 
												if(it_1 != firstParent.end() && it_1->second == node_index+1){
												//	std::cout << "i m testing! "<<std::endl;
													it_1->second = node_index+2;
													it1->second = node_index+2;
												//	std::set<size_t>::iterator it2 = IndexOfFirstNodesAfterRoot.find(node_index+1);
												//	assert(it2 != IndexOfFirstNodesAfterRoot.end());
												//	IndexOfFirstNodesAfterRoot.erase(it2);
												//	IndexOfFirstNodesAfterRoot.insert(node_index+2);
												}
											}else{
												making_tree = false;
												break;
											}
										}
										nodes.push_back(last_node);
									}
									else{//Indices of all the nodes after that should be shifted 
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										std::vector<size_t> second_last_node = nodes.at(nodes.size()-2);
										if(other.size() != 0){
											for(size_t j =nodes.size()-1; j > node_index+1; j--){
												nodes.at(j) = nodes.at(j-2);
											}
											nodes.at(node_index +1) =other;
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.at(node_index+2) = new_suffix;
											size_t upper_bound = nodes.size()-1;
											size_t shifting_value = 2;
											shift_node_relation(node_index,  upper_bound , shifting_value);
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
										//	std::cout << " node relation for node  " << node_index <<std::endl;
										//	pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(node_index);
										//	for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
										//		std::cout << it1->second << " " ;
										//	}
										//	std::cout << " " <<std::endl;
											nodes.push_back(second_last_node);
											nodes.push_back(last_node);
											size_t shift_value = 2;
											shift_first_parent(node_index,shift_value,common_part);
										//	std::cout<<"nodes relation1: "<<std::endl;
										//	for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
										//		std::cout << it->first <<" "<< it->second << std::endl;
										//	}
										}else{//when the entire current string is completely on the parent node and parent has no child node 
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											if(current_parent.at(current_parent.size()-1) != powerOfTwo.at(31)){
												for(size_t j =nodes.size()-1; j > node_index+1; j--){
													nodes.at(j) = nodes.at(j-1);
												}
										//		std::cout << "here!"<<std::endl;
												nodes.at(node_index +1) =new_suffix;
												//The rest of node relation should be updated
												size_t upper_bound = nodes.size()-1;
												size_t shifting_value = 1;
												shift_node_relation(node_index,  upper_bound , shifting_value);
										//		for(size_t shift = nodes.size()-1; shift > node_index;shift--){
										//			pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
										//			std::vector<size_t> counter;
										//			for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
										//				if(it1 != nodes_relation.end()){
										//					counter.push_back(it1->second);
										//				}
										//			}
										//	//		std::multimap<size_t,size_t>::iterator p2 = nodes_relation.find(shift);
										//			if(counter.size() > 0){
										//				nodes_relation.erase(shift);
										//				for(size_t in = 0; in < counter.size(); in++){
										//					intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
										//				}
										//			}
										//		}//NEW CODES: checked on ecoli3, seems it works! 
										//		size_t ItsParent = nodes.size();
										//		size_t ItsfirstPar = nodes.size();
										//		find_parent(node_index,ItsParent);
										//		first_parent_index(node_index,ItsfirstPar);
										//		if(ItsParent != nodes.size()){
										//			while(ItsParent > ItsfirstPar){
										//				std::vector<size_t> children;
										//				find_child_nodes(ItsParent,children);
										//				for(size_t in = 0; in < children.size();in++){
										//					if(children.at(in) > node_index){
										//						delete_relation(ItsParent,children.at(in));
										//						intermediate.insert(make_pair(ItsParent,children.at(in)+1));	
										//						std::cout << "ItsParent " << ItsParent << " new children index " << children.at(in)+1 <<std::endl;
										//					}
										//				}
										//				size_t par = ItsParent;
										//				find_parent(par,ItsParent);
										//			}
										//			if(ItsParent == ItsfirstPar){
										//				std::vector<size_t> children;
										//				find_child_nodes(ItsParent,children);
										//				for(size_t in = 0; in < children.size();in++){
										//					if(children.at(in) > node_index){
										//						delete_relation(ItsParent,children.at(in));
										//						intermediate.insert(make_pair(ItsParent,children.at(in)+1));
										//						std::cout << "ItsParent " << ItsParent << " new children index " << children.at(in)+1 <<std::endl;	
										//					}
										//				}
										//			}
										//		}
										//		for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
										//			nodes_relation.insert(make_pair(it->first, it->second));
										//		}
												nodes_relation.insert(make_pair(node_index,node_index+1)); 
												nodes.push_back(last_node);
												size_t shift_value = 1;
												shift_first_parent(node_index, shift_value,common_part);
										//		for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){//size -1
										//			std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+1));
										//			if(it != firstParent.end() && it->second == shift){
										//				it->second = it->second+1;
										//			}
										//		}
										//		std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
										//		if(it != firstParent.end() && it->second == node_index){					
										//			firstParent.erase(it);
										//			firstParent.insert(make_pair(common_part,node_index));
										//		}	
												nodes.at(node_index) = common_part;// Seems like an extra work since we already knew the entire parent is covered
											//	std::cout << "# is pushed back!" <<std::endl;
											//	for(size_t i = 0; i < nodes.size(); i++){
											//		std::vector<size_t> node = nodes.at(i);
											//		std::cout << "node at " << i << " is ";
											//		for(size_t j =0; j < node.size();j++){
											//			std::cout << node.at(j) << " " ;
											//		}
											//		std::cout << " " <<std::endl;
											//	}

											}else{
												std::cout<< "we are here 35" << std::endl;
												making_tree = false;
												break;
											}
										}

									}
									making_tree = false;
								}else{//if it already has some children, here i need to check all the child nodes for the common context
									bool ItIsNotOnAnyOfChildren = false;
									size_t new_parent_index;
									for(size_t j = 0; j < childs.size();j++){
										if(other.size() != 0){
											if(nodes.at(childs.at(j)).at(0)== other.at(0)){
												ItIsNotOnAnyOfChildren = true;
												current_parent = nodes.at(childs.at(j));
												new_parent_index = childs.at(j);
												break;
											}
										}
									}
									if(ItIsNotOnAnyOfChildren == false){ //In this case current string is added as a new child node.
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										if(other.size() != 0){
											for(size_t j =nodes.size()-1; j > node_index+1; j--){
												nodes.at(j) = nodes.at(j-1);
											}
											nodes.at(node_index+1)=other;
										}else{
											bool AddNoMore = false;
											for(size_t j = 0; j < childs.size(); j++){
												std::vector<size_t> new_suffix;
												new_suffix.push_back(powerOfTwo.at(31));
												if(nodes.at(childs.at(j))== new_suffix){
													AddNoMore = true;
													break;
												}
											}
											if(AddNoMore == true){
												making_tree = false;
												break;
											}else{
												for(size_t j =nodes.size()-1; j > node_index+1; j--){
													nodes.at(j) = nodes.at(j-1);
												}
												std::vector<size_t> new_suffix;
												new_suffix.push_back(powerOfTwo.at(31));
												nodes.at(node_index+1)=new_suffix;
											}
										}
										nodes.push_back(last_node);
									//	std::multimap<size_t,size_t> intermediate;
										size_t shifting_value = 1;
										size_t upper_bound = nodes.size()-2;
										shift_node_relation(node_index,  upper_bound , shifting_value);
									//	for(size_t shift = nodes.size()-2; shift > node_index;shift--){
									//		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
									//		std::vector<size_t> counter;
									//		for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
									//			if(it1 != nodes_relation.end()){
									//				counter.push_back(it1->second);
									//			}
									//		}
									//		if(counter.size() > 0){
									//			nodes_relation.erase(shift);
									//			for(size_t in = 0; in < counter.size(); in++){
									//				std::cout <<"shift+1 " << shift + 1<<std::endl;
									//				intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
									//			}
									//		}
									//	}//NEW CODE:
									//	size_t ItsParent = nodes.size();
									//	size_t ItsfirstPar = nodes.size();
									//	find_parent(node_index,ItsParent);
									//	first_parent_index(node_index,ItsfirstPar);
									//	if(ItsParent != nodes.size()){
									//		while(ItsParent > ItsfirstPar){
									//			std::vector<size_t> children;
									//			find_child_nodes(ItsParent,children);
									//			for(size_t in = 0; in < children.size();in++){
									//				if(children.at(in) > node_index){
									//					delete_relation(ItsParent,children.at(in));
									//					intermediate.insert(make_pair(ItsParent,children.at(in)+1));	
									//					std::cout << "ItsParent1 " << ItsParent << " children index + 1 " << children.at(in)+1 <<std::endl;
									//				}
									//			}
									//			size_t par = ItsParent;
									//			find_parent(par,ItsParent);
									//		}
									//		if(ItsParent == ItsfirstPar){
									//			std::vector<size_t> children;
									//			find_child_nodes(ItsParent,children);
									//			for(size_t in = 0; in < children.size();in++){
									//				if(children.at(in) > node_index){
									//					delete_relation(ItsParent,children.at(in));
									//					intermediate.insert(make_pair(ItsParent,children.at(in)+1));
									//					std::cout << "ItsParent " << ItsParent << " children index + 1 " << children.at(in)+1 <<std::endl;	
									//				}
									//			}
									//		}
									//	}
										for(size_t j = 0 ; j < childs.size();j++){//Here we add new relation between node_index and its children.
											delete_relation(node_index,childs.at(j));
											nodes_relation.insert(make_pair(node_index, childs.at(j)+1));
										}
									//	for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
									//		nodes_relation.insert(make_pair(it->first, it->second));
									//	}
									//	std::cout<< "here!"<<std::endl;
										nodes_relation.insert(make_pair(node_index,node_index+1));
									//	std::cout <<"node index + 1 " << node_index + 1<<std::endl;
										size_t shift_value = 1;
										shift_first_parent(node_index, shift_value,common_part);
									//	for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){
									//		std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+1));		
									//		if(it != firstParent.end() && it->second == shift){
									//			std::cout << "shift for first parent: " << shift <<std::endl;
									//			it->second = it->second+1;
									//		}
									//	}
									//	std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
									//	if(it != firstParent.end() && it->second == node_index){
									//		firstParent.erase(it);
									//		firstParent.insert(make_pair(common_part,node_index));	
									//	}
										nodes.at(node_index)= common_part;
									//	std::cout<<"nodes relation2: "<<std::endl;
									//	for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
									//		std::cout << it->first <<" "<< it->second << std::endl;
									//	}
										making_tree = false;
									}else{//calculate the common_part (This is the only continious condition)
										std::cout<<"continuous one! "<<std::endl;
										nodes.at(node_index) = common_part;//This is just before changing the index to the new one
										std::cout << node_index << std::endl;
										node_index = new_parent_index;
										std::cout << node_index << std::endl;
										current_string = other;
										shorter_length = current_string.size();
										if(current_parent.size()<= shorter_length){
											shorter_length = current_parent.size();
										}
										common_part.clear();
										for(size_t j = 0; j < shorter_length ;  j++){
											if(current_parent.at(j)==current_string.at(j)){
												common_part.push_back(current_parent.at(j));
											}else break;
										}
									//	std::cout<<"common_part:"<<std::endl;
									//	for(size_t j =0; j < common_part.size(); j++){
									//		std::cout << common_part.at(j)<< " ";
									//	}
									//	std::cout << " " <<std::endl;
									}
								}
							}
							//The second case://The non common part is added as a child node and all old child nodes become children of the non-common part (The third and the fourth one was added here)
							else if((shorter_length == current_string.size() && common_part.size() == shorter_length) || common_part.size() < shorter_length ){
								std::vector<size_t> non_common;
								std::vector<size_t> non_common_2;
								for(size_t j =common_part.size(); j < current_parent.size(); j++){
									non_common.push_back(current_parent.at(j));
								}
								if(common_part.size() < shorter_length){
									for(size_t j =common_part.size(); j < current_string.size(); j++){
										non_common_2.push_back(current_string.at(j));
									}	
								}else { std::cout <<" there is no non common 2" << std::endl;}
							//	std::cout << "non_common "<<std::endl;
								assert(non_common.size() != 0);
							//	for( size_t j =0; j < non_common.size() ; j ++){
							//		std::cout << non_common.at(j)<< " ";
							//	}
							//	std::cout<< " " << std::endl;
							//	std::cout << "non_common_2"<<std::endl;
							//	for( size_t j =0; j < non_common_2.size() ; j ++){
							//		std::cout << non_common_2.at(j)<< " ";
							//	}
							//	std::cout<< " " << std::endl;
								vector<size_t> childs;
								find_child_nodes(node_index,childs);
								if(childs.size() == 0){
									std::cout << "node_index:: "<< node_index <<std::endl;
									if(node_index == nodes.size()-1){//If it is the last node on the tree
										if(common_part.size() < shorter_length){
											nodes.push_back(non_common_2);
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.push_back(new_suffix);
										}
										nodes.push_back(non_common);//Node relation is modified later on.
									}else if(node_index == nodes.size()-2){
									//	std::cout<<"node size - 2" <<std::endl;
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										if(common_part.size() < shorter_length){
											nodes.at(node_index+1)=non_common_2;
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.at(node_index+1)= new_suffix;
										}
										nodes.push_back(non_common);
										nodes.push_back(last_node);
										std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(last_node);
										std::map<size_t,size_t>::iterator it1 = FirstCenterOnFirstNodes.find(last_node.at(0));//XXX HERE! 
										if(it != firstParent.end() && it->second == node_index+1){
											firstParent.insert(make_pair(last_node,node_index+3));
											it1->second = node_index+3;
										//	std::set<size_t>::iterator it2 = IndexOfFirstNodesAfterRoot.find(node_index+1);
										//	assert(it2 != IndexOfFirstNodesAfterRoot.end());
										//	IndexOfFirstNodesAfterRoot.erase(it2);
										//	IndexOfFirstNodesAfterRoot.insert(node_index+3);
										}
									}else{//Indices of all the nodes after that should be shifted 
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										std::vector<size_t> second_last_node = nodes.at(nodes.size()-2);
									//	std::cout<< "node index "<< node_index<<std::endl;
										for(size_t j =nodes.size()-1; j > node_index+2; j--){
											nodes.at(j) = nodes.at(j-2);
										}
									//	nodes.at(node_index) = common_part;
										if(common_part.size() < shorter_length){
											nodes.at(node_index+1)=non_common_2;
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.at(node_index+1)= new_suffix;
										}
										nodes.at(node_index+2)= non_common;
										nodes.push_back(second_last_node);
										nodes.push_back(last_node);
										// Shifting the rest of relations as well:
									//	std::multimap<size_t,size_t> intermediate;
										size_t upper_bound = nodes.size()-3;
										size_t shifting_value = 2;
										shift_node_relation(node_index,  upper_bound , shifting_value);
									//	for(size_t shift = nodes.size()-3; shift > node_index;shift--){
									//		std::cout << "shift "<< shift << std::endl;
									//		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
									//		std::vector<size_t> counter;
									//		for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
									//			if(it1 != nodes_relation.end()){
									//				counter.push_back(it1->second);
									//			}
									//		}
									//		if(counter.size() > 0){
									//			nodes_relation.erase(shift);
									//			for(size_t in = 0; in < counter.size(); in++){
									//				intermediate.insert(make_pair(shift+2, counter.at(in)+2));					
									//			}
									//		}
									//	}
										//NEW_CODE: 
								//		size_t ItsParent = nodes.size();
								//		size_t ItsfirstPar = nodes.size();
								//		find_parent(node_index,ItsParent);
								//		first_parent_index(node_index,ItsfirstPar);
								//		if(ItsParent != nodes.size()){
								//			while(ItsParent > ItsfirstPar){
								//				std::vector<size_t> children;
								//				find_child_nodes(ItsParent,children);
								//				for(size_t in = 0; in < children.size();in++){
								//					if(children.at(in) > node_index){
								//						delete_relation(ItsParent,children.at(in));
								//						intermediate.insert(make_pair(ItsParent,children.at(in)+2));	
								//						std::cout << "ItsParent1 " << ItsParent << " new children index1 " << children.at(in)+2 <<std::endl;
								//					}
								//				}
								//				size_t par = ItsParent;
								//				find_parent(par,ItsParent);
								//			}
								//			if(ItsParent == ItsfirstPar){
								//				std::vector<size_t> children;
								//				find_child_nodes(ItsParent,children);
								//				for(size_t in = 0; in < children.size();in++){
								//					if(children.at(in) > node_index){
								//						delete_relation(ItsParent,children.at(in));
								//						intermediate.insert(make_pair(ItsParent,children.at(in)+2));
								//						std::cout << "ItsParent " << ItsParent << " new children index2 " << children.at(in)+2 <<std::endl;	
								//					}
								//				}
								//			}
								//		}
									//	for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
									//		nodes_relation.insert(make_pair(it->first, it->second));
									//	}
									//	std::cout << "test!"<<std::endl;
									}
							//		nodes.at(node_index) = common_part;
									nodes_relation.insert(make_pair(node_index,node_index+1));
									nodes_relation.insert(make_pair(node_index,node_index+2));
							}else{//It has children
								//This case is seperated into two different cases, if node_index = nodes.size -2 and else.
								std::vector<size_t> last_node = nodes.at(nodes.size()-1);
								std::vector<size_t> new_suffix;
								new_suffix.push_back(powerOfTwo.at(31));
								if(node_index == nodes.size()-2){//Note that since it has an at least one child node, it can not be the last node, also that last node can not be a first parent.
									if(common_part.size() < shorter_length){
										nodes.at(node_index+1)=non_common_2;
									}else{
										nodes.at(node_index+1)= new_suffix;
									}
									nodes.push_back(non_common);
									nodes.push_back(last_node);
									assert(childs.size()==1);
								}else{
									std::vector<size_t> second_last_node;
									second_last_node = nodes.at(nodes.size()-2);
									for(size_t j =nodes.size()-1; j > node_index+2; j--){
										nodes.at(j) = nodes.at(j-2);
									}
									if(common_part.size() < shorter_length){
										nodes.at(node_index+1)=non_common_2;
									//	std::cout<< "HEya!" <<std::endl;
									}else{
										nodes.at(node_index+1)= new_suffix;
									//	std::cout<< "HEYA!" <<std::endl;
									}
									nodes.at(node_index+2)= non_common;
									nodes.push_back(second_last_node);
									nodes.push_back(last_node);
									size_t upper_bound = nodes.size()-2;
									size_t shifting_value = 2;
									shift_node_relation(node_index,  upper_bound , shifting_value);
								}
								for(size_t j = 0 ; j < childs.size();j++){//Here we add the splitted part of the parent node and all the previous child nodes.
									delete_relation(node_index,childs.at(j));
									nodes_relation.insert(make_pair(node_index+2, childs.at(j)+2));
								}
							//	std::cout<<"nodes relation4: "<<std::endl;
							//	for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
							//		std::cout << it->first <<" "<< it->second << std::endl;
							//	}
							//	nodes.at(node_index) = common_part;
								nodes_relation.insert(make_pair(node_index,node_index+1));
								nodes_relation.insert(make_pair(node_index,node_index+2));
							}
							size_t shift_value = 2;
							shift_first_parent(node_index,shift_value, common_part);
						//	std::cout<<"nodes relation here is : "<<std::endl;
						//	for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
						//		std::cout << it->first <<" "<< it->second << std::endl;
						//	}
							nodes.at(node_index) = common_part;
							making_tree = false;					
						}
						//The third case
						//else if(length == current_parent.size() && common_part.size() < length){ current parent is broken(similar to case two)}
						//The forth case:
						//else if(length == current_string.size() && common_part.size() < length){ Similar to case two }
					}
				}
			}
			}
		}
		std::cout<< "all the nodes:"<<std::endl;
		for(size_t i = 0; i < nodes.size(); i++){
			std::vector<size_t> node = nodes.at(i);
			std::cout << "node at " << i << " is ";
			for(size_t j =0; j < node.size();j++){
				std::cout << node.at(j) << " " ;
			}
			std::cout << " " <<std::endl;
		}
		std::cout<<"nodes relation: "<<std::endl;
		for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
			std::cout << it->first <<" "<< it->second << std::endl;
		}
	}
		count_branches();
	}
	void suffix_tree::count_branches(){
		std::cout << "counting branches "<<std::endl;
		branch_counter.clear();
		for(size_t seq =0; seq < data.numSequences(); seq++){
			std::cout << "seq: "<<seq<< "suffix size " << suffixes.at(seq).size() << std::endl;
			for(size_t i = 0; i < suffixes.at(seq).size(); i++){
				std::cout<< "i: "<<i << std::endl;
				for(size_t q = 0; q < suffixes.at(seq).at(i).size();q++){
					std::vector<size_t> current = suffixes.at(seq).at(i).at(q);
					for(size_t j =0; j < current.size(); j++){
						std::cout << current.at(j) << " ";
					}
					std::cout << " " << std::endl;
					vector<size_t> branch;//insert this vector to the branch _counter map
					for(map<std::vector<size_t>,size_t>::iterator it = firstParent.begin();it!=firstParent.end();it++){
						std::vector<size_t> first_parent = it->first;
						if(first_parent.at(0)==current.at(0)){
							std::vector<size_t> updated_current;
							if(first_parent.size()==current.size()){
								std::cout<<"All the current is on the first parent!"<<std::endl;
							}else{//Note that first parent can not be longer than current.
								std::cout << "first updated current"<<std::endl;
								for(size_t j =first_parent.size() ; j < current.size(); j++){
									updated_current.push_back(current.at(j));
									std::cout << current.at(j)<<std::endl;
								}
							}
							size_t parent_index = it->second;
							branch.push_back(it->second);
							current = updated_current;
							while(current.size()>0){
								std::cout << "parent_index "<< parent_index <<std::endl;
								multimap<size_t,size_t>::iterator check = nodes_relation.find(parent_index);
								if(check != nodes_relation.end()){
									pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(parent_index);
									bool ItIsNotOnAnyChildNode = true;
									for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
										size_t CurrentParentIndex = it1->second;
										std::cout << "CurrentParentIndex "<< CurrentParentIndex <<std::endl;
										if(current.at(0)==nodes.at(CurrentParentIndex).at(0)){
											ItIsNotOnAnyChildNode = false;										
											branch.push_back(CurrentParentIndex);
											std::vector<size_t> updatedCurrent;
											if(nodes.at(CurrentParentIndex).size()==current.size()){
												std::cout<<"All the current is on the parent!"<<std::endl;
											}else{
												std::cout<< "next updated current: "<<std::endl;
												for(size_t j = nodes.at(CurrentParentIndex).size(); j<current.size(); j++){
													updatedCurrent.push_back(current.at(j));
													std::cout << current.at(j)<<std::endl;
												}
											}
											current = updatedCurrent;
											parent_index = CurrentParentIndex;
											break;
										}
									}
									if(ItIsNotOnAnyChildNode == true){
										std::cout<< "There is something wrong! "<<std::endl;
									}
								} else{ 
									std::cout<< "Parent has no kid! It shouldn't happen"<<std::endl;
									break;
								}
							}
							break;
						}
					}
					std::vector<size_t> sub_branch;
					cout<< "branch size: "<< branch.size() << endl;
					for(size_t i = 0 ; i < branch.size(); i++){
						sub_branch.push_back(branch.at(i));
						std::map<std::vector<size_t>,size_t>::iterator it1 = branch_counter.find(sub_branch);
						if(it1 == branch_counter.end()){
							branch_counter.insert(make_pair(sub_branch,0));
							it1 = branch_counter.find(sub_branch);
						}
					//	for(size_t j = 0; j < sub_branch.size();j++){
					//		std::cout << sub_branch.at(j)<< " ";
					//	}					
						it1->second = it1->second +1;
					//	std::cout<< "number of happening "<< it1->second << std::endl;
					}
				}			
			}

		}
		cout<< "branch counter: "<<endl;
	//	for(std::map<std::vector<size_t> , size_t>::iterator it = branch_counter.begin(); it != branch_counter.end(); it++){
	//		vector<size_t> br = it->first;
	//		for(size_t i =0; i < br.size();i++){
	//			cout<< br.at(i)<< " ";
	//		}
	//		cout<< " number "<<it->second;
	//		cout << " " <<endl;
	//		
	//	}
		std::cout << "counting the branches is done! "<<std::endl;
	}
	std::vector<std::vector<size_t> > suffix_tree::get_nodes()const{
		return nodes;
	}
	std::map<std::vector<size_t>, size_t> suffix_tree::get_count()const{
		return branch_counter;
	}
	std::vector<size_t> suffix_tree::get_first_parent()const{
		std::vector<size_t> first_parent;
		for(std::map<std::vector<size_t>,size_t>::const_iterator it = firstParent.begin(); it != firstParent.end(); it++){
			first_parent.push_back(it->second);
		}
		sort(first_parent.begin(),first_parent.end());	
	//	std::cout<<"first_parent: "<< std::endl;
	//	for(size_t i =0 ; i < first_parent.size(); i++){
	//		std::cout<< first_parent.at(i)<< " " ;
	//	}
	//	std::cout << " " <<std::endl;
		return first_parent;
	}
	size_t suffix_tree::get_power_of_two(size_t & power)const{
		size_t power_of_two = powerOfTwo.at(power);
		return power_of_two;
	}
	void suffix_tree::shift_node_relation( size_t & node_index, size_t & upper_bound , size_t & shift){
	//	std::cout << "node_index "<< node_index << " upper_bound "<< upper_bound << " shift "<< shift <<std::endl;
		std::multimap<size_t,size_t> intermediate;
		for(size_t i = upper_bound; i > node_index;i--){
			pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(i);
			std::vector<size_t> counter;
			for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
				if(it1 != nodes_relation.end()){
					counter.push_back(it1->second);
				}
			}
			if(counter.size() > 0){
				nodes_relation.erase(i);
				for(size_t in = 0; in < counter.size(); in++){
			//		std::cout <<"i + shift " << i + shift <<std::endl;
					intermediate.insert(make_pair(i + shift, counter.at(in)+shift));					
				}
			}
		}
		size_t ItsParent = nodes.size();
		size_t ItsfirstPar = nodes.size();
		find_parent(node_index,ItsParent);
		first_parent_index(node_index,ItsfirstPar);
	//	std::cout <<"node index "<< node_index<< "ItsParent " << ItsParent << " ItsFirstPar "<< ItsfirstPar << " node size "<< nodes.size() << std::endl;
		if(ItsParent != nodes.size()){
			while(ItsParent > ItsfirstPar){
				std::vector<size_t> children;
				find_child_nodes(ItsParent,children);
			//	std::cout << "children size "<< children.size() <<std::endl;
				for(size_t in = 0; in < children.size();in++){
					if(children.at(in) > node_index){
						delete_relation(ItsParent,children.at(in));
						intermediate.insert(make_pair(ItsParent,children.at(in)+shift));	
					//	std::cout << "ItsParent1 " << ItsParent << " children index + shift " << children.at(in)+ shift <<std::endl;
					}
				}
				size_t par = ItsParent;
				find_parent(par,ItsParent);
			}
			if(ItsParent == ItsfirstPar){
				std::vector<size_t> children;
				find_child_nodes(ItsParent,children);
				for(size_t in = 0; in < children.size();in++){
					if(children.at(in) > node_index){
						delete_relation(ItsParent,children.at(in));
						intermediate.insert(make_pair(ItsParent,children.at(in)+shift));
					//	std::cout << "ItsParent " << ItsParent << " children index + shift " << children.at(in)+shift <<std::endl;	
					}
				}
			}
		}
		for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
			nodes_relation.insert(make_pair(it->first, it->second));
		}
	}
	void suffix_tree::shift_first_parent(size_t & node_index , size_t & shifting_value, std::vector<size_t> & common_part){
		for(size_t i = node_index+1; i <nodes.size()-shifting_value; i++){
			std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(i+shifting_value));
			std::map<size_t,size_t>::iterator it1 = FirstCenterOnFirstNodes.find(nodes.at(i+shifting_value).at(0));
		//	std::cout<< "nodes.at(i+shifting_value).at(0) "<< nodes.at(i+shifting_value).at(0) <<std::endl;
			if(it != firstParent.end() && it->second == i){
		//		std::cout << it->second << " " << i << " " << shifting_value <<std::endl;

			//	assert(it1 != FirstCenterOnFirstNodes.end() );
			//	for(std::set<size_t>::iterator it2 = IndexOfFirstNodesAfterRoot.begin(); it2 != IndexOfFirstNodesAfterRoot.end(); it2++){//TODO
			//		std::cout << *it2 <<std::endl;
			//	}
			//	std::set<size_t>::iterator it2 = IndexOfFirstNodesAfterRoot.find(it->second);
			//	assert(it2 != IndexOfFirstNodesAfterRoot.end());
			//	IndexOfFirstNodesAfterRoot.erase(it2);
			//	IndexOfFirstNodesAfterRoot.insert(it->second+shifting_value);
				it->second = it->second+shifting_value;
				it1->second = it1->second+shifting_value;
			 
			}
		}
	//	std::cout<< "first parent shift function: " << " node index " << node_index << " its content size " << nodes.at(node_index).size()<<std::endl;
		std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
		if(it != firstParent.end() && it->second == node_index){
			firstParent.erase(it);
			firstParent.insert(make_pair(common_part,node_index));
		}
	}
	std::map<size_t, std::vector<size_t> > suffix_tree::get_center_on_a_sequence(size_t & seq_id)const{
		return successive_centers.at(seq_id);
	}
	std::vector<std::vector<std::vector<size_t> > > suffix_tree::get_suffix_of_sequence(size_t & seq_id)const{
		return suffixes.at(seq_id);
	}
*/

#include "overlap.cpp"

#endif
