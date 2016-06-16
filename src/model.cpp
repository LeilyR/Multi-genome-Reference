#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP
#define PRINT 0
template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::compute(overlap_type & o) {

//	compute_simple_lazy_splits(o);
//	compute_simple(o);
	compute_vcover_clarkson(o);

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
	std::cout<< "al in lazy split recurssion level " << level << std::endl;
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
		al->print();
		cout << endl;
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
// TODO
			o.test_partial_overlap();
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

		al->get_lr1(left1, right1);//TODO Parallel
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
	pair<std::multimap<double, size_t>::iterator, std::multimap<double, size_t>::iterator > eqr = weight_to_index.equal_range(weight);
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


	set<size_t> removed; // nodes removed = vertex cover 
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
	compute_simple_lazy_splits(o);



}

template<typename T, typename overlap_type>
void initial_alignment_set<T,overlap_type>::find_als_weight(std::set<const pw_alignment*>& independent_set, std::set<size_t>& removed_set, std::vector<const pw_alignment*>& backup, size_t & at){
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
		sorter.insert(make_pair(weight, *it));
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
							size_t overlap = (left1-r1)/(length+al->alignment_length());
							if((overlap*100)>0.95){
								highly_overlapped.insert(al);
								seen.insert(al);
							}


						}
						else if(l2<right1 && right1< r2){
							size_t overlap = (right1-l2)/(length+al->alignment_length());
							if((overlap*100)>0.95){
								highly_overlapped.insert(al);
								seen.insert(al);
							}

						}
						else if(left1<r2 && right1 > r2){
							size_t overlap = (left1-r2)/(length+al->alignment_length());
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
#pragma omp critical(seen)
{
	alind.super_search_overlap_and_remove(current, left, right, seen1, touched_intervals);
//	std::cout << "seen1 "<<seen1.size() << " " << touched_intervals.size()<<std::endl;
	for(size_t i=0; i<seen1.size(); ++i) {
		seen.insert(seen1.at(i)); // TODO for now we keep seen, but it should not be necessary as all alignments are deleted from the index
		cc.insert(seen1.at(i));
	}
}
//	for(std::multimap<size_t, std::pair<size_t, size_t> >::iterator it= touched_intervals.begin(); it!=touched_intervals.end(); ++it) {
//		size_t ref = it->first;
//		size_t l = it->second.first;
//		size_t r = it->second.second;
//		std::cout<< "ref " << ref << " l " << l << " r "<< r <<std::endl;
//	}

	for(std::multimap<size_t, std::pair<size_t, size_t> >::iterator it= touched_intervals.begin(); it!=touched_intervals.end(); ++it) {
		size_t ref = it->first;
		size_t l = it->second.first;
		size_t r = it->second.second;
	//	std::cout<< " touched ref " << ref << " l " << l << " r "<< r <<std::endl;

		cc_step(ref, l, r, cc, seen);
	}
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
	if(al.getbegin1()<al.getend1()){
		sstr1 << 0 <<":"<< al.getreference1()<<":"<<left1;
	}else{
		sstr1 << 1 <<":"<< al.getreference1()<<":"<<left1;		
	}
	std::stringstream sstr2;
	size_t left2, right2;
	al.get_lr2(left2, right2);
	if(al.getbegin2()<al.getend2()){
		sstr2 << 0 <<":"<< al.getreference2()<<":"<<left2;
	}else{
		sstr2 << 1 <<":"<< al.getreference2()<<":"<<left2;		
	}
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
//Finds all the centers the happen on a sequence:
	finding_centers::finding_centers(all_data & d):data(d),AlignmentsFromClustering(data.numSequences()),centersOnSequence(data.numSequences()), centersOfASequence(data.numSequences()),initial_suffixes(data.numSequences()){
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				if(it->second.size() != 0){
					string center = it->first;
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int center_dir = atoi(center_parts.at(0).c_str());
					unsigned int center_ref = atoi(center_parts.at(1).c_str());
					unsigned int center_left = atoi(center_parts.at(2).c_str());
					for(size_t j = 0; j < it->second.size();j++){
						pw_alignment * p = & it->second.at(j);
						size_t left1; 
						size_t left2;
						size_t right1;
						size_t right2;
						p->get_lr1(left1,right1);
						p->get_lr2(left2,right2);
						if(i != center_ref){
							if(p->getreference1()== i && p->getreference2()== center_ref &&left2== center_left){
								std::map<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
								if(it1 == AlignmentsFromClustering.at(i).end()){
									AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
								}
							} 
							if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
								std::map<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
								if(it2 == AlignmentsFromClustering.at(i).end()){
									AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
								}
							}
						}
						if(i == center_ref){
							if(p->getreference1()== i && p->getreference2()== i){
								std::map<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(left1);				
								if(it == AlignmentsFromClustering.at(i).end()){
									AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
								}
								std::map<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left2);				
								if(it1 == AlignmentsFromClustering.at(i).end()){
									AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
								}
							}else{
								std::map<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(center_left);
								if(it == AlignmentsFromClustering.at(i).end()){
									AlignmentsFromClustering.at(i).insert(make_pair(center_left,p));
								}else continue;
							}
						}
					}
				}
			}
		/*	std::cout << "on sequence " << i << ":" << std::endl;
			for(std::map<std::string, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).begin(); it1 != AlignmentsFromClustering.at(i).end(); it1++){
				pw_alignment * p = it1->second;
			//	p->print();
				std::cout << "position is  "<< it1->first << std::endl;
			}*/
		}
	}
	void finding_centers::findMemberOfClusters(map<string,vector<pw_alignment> > & alignmentsOfClusters){//fills in a map with centers and their members and removes the direction from members
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
			assert(it->second.size() != 0);
			if(it->second.size() != 0){
				std::cout << "no of als " << it->second.size()  << " " << it->first << std::endl;
				for(size_t i = 0; i < it->second.size(); i++){
					size_t ref1;
					size_t ref2;
					size_t left1;
					size_t left2;
					size_t right1;
					size_t right2;
					pw_alignment p =it->second.at(i);
					p.print();
					p.get_lr1(left1,right1);
					p.get_lr2(left2,right2);
					ref1 = p.getreference1();
					ref2 = p.getreference2();
					unsigned int al_dir1;
					unsigned int al_dir2;
					if(p.getbegin1()<p.getend1()){
						al_dir1 = 0;
					}else al_dir1 = 1;
					if(p.getbegin2()<p.getend2()){
						al_dir2 = 0;
					}else al_dir2 = 1;
					stringstream sample1;
					stringstream sample2;
					sample1<< al_dir1<< ":" << ref1 << ":" << left1;
					sample2<< al_dir2 << ":" << ref2 << ":" << left2;
					stringstream member1;
					stringstream member2;
					member1 << ref1 << ":"<<left1;
					member2 << ref2 << ":" << left2;
					if(sample1.str()==it->first){
						memberOfCluster.insert(make_pair(member2.str(),it->first));
					}else{
						memberOfCluster.insert(make_pair(member1.str(),it->first));
					}
					std::vector<std::string> center_parts;
					strsep(it->first, ":" , center_parts);
					unsigned int ref = atoi(center_parts.at(1).c_str());
					unsigned int left = atoi(center_parts.at(2).c_str());
					stringstream cent;
					cent << ref <<":"<<left;
					memberOfCluster.insert(make_pair(cent.str(),it->first));
				}
			}
		}
	//	std::cout << "memberOfCluster size "<< memberOfCluster.size() << std::endl;
	}
	void finding_centers::center_frequency(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::vector<std::map<size_t, std::string> > & centerOnseq){//it basically returns indices of centers on each sequence
		setOfAlignments(alignmentsOfClusters);//set all the alignments on each sequence //TODO change the potential long center map (centersOfASequence) remove those which have less than allowed gap but dont go through the same direction!
		findMemberOfClusters(alignmentsOfClusters);// returns strings that are member of a cluster in memberOfCluster
		//First: fill in the "center_index" vector
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it2=alignmentsOfClusters.begin(); it2 != alignmentsOfClusters.end();it2++){ //This is 'al_of_a_ccs' map in the main function. It contains centers with both directions
			std::string center = it2->first;
			center_index.push_back(center);
		}
		cout<< "center index size" << center_index.size()<<endl;
		for(size_t i = 0; i< data.numSequences(); i++){
			cout << "sequence: " << i << endl;
			const dnastring & sequence = data.getSequence(i);
			// Second: Find all the centers on each sequence
			for(size_t n= 0; n < sequence.length(); n++){
				string center;
				std::map<size_t, pw_alignment*>::iterator it=AlignmentsFromClustering.at(i).find(n);
				if(it != AlignmentsFromClustering.at(i).end()){
					stringstream member;
					member<< i << ":" << n;
					std::cout << "mem" << member.str() <<std::endl;
					std::map<std::string,std::string>::iterator it1 = memberOfCluster.find(member.str());
					if(it1 != memberOfCluster.end()){
						center = it1->second;					
						cout<< "center: "<< center<<endl;
						size_t cent_index = 0;
						for(size_t j =0; j < center_index.size(); j ++){
							if(center_index.at(j)==center){
								cent_index = j;
								std::cout<< " n "<< n << " cent_index " << cent_index<<endl;
								centersOnSequence.at(i).insert(make_pair(n,center));
								centerOnseq.at(i).insert(make_pair(n,center));//This is the global one. It is initialized in main function
								break;
							}
						}
					/*	int cent_index = 0;
						std::map<std::string , int>::iterator it2 = oriented_index.find(center);
						assert(it2 != oriented_index.end());
						cent_index = it2->second;
						cout<< "cent_index " << cent_index<<endl;
						centersOnSequence.at(i).insert(make_pair(n,center));
						centerOnseq.at(i).insert(make_pair(n,center));//This is the global one	*/						
					}
				}
			}
			//Third: Find potential candidates for being a long center. 
			//We are going to find those centers that have less than ALLOWED_GAP base pairs  distances. Their indices are saved in a vector.
			size_t last_position;
			size_t first_position;
			std::map<size_t, std::vector<size_t> > AllConnectedOnes;
			vector<size_t> vectorOfcenters;
			std::vector<size_t> position_of_each_suffix; //First one is position, second of is center
			bool direction = true;
			size_t dir;
			size_t dir1;
			for(std::map<size_t, std::string>::iterator it = centersOnSequence.at(i).begin(); it != centersOnSequence.at(i).end();it++){
				std::map<size_t, pw_alignment*>::iterator it1=AlignmentsFromClustering.at(i).find(it->first);
				pw_alignment * al = it1->second;
				size_t l1,l2,r1,r2;
				al->get_lr1(l1,r1);
				al->get_lr2(l2,r2);
				if(l1 == it->first && al->getreference1()== i){
					if(al->getbegin1()<al->getend1()){
						dir1 = 1;
					}else{
						dir1 = 0;
					}
				}
				if(l2 == it->first && al->getreference2() == i){
					if(al->getbegin2()< al->getend2()){
						dir1 = 1;
					}else{
						dir1 = 0;
					}
				}
				if(vectorOfcenters.size()==0){
					last_position = it->first;
					first_position = it->first;
					dir = dir1;
				}
				if(dir ==dir1){
					direction = true;
				}else{
					direction = false;
				}
				std::cout << "position " << it->first <<std::endl;
				if((it->first - last_position) < ALLOWED_GAP && direction == true){
					if(vectorOfcenters.size()!=0){
						std::cout << "smaller than ALLOWED_GAP! "<<std::endl;
					}
					size_t index;
					for(size_t k= 0; k < center_index.size(); k++){
						if( it->second == center_index.at(k)){
							index = k;
							break;
						}
					}
					assert(index >=0 && index < center_index.size());
				/*	int index;
					std::map<std::string, int>::iterator it2 = oriented_index.find(it->second);
					assert(it2 != oriented_index.end());
					index = it2->second;*/
					vectorOfcenters.push_back(index);
					position_of_each_suffix.push_back(it->first);
					std::cout << "index "<< index << std::endl;
					
				}else{
					std::cout<< "bigger than ALLOWED_GAP: " <<std::endl;
					for(size_t i =0; i < vectorOfcenters.size();i++){
						std::cout << vectorOfcenters.at(i)<< " ";
					}
					std::cout << "-------"<<std::endl;
					AllConnectedOnes.insert(make_pair(first_position, vectorOfcenters));
					//Creates suffix and their related position , then clear the 'position_of_each_suffix' 
					if(vectorOfcenters.size()>1){
						for(size_t j =0; j < vectorOfcenters.size(); j++){
							std::vector<size_t> suffix;
							for(size_t k =j; k < vectorOfcenters.size();k++){
								suffix.push_back(vectorOfcenters.at(k));
							}
							initial_suffixes.at(i).insert(make_pair(suffix,position_of_each_suffix.at(j)));
						}
					}
					position_of_each_suffix.clear();
					vectorOfcenters.clear();
				}
				std::map<size_t, pw_alignment*>::iterator it2=AlignmentsFromClustering.at(i).find(it->first);
				pw_alignment * p = it2->second;
				size_t left1,left2,right1,right2;
				p->get_lr1(left1,right1);
				p->get_lr2(left2,right2);
				if(left1 == it->first && p->getreference1()== i){
					last_position = right1;
					if(p->getbegin1()<p->getend1()){
						dir = 1;
					}else{
						dir = 0;
					}
				}
				if(left2 == it->first && p->getreference2() == i){
					last_position = right2;
					if(p->getbegin2()< p->getend2()){
						dir = 1;
					}else{
						dir = 0;
					}
				}
				std::cout << "last position "<< last_position << std::endl;
			}
			//If the last center had ALLOWED_GAP with its previous one and there is no more left to go through the else
			AllConnectedOnes.insert(make_pair(first_position, vectorOfcenters));
			if(vectorOfcenters.size()>1){
				for(size_t j =0; j < vectorOfcenters.size(); j++){
					std::vector<size_t> suffix;
					for(size_t k =j; k < vectorOfcenters.size();k++){
						suffix.push_back(vectorOfcenters.at(k));
					}
					initial_suffixes.at(i).insert(make_pair(suffix,position_of_each_suffix.at(j)));
				}
			}

			for(std::map<size_t, std::vector<size_t> >::iterator it = AllConnectedOnes.begin(); it != AllConnectedOnes.end(); it++){
				if(it->second.size() != 1 && it->second.size() != 0){
					std::cout<< "pos "<< it->first <<std::endl;
					 centersOfASequence.at(i).insert(make_pair(it->first, it->second));//connected centers that can be potential long centers
					for(size_t i =0; i < it->second.size();i++){
						std::cout << it->second.at(i)<< " ";
					}
					std::cout << "-------"<<std::endl;
				}
			}
			//At this stage i need to check for fully reversed centers and merge two clusters in to one.
			std::map<size_t , std::vector<size_t> > intermediate;
			for(std::map<size_t,std::vector<size_t> >::iterator it = centersOfASequence.at(i).begin() ; it != centersOfASequence.at(i).end();it++){
				if(it == intermediate.end()){
					std::string center = center_index.at(it->second.at(0));
					std::vector<std::string> center_parts;
					strsep(center, ":" , center_parts);
					unsigned int dir = atoi(center_parts.at(0).c_str());
					unsigned int ref = atoi(center_parts.at(1).c_str());
					unsigned int left = atoi(center_parts.at(2).c_str());
					for(std::map<size_t , std::vector<size_t> >::iterator it1 = centersOfASequence.at(i).begin() ; it1 != centersOfASequence.at(i).end();it1++){
						std::string cent = center_index.at(it1->second.at(it1->second.size()-1));
						std::vector<std::string> cent_parts;
						strsep(cent, ":" , cent_parts);
						unsigned int dir1 = atoi(center_parts.at(0).c_str());
						unsigned int ref1 = atoi(center_parts.at(1).c_str());
						unsigned int left1 = atoi(center_parts.at(2).c_str());
						bool Fullyreverse = true;
						if(ref == ref1 && left == left1 && dir != dir1 && it1->second.size() == it->second.size()){
							for(size_t i = 0; i < it1->second.size();i++){//We are sure that the length is higher than 1
							// All the next left and ref should be checked!
								std::string center1 = center_index.at(it1->second.at(i));
								std::vector<std::string> center1_parts;
								strsep(cent, ":" , center1_parts);
								unsigned int dir_c1 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c1 = atoi(center_parts.at(1).c_str());
								unsigned int left_c1 = atoi(center_parts.at(2).c_str());
								std::string center2 = center_index.at(it->second.at(it1->second.size()-1+i));
								std::vector<std::string> center2_parts;
								strsep(cent, ":" , center2_parts);
								unsigned int dir_c2 = atoi(center_parts.at(0).c_str());
								unsigned int ref_c2 = atoi(center_parts.at(1).c_str());
								unsigned int left_c2 = atoi(center_parts.at(2).c_str());
								if(ref_c1 == ref_c2 && left_c1 == left_c2 && dir_c1 != dir_c2){
								
								}else{
									Fullyreverse = false;
									break;
								}
							}
						}else{
							Fullyreverse = false;
						}
						if(Fullyreverse == true){//The fully reverse exists and we keep one of them
							std::cout << "fully reversed! "<<std::endl;
							intermediate.insert(make_pair(it1->first,it1->second));
							break;
						}
					}
				}

			}
			//Remove the intermediate from centersOfASequence
			for(std::map< size_t , std::vector<size_t> >::iterator it = intermediate.begin();it != intermediate.end();it++){
				centersOfASequence.at(i).erase(it);
			}
			//Print the potential long centers:
			for(std::map< size_t , std::vector<size_t> >::iterator it = centersOfASequence.at(i).begin();it != centersOfASequence.at(i).end();it++){
				std::cout << "position "<< it->first <<std::endl;
				for(size_t i =0; i < it->second.size();i++){
					std::cout << it->second.at(i)<<  " ";
				}
				std::cout << ""<<std::endl;
			}
		}
	}
	std::map< size_t, std::string> finding_centers::get_sequence_centers(size_t& id)const{
		return centersOnSequence.at(id);
	}
	std::string finding_centers::find_center_name(size_t & centerIndex)const{
	/*	std::string center;
		for(std::map<std::string, int>::iterator it = oriented_index.begin(); it != oriented_index.end(); it++){
			if(it->second == centerIndex){
				center = it->first;
				break;
			}
		}
		assert(center.size()>0);
		return center;*/
		return center_index.at(centerIndex);
	}
	std::vector<size_t> finding_centers::get_long_center_position(size_t & seq_id , std::vector<std::string> & long_center){
		std::vector<size_t> position;
		std::vector<size_t> centers;//TODO it is too bad should be changed to another data container like a map.
		for(size_t i =0; i < long_center.size();i++){
			for(size_t j =0; j < center_index.size();j++){
				if(center_index.at(j)==long_center.at(i)){
					centers.push_back(j);
					break;
				}
			}
		}
		pair<std::multimap<std::vector<size_t>,size_t>::iterator , std::multimap<std::vector<size_t>,size_t>::iterator > it = initial_suffixes.at(seq_id).equal_range(centers);
		for(std::multimap<std::vector<size_t>,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			position.push_back(it1->second);
		}	
		return position;
	}
	std::map<size_t, std::vector<size_t> >finding_centers::get_center(size_t & seq_id)const{
		return centersOfASequence.at(seq_id);
	}
	size_t finding_centers::get_number_of_centers()const{
		return center_index.size();
	}
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
/*		if(successive_centers.at(seq_id).size() > 1){//Adding an extar char at the end of each suffix if the last center happens more than once!
			size_t last_center = successive_centers.at(seq_id).at(successive_centers.at(seq_id).size()-1);
			for(size_t i =0; i < successive_centers.at(seq_id).size()-1; i ++){
				size_t current_center = successive_centers.at(seq_id).at(i);
			//	std::cout<< " current " << current_center << " last center "<< last_center<<std::endl;
				if(current_center == last_center){
					for(size_t j =0;j < suffixes.at(seq_id).size();j++){
						std::vector<size_t> new_suffix = suffixes.at(seq_id).at(j);
						new_suffix.push_back(powerOfTwo.at(31));
						suffixes.at(seq_id).at(j) = new_suffix;
					}
					std::vector<size_t> new_suffix;
					new_suffix.push_back(powerOfTwo.at(31));
					suffixes.at(seq_id).push_back(new_suffix);
					break;
				}
			}
		}
		if(successive_centers.size()==1){
			for(size_t j =0;j < suffixes.at(seq_id).size();j++){
				std::vector<size_t> new_suffix = suffixes.at(seq_id).at(j);
						new_suffix.push_back(powerOfTwo.at(31));
						suffixes.at(seq_id).at(j) = new_suffix;
			}
			std::vector<size_t> new_suffix;
			new_suffix.push_back(powerOfTwo.at(31));
			suffixes.at(seq_id).push_back(new_suffix);
		}
		std::cout<< "seq id "<< seq_id <<std::endl;
		if(suffixes.at(seq_id).size() > 0){
			std::cout << " suffixes are " << std::endl;
			for(size_t i =0; i < suffixes.at(seq_id).size();i++){
				std::vector<std::size_t> suf = suffixes.at(seq_id).at(i);
				for( size_t j =0; j < suf.size() ; j ++){
					std::cout << suf.at(j)<< " ";
				}
					std::cout<< " " << std::endl;
			}
		}*/
	/*	std::cout<< "all the suffixes: "<<std::endl;
		for(size_t i =0; i < suffixes.at(seq_id).size(); i++){
			for(size_t j =0; j < suffixes.at(seq_id).at(i).size(); j++){
				for(size_t  k =0; k < suffixes.at(seq_id).at(i).at(j).size();k++){
					std::cout << suffixes.at(seq_id).at(i).at(j).at(k) << " ";
				}
				std::cout<< " " <<std::endl;
			}
		}*/
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
										/*		for(size_t shift = nodes.size()-1; shift > node_index;shift--){
													pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
													std::vector<size_t> counter;
													for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
														if(it1 != nodes_relation.end()){
															counter.push_back(it1->second);
														}
													}
											//		std::multimap<size_t,size_t>::iterator p2 = nodes_relation.find(shift);
													if(counter.size() > 0){
														nodes_relation.erase(shift);
														for(size_t in = 0; in < counter.size(); in++){
															intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
														}
													}
												}//NEW CODES: checked on ecoli3, seems it works! 
												size_t ItsParent = nodes.size();
												size_t ItsfirstPar = nodes.size();
												find_parent(node_index,ItsParent);
												first_parent_index(node_index,ItsfirstPar);
												if(ItsParent != nodes.size()){
													while(ItsParent > ItsfirstPar){
														std::vector<size_t> children;
														find_child_nodes(ItsParent,children);
														for(size_t in = 0; in < children.size();in++){
															if(children.at(in) > node_index){
																delete_relation(ItsParent,children.at(in));
																intermediate.insert(make_pair(ItsParent,children.at(in)+1));	
																std::cout << "ItsParent " << ItsParent << " new children index " << children.at(in)+1 <<std::endl;
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
																intermediate.insert(make_pair(ItsParent,children.at(in)+1));
																std::cout << "ItsParent " << ItsParent << " new children index " << children.at(in)+1 <<std::endl;	
															}
														}
													}
												}*/
										//		for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
										//			nodes_relation.insert(make_pair(it->first, it->second));
										//		}
												nodes_relation.insert(make_pair(node_index,node_index+1)); 
												nodes.push_back(last_node);
												size_t shift_value = 1;
												shift_first_parent(node_index, shift_value,common_part);
										/*		for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){//size -1
													std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+1));
													if(it != firstParent.end() && it->second == shift){
														it->second = it->second+1;
													}
												}
												std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
												if(it != firstParent.end() && it->second == node_index){					
													firstParent.erase(it);
													firstParent.insert(make_pair(common_part,node_index));
												}*/	
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
									/*	size_t ItsParent = nodes.size();
										size_t ItsfirstPar = nodes.size();
										find_parent(node_index,ItsParent);
										first_parent_index(node_index,ItsfirstPar);
										if(ItsParent != nodes.size()){
											while(ItsParent > ItsfirstPar){
												std::vector<size_t> children;
												find_child_nodes(ItsParent,children);
												for(size_t in = 0; in < children.size();in++){
													if(children.at(in) > node_index){
														delete_relation(ItsParent,children.at(in));
														intermediate.insert(make_pair(ItsParent,children.at(in)+1));	
														std::cout << "ItsParent1 " << ItsParent << " children index + 1 " << children.at(in)+1 <<std::endl;
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
														intermediate.insert(make_pair(ItsParent,children.at(in)+1));
														std::cout << "ItsParent " << ItsParent << " children index + 1 " << children.at(in)+1 <<std::endl;	
													}
												}
											}
										}*/
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
									/*	for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){
											std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+1));		
											if(it != firstParent.end() && it->second == shift){
												std::cout << "shift for first parent: " << shift <<std::endl;
												it->second = it->second+1;
											}
										}
										std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
										if(it != firstParent.end() && it->second == node_index){
											firstParent.erase(it);
											firstParent.insert(make_pair(common_part,node_index));	
										}*/
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
									/*	for(size_t shift = nodes.size()-3; shift > node_index;shift--){
											std::cout << "shift "<< shift << std::endl;
											pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
											std::vector<size_t> counter;
											for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												if(it1 != nodes_relation.end()){
													counter.push_back(it1->second);
												}
											}
											if(counter.size() > 0){
												nodes_relation.erase(shift);
												for(size_t in = 0; in < counter.size(); in++){
													intermediate.insert(make_pair(shift+2, counter.at(in)+2));					
												}
											}
										}*/
										//NEW_CODE: 
								/*		size_t ItsParent = nodes.size();
										size_t ItsfirstPar = nodes.size();
										find_parent(node_index,ItsParent);
										first_parent_index(node_index,ItsfirstPar);
										if(ItsParent != nodes.size()){
											while(ItsParent > ItsfirstPar){
												std::vector<size_t> children;
												find_child_nodes(ItsParent,children);
												for(size_t in = 0; in < children.size();in++){
													if(children.at(in) > node_index){
														delete_relation(ItsParent,children.at(in));
														intermediate.insert(make_pair(ItsParent,children.at(in)+2));	
														std::cout << "ItsParent1 " << ItsParent << " new children index1 " << children.at(in)+2 <<std::endl;
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
														intermediate.insert(make_pair(ItsParent,children.at(in)+2));
														std::cout << "ItsParent " << ItsParent << " new children index2 " << children.at(in)+2 <<std::endl;	
													}
												}
											}
										}*/
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
	/*	for(std::map<std::vector<size_t> , size_t>::iterator it = branch_counter.begin(); it != branch_counter.end(); it++){
			vector<size_t> br = it->first;
			for(size_t i =0; i < br.size();i++){
				cout<< br.at(i)<< " ";
			}
			cout<< " number "<<it->second;
			cout << " " <<endl;
			
		}*/
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
	merging_centers::merging_centers(all_data & d, finding_centers & cent , suffix_tree & t):data(d), centers(cent),tree(t){}
	merging_centers::~merging_centers(){}
	void merging_centers::updating_centers(std::vector<size_t> & center_string, size_t & index){//replace 'center_string' with its new 'index' where ever 'center_sting' occurs on the tree and updates the number of happening.
		//First we make a new tree ! 
		tree.create_tree(center_string,index);
		std::vector<std::vector<size_t> > nodes = tree.get_nodes();
		std::map<std::vector<size_t>, size_t > counts = tree.get_count(); // Updated number of happening for each updated string of centers(Notice that it has node number not the centers indices)
	//	std::cout << "size of the tree: " << nodes.size() << " number of branches: "<< counts.size() <<std::endl;
		std::map<std::vector<size_t>, int> intermediate;
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
		//	std::cout << "number of nodes on a path : " << br.size() <<std::endl;
		//	for(size_t i =0; i < br.size();i++){
		//		std::cout << br.at(i)<< " ";
		//	}
		//	std::cout << " "<< std::endl;
			size_t number = it->second;//number of the path happening
			std::vector<size_t> seriesOfCenters;
		//	std::cout<< "centers on these nodes are "<<std::endl;
			for(size_t j = 0 ; j < br.size(); j++){
				for(size_t k =0; k < nodes.at(br.at(j)).size();k++){
					seriesOfCenters.push_back(nodes.at(br.at(j)).at(k));
				//	std::cout << nodes.at(br.at(j)).at(k) <<std::endl;
				}
			}
			size_t power = 31 ;
			if(seriesOfCenters.at(seriesOfCenters.size()-1) == tree.get_power_of_two(power)){
				seriesOfCenters.pop_back();	
			}
		//	std::cout << "current number of centers " << seriesOfCenters.size() << std::endl;  
			std::map<std::vector<size_t>, int>::iterator it1 = gains.find(seriesOfCenters);
			if(it1 != gains.end()){//It is kept as it is
				intermediate.insert(make_pair(it1->first,it1->second));
			}else{
				size_t number_of_new_center = 0;
				for(size_t i = 0; i < seriesOfCenters.size(); i++){
					if(seriesOfCenters.at(i)== index){//??
						number_of_new_center = number_of_new_center + 1;
					}
				}
		//		size_t old_length = seriesOfCenters.size()+ (number_of_new_center*center_string.size());
		//		std::cout << "old_length " << old_length << std::endl;
//				size_t gain = number*old_length - (number + seriesOfCenters.size());
				size_t gain = number*seriesOfCenters.size() - (number + seriesOfCenters.size());
				intermediate.insert(make_pair(seriesOfCenters,gain));
			}
		}
		gains.clear();
		for(std::map<std::vector<size_t>,int>::iterator it = intermediate.begin(); it != intermediate.end(); it++){
			gains.insert(make_pair(it->first,it->second));
		}

	}

	void merging_centers::merg_gain_value(){//calculates the gain value and builds second tree and so on iteratively
		std::vector<std::vector<size_t> > nodes = tree.get_nodes();//Includes the centers of each node
		cout<< "size of the tree " << nodes.size()<<endl;
		for(size_t i =0; i < nodes.size();i++){
			for(size_t j =0; j < nodes.at(i).size();j++){
				std::cout<< nodes.at(i).at(j)<< " ";
			}
			std::cout << "" <<std::endl;
		}
		size_t original_center_numbers = centers.get_number_of_centers();
	//	std::cout << "original center number: "<< original_center_numbers << std::endl;
		std::map<std::vector<size_t> , size_t> counts = tree.get_count();//vector<size_t> shows a path(size_t s are node indices) and size_t is its number of happening.
		//calculating the initial gain values:
		for(map<vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;
			std::vector<size_t> centers;
			int gain = 0;
			for(size_t j = 0 ; j < br.size(); j++){
				for(size_t k =0; k < nodes.at(br.at(j)).size();k++){
					centers.push_back(nodes.at(br.at(j)).at(k));//push back all the centers of all the nodes on the path br
				}
			}
			size_t power = 31 ;
			if(centers.at(centers.size()-1) == tree.get_power_of_two(power)){//The extra ending character is removed
				centers.pop_back();	
			}
		//	std::cout<< "centers " ;
		//	for(size_t j =0; j < centers.size();j++){
		//		cout<< centers.at(j) << " " ;
		//	}
		//	cout << " " <<endl;
		//	std::cout<< "center size: " << centers.size() << "  number of happening: "<< number << std::endl;
			gain = number*(centers.size())-(number+centers.size());
			if(centers.size()!= 0){// We may get zero when the original path only had the extra endign char.
				gains.insert(make_pair(centers,gain));//centers-->index of centers, gain of the long center
				cout<< "gain "<<gain <<endl;
			}
		}
		int highest_gain=0;
		std::vector<size_t> highest_path;
		for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
			if(it->second > highest_gain){
					highest_gain = it->second;
					highest_path = it->first;
			}else continue;
		}
		if(highest_gain > 0){
			merged_centers.insert(make_pair(highest_path, original_center_numbers+1));//Making the first new center with a new index which is 'original_center_numbers+1'
		}
		size_t center_numbers;
		center_numbers = original_center_numbers + 1;//it will be used when we are going to insert the next megerd center to the merged_centers map.
	//	std::cout << "cent number: "<< center_numbers << std::endl;
		while(highest_gain > 0){
			std::cout << "highest gain: "<< highest_gain << std::endl; 
			std::cout <<" highest gain path " ;
			for(size_t i = 0 ; i < highest_path.size(); i++){
				std::cout << highest_path.at(i)<< " ";
			}
			std::cout << "" << std::endl;
			updating_centers(highest_path, center_numbers);//update all the strings, their number of happpening and gains! Notice that new trees are made in this function!
			highest_gain=0;
			std::map<std::vector<size_t>, size_t>::iterator hi = merged_centers.find(highest_path);
			std::vector<size_t> hi_from_map;
			hi_from_map.push_back(hi->second);//Because we dont want to use the same center again as the one with the highest gain
			assert(hi != merged_centers.end());
			for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
				if(it->second > highest_gain && it->first != hi_from_map && it->first.size() != 1){
			//		std::cout<< "gain "<< it->second <<std::endl;
					highest_gain = it->second;
					highest_path = it->first;
				}else continue;
			}
			if(gains.size() != 0){
				std::cout << "check highest path: " <<std::endl;
				for(size_t j = 0; j < highest_path.size(); j++){
					std::cout << highest_path.at(j)<< " ";
				}
			//	std::cout << " " <<std::endl;
			}
			center_numbers = center_numbers +1;
			merged_centers.insert(make_pair(highest_path, center_numbers));
		}
//		std::cout << "final result: " << std::endl;
//		for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
//			for(size_t j =0; j < it->first.size(); j++){
//				std::cout<< it->first.at(j)<< " ";
//			}
//			std::cout<< " gain is " << it->second << std::endl;
//		}
//		std::cout << "new centers: "<<std::endl;
//		for(std::map<std::vector<size_t>, size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){
//			if(it != merged_centers.end()){
//				for(size_t i = 0; i < it->first.size(); i ++){
//					std::cout << it->first.at(i)<< " ";
//				}
//				std::cout << " its index is " << it->second <<std::endl;
//			}else {std::cout << "there is no merged center! " <<std::endl;}
//			
//		}		
	}
	void merging_centers::adding_new_centers(std::vector<std::vector<std::string> > & long_centers, std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){//Filling in the long centers vector and centersPositionOnASeq which contains all the long centers 
		//First the initial tree is built:
		std::vector<size_t> h_gain;//it is only used for the first time of making tree. for the next times the center with the highest gain is used.
		size_t index = 0; // The index of the center with the highest gain is used frome the next rounds.
		tree.create_tree(h_gain,index);//First suffix tree
		merg_gain_value();
		size_t biggest_index = 0;
		std::vector<std::string> sequence_of_centers;
		std::cout << "merged centers size: "<< merged_centers.size()<<std::endl;
		for(std::map<std::vector<size_t>,size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){//The new indices are converted to their name
			std::cout << " merged_center: " << it->second << std::endl;
		//	std::vector<size_t> updated_center = it->first;
			std::vector<size_t> list;
			list = it->first;
			size_t id = list.at(0);	
			bool ThereIsStillABigID = true;
			while (ThereIsStillABigID == true){
			//	std::cout << "here!" << std::endl;
			//	std::cout << "number of original centers: "<< centers.get_number_of_centers()<<std::endl;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);
				//	std::cout << "id " << id <<std::endl;
					if(id > centers.get_number_of_centers()){
						for(std::map<std::vector<size_t>,size_t>::iterator it1 = merged_centers.begin(); it1 != merged_centers.end(); it1++){
							if(it1->second == id){
								std::vector<size_t> temp;
								std::cout << "j "<< j << std::endl;
								if(j != 0){
									for(size_t i = 0; i < j;i++){
										temp.push_back(list.at(i));	
									}
									std::cout<< "temp0: "<<std::endl;
									for(size_t i = 0; i < temp.size(); i ++){
										std::cout << temp.at(i)<< " ";
									}
									std::cout << " " <<std::endl;
								}
								for(size_t i = 0; i < it1->first.size();i++){
									temp.push_back(it1->first.at(i));
								}
								std::cout<< "temp1: "<<std::endl;
								for(size_t i = 0; i < temp.size(); i ++){
									std::cout << temp.at(i)<< " ";
								}
								std::cout << " " <<std::endl;
								for(size_t i = j+1; i < list.size();i++){
									temp.push_back(list.at(i));
								}
								std::cout<< "temp: "<<std::endl;
								for(size_t i = 0; i < temp.size(); i ++){
									std::cout << temp.at(i)<< " ";
								}
								list = temp;
								j = j + it1->first.size()-1;
								std::cout << "j1 "<< j <<std::endl;
								break;
							} else continue;
						}
					}
				}
				ThereIsStillABigID = false;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);
					std::cout << "j2 " << j <<std::endl;
					if(id > centers.get_number_of_centers()){
						ThereIsStillABigID = true;
						std::cout << "Still a new center! "<<std::endl;
						break;
					}else continue;
				}
			}
			std::cout << "list size is "<< list.size() <<std::endl;
			for(size_t j =0; j < list.size(); j++){
				size_t id = list.at(j);
				std::string center = centers.find_center_name(id);
				sequence_of_centers.push_back(center);
				std::cout << center << std::endl;
			}
			long_centers.push_back(sequence_of_centers);
			for(size_t i = 0; i < data.numSequences(); i++){
				find_new_centers(it->second,sequence_of_centers,i, centersPositionOnASeq);
			}
			sequence_of_centers.clear();
		}
	}
	void merging_centers::index_centers(std::map<std::string , std::vector<pw_alignment> > & al_of_a_ccs){
		int id = 1;
		center_index.clear();
		for(std::map<std::string , std::vector<pw_alignment> >::iterator it = al_of_a_ccs.begin();it != al_of_a_ccs.end();it++){
			std::string center = it->first;
			std::vector<std::string> cent_parts;
			strsep(center, ":" , cent_parts);
			unsigned int cent_dir = atoi(cent_parts.at(0).c_str());
			unsigned int cent_ref = atoi(cent_parts.at(1).c_str());
			unsigned int cent_left = atoi(cent_parts.at(2).c_str());
			stringstream rev_cent;
			if(cent_dir==0){
				rev_cent <<1<<":" << cent_ref <<":"<< cent_left;
			}else{
				rev_cent <<0<<":" << cent_ref <<":"<< cent_left;
			}
			std::string reverse = rev_cent.str();
			std::map<std::string,int>::iterator it1 = center_index.find(reverse);
			if(it1 == center_index.end()){
				std::map<std::string, int>::iterator it2 = center_index.find(center);
				assert(it2 == center_index.end());
				if(cent_dir == 0){
					center_index.insert(make_pair(center,id));
				}else{
					center_index.insert(make_pair(center, (-1*(id))));
				}
				id ++;
			}else{
				center_index.insert(make_pair(center,-1*(it1->second)));
			}			
		}
		std::cout << "indices"<<std::endl;
		for(std::map<std::string, int>::iterator it = center_index.begin(); it!=center_index.end();it++){
			std::cout << it->first << " " << it->second <<std::endl;
		}		
	}
	void merging_centers::add_long_centers_to_map(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers){
		for(size_t i =0; i < long_centers.size(); i ++){
			vector<std::string> long_center = long_centers.at(i);
			std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(long_center);
			if(new_cent == new_centers.end()){
				new_centers.insert(make_pair(long_center,std::vector<pw_alignment>()));
			}
		}
	}
/*	void merging_centers::set_samples(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster,std::vector<std::string> & long_center, size_t & seq_id, std::vector<bool>& smaple1, std::vector<bool>& sample2){
		for(size_t i = 0; i < long_center.size()){
			std::cout<< " center "<< long_center.at(i)<<std::endl;
			std::string center = long_center.at(i);
			std::vector<std::string> cent_parts;
			strsep(center, ":" , cent_parts);
			unsigned int cent_dir = atoi(cent_parts.at(0).c_str());
			unsigned int cent_ref = atoi(cent_parts.at(1).c_str());
			unsigned int cent_left = atoi(cent_parts.at(2).c_str());
			std::map<std::string , std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(center);
			assert(it3 != alignments_in_a_cluster.end());
			size_t left_of_a_sample;
			for(std::map<size_t , std::string>::iterator leftOfSample = centerOnSequence.at(seq_id).begin(); leftOfSample !=  centerOnSequence.at(seq_id).end(); leftOfSample++){
				if (center == leftOfSample->second && leftOfSample->first >= it->first && leftOfSample->first >= right){
					left_of_a_sample = leftOfSample->first;
					std::cout<< "left is set! " << left_of_a_sample <<std::endl;
					break;
				}
			}
			std::cout << "number of als " << it3->second.size()<<std::endl;
			//If center itself is on the ref
			for(size_t k = 0; k < it3->second.size();k++){
				pw_alignment p = it3->second.at(k);
				size_t r1,r2,l1,l2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				std::vector<bool> sample1_p = p.getsample1();
				std::vector<bool> sample2_p = p.getsample2();
				std::cout << " sample1_p.size() " << sample1_p.size()<< " " << sample2_p.size()<< " " << p.alignment_length() <<std::endl;
				std::cout<< "ref 1 "<< p.getreference1() <<" ref2 "<< p.getreference2()<<" seq_id " << seq_id << " l1 " << l1 << " l2 " << l2 << " left_of_a_sample "<< left_of_a_sample << std::endl;
				std::vector< std::vector<bool> >sample;
				p.get_reverse_complement_sample(sample);
				std::cout << "r1 " << r1 << " r2 " << r2 <<std::endl;
				if(p.getreference1()== seq_id && l1 == left_of_a_sample && p.getreference2()== cent_ref && l2 == cent_left){
					std::cout << "center on the second ref "<<std::endl;
					right = r1;
					if(p.getbegin1() < p.getend1()){
						for(size_t m =0; m < sample1_p.size(); m++){
							sample1.push_back(sample1_p.at(m));
							sample2.push_back(sample2_p.at(m));
						}
					}else{//I should always think about turning the corresponding center!
						for(size_t m = 0; m < sample.at(0).size();m++){
							sample1.push_back(sample.at(0).at(m));
							sample2.push_back(sample.at(1).at(m));
						}
					}
									break;
								}
								else if(p.getreference2() == j && l2 == left_of_a_sample && p.getreference1() == cent_ref && l1 == cent_left){
									std::cout << "center on the first ref" << std::endl;
									right = r2;
									if(p.getbegin2()<p.getend2()){
										std::cout << "forward "<<std::endl;
										for(size_t m =0; m < sample2_p.size();m++){
											sample1.push_back(sample2_p.at(m));
											sample2.push_back(sample1_p.at(m));
										}
									}else{
										std::cout << "reverse "<<std::endl;
										for(size_t m = 0; m < sample.at(1).size();m++){
											sample1.push_back(sample.at(1).at(m));
											sample2.push_back(sample.at(0).at(m));
										}
									}
									break;
								}else if(p.getreference2() == j && cent_ref == j && cent_left == l2){//Instead of using samples sequence bases should be used! Samples are including gaps!
									std::cout << "center on the ref - second sample "<<std::endl;
									right = r2;
								//	if(cent_dir == 0){ //seq from l2 to r2 for both refs
										std::vector<bool> intermediate;
										for(size_t m =l2; m <= r2; m++){
											char base = data.getSequence(j).at(m);
										//	std::cout << base <<std::endl;
											std::vector<bool> bits;
											pw_alignment::get_bits(base,bits);
										//	std::cout << "bits " << bits.at(0) << " " << bits.at(1) << " " <<bits.at(2) <<std::endl;
											for(size_t n = 0; n < 3; n++){
												intermediate.push_back(bits.at(n));
											}
										}
										std::cout << intermediate.size() << std::endl;
										for(size_t m =0 ; m < intermediate.size();m++){
											sample1.push_back(intermediate.at(m));
											sample2.push_back(intermediate.at(m));
										}
								//	}
							//already cmed		else{ //rev_comp of seq from l2 to r2 
							//			for(size_t m = 0; m < sample.at(1).size();m++){
							//				sample1.push_back(sample.at(1).at(m));
							//				sample2.push_back(sample.at(1).at(m));
							//			}
							//				std::cout << "dunno! " <<std::endl;
							//		}
									break;
								}else if(p.getreference1()== j && cent_ref == j && cent_left == l1){
									std::cout << "center on the ref - first sample "<<std::endl;
									right = r1;
							//		if(cent_dir == 0){
										std::vector<bool> intermediate;
										for(size_t m =l1; m <= r1; m++){
											char base = data.getSequence(j).at(m);
											vector<bool> bits;
											pw_alignment::get_bits(base,bits);
											for(size_t n = 0; n < 3; n++){
												intermediate.push_back(bits.at(n));
											}
										}
										for(size_t m =0 ; m < intermediate.size();m++){
											sample1.push_back(intermediate.at(m));
											sample2.push_back(intermediate.at(m));
										}
							//already cmed	}else{
								//		for(size_t m = 0; m < sample.at(0).size();m++){
								//			sample1.push_back(sample.at(0).at(m));
								//			sample2.push_back(sample.at(0).at(m));
								//		}
								//			std::cout << "dunno! " <<std::endl;
								//
								//	}
									break;
								}
								std::cout << sample1.size() << " " << sample2.size() << std::endl;
							}
							if(i != it->second.size()-1){//gap between two centers
								size_t left_of_next_center=0;
								for(std::map<size_t , std::string>::iterator leftOfNextCenter = centerOnSequence.at(j).begin(); leftOfNextCenter != centerOnSequence.at(j).end(); leftOfNextCenter++){
									if (it->second.at(i+1)== leftOfNextCenter->second){
										left_of_next_center = leftOfNextCenter->first;
										std::cout<< "left of next center is set! " << it->second.at(i+1) <<" "<< left_of_next_center <<std::endl;
										break;
									}
								}
								vector<bool> middle_part_of_sample;
								vector<bool> gap_sample;
								std::cout << "right+1 "<< right+1 << " left "<< left_of_next_center <<std::endl;
							//	size_t counter = 0;
								for(size_t m = right+1 ; m < left_of_next_center ; m++){
							//		counter +=1;
									char base = data.getSequence(j).at(m);
									char gap = '-';
									vector<bool> bits;
									vector<bool> gap_bit;
									pw_alignment::get_bits(base,bits);
									pw_alignment::get_bits(gap,gap_bit);
									for(size_t n = 0; n < 3; n++){
										middle_part_of_sample.push_back(bits.at(n));
										gap_sample.push_back(gap_bit.at(n));
									}
								}
							//	std::cout << "counter "<< counter <<std::endl;
								std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
								for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
									sample1.push_back(middle_part_of_sample.at(m));
									sample2.push_back(gap_sample.at(m));
								}
							}
						}
	}*/
	void merging_centers::create_alignment(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::vector<std::string> > > & centersPositionOnASeq, std::vector<std::map<size_t , std::string> > & centerOnSequence){
		index_centers(alignments_in_a_cluster);
		size_t artificial_ref = data.numSequences()+new_centers.size();
		add_long_centers_to_map(long_centers,new_centers);
		for(size_t i =0; i < long_centers.size();i++){
			std::cout << "current long center " << std::endl;
			for(size_t k =0; k < long_centers.at(i).size();k++){
				std::cout << long_centers.at(i).at(k)<< " ";
			}
			std::cout << " " <<std::endl;
			for(size_t j = 0 ; j < data.numSequences(); j++){
				std::cout << "num seq "<< j <<std::endl;
				for(std::map<size_t , std::vector<std::string> >::iterator it = centersPositionOnASeq.at(j).begin(); it != centersPositionOnASeq.at(j).end(); it++){//It includes all the long centers on a sequence //TODO make a better container that one doesn't need to loop over!
					size_t reverse = 0;
					if(it->second == long_centers.at(i)){//If the long center is on that sequence
						//first ref is considered on the forward strand of a sequence 
						for(size_t k =0; k < it->second.size(); k++){
							std::cout<< it->second.at(k)<<std::endl;
						}
						std::vector<bool> sample1;
						std::vector<bool> sample2;
						pw_alignment al;
						pw_alignment self_al;
						al.setreference1(j);
						al.setbegin1(it->first);
						size_t left_one = it->first; 
						size_t right_of_last_piece=0;
						size_t length = 0;
						size_t position = it->first;
						std::cout << "position " << position << std::endl;
						find_long_center_length(it->second, alignments_in_a_cluster,j,position,length,right_of_last_piece,centerOnSequence);//It returns length of the long center and its right position on the sequence
						std::cout << "end of last piece "<< right_of_last_piece << " length "<< length << std::endl;
						al.setend1(right_of_last_piece);
						size_t right_one = right_of_last_piece;
						//Up to now first reference of an pw-alignment is set. Center is always considered as the second reference.
						std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(it->second);
						assert(new_cent != new_centers.end());
						if(new_cent->second.size() == 0){//when we are at the first sequence that this long center is located	
							al.setreference2(artificial_ref);
							artificial_ref = artificial_ref+1;
							al.setbegin2(0);//Long center is considered as an artificial sequence
							al.setend2(length-1);
						}else{
							pw_alignment p1 = new_cent->second.at(0);
							al.setreference2(p1.getreference2());
							assert(p1.getreference2() >=data.numSequences());
							al.setbegin2(p1.getbegin2());//Long center is considered as an artificial sequence
							al.setend2(p1.getend2());
						}
						std::cout<< "begin2 "<< al.getbegin2() << " end2 " << al.getend2() << std::endl;
						//Second reference is already set.
						//We are setting samples here:
						size_t right = 0; //end of the previous center on a long center
						std::cout<< "it second size " << it->second.size() <<std::endl;
					//	set_smaples(it->second,j);
						for(size_t i =0; i < it->second.size();i++){
							std::cout<< " center "<< it->second.at(i)<<std::endl;
							std::string center = it->second.at(i);
							std::vector<std::string> cent_parts;
							strsep(center, ":" , cent_parts);
							unsigned int cent_dir = atoi(cent_parts.at(0).c_str());
							unsigned int cent_ref = atoi(cent_parts.at(1).c_str());
							unsigned int cent_left = atoi(cent_parts.at(2).c_str());
							std::map<std::string , std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(center);
							size_t left_of_a_sample;
							for(std::map<size_t , std::string>::iterator leftOfSample = centerOnSequence.at(j).begin(); leftOfSample !=  centerOnSequence.at(j).end(); leftOfSample++){
								if (center == leftOfSample->second && leftOfSample->first >= it->first && leftOfSample->first >= right){
									left_of_a_sample = leftOfSample->first;
									std::cout<< "left is set! " << left_of_a_sample <<std::endl;
									break;
								}
							}
							assert(it3 != alignments_in_a_cluster.end());
							std::cout << "number of als " << it3->second.size()<<std::endl;
							for(size_t k = 0; k < it3->second.size();k++){
								pw_alignment p = it3->second.at(k);
								size_t r1,r2,l1,l2;
								p.get_lr1(l1,r1);
								p.get_lr2(l2,r2);
								std::vector<bool> sample1_p = p.getsample1();
								std::vector<bool> sample2_p = p.getsample2();
								std::cout << " sample1_p.size() " << sample1_p.size()<< " " << sample2_p.size()<< " " << p.alignment_length() <<std::endl;
								std::cout<< "ref 1 "<< p.getreference1() <<" ref2 "<< p.getreference2()<<" j " << j << " l1 " << l1 << " l2 " << l2 << " left_of_a_sample "<< left_of_a_sample << std::endl;
								std::vector< std::vector<bool> >sample;
								p.get_reverse_complement_sample(sample);
								std::cout << "r1 " << r1 << " r2 " << r2 <<std::endl;
								if(p.getreference1()== j && l1 == left_of_a_sample && p.getreference2()== cent_ref && l2 == cent_left){
									std::cout << "center on the second ref "<<std::endl;
									right = r1;
									if(p.getbegin1() < p.getend1()){
										for(size_t m =0; m < sample1_p.size(); m++){
											sample1.push_back(sample1_p.at(m));
											sample2.push_back(sample2_p.at(m));
										}
									}else{//I should always think about turning the corresponding center!
										for(size_t m = 0; m < sample.at(0).size();m++){
											sample1.push_back(sample.at(0).at(m));
											sample2.push_back(sample.at(1).at(m));
											reverse++;
										}
									}
									break;
								}
								else if(p.getreference2() == j && l2 == left_of_a_sample && p.getreference1() == cent_ref && l1 == cent_left){
									std::cout << "center on the first ref" << std::endl;
									right = r2;
									if(p.getbegin2()<p.getend2()){
										std::cout << "forward "<<std::endl;
										for(size_t m =0; m < sample2_p.size();m++){
											sample1.push_back(sample2_p.at(m));
											sample2.push_back(sample1_p.at(m));
										}
									}else{
										std::cout << "reverse "<<std::endl;
										for(size_t m = 0; m < sample.at(1).size();m++){
											sample1.push_back(sample.at(1).at(m));
											sample2.push_back(sample.at(0).at(m));
											reverse ++;
										}
									}
									break;
								}else if(p.getreference2() == j && cent_ref == j && cent_left == l2 && cent_left == left_of_a_sample){//Instead of using samples sequence bases should be used! Samples are including gaps! 
									std::cout << "center on the ref - second sample "<<std::endl;
									right = r2;
								//	if(cent_dir == 0){ //seq from l2 to r2 for both refs
										std::vector<bool> intermediate;
										for(size_t m =l2; m <= r2; m++){
											char base = data.getSequence(j).at(m);
										//	std::cout << base <<std::endl;
											std::vector<bool> bits;
											pw_alignment::get_bits(base,bits);
										//	std::cout << "bits " << bits.at(0) << " " << bits.at(1) << " " <<bits.at(2) <<std::endl;
											for(size_t n = 0; n < 3; n++){
												intermediate.push_back(bits.at(n));
											}
										}
										std::cout << intermediate.size() << std::endl;
										for(size_t m =0 ; m < intermediate.size();m++){
											sample1.push_back(intermediate.at(m));
											sample2.push_back(intermediate.at(m));
										}
								//	}
								/*	else{ //rev_comp of seq from l2 to r2 
										for(size_t m = 0; m < sample.at(1).size();m++){
											sample1.push_back(sample.at(1).at(m));
											sample2.push_back(sample.at(1).at(m));
										}
											std::cout << "dunno! " <<std::endl;
									}*/
									break;
								}else if(p.getreference1()== j && cent_ref == j && cent_left == l1 && cent_left == left_of_a_sample ){ //XXX Just added cent_left == left_of_a_sample
									std::cout << "center on the ref - first sample "<<std::endl;
									right = r1;
							//		if(cent_dir == 0){
										std::vector<bool> intermediate;
										for(size_t m =l1; m <= r1; m++){
											char base = data.getSequence(j).at(m);
											vector<bool> bits;
											pw_alignment::get_bits(base,bits);
											for(size_t n = 0; n < 3; n++){
												intermediate.push_back(bits.at(n));
											}
										}
										for(size_t m =0 ; m < intermediate.size();m++){
											sample1.push_back(intermediate.at(m));
											sample2.push_back(intermediate.at(m));
										}
								/*	}else{
										for(size_t m = 0; m < sample.at(0).size();m++){
											sample1.push_back(sample.at(0).at(m));
											sample2.push_back(sample.at(0).at(m));
										}
											std::cout << "dunno! " <<std::endl;

									}*/
									break;
								}
								std::cout << sample1.size() << " " << sample2.size() << std::endl;
							}
							if(i != it->second.size()-1){//gap between two centers
								size_t left_of_next_center=0;
								for(std::map<size_t , std::string>::iterator leftOfNextCenter = centerOnSequence.at(j).begin(); leftOfNextCenter != centerOnSequence.at(j).end(); leftOfNextCenter++){
									size_t left_1,right_1;
									al.get_lr1(left_1,right_1);
									if (it->second.at(i+1)== leftOfNextCenter->second && leftOfNextCenter->first >= left_1 && leftOfNextCenter->first <= right_1){//XXX a center can happen couple of times on a sequence, the right one should be chosen!
										left_of_next_center = leftOfNextCenter->first;
										std::cout<< "left of next center is set! " << it->second.at(i+1) <<" "<< left_of_next_center <<std::endl;
										break;
									}
								}
								vector<bool> middle_part_of_sample;
								vector<bool> gap_sample;
								std::cout << "right+1 "<< right+1 << " left "<< left_of_next_center <<std::endl;
							//	size_t counter = 0;
								for(size_t m = right+1 ; m < left_of_next_center ; m++){
							//		counter +=1;
									char base = data.getSequence(j).at(m);
									char gap = '-';
									vector<bool> bits;
									vector<bool> gap_bit;
									pw_alignment::get_bits(base,bits);
									pw_alignment::get_bits(gap,gap_bit);
									for(size_t n = 0; n < 3; n++){
										middle_part_of_sample.push_back(bits.at(n));
										gap_sample.push_back(gap_bit.at(n));
									}
								}
							//	std::cout << "counter "<< counter <<std::endl;
								std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
								for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
									sample1.push_back(middle_part_of_sample.at(m));
									sample2.push_back(gap_sample.at(m));
								}
							}
						}
						std::cout << "sam1 "<< sample1.size() << " sam2 " << sample2.size() << std::endl;
						al.set_alignment_bits(sample1,sample2);
						std::cout << "alignment is: " <<std::endl;
					//	al.print();
					//	std::vector<char> seq_base;
					//	std::vector<char>al_base;
					//	std::cout << al.alignment_length()<< " " << al.getbegin1() << " " << al.getend1() << std::endl;
					//	for(size_t m = al.getbegin1(); m<=al.getend1();m++){
					//		char base = data.getSequence(j).at(m);
					//		seq_base.push_back(base);
					//	}
					//	for(size_t m =0; m < al.alignment_length();m++){
					//		if(al.get_al_ref1().at(m) != '-'){
					//	//		al_base.push_back(al.get_al_ref1().at(m));
					//		}else{
					//			std::cout<< "there is a gap!"<<std::endl;
					//		}
					//	}
					//	std::cout << "al base size "<< al_base.size() << " seq base size "<< seq_base.size() << std::endl;
					//	assert(seq_base.size() == al_base.size());
					//	for(size_t m =0; m < al_base.size();m++){
					//		assert(al_base.at(m)==seq_base.at(m));
					//	}
						std::cout << "reverse is " << reverse << " " << it->second.size()<<std::endl;
						if(reverse==it->second.size()){//TODO it is not the most thoughtful way, the better way is making it in a right direction from the scratch, it is just to see if it works this way
							pw_alignment p;
							al.get_reverse(p);
							assert((p.getbegin2()== length -1) && (p.getend2()==0));
							p.setbegin2(0);
							p.setend2(length-1);
							new_cent->second.push_back(p);
							std::cout<< "p is added "<<std::endl;
							p.print();

						}else{
							new_cent->second.push_back(al);
							std::cout << "al is added! "<<std::endl;
							al.print();
						}
					}// no break is needed at the end of this if loop becasue a long center can occur more than one time on a sequence
				}
			}
		}
		remove_fully_reverse_refs(long_centers, new_centers);//Edits new_centers
	}
	void merging_centers::remove_fully_reverse_refs(vector<vector<std::string> > & local_long_centers, std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers){
	//Translate centers to their indices and check for negetive reverse of each center. If it exists one of them is chosen and the other one is considered as the reverse of that reference.	
		std::cout << "remove reverse centers "<< std::endl;
		std::map<std::vector<int> , std::vector<string> >long_center_index;
		for(size_t j =0;j< local_long_centers.size();j++){
			std::vector<std::string> center = local_long_centers.at(j);
			std::vector<int> indices;
			for(size_t i =0; i < center.size();i++){
				std::map<std::string,int>::iterator it1 = center_index.find(center.at(i));
			//	std::cout << "cent " << it1->first << " index " << it1->second << std::endl;
				assert(it1!=center_index.end());
				indices.push_back(it1->second);
			}
			long_center_index.insert(make_pair(indices, center));
		}
		for(std::map<std::vector<int> , std::vector<std::string> >::iterator it = long_center_index.begin(); it != long_center_index.end();it++){
			for(size_t i =0 ; i < it->first.size();i++){
				std::cout << it->first.at(i) << " ";
			}
			std::cout << " " << std::endl;

		}
		for(std::map<std::vector<int> , std::vector<std::string> >::iterator it = long_center_index.begin(); it != long_center_index.end();it++){
			//Pick a vector make its negetive reverse and find it in "long_center_index" if it exists check for the one who has more als since it helps us use less reverse_flag in encoding for the others use the same virtual ref in the opposite direction
			std::vector<int> rev_indices;
			for(size_t i = it->first.size() ; i >0 ;i--){
				int rev = -1 * it->first.at(i-1);
				std::cout << rev << std::endl;
				rev_indices.push_back(rev);
			}
			for(size_t i =0; i < rev_indices.size();i++){
				std::cout << rev_indices.at(i) << " ";
			}
			std::cout << " " << std::endl;
			std::map<std::vector<int>, std::vector<std::string> >::iterator it1 = long_center_index.find(rev_indices);
			if(it1 != long_center_index.end()){//if the negative reverse exists
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it2 = new_centers.find(it1 ->second);
				assert(it2 != new_centers.end());
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it3 = new_centers.find(it->second);
				assert(it3 != new_centers.end());
				if(it2->second.size() > it3->second.size()){
					pw_alignment al = it2->second.at(0);
					for(size_t i =0; i < it3->second.size(); i++){
						pw_alignment p = it3->second.at(i);
						p.setreference2(al.getreference2());
						size_t begin = p.getbegin2();
						size_t end = p.getend2();
						p.setbegin2(end);
						p.setend2(begin);
						std::cout<< "negetive reverse exists: "<<std::endl;
						p.print();
					}
				}else{
					pw_alignment al = it3->second.at(0);
					for(size_t i =0; i < it2->second.size(); i++){
						pw_alignment p = it2->second.at(i);
						p.setreference2(al.getreference2());
						size_t begin = p.getbegin2();
						size_t end = p.getend2();
						p.setbegin2(end);
						p.setend2(begin);
						std::cout<< "negetive reverse exists: "<<std::endl;
						p.print();
					}
				}
			}
		}
	}
//The following function finds long centers on dna sequences
	void merging_centers::find_new_centers(size_t & center_indices, std::vector<std::string > & current_long_center, size_t & seq_id , std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){
	//	std::map<size_t,std::vector<size_t> > all_centers = tree.get_center_on_a_sequence(seq_id);
		std::vector<size_t> positions = centers.get_long_center_position(seq_id , current_long_center);
	//	std::vector<std::vector<std::vector<size_t> > > all_suffix = tree.get_suffix_of_sequence(seq_id);
	//	std::cout << "all the suffixes on "<< seq_id << std::endl;
	//	for(size_t i =0;i <all_suffix.size();i++){
	//		for(size_t j = 0; j < all_suffix.at(i).size();j++){
	//			for(size_t k =0; k < all_suffix.at(i).at(j).size();k++){
	//				std::cout <<all_suffix.at(i).at(j).at(k)<< " ";
	//			}
	//			std::cout << " "<<std::endl;
	//		}
	//		std::cout << " "<<std::endl;
	//	}
		std::cout<<"sequence is "<< seq_id << std::endl;
		std::cout << "current index "<< center_indices<< std::endl;
		std::cout << "long center is "<<std::endl;
		for(size_t i =0; i < current_long_center.size();i++){
			std::cout << current_long_center.at(i) << " ";
		}
		std::cout << " " <<std::endl;
		for(size_t  i =0; i < positions.size(); i++){
			std::cout << "at position "<< positions.at(i) << std::endl;
			centersPositionOnASeq.at(seq_id).insert(make_pair(positions.at(i),current_long_center));			
		}
	//	for(std::map<size_t,std::vector<size_t> >::iterator it1 = all_centers.begin(); it1!= all_centers.end();it1++){
	//			std::cout << " " <<std::endl;
	//			std::cout << "at position "<< it1->first << std::endl;
	//			for(size_t i =0; i < it1->second.size(); i++){
	//				std::cout << it1->second.at(i)<< " ";
	//			}
	//			std::cout << " " <<std::endl;
	//			if(it1->second.at(0) == center_indices){//XXX Why at(0)?
	//				centersPositionOnASeq.at(seq_id).insert(make_pair(it1->first,current_long_center));
	//				std::cout << "seq "<<seq_id << "position "<<it1->first << std::endl;
	//			}
	//			for(size_t i =0; i < it1->second.size();i++){
	//				if(it1->second.at(i)==center_indices){
//
//					}
//				}
//		}
	}
	void merging_centers::find_long_center_length(std::vector<std::string> & centers, std::map<std::string,std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & seq_id,size_t & position ,size_t & center_length, size_t & end_of_last_piece, std::vector<std::map<size_t , std::string> > & centerOnSequence){
//		end_of_last_piece = 0;
		size_t current_position  = position;
		center_length = 0;
		for(size_t i =0; i < centers.size();i++){
			std::string center = centers.at(i);
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
		//	unsigned int center_dir = atoi(center_parts.at(0).c_str());
			unsigned int center_ref = atoi(center_parts.at(1).c_str());
			unsigned int center_left = atoi(center_parts.at(2).c_str());
			std::map<std::string , std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(centers.at(i));
			std::cout << "center is "<< centers.at(i) <<std::endl;
			assert(it != alignments_in_a_cluster.end());
		//	std::cout<< "al size "<<it->second.size() << std::endl;
		//	bool selfAligned = true;
			for(size_t j =0; j < it->second.size();j++){
				pw_alignment p = it->second.at(j);
				size_t r1,r2,l1,l2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
		//		size_t centerPosition = data.getSequence(seq_id).length();
				std::cout << "l1 " << l1 << " r1 "<< r1 << " l2 "<< l2 << " r2 " << r2 << " ref1 "<<ref1 << " ref2 "<< ref2<<std::endl;
				std::cout << "seq id " << seq_id << " cent ref "<< center_ref << " current pos " << current_position << " cent lef "<< center_left << std::endl;
				//If center_ref is not the seq_id 
		/*		for ( std::map<size_t , std::string>::iterator it1 = centerOnSequence.at(seq_id).begin() ; it1 != centerOnSequence.at(seq_id).end();it1++){
					if(it1->second == centers.at(i)){//It is wrong cus a single center may happen more than once on a seq
						centerPosition = it1->first;
						std::cout << "center position " << centerPosition << std::endl;
						break;
					}
				}
				assert(centerPosition != data.getSequence(seq_id).length());*/
				if((ref1== seq_id && ref2== center_ref && l1 >= current_position && l1 <= current_position + ALLOWED_GAP && l2 ==center_left)){
					std::cout<< "on ref 1"<<std::endl;
					center_length += r2-l2+1;
					current_position = r1;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
						std::cout << "end point0 "<< end_of_last_piece << " r1 "<< r1 <<std::endl;
					}
					break;
				}
				else if(ref2== seq_id && ref1== center_ref && l2 >= current_position && l2<= current_position + ALLOWED_GAP && l1 ==center_left){
					center_length +=r1-l1+1;
					current_position = r2;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
						std::cout << "end point1 "<< end_of_last_piece << " r2 "<< r2 <<std::endl;

					}
					break;
				}
				//If we need to align the same piece against itself
				else if(ref1 == seq_id && seq_id == center_ref && center_left == l1 && center_left >= current_position && center_left <= current_position + ALLOWED_GAP){
					center_length +=r1-l1+1;
					current_position = r1;
					std::cout << "center length here! "<< center_length  <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
						std::cout << "end point2 "<< end_of_last_piece << " r1 "<< r1 <<std::endl;
					}
					break;
				}
				else if(ref2 == seq_id && seq_id == center_ref && center_left == l2 && center_left >= current_position && center_left <= current_position + ALLOWED_GAP){
					center_length +=r2-l2+1;
					current_position = r2;
					std::cout << "center length there! " << center_length <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
						std::cout << "end point3 "<< end_of_last_piece << " r2 "<< r2 <<std::endl;

					}
					break;

				}


			}
		}
	}


	

#endif


