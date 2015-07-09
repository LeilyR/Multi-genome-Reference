#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP


template<typename T>
void initial_alignment_set<T>::compute(overlap & o) {

//	compute_simple_lazy_splits(o);
//	compute_simple(o);
	compute_vcover_clarkson(o);

}




/*
	Insert remove dynamics:
	we never want to insert an alignment that was just removed (in same function level)
	alignments that were just inserted, might be removed by the next call. In that case we can permanentely delete

*/
template<typename T>
void initial_alignment_set<T>::insert_alignment_sets(overlap & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_ins, std::set<pw_alignment, compare_pw_alignment> & all_rem, std::vector<pw_alignment> & this_ins, std::vector<pw_alignment> & this_rem) {

	for(size_t i=0; i<this_ins.size(); ++i) {
		all_ins.insert(this_ins.at(i));
	}

	for(size_t i=0; i<this_rem.size(); ++i) {
		std::set<pw_alignment>::iterator findr = all_ins.find(this_rem.at(i));
		if(findr!=all_ins.end()) {
			pw_alignment al = this_rem.at(i);
			all_ins.erase(al);
			all_rem.erase(al);
			ovrlp.remove_alignment(al); 
		} else {
			all_rem.insert(this_rem.at(i));
		}

	}
}


template<typename T>
void initial_alignment_set<T>::local_undo(overlap & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_inserted, std::set<pw_alignment , compare_pw_alignment> & all_removed) {
	for(std::set<pw_alignment, compare_pw_alignment >::iterator it = all_inserted.begin(); it!=all_inserted.end(); ++it) {
		ovrlp.remove_alignment(*it);
	}

	for(std::set<pw_alignment, compare_pw_alignment >::iterator it = all_removed.begin(); it!= all_removed.end(); ++it ) {
		const pw_alignment & ral = *it;
		ovrlp.insert_without_partial_overlap(ral);
	}

}

template<typename T>
void initial_alignment_set<T>::all_push_back(std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, std::set<pw_alignment, compare_pw_alignment> & all_inserted, std::set<pw_alignment, compare_pw_alignment> & all_removed ) {
	for(std::set<pw_alignment>::iterator it = all_inserted.begin(); it!=all_inserted.end(); ++it) {
		inserted_alignments.push_back(*it);
	}

	for(std::set<pw_alignment>::iterator it = all_removed.begin(); it!= all_removed.end(); ++it ) {
		removed_alignments.push_back(*it);
	}

}

template<typename T>
void initial_alignment_set<T>::lazy_split_full_insert_step(overlap & ovrlp, size_t level, size_t & rec_calls, pw_alignment alin, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {

	std::cout << " full insert step on level " << level  << std::endl;

	std::set<pw_alignment, compare_pw_alignment> remove_als; 
	std::vector<pw_alignment> insert_als;
	splitpoints spl2(alin, ovrlp, data);
	spl2.recursive_splits();
	spl2.split_all(remove_als, insert_als);


	size_t inserted_counter = 0;
	if(remove_als.empty()) {
		local_gain = 0;
		for(size_t i=0; i<insert_als.size(); ++i) {
			double g1;
			double g2;
			common_model.gain_function(insert_als.at(i), g1, g2);
			double avg = (g1 + g2) / 2.0 - base_cost;
			if(avg > 0) {
				ovrlp.insert_without_partial_overlap(insert_als.at(i));
				inserted_alignments.push_back(insert_als.at(i));
				local_gain+=avg;
				inserted_counter++;
			}
		}		
		std::cout << "full insert: " << level << " on al length " << alin.alignment_length() << " finally inserted pieces: " << inserted_counter << " real local gain " << local_gain << std::endl;		
	// insert pieces for recursion:
	} else {
		std::set<pw_alignment, compare_pw_alignment> all_removed;
		double lost_gain = 0;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			const pw_alignment & ral = *it;
			double gain1, gain2;
			common_model.gain_function(ral, gain1, gain2);
			double rav_gain = (gain1 + gain2)/2 - base_cost;
			lost_gain += rav_gain;

			ovrlp.remove_alignment(ral);
			all_removed.insert(ral);
		//			std::cout << "REM2 " << std::endl;
		//	       		ral->print();
		//			std::cout << std::endl;	
		}

		double max_parts_gain = 0;
		std::multimap<double, pw_alignment> ordered_parts;
		for(size_t i=0; i<insert_als.size(); ++i) {
			double gain1, gain2;
			common_model.gain_function(insert_als.at(i), gain1, gain2);
			double iav_gain = (gain1 + gain2)/2 - base_cost;
			if(iav_gain > 0) {
				ordered_parts.insert(std::make_pair(iav_gain, insert_als.at(i)));
				max_parts_gain+=iav_gain;
			}

		}
		std::cout << "2nd split: " << level << " on al length " << alin.alignment_length() <<  " new remove " << remove_als.size() << " lost gain " << lost_gain <<
		       	" new insert " << ordered_parts.size() << " gain "<<  max_parts_gain << std::endl;

		double sum_of_gain = 0;
		size_t alnum = 0;

		std::set<pw_alignment, compare_pw_alignment> all_inserted;

		for(std::multimap<double, pw_alignment>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); ++it) {
			double thisgain;
			std::vector<pw_alignment> this_inserted;
			std::vector<pw_alignment> this_removed;
			pw_alignment thisal = it->second;

			std::cout <<"full insert function level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << std::endl; 
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
template<typename T>
void initial_alignment_set<T>::lazy_split_insert_step(overlap & ovrlp, size_t level, size_t & rec_calls, pw_alignment al, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain) {
// start with computing gain of al
	rec_calls++; // count number of calls
	double gain1;
	double gain2;
	common_model.gain_function(al, gain1, gain2);
	double av_al_gain = (gain1 + gain2) / 2 - base_cost;
	// we continue if information gain is possible from the current alignment
	local_gain = 0;
	if(av_al_gain > 0) {
		splitpoints spl(al, ovrlp, data);
		spl.nonrecursive_splits();
		// sets of alignments that need to be removed and inserted if we want the current alignment 
		std::set<pw_alignment, compare_pw_alignment> remove_als; // remove alignments are pointers to objects contained in the overlap structure
		std::vector<pw_alignment> insert_als;
		spl.split_all(remove_als, insert_als);

		std::cout <<"initial: level " << level << " positive gain: " << av_al_gain << " split res rem " << remove_als.size() << " ins " << insert_als.size() << std::endl;
	
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
		std::multimap<double, pw_alignment> ordered_parts;
		double max_parts_gain = 0;
		for(size_t i=0; i<insert_als.size(); ++i) {
			common_model.gain_function(insert_als.at(i), gain1, gain2);
			double iav_gain = (gain1 + gain2)/2 - base_cost;
			if(iav_gain > 0) {
				max_parts_gain += iav_gain;
				ordered_parts.insert(std::make_pair(iav_gain, insert_als.at(i)));
			}
//			std::cout << "INS " << std::endl;
//			insert_als.at(i).print();
//			std::cout << std::endl;
		}

		std::cout << "splits: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " causes to remove " << remove_als.size() << " alignments with " << lost_information_gain << 
			" gain. We could insert " << ordered_parts.size() << " small pieces. Upper bound of gain is " << max_parts_gain << std::endl;

		// here: we actually decide to insert a small part (in the next function we check more carefully for indirectly induced splits)
		if(remove_als.empty() && ordered_parts.size() == 1) {
			pw_alignment nal = ordered_parts.begin()->second;
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

			// remove before recursion
			size_t remove_count = 0;
			for(std::set<pw_alignment, compare_pw_alignment>::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
				const pw_alignment & pwa = *it;
				all_removed.insert(pwa); 			
	//			std::cout <<"level " << level << "REMOVE NOW: " << std::endl;
	//			pwa->print();
				std::cout << " ovrlp size " << ovrlp.size() << std::endl;
				
				ovrlp.remove_alignment(pwa);


				remove_count++;
			}

			std::cout << " before recursion, we have removed " << remove_count << " alignments with gain " << lost_information_gain <<  " upper bound of gain is " << upper_bound_of_gain << std::endl;
			size_t alnum = 0;
			for(std::map<double, pw_alignment>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); ++it) {
				double ub_of_thisgain = it->first;
				double thisgain;
				std::vector<pw_alignment> this_inserted;
				std::vector<pw_alignment> this_removed;
				pw_alignment thisal = it->second;

			//	thisal.print();
			//	std::cout << std::endl;

				std::cout <<"start recursion: level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << std::endl; 
				// this function only changes ovrlp if this was locally good. We may need to undo the changes if they were globally bad
				lazy_split_insert_step(ovrlp, level + 1, rec_calls,  thisal, this_inserted, this_removed, thisgain);
				double lost_gain = ub_of_thisgain - thisgain; // we have stepped down that much from the previous upper bound of gain
				all_lost_gain += lost_gain;
				total_recursive_gain += thisgain;


				std::cout <<"end recursion: level "<< level << " end rec "<< alnum << " on al length " << thisal.alignment_length() << " subtree local gain " << thisgain << " ub was " << ub_of_thisgain << " lost gain is " << lost_gain << " total lost gain is " << all_lost_gain << " total gain of this level until now: " << total_recursive_gain <<  " ub is " << upper_bound_of_gain << std::endl; 
			

				insert_alignment_sets(ovrlp, all_inserted, all_removed, this_inserted, this_removed);

				// if we have lost more than we may gain, we undo the entire recursive insertion subtree from here
				if(all_lost_gain > upper_bound_of_gain) {
					
					local_undo(ovrlp, all_inserted, all_removed);


					std::cout << "undo subtree: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " ABORT because of low upper bound of gain" << std::endl;
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

	std::cout << "take subtree: level " << level << " on al length " << al.alignment_length() << " gain " << av_al_gain << " causes to remove " << all_removed.size() << " alignments with " << lost_information_gain << 
			" gain. We want to insert " << all_inserted.size() << " small pieces. TOTAL gain is " << local_gain << std::endl;

				return;
			
			} else {
				local_gain = 0;
				local_undo(ovrlp, all_inserted, all_removed);
			}
		} // upper bound gain pos

	


	} // pos gain

}



template<typename T>
void initial_alignment_set<T>::compute_simple_lazy_splits(overlap & o){
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;

	for (size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment & al = sorted_original_als.at(i);
		double gain_of_al = 0;

		cout << " start recursive lazy insert tree on " << endl;
		al.print();
		cout << endl;

		std::vector<pw_alignment> ins_als;
		std::vector<pw_alignment> rem_als;
		size_t count_rec_calls = 0;
		lazy_split_insert_step(o,0, count_rec_calls,  al, ins_als, rem_als, gain_of_al);
		
		double gain1;
		double gain2;
		common_model.gain_function(al, gain1, gain2);
		double vgain = (gain1 + gain2) / 2 - base_cost;


		double rgain = 0;
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
		

		double checkgain = igain - rgain;

		if(!ins_als.empty()) {
			used++;
			total_gain+=gain_of_al;
			pcs_ins+=ins_als.size();
			pcs_rem+=rem_als.size();
			std::cout << " alignment taken "<< std::endl;
		} else {
			not_used++;
			std::cout << " alignment not taken " << std::endl;
		}
		
		std::cout << " on alignment " << i << " length " << al.alignment_length() << " gain " << vgain << " recursive lazy split insert results:" << std::endl;
		std::cout << " al " << i << " gain in " << vgain << " gain out " << gain_of_al << " check gain (ins - rem) " << checkgain << std::endl;
		std::cout << " out: " << rem_als.size() << " gain " << rgain << " in: " << ins_als.size() << " gain " << igain << std::endl;
		std::cout << " total gain until here: " << total_gain << " needed " << count_rec_calls<< " recursive lazy split insert calls " << std::endl;
		
		
	}
	result_gain = total_gain;
	used_alignments = used;
	std::cout << "Input size " << sorted_original_als.size() << " max gain " << max_gain << std::endl;
	std::cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << std::endl;
	std::cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << std::endl;



}



/**
	greedy test function 
**/
template<typename T>
void initial_alignment_set<T>::compute_simple(overlap & o) {
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;
//	std::cout<< "sorted alignment size" << sorted_original_als.size()<<std::endl;	
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;

		// TODO remove
		double gain1, gain2;
		common_model.gain_function(*(al), gain1, gain2);
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


		splitpoints spl(*al, o, data);
		spl.recursive_splits();
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
				o.remove_alignment(*it);
				pcs_rem++;
			}	
			for(size_t j=0; j<insert_als.size(); ++j) {
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



template<typename T>
void initial_alignment_set<T>::search_area(size_t from, size_t to, const std::multimap<size_t, size_t> & data, std::set<size_t> & res) {
	std::multimap<size_t, size_t>::const_iterator search_begin = data.lower_bound(from);
	std::multimap<size_t, size_t>::const_iterator search_end = data.upper_bound(to);

//	cout << " search fr " << from << " to " << to << endl;
	for(std::multimap<size_t, size_t>::const_iterator it = search_begin; it!= search_end; ++it) {
		res.insert(it->second);
//		cout << " found " << it->first << " al " << it->second << endl;
	}
}


template<typename T>
void initial_alignment_set<T>::overlap_graph(std::vector<std::set<size_t> > & edges) {

	std::vector<std::multimap<size_t, size_t> > all_left_pos(data.numSequences()); // reference -> alignment left pos, alignment index

	// all alignment left pos into a vector of multimaps
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment & al = sorted_original_als.at(i);
		size_t r1 = al.getreference1();
		size_t r2 = al.getreference2();
		size_t left, right;
		al.get_lr1(left, right);
		all_left_pos.at(r1).insert(std::make_pair(left, i));
		al.get_lr2(left, right);
		all_left_pos.at(r2).insert(std::make_pair(left, i));
	}

	edges = std::vector<std::set<size_t> >(sorted_original_als.size()); // has to be completed with other direction

	// for each alignment: if there is overlap, we will find a left pos of another alignment between it left and right (or the other way around)
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment & al = sorted_original_als.at(i);
		size_t r1 = al.getreference1();
		size_t r2 = al.getreference2();
		size_t left, right;
	//	cout << " on al " << i << endl;
		al.get_lr1(left, right);
		search_area(left, right, all_left_pos.at(r1), edges.at(i));
		al.get_lr2(left, right);
		search_area(left, right, all_left_pos.at(r2), edges.at(i));
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

	cout << " on a graph with " << edges.size() << " nodes and " << edge_count/2 << " edges" << endl;


}

template<typename T>
void initial_alignment_set<T>::remove_val(size_t index, std::map<size_t, double> & index_to_weight,  std::multimap<double, size_t> & weight_to_index) {
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
template<typename T>
void initial_alignment_set<T>::update_remove(size_t index, std::vector<std::set<size_t> > & edges, std::map<size_t, double> & index_to_weight, std::multimap<double, size_t> & weight_to_index) {
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

template<typename T>
void initial_alignment_set<T>::compute_vcover_clarkson(overlap & o) {
	std::vector<std::set<size_t> > edges;
	overlap_graph(edges);


	// initialize the clarkson weights (wc(v)/dc(v) in terminology of the original paper)
	std::map<size_t, double> index_to_weight;
	std::multimap<double, size_t> weight_to_index;
	vector<double> orig_weights(sorted_original_als.size());
	double total_in_weight = 0;
	for(size_t i=0; i<edges.size(); ++i) {
		pw_alignment  al = sorted_original_als.at(i);
		double gain1, gain2;
		common_model.gain_function(al, gain1, gain2);
		double vgain = (gain1+gain2)/2 - base_cost;
	//	cout << " al " << i << " gain " << vgain <<  " at " << al << endl;
		assert(vgain >= 0);
		total_in_weight += vgain;
		orig_weights.at(i) = vgain;
		double weight = vgain / edges.at(i).size();
		index_to_weight.insert(std::make_pair(i, weight));
		weight_to_index.insert(std::make_pair(weight, i));	
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
		cout << " remove: " << minind << " now removed: " << removed.size() << endl;
		cout << "orig weight " << orig_weights.at(minind) << " clarkson weight " << minit->first << endl;
		update_remove(minind, edges, index_to_weight, weight_to_index);
		cout << " size " << weight_to_index.size() << endl;


		size_t numedges = 0;
		for(size_t i=0; i<edges.size(); ++i) {
			numedges+=edges.at(i).size();
		}
		cout << " num edges now " << numedges << endl;

		if(numedges==0) {
			break;
		}

	
	}

	
	std::vector<pw_alignment > backup_soa(sorted_original_als);

	// add independent set
	size_t at = 0;
	double weight_in_independent_set = 0;
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		std::set<size_t>::iterator findi = removed.find(i);
		if(findi == removed.end()) {
			sorted_original_als.at(at) = backup_soa.at(i);
			weight_in_independent_set += orig_weights.at(i);
			at++;
			cout << " keep " << i << " weight " << orig_weights.at(i) << " degree " << edges.at(i).size() << endl;
		}
	}
	size_t indep_set_size = at;
	double remainder_weight = 0;
	// removed alignments after the others
	for(std::set<size_t>::iterator it = removed.begin(); it!=removed.end(); ++it) {
		sorted_original_als.at(at) = backup_soa.at(*it);
		remainder_weight += orig_weights.at(*it);
		at++;
	}
	assert(at == sorted_original_als.size());

	cout << " on " << backup_soa.size() << " alignments with total gain of " << total_in_weight << endl;
	cout << " found an INDEPENDENT SET of " << indep_set_size << " alignments with total gain " << weight_in_independent_set <<endl;
	cout << " remainder contains " << removed.size() << " alignments with total gain of " << remainder_weight << endl;

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

compute_cc::compute_cc(const all_data & dat): als_on_reference(dat.numSequences()), last_pos(dat.numSequences(), 0) {
	max_al_ref_length = 0;
	for(size_t i=0; i<dat.numAlignments(); ++i) {
		const pw_alignment & a = dat.getAlignment(i);
		alignments.insert(a);
		add_on_mmaps(a);
	}

}


void compute_cc::add_on_mmaps(const pw_alignment & pwa) {
	size_t ref1 = pwa.getreference1();
	size_t ref2 = pwa.getreference2();
/*	
	std::cout << "INS " <<std::endl;
	pwa->print();
	std::cout << std::endl;	
*/
	als_on_reference.at(ref1).insert(std::make_pair(pwa.getbegin1(), pwa));
	als_on_reference.at(ref1).insert(std::make_pair(pwa.getend1(), pwa));
	als_on_reference.at(ref2).insert(std::make_pair(pwa.getbegin2(), pwa));
	als_on_reference.at(ref2).insert(std::make_pair(pwa.getend2(), pwa));

	size_t left, right;
	pwa.get_lr1(left, right);
	size_t length = right - left + 1;
	if(length > max_al_ref_length) max_al_ref_length = length;
	if(right > last_pos.at(ref1)) last_pos.at(ref1) = right;
	pwa.get_lr2(left, right);
	length = right - left + 1;
	if(length > max_al_ref_length) max_al_ref_length = length;
	if(right > last_pos.at(ref2)) last_pos.at(ref2) = right;

}


void compute_cc::remove_on_mmaps(const pw_alignment & al) {
	size_t ref1 = al.getreference1();
	size_t left, right;
	al.get_lr1(left, right);

	std::pair<std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator > l1 = 
		als_on_reference.at(ref1).equal_range(left);
	for(std::multimap<size_t, pw_alignment>::iterator it = l1.first; it!=l1.second; ++it) {
		if(al.equals(it->second)) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}
	std::pair<std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator  > r1 = 
		als_on_reference.at(ref1).equal_range(left);
	for(std::multimap<size_t, pw_alignment>::iterator it = r1.first; it!=r1.second; ++it) {
		if(al.equals(it->second)) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}

	size_t ref2 = al.getreference2();
	al.get_lr2(left, right);
	std::pair<std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator  > l2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(std::multimap<size_t, pw_alignment>::iterator it = l2.first; it!=l2.second; ++it) {
		if(al.equals(it->second)) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}
	std::pair<std::multimap<size_t, pw_alignment>::iterator, std::multimap<size_t, pw_alignment>::iterator  > r2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(std::multimap<size_t, pw_alignment>::iterator it = r2.first; it!=r2.second; ++it) {
		if(al.equals(it->second)) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}

}


void compute_cc::compute(std::vector<std::set< pw_alignment , compare_pw_alignment> > & ccs) {
	std::set <pw_alignment, compare_pw_alignment> seen;
	std::multimap<size_t , std::set<pw_alignment, compare_pw_alignment> > sorter; // sorts ccs to give largest first
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment & al = *it;
		std::set<pw_alignment, compare_pw_alignment>::iterator seenal = seen.find(al);
	//	std::cout << " seen " << seen.size() << std::endl;
		if(seenal == seen.end()) {
	//		std::cout << " getcc" << std::endl;
			std::set<pw_alignment, compare_pw_alignment> cc;
			get_cc(al, cc, seen);
	//		std::cout << "FOUND CC size " << cc.size() << std::endl;
			sorter.insert(std::make_pair(cc.size(), cc));
		}
	
	}
	for(std::multimap<size_t , std::set<pw_alignment, compare_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); ++it) {
		ccs.push_back(it->second);
	}



}

void compute_cc::get_cc(const pw_alignment & al, std::set <pw_alignment, compare_pw_alignment> & cc, std::set <pw_alignment, compare_pw_alignment> & seen) {
	
	size_t left, right;
	al.get_lr1(left, right);
	cc_step(al.getreference1(), left, right, cc, seen);
	al.get_lr2(left, right);
	cc_step(al.getreference2(), left, right, cc, seen);
}




// TODO further improvements to this function are possible if we store intervals on the references in which all alignments were already processed
void compute_cc::cc_step(size_t ref, size_t left, size_t right, std::set <pw_alignment, compare_pw_alignment> & cc, std::set <pw_alignment , compare_pw_alignment>  & seen ) {
	// search bounds (where could other alignments which overlap with the current one start or end)
	// leftbound: all alignments starting at leftbound or ealier either have the end in the search interval or no overlap with the search interval
//	std::cout << " cc step " << ref << " fr " << left << " to " << right << " seen is " << seen.size() << " we are on " << alignments.size() << " alignments" <<  std::endl; 
	std::multimap<size_t, pw_alignment>::iterator searchbegin;
	std::multimap<size_t, pw_alignment>::iterator searchend;
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

	std::set <pw_alignment, compare_pw_alignment> seen1;
	std::set <pw_alignment, compare_pw_alignment> seen2;


	// search for overlap first, then do all recursive calls after overlapping alignments have been put to seen 
	// this reduces the maximal recursion level
	size_t numseen = 0;
	for(std::multimap<size_t, pw_alignment >::iterator it = searchbegin; it!=searchend; ++it) {
		const pw_alignment & al = it->second;


//		std::cout << " at " << it->first << std::endl;
		std::set <pw_alignment , compare_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) { // if current al not contained in any connected component
//			std::cout << " not seen" << std::endl;
			size_t aleft, aright;	
	//		size_t leftmost_point_of_al_on_ref = numeric_limits<size_t>::max(); 
			if(al.getreference1()==ref) {
				al.get_lr1(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					seen1.insert(al);
					cc.insert(al);
			//		std::cout << "ovlr " << cc.size() << " "  << seen.size() <<  " ref "<< ref << " : " << left << " " << right << " ovrlaps " << std::endl;
			//		al->print();
				//	al->get_lr2(aleft, aright);
				//	cc_step(al->getreference2(), aleft, aright, cc, seen);
				}
		//		if(aleft < leftmost_point_of_al_on_ref ) {
		//			 leftmost_point_of_al_on_ref = aleft;
		//		}
			}
			if(al.getreference2()==ref) {
				al.get_lr2(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					seen2.insert(al);
					cc.insert(al);
				//	std::cout << "ovlr " << cc.size() << " "  << seen.size() << " ref "<< ref << " : " << left << " " << right << " ovrlaps " << std::endl;
				//	al->print();
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
		}
	}
//	std::cout << " found overlap with " << seen1.size() << " and " << seen2.size() << " alignments, already seen before: " << numseen << std::endl;


	// now remove all seen alignments to be faster

	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		remove_on_mmaps(*it);
	}
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		remove_on_mmaps(*it);
	}

	size_t debugsum = 0;
	for(size_t i=0; i<als_on_reference.size(); ++i) {
		debugsum+=als_on_reference.at(i).size();
	}
//	std::cout << " mmaps length " <<debugsum << std::endl;

	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		const pw_alignment & al = *it;
		size_t aleft, aright;	
		al.get_lr2(aleft, aright);
		cc_step(al.getreference2(), aleft, aright, cc, seen);
	} 	
	for(std::set<pw_alignment, compare_pw_alignment>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		const pw_alignment & al = *it;
		size_t aleft, aright;	
		al.get_lr1(aleft, aright);
		cc_step(al.getreference1(), aleft, aright, cc, seen);
	} 	
	
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



template<typename tmodel>
void affpro_clusters<tmodel>::add_alignment(const pw_alignment & al) {
	// Get identifiers for both parts of the pairwise alignment
	std::stringstream sstr1;
//	std::cout<<"data1 ad in add_al: "<< & dat << std::endl;	
	size_t left1, right1;
	al.get_lr1(left1, right1);
	sstr1 << al.getreference1()<<":"<<left1;
	std::stringstream sstr2;
	size_t left2, right2;
	al.get_lr2(left2, right2);
	sstr2 << al.getreference2()<<":"<<left2;
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


}

	template<typename tmodel>
	size_t affpro_clusters<tmodel>::get_sequence_length(size_t ref_idx)const{	
		return 	sequence_lengths.at(ref_idx);
	}
	finding_centers::finding_centers(all_data & d):data(d),AlignmentsFromClustering(data.numSequences()),centersOfASequence(data.numSequences(),std::vector<size_t>()){
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				string center = it->first;
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_ref = atoi(center_parts.at(0).c_str());
				unsigned int center_left = atoi(center_parts.at(1).c_str());
				for(size_t j = 0; j < it->second.size();j++){
					pw_alignment * p = & it->second.at(j);
					size_t left1; 
					size_t left2;
					size_t right1;
					size_t right2;
					p->get_lr1(left1,right1);
					p->get_lr2(left2,right2);
					if(i != center_ref){
						if(p->getreference1()== i && p->getreference2()== center_ref&&left2==center_left){
							std::multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							}else continue;
						}else if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
							std::multimap<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
							if(it2 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							}else continue;
						}else continue;
					}
					if(i == center_ref){
						if(p->getreference1()== i && p->getreference2()== i){
							std::multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(left1);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							}
							std::multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left2);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							}
						}else{
							std::multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(center_left);
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(center_left,p));
							}else continue;
						}
					}
				}
			}
		}
	}
	void finding_centers::findMemberOfClusters(map<string,vector<pw_alignment> > & alignmentsOfClusters){
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
			for(size_t i = 0; i < it->second.size(); i++){
				size_t ref1;
				size_t ref2;
				size_t left1;
				size_t left2;
				size_t right1;
				size_t right2;
				pw_alignment p =it->second.at(i);
				p.get_lr1(left1,right1);
				p.get_lr2(left2,right2);
				ref1 = p.getreference1();
				ref2 = p.getreference2();
				stringstream sample1;
				stringstream sample2;
				sample1 << ref1 << ":" << left1;
				sample2 << ref2 << ":" << left2;
				if(sample1.str()==it->first){
					memberOfCluster.insert(make_pair(sample2.str(),it->first));
				}else{
					memberOfCluster.insert(make_pair(sample1.str(),it->first));
				}
					memberOfCluster.insert(make_pair(it->first,it->first));

				
			}
		}
	}
	void finding_centers::center_frequency(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//it basically returns indices of centers on each sequence.
		setOfAlignments(alignmentsOfClusters);
		findMemberOfClusters(alignmentsOfClusters);	
		std::vector<std::string> center_index;
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it2=alignmentsOfClusters.begin(); it2 != alignmentsOfClusters.end();it2++){
			string cent = it2->first;
			center_index.push_back(cent);
		}
		cout<< "center index size" << center_index.size()<<endl;
		for(size_t i = 0; i< data.numSequences(); i++){
			cout << "sequence: " << i << endl;
			const dnastring & sequence = data.getSequence(i);
			for(size_t n= 0; n < sequence.length(); n++){
				std::multimap<size_t, pw_alignment*>::iterator it=AlignmentsFromClustering.at(i).find(n);
				if(it != AlignmentsFromClustering.at(i).end()){
					stringstream member;
					member << i << ":" << n;
					std::map<std::string,std::string>::iterator it1 = memberOfCluster.find(member.str());
					assert(it1 != memberOfCluster.end());
					string center = it1->second;
					cout<< "center: "<< center<<endl;
					size_t cent_index = 0;
					for(size_t j =0; j < center_index.size(); j ++){
						if(center_index.at(j)==center){
							cent_index = j;
							centersOfASequence.at(i).push_back(cent_index);
							cout<< "cent_index " << cent_index<<endl;
							break;
						}
					}
				}
			}
		} 
		/*
		for(size_t i = 0; i < centersOfASequence.size(); i ++){
			cout<< " centers on sequence " << i << " are " <<endl;
			for(size_t j =0 ; j < centersOfASequence.at(i).size(); j++){
				cout<< centersOfASequence.at(i).at(j) <<endl;
			}
		}
		*/

	}
	std::vector<size_t>  finding_centers::get_center(size_t seq_id)const{
		return centersOfASequence.at(seq_id);
	}
	suffix_tree::suffix_tree(all_data & d, finding_centers & c):data(d),centers(c){
	
	}
	suffix_tree::~suffix_tree(){
	
	}
	void suffix_tree::create_suffix(size_t seq_id){
		suffixes.clear();
		std::vector<size_t> successive_centers = centers.get_center(seq_id);
		for(size_t j =0; j < successive_centers.size(); j++){
			string suffix;
			for(size_t k =j; k < successive_centers.size();k++){
				suffix +=successive_centers.at(k);
			}
			suffixes.push_back(suffix);
		}
		if(successive_centers.size() > 1){
			size_t last_center = successive_centers.at(successive_centers.size()-1);
			for(size_t i =0; i < successive_centers.size()-1; i ++){
				size_t current_center = successive_centers.at(i);
				std::cout<< " current " << current_center << " last center "<< last_center<<std::endl;
				if(current_center == last_center){
					for(size_t j =0;j < suffixes.size();j++){
						string new_suffix = suffixes.at(j) + '#';
						suffixes.at(j) = new_suffix;
					}
					suffixes.push_back("#");
					break;
				}
			}
		}
		if(successive_centers.size()==1){
			for(size_t j =0;j < suffixes.size();j++){
				string new_suffix = suffixes.at(j) + '#';
				suffixes.at(j) = new_suffix;
			}
			suffixes.push_back("#");
		}
		if(suffixes.size() > 0){
			std::cout << " suffixes are " << std::endl;
			for(size_t i =0; i < suffixes.size();i++){
				string suf = suffixes.at(i);
				for( size_t j =0; j < suf.size() ; j ++){
					std::cout << int(suf.at(j))<< " ";
				}
					std::cout<< " " << std::endl;
			}
		}
	}
	void suffix_tree::find_a_node(size_t& node_number,size_t& parent_node, std::string& node){
		for(size_t i = parent_node; i < nodes.size(); i ++){
		//	std::cout << "size of nodes: " << nodes.size() <<std::endl;
			std::string current_node = nodes.at(i);
			if(node == current_node){
				std::cout << "i "<< i << std::endl;
				node_number = i;
				break;
			}
		}
	}
	void suffix_tree::find_sibling(size_t& current_node, vector<size_t>& siblings){
		size_t CommonPar;
		for(std::multimap<size_t,size_t>::iterator par =nodes_relation.begin();par!=nodes_relation.end();par++){
			size_t common_par = par->first;
			pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(common_par);
			for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				if(it1->second == current_node){
					CommonPar = common_par;	
					break;
				}
			}
		}
		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(CommonPar);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			siblings.push_back(it1->second)	;
		}
	}
	void suffix_tree::find_parent(size_t & node_index, size_t & parent){
		for(std::multimap<size_t,size_t>::iterator par =nodes_relation.begin();par!=nodes_relation.end();par++){
			size_t common_par = par->first;
			pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(common_par);
			for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
				if(it1->second == node_index){
					parent = common_par;	
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
					std::cout<<"child_node:"<< child_node << it1->second <<std::endl;
					break;
				}
			}
			std::cout<<"nodes relation: "<<std::endl;
				for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
					std::cout << it->first <<" "<< it->second << std::endl;
			}
	}
	void suffix_tree::make_a_tree(){
		size_t active_length = 0;
		size_t active_node = 0;
		std::vector<std::string> first_parent;//To Do: remove it and try to use firstParent map instead of this vector, because they are basically the same thing
		for(size_t seq_id =0; seq_id < data.numSequences(); seq_id++){
			create_suffix(seq_id);
			std::cout<< "seq id"<<seq_id<<std::endl;
			for(size_t i = 0; i < suffixes.size(); i++){
				std::cout<< "suffix i " << i << std::endl;
				std::string current = suffixes.at(i);
				std::cout<<"nodes relation: "<<std::endl;
				for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
					std::cout << it->first <<" "<< it->second << std::endl;
				}
std::cout << " first parent map: "<<std::endl;
		for(std::map<std::string, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::string parent = it->first;
			for(size_t i =0; i< parent.size();i++){
				std::cout<< int(parent.at(i))<< " ";
			}
			std::cout << " , " << it->second << std::endl;
		}
				size_t temp = active_node;
				size_t last_common_index = 0;
				if(first_parent.size()!= 0){
					for(size_t j= 0; j< first_parent.size(); j++){
						if (current.at(0) == first_parent.at(j).at(0)){
							string current_parent = first_parent.at(j);
							size_t node_index;//current parent node index
						//	vector<size_t> parent_path;
							std::map<std::string, size_t>::iterator it = firstParent.find(current_parent);
							find_a_node(node_index,it->second,current_parent);
							cout<< " node_index " << node_index << " from map " << it-> second <<endl;
							active_node += 1; //It is used to show that the current node have another parent rather than root
							string common_part;
							size_t commonIndex = 0;
							cout<<int(current.at(0)) << " " <<int(first_parent.at(j).at(0))<<endl;
							cout <<"currentsize "<<current.size() << "currentparentsize "<< current_parent.size()<<endl;
							size_t num_kid =0;
							while(current.size()>current_parent.size()&& (current.at(0)==current_parent.at(0))){
								size_t common_index = 0;
								size_t common_piece;
								size_t first_index_after_parent;
								const size_t n_index = node_index;
								cout << "n_index: "<< n_index<<endl;
								pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator > p1 = nodes_relation.equal_range(n_index);
								size_t biggest_child =0;
								vector<size_t> all_children;
								for(std::multimap<size_t, size_t>::iterator it = p1.first; it!=p1.second; ++it){
									if(it != nodes_relation.end()){
										all_children.push_back(it->second);
										if(it->second > biggest_child){
											biggest_child = it->second;
										}else continue;
										cout << "biggest child: " << biggest_child <<endl;
									}else{cout<<"it has no kid!"<<endl;}
								}
								for(size_t k = 0 ; k < current_parent.size(); k ++){// it may only makes sense for the case that current suffix is longer than its parent and it doesnt completely fit on it or any of its child nodes.
									cout<< "parent at " << k << " is " << int(current_parent.at(k))<< " current at " << k << " is " <<int(current.at(k))<<endl;
									if(current.at(k)==current_parent.at(k)){
										last_common_index += 1;
										common_part += current.at(k);
										first_index_after_parent = k +1;
									}else break;
								}
								commonIndex = last_common_index;
								size_t adding_a_child =0;
								if(first_index_after_parent == current_parent.size()){//Parent doesn't break
									//It checkes if the current parent already has a kid!
									if(biggest_child > 0){
										for(size_t child_node = 0; child_node < all_children.size();child_node ++){
											size_t ChildNode = all_children.at(child_node);
											if(nodes.at(ChildNode).at(0)== current.at(first_index_after_parent)){
												adding_a_child ++;
												string updated_current;
												for(size_t update =first_index_after_parent; update< current.size(); update ++){
													updated_current += current.at(update);
												}
												current = updated_current;
												cout << "ChildNode: "<<ChildNode << endl;
												current_parent = nodes.at(ChildNode);
												break;
											}else continue;
										}
									}
									if(adding_a_child == 0){//Creating a new child node
										cout << "here!"<<endl;
										for(size_t shift = node_index+1; shift <nodes.size(); shift++){
											std::map<std::string, size_t>::iterator it = firstParent.find(nodes.at(shift));
											if(it != firstParent.end() && it->second == shift){
												it->second = it->second+1;
											}
										}
										std::string last_node = nodes.at(nodes.size()-1);
										for(size_t shift =nodes.size()-1 ; shift > node_index+1; shift --){
											nodes.at(shift)=nodes.at(shift-1);
										}
										nodes.push_back(last_node);
										std::string updated_current;
										for(size_t update =first_index_after_parent; update< current.size(); update ++){
											updated_current += current.at(update);
										}
										nodes.at(node_index+1)=updated_current;
										if(biggest_child > 0){
											for(size_t shift = nodes.size()-2; shift > node_index;shift--){
												pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												std::vector<size_t> counter;
												std::multimap<size_t,size_t> intermediate;
												for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														std::cout<<"child in map: "<<it1->second<<std::endl;
													}
												}
												std::cout<< "count: "<<counter.size() <<std::endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){																									intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
														std::cout << "relation in problematic part"<<std::endl;
	 													std::cout<< shift+1 << " " << counter.at(in)+1 << std::endl;
													}
													for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
														nodes_relation.insert(make_pair(it->first, it->second));
													}
												}
												size_t ItsParent;
												find_parent(shift,ItsParent);
												if(ItsParent == node_index){
													delete_relation(ItsParent,shift);
													nodes_relation.insert(make_pair(ItsParent,shift+1));
													std::cout<< "shifting children"<<endl;
													std::cout << ItsParent<< " " <<shift+1 <<std::endl;
												}
											}
											nodes_relation.insert(make_pair(node_index,node_index+1));
											for(size_t m = 0; m < nodes.size(); m++){
												string node = nodes.at(m);
												std::cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													std::cout << int(node.at(l)) << " " ;
												}
											}
											std::cout << " " << std::endl;	
										}else{//get sure you creat a node with # later!
											for(size_t shift = node_index+1; shift <nodes.size(); shift++){
												std::map<std::string, size_t>::iterator it = firstParent.find(nodes.at(shift));
												if(it != firstParent.end() && it->second == shift){
													it->second = it->second+1;
												}
											}
											std::string last_node = nodes.at(nodes.size()-1);
											for(size_t shift =nodes.size()-1 ; shift > node_index+1; shift --){
												nodes.at(shift)=nodes.at(shift-1);
											}
											nodes.push_back(last_node);
											nodes.at(node_index+1)="#";
											for(size_t shift = nodes.size()-2; shift >=node_index+1;shift--){
												std::pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												std::vector<size_t> counter;
												std::multimap<size_t,size_t> intermediate;
												for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														std::cout<<"child in map: "<<it1->second<<std::endl;
													}
												}
												std::cout<< "count: "<<counter.size() <<std::endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){																									intermediate.insert(make_pair(shift+2, counter.at(in)+2));					
														std::cout << "relation in problematic part"<<std::endl;
	 													std::cout<< shift+2 << " " << counter.at(in)+2 << std::endl;
													}
													for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
														nodes_relation.insert(std::make_pair(it->first, it->second));
													}
												}
											}
											nodes_relation.insert(std::make_pair(node_index,node_index+1));
											nodes_relation.insert(std::make_pair(node_index,node_index+2));
											std::cout<< "node index: "<< node_index << " node index+1 "<<node_index+1<< " node index +2  " << node_index +2 << std::endl;
											current = updated_current;
											nodes.at(node_index+2)=current;
											num_kid = 1;
											std::cout<< "here num kid is changed to 1"<<endl;
										}
										current = updated_current;
										break;
										//need to break while loop here!inam check kon!How?
									}
								}else{//when we need to break the current parent and make a new branch with extra context, make current kids kids of this new child(**)
									num_kid = 1;
									std::cout<< "problem in else"<<std::endl;
									std::string updated_current;
									std::cout<< "updated_current: ";
									for(size_t update =first_index_after_parent; update< current.size(); update ++){
										updated_current += current.at(update);
										std::cout<< int(current.at(update));
									}
									std::cout<< " " <<std::endl;
									std::string extra_part_of_parent;
									std::cout<< "extra part of parent: ";
									for(size_t update = first_index_after_parent; update<current_parent.size(); update ++){
										extra_part_of_parent += current_parent.at(update);
										std::cout<< int(current_parent.at(update));
									}
									std::cout << " " << std::endl;
									std::cout<< "biggest_child: "<<biggest_child<<std::endl;
									for(size_t shift = node_index+1; shift <nodes.size(); shift++){
										std::map<std::string, size_t>::iterator it = firstParent.find(nodes.at(shift));
										if(it != firstParent.end() && it->second == shift){
											it->second = shift+2;
										}
									}
									if(biggest_child > 0){
									//	for(size_t child_node = node_index+1; child_node <= biggest_child;child_node ++)
									//		nodes.at(child_node) = extra_part_of_parent + nodes.at(child_node);//in ghalate!
									//	
										std::string last_node = nodes.at(nodes.size()-1);
										std::string second_last_node = nodes.at(nodes.size()-2);
										for(size_t shift =nodes.size()-1 ; shift >= node_index+1; shift --){
											nodes.at(shift)=nodes.at(shift-2);
										}
										if(current_parent == first_parent.at(j)){
											first_parent.at(j) = common_part;
											std::map<std::string, size_t>::iterator it1 = firstParent.find(current_parent);
											if(it1 != firstParent.end()){
												firstParent.erase(it1);
											}
												std::map<std::string, size_t>::iterator it = firstParent.find(common_part);
												if(it == firstParent.end()){
													firstParent.insert(std::make_pair(common_part,node_index));
												}else{
													it->second = node_index;
												}
										}
										for(size_t shift = nodes.size()-1; shift >= node_index;shift--){
											pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
											std::vector<size_t> counter;
											std::multimap<size_t,size_t> intermediate;
											for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												if(it1 != nodes_relation.end()){
													counter.push_back(it1->second);
													std::cout<<"child in map: "<<it1->second<<std::endl;
												}
											}
											std::cout<< "count: "<<counter.size() <<std::endl;
											if(counter.size() > 0){
												nodes_relation.erase(shift);
												for(size_t in = 0; in < counter.size(); in++){
													intermediate.insert(make_pair(shift+2, counter.at(in)+2));//First insert them to an intermadiate map
													std::cout << "relation in problematic part1"<<std::endl;
	 												std::cout<< shift+2 << " " << counter.at(in)+2 << std::endl;
												}
												for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
													nodes_relation.insert(make_pair(it->first, it->second));
												}
											}
										}
										nodes_relation.insert(make_pair(node_index,node_index+1));
										nodes_relation.insert(make_pair(node_index,node_index+2));
										std::cout<< "node index: "<< node_index << " node index+1 "<<node_index+1<< " node index +2  " << node_index +2 << std::endl;
										nodes.push_back(second_last_node);
										nodes.push_back(last_node);
										nodes.at(node_index+2)= extra_part_of_parent;
										nodes.at(node_index+1)=updated_current;
										nodes.at(node_index)=common_part;
				
									}else{
										std::string last_node = nodes.at(nodes.size()-1);
										std::string second_last_node = nodes.at(nodes.size()-2);
										for(size_t shift =nodes.size()-1 ; shift > node_index+2; shift --){
											nodes.at(shift)=nodes.at(shift-2);
										}
										nodes.push_back(second_last_node);
										nodes.push_back(last_node);
										nodes.at(node_index+2) = extra_part_of_parent;
										nodes.at(node_index+1) = updated_current;
										nodes.at(node_index) = common_part;
										if(current_parent == first_parent.at(j)){
											first_parent.at(j) = common_part;
											std::map<std::string, size_t>::iterator it1 = firstParent.find(current_parent);
											if(it1 != firstParent.end()){
												firstParent.erase(it1);
											}
											std::map<std::string, size_t>::iterator it = firstParent.find(common_part);
											if(it == firstParent.end()){
												firstParent.insert(make_pair(common_part,node_index));
											}else{
												it->second = node_index;
											}
										}
										for(size_t shift = nodes.size()-3; shift > node_index;shift--){
											pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
											std::vector<size_t> counter;
											std::multimap<size_t, size_t> intermediate;
											for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												if(it1 != nodes_relation.end()){
													counter.push_back(it1->second);
												}
											}
											std::cout<< "count: "<<counter.size() <<std::endl;
											if(counter.size() > 0){
												nodes_relation.erase(shift);
												for(size_t in = 0; in < counter.size(); in++){
													intermediate.insert(make_pair(shift+2, counter.at(in)+2));							
													std::cout << "relation in problematic part2"<<std::endl;
	 												std::cout<< shift+2 << " " << counter.at(in)+2 << std::endl;
												}
												for(std::multimap<size_t, size_t>::iterator it = intermediate.begin(); it != intermediate.end();it++){
													nodes_relation.insert(make_pair(it->first,it->second));
												}
											}
										}
										nodes_relation.insert(make_pair(node_index,node_index+1));
										nodes_relation.insert(make_pair(node_index,node_index+2));
										for(size_t m = 0; m < nodes.size(); m++){
											std::string node = nodes.at(m);
											std::cout << "nodes at " << m << " : ";
											for(size_t l =0; l < node.size();l++){
												std::cout << int(node.at(l)) << " " ;
											}
										}
										std::cout << " " << std::endl;
									//	current = updated_current;
									//	current_parent = common_part;
									/*	for(size_t shift = node_index+2; shift <nodes.size(); shift++){
											map<string, size_t>::iterator it = firstParent.find(nodes.at(shift));
											if(it != firstParent.end() && it->second == shift){
												it->second = it->second+2;
											}
										}*/
									}
									current = updated_current;
									current_parent = common_part;
									break;	//break the while loop					
								}
								size_t parentInd = node_index;
								find_a_node(node_index, parentInd, current_parent);
								std::cout << "node index" << node_index << " , "<< parentInd<<std::endl;// are they the same? No, node index is from updated parent!
							}//Here is end of while loop!
							std::cout << "currentsize "<<current.size()<<std::endl;
							std::cout<< "current parent size: "<<current_parent.size()<<std::endl;
							std::cout<<"num kid: "<<num_kid<<std::endl;
							last_common_index =0;
							std::string commonPart;
							if(num_kid ==0){
								for(size_t k = 0 ; k < current.size(); k ++){
									std::cout<< "parent at " << k << " is " << int(current_parent.at(k))<< " current at " << k << " is " <<int(current.at(k))<<std::endl;
									if(current.at(k)==current_parent.at(k)){
										last_common_index += 1;
										commonPart += current.at(k);
									}else break;
								}
								std::cout<<"last common index "<<last_common_index<<std::endl;
								if(last_common_index != 0){
									std::cout << "nodes size: " << nodes.size() << "node index" << node_index << std::endl;
									nodes.at(node_index)= commonPart;
									std::cout<< "common: " << std::endl;
									for(size_t c=0; c< commonPart.size(); c++){
										std::cout<< int(commonPart.at(c))<< " " ;
									}
									std::cout << " " <<std::endl;
									if(current_parent == first_parent.at(j)){
										first_parent.at(j) = commonPart;
										std::map<std::string, size_t>::iterator it1 = firstParent.find(current_parent);
										if(it1 != firstParent.end()){
											firstParent.erase(it1);
										}
										std::map<std::string, size_t>::iterator it = firstParent.find(commonPart);
										if(it == firstParent.end()){
											firstParent.insert(make_pair(commonPart,node_index));
										}else{
											it->second = node_index;
										}
									}
									std::string other_branch;
									if((commonPart == current&&current_parent.size()>current.size())||commonPart.size() < current.size()){
									//	bool new_relation =false;
										if(commonPart == current&&current_parent.size()>current.size()){
											std::cout << "all the current is on the parent!"<<std::endl;
											for(size_t k = last_common_index; k < current_parent.size();k++){
												other_branch += current_parent.at(k);
											}
											std::cout<< "other1: " << std::endl;
											for(size_t c=0; c< other_branch.size(); c++){
												std::cout<< int(other_branch.at(c))<< " " ;
											}
											std::cout << " " <<std::endl;	
										//	new_relation = true;
										}
										if(commonPart.size() < current.size()){
											for(size_t k = last_common_index; k < current_parent.size();k++){
												other_branch += current_parent.at(k);
											}
											std::cout<< "other2: " << std::endl;
											for(size_t c=0; c< other_branch.size(); c++){
												std::cout<< int(other_branch.at(c))<< " " ;
											}
											std::cout << " " <<std::endl;
										//	new_relation = true;
										}
std::cout << " first parent map: "<<std::endl;
		for(std::map<std::string, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::string parent = it->first;
			for(size_t i =0; i< parent.size();i++){
				std::cout<< int(parent.at(i))<< " ";
			}
			std::cout << " , " << it->second << std::endl;
		}
	
										if(node_index == nodes.size()-1){
											std::cout << "when we are at if"<<std::endl;
											nodes.push_back(other_branch);
											if(i+last_common_index <= suffixes.size()-1){
												nodes.push_back(suffixes.at(i+last_common_index));
											}else{
												nodes.push_back("#");
											}
											nodes_relation.insert(make_pair(node_index,nodes.size()-1));
											nodes_relation.insert(make_pair(node_index,nodes.size()-2));
											std::cout<< node_index << " " << nodes.size()-1 << " " << nodes.size()-2 << " " << std::endl;
											for(size_t m = 0; m < nodes.size(); m++){
												std::string node = nodes.at(m);
												std::cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													std::cout << int(node.at(l)) << " " ;
												}
											}
											std::cout << " " << std::endl;	
										}else{// all the nodes in branch counter should be shifted by two, worth it to creat a new function that does that
											std::cout<< "if we are at else"<<std::endl;
											std::string last_node = nodes.at(nodes.size()-1);
											std::string second_last_node = nodes.at(nodes.size()-2);
											std::map<std::string,size_t> InterMediate;
											for(size_t shift = node_index+1; shift <nodes.size(); shift++){
												std::map<std::string, size_t>::iterator it = firstParent.find(nodes.at(shift));
												if(it != firstParent.end() && it->second == shift){
													InterMediate.insert(make_pair(nodes.at(shift),it->second+2));
												}
											}
											for(std::map<std::string , size_t>::iterator p1 = InterMediate.begin(); p1 != InterMediate.end(); p1++){
												std::map<std::string , size_t>::iterator p2 = firstParent.find(p1->first);
												firstParent.erase(p2);
											}
											for(std::map<std::string , size_t>::iterator p1 = InterMediate.begin(); p1 != InterMediate.end(); p1++){
												firstParent.insert(make_pair(p1->first,p1->second));
											}

											pair < std::multimap<size_t,size_t>::iterator, std::multimap<size_t,size_t>::iterator > check_kids = nodes_relation.equal_range(node_index);
											size_t kid_count = 0;
										//	if(new_relation==true){
												for(std::multimap<size_t,size_t>::iterator it1 = check_kids.first ; it1 != check_kids.second; it1++){
													if(it1 != nodes_relation.end()){
														kid_count++;
													}
												}
										//	}
											for(size_t shift =nodes.size()-1 ; shift > node_index+2; shift --){
												nodes.at(shift)=nodes.at(shift-2);
											}
											nodes.push_back(second_last_node);
											nodes.push_back(last_node);
											nodes.at(node_index+2) = other_branch;
											if(i+last_common_index+commonIndex <= suffixes.size()-1){
												std::cout<< " i "<< i << " i+ last_common_index "<< i+last_common_index << "commonIndex "<< commonIndex <<std::endl;
												nodes.at(node_index+1) = suffixes.at(i+last_common_index+commonIndex);
											}else{
												nodes.at(node_index+1) = '#';
												std::cout<< "node at " << node_index + 1<< " is "  << nodes.at(node_index+1) <<std::endl;
											}
											for(size_t shift = nodes.size()-2; shift >= node_index+1;shift--){
												std::cout << "shift is: "<<shift<<std::endl;
												pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												std::vector<size_t> counter;
												std::multimap<size_t,size_t>intermediate;
												for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														std::cout<< "parent "<< shift << " , "<< it1->first<<"child in map: "<<it1->second<<std::endl;
													}
												}
												std::cout<< "count: "<<counter.size() <<std::endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){
														intermediate.insert(make_pair(shift+2, counter.at(in)+2));							
														std::cout << "er"<<std::endl;
	 													std::cout<< shift+2 << " " << counter.at(in)+2 << std::endl;
													}
													for(std::multimap<size_t, size_t>::iterator it =intermediate.begin();it !=intermediate.end();it++){
														nodes_relation.insert(make_pair(it->first,it->second));
													}
												}
												if(counter.size()==0){
													vector<size_t> siblings;
													size_t ItsParent;
													find_sibling(node_index,siblings);
													find_parent(node_index,ItsParent);
													for(size_t s = 0; s < siblings.size(); s++){
														if(siblings.at(s)==shift && siblings.at(s)> node_index){
															delete_relation(ItsParent,siblings.at(s));
															nodes_relation.insert(make_pair(ItsParent,siblings.at(s)+2));
															std::cout<< "shifting bigger siblings"<<endl;
															std::cout << ItsParent<< " " <<siblings.at(s)+2 <<std::endl;
														}
													}
												}
											}
											std::cout<< "nodeindex: "<< node_index << std::endl;//always check before inserting
											nodes_relation.insert(make_pair(node_index, node_index+1));
											nodes_relation.insert(make_pair(node_index, node_index+2));
											std::cout << "node index +1 " << node_index+1 << " node index +2 "<< node_index+2 <<std::endl;
											if(kid_count> 0){
												std::cout<< "has sub kids! "<<std::endl;
												for(size_t i = 0;  i < kid_count; i++){
													nodes_relation.insert(make_pair(node_index+2, node_index+3+i));
													std::cout << "node index +2 " << node_index+2 << " node index +3 +i"<< node_index+3 + i <<std::endl;
												}
											}
											for(size_t m = 0; m < nodes.size(); m++){
												std::string node = nodes.at(m);
												std::cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													std::cout << int(node.at(l)) << " " ;
													}
											}
											std::cout << " " << std::endl;		
										}
									}	
								}
							}
							break;
						}
					}
				}
				if(temp==active_node){//it is used only for adding a new branch to the root.
					std::cout<< " node size: " << nodes.size() << " current: " << std::endl;
					for(size_t l =0; l < current.size(); l++){
						std::cout<< int(current.at(l))<< " ";
					}
					std::cout << " " << std::endl;
					first_parent.push_back(current);
					nodes.push_back(current);
				//	vector<size_t> path;
				//	path.push_back(nodes.size()-1);
				//	branch_counter.insert(make_pair(path,1));
					std::map<std::string, size_t>::iterator it = firstParent.find(current);
					if(it == firstParent.end()){
						firstParent.insert(make_pair(current,nodes.size()-1));
					}else{
						it->second = nodes.size()-1;
					}
				}
			}
		}
		std::cout<<"first parents: "<<std::endl;
		for(size_t k=0; k < first_parent.size(); k++){
			for(size_t i =0; i< first_parent.at(k).size();i++){
				std::cout<< int(first_parent.at(k).at(i))<< " ";
			}
			std::cout << " , " ;
		}
		std::cout<< " "<<std::endl;
		std::cout << " first parent map: "<<std::endl;
		for(std::map<std::string, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::string parent = it->first;
			for(size_t i =0; i< parent.size();i++){
				std::cout<< int(parent.at(i))<< " ";
			}
			std::cout << " , " << it->second << std::endl;
		}
		std::cout<< "all the nodes:"<<std::endl;
		for(size_t i = 0; i < nodes.size(); i++){
			std::string node = nodes.at(i);
			std::cout << "node at " << i << " is ";
			for(size_t j =0; j < node.size();j++){
				std::cout << int(node.at(j)) << " " ;
			}
			std::cout << " " <<std::endl;
		}
		std::cout<<"nodes relation: "<<std::endl;
		for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
				std::cout << it->first <<" "<< it->second << std::endl;
		}
	}
	void suffix_tree::count_paths(){
		for(size_t seq =0; seq < data.numSequences(); seq++){
			create_suffix(seq);
			std::cout << "seq: "<<seq<<std::endl;
			for(size_t i = 0; i < suffixes.size(); i++){
				std::cout<< "i: "<<i << std::endl;
				string current = suffixes.at(i);
				vector<size_t> branch;//insert this vector to the branch _counter map
				for(map<string,size_t>::iterator it = firstParent.begin();it!=firstParent.end();it++){
					string first_parent = it->first;
					if(first_parent.at(0)==current.at(0)){
						string updated_current;
						if(first_parent.size()==current.size()){
							std::cout<<"No kid!"<<std::endl;
						}else{
							for(size_t j =first_parent.size() ; j < current.size(); j++){
								updated_current+=current.at(j);
							}
						}
							size_t parent_index = it->second;
							branch.push_back(it->second);
							current = updated_current;
						while(current.size()>0){
							multimap<size_t,size_t>::iterator check = nodes_relation.find(parent_index);
							if(check != nodes_relation.end()){
								pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(parent_index);
								for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
									size_t parentIndex = it1->second;
									if(current.at(0)==nodes.at(parentIndex).at(0)){
										branch.push_back(parentIndex);
										string updatedCurrent;
										if(nodes.at(parentIndex).size()==current.size()){
										}else{
											for(size_t j = nodes.at(parentIndex).size(); j<current.size(); j++){
												updatedCurrent+=current.at(j);
											}
										}
										current = updatedCurrent;
										parent_index = parentIndex;
										break;
									}
								}
							} else{ 
								std::cout<< "has no kid"<<std::endl;
								break;
								}
						}
						break;
					}
				}
				cout<< "branch size: "<< branch.size() << endl;
				for(size_t i = 0; i < branch.size(); i++){
					vector<size_t>sub_branch;
					std::cout<< " push back to sub branch: " <<std::endl;
					for(size_t j = i; j< branch.size();j++){
						sub_branch.push_back(branch.at(j));
						std::cout<< branch.at(j)<<std::endl;
						map<vector<size_t>,size_t>::iterator it1 = branch_counter.find(sub_branch);
						if(it1 == branch_counter.end()){
							branch_counter.insert(make_pair(sub_branch,0));
							it1 = branch_counter.find(sub_branch);
						}
						it1->second = it1->second +1;
					}
				}
				
			}

		}
		cout<< "branch counter: "<<endl;
		for(map<vector<size_t> , size_t>::iterator it = branch_counter.begin(); it != branch_counter.end(); it++){
			vector<size_t> br = it->first;
			for(size_t i =0; i < br.size();i++){
				cout<< br.at(i)<< " ";
			}
			cout<< " number "<<it->second;
			cout << " " <<endl;
			
		}
	}
	std::vector<std::string> suffix_tree::get_nodes()const{
		return nodes;
	}
	std::map<std::vector<size_t>, size_t> suffix_tree::get_count()const{
		return branch_counter;
	}
	std::vector<size_t> suffix_tree::get_first_parent()const{
		std::vector<size_t> first_parent;
		for(std::map<std::string,size_t>::const_iterator it = firstParent.begin(); it != firstParent.end(); it++){
			first_parent.push_back(it->second);
		}
		sort(first_parent.begin(),first_parent.end());	
		std::cout<<"first_parent: "<< std::endl;
		for(size_t i =0 ; i < first_parent.size(); i++){
			std::cout<< first_parent.at(i)<< " " ;
		}
		std::cout << " " <<std::endl;
		return first_parent;
	}
	merging_centers::merging_centers(finding_centers & cent , suffix_tree & t):centers(cent),tree(t){}
	merging_centers::~merging_centers(){}
	void merging_centers::merg_gain_value(){
	//	tree.get_first_parent();
		vector<string> nodes = tree.get_nodes();
		cout<< "size: " << nodes.size()<<endl;
		map<vector<size_t> , size_t> counts = tree.get_count();//vector<size_t> shows a path and size_t is its number of happening.
	//	map<vector<size_t>, size_t> tree_gains;	// Just final gain values are added, vector<size_t> a path (long center) and size_t its gain.
		cout<< "first parent size: "<< tree.get_first_parent().size()<<endl;
		//for(size_t i =0; i < tree.get_first_parent().size();i++){
		//	cout << "i "<< i <<endl;
		//	size_t firstparent = tree.get_first_parent().at(i);
			map<vector<size_t>, size_t> subtree_gains;//gain values
			for(map<vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
				vector<size_t> br = it->first;
			//	if(br.at(0)== firstparent){
					size_t number = it->second;
					string centers;
					int gain = 0;
					for(size_t j = 0 ; j < br.size(); j++){
						centers += nodes.at(br.at(j));
					}
					cout<< "number "<< number << "centers " ;
					for(size_t j =0; j < centers.size();j++){
						cout<< int(centers.at(j)) << " " ;
					}
					cout << " " <<endl;
					gain = number*(centers.size())-(number+centers.size());
						subtree_gains.insert(make_pair(br,gain));
						cout<< "gain "<<gain <<endl;
			//	}
			}
			size_t highest_gain=0;
			vector<size_t>highest_path;
			for(map<vector<size_t>, size_t>::iterator it = subtree_gains.begin(); it != subtree_gains.end(); it++){
				if(it->second > highest_gain){
						highest_gain = it->second;
						highest_path = it->first;
				}else continue;
			}
			while(highest_gain > 0){
				for(map<vector<size_t>, size_t>::iterator it = subtree_gains.begin(); it != subtree_gains.end(); it++){
				//	vector<size_t>common_path;
					string center_index;
					string highest_index;
					for(size_t i =0; i < highest_path.size(); i++){
						highest_index += nodes.at(highest_path.at(i));
					}
					for(size_t i =0; i < it->first.size(); i++){
						center_index += nodes.at(it->first.at(i));
					}
					std::cout << "center_index: "<< endl;
					for(size_t i = 0; i < center_index.size(); i++){
						std::cout << center_index.at(i);
					}
					std::cout << "" <<std::endl;
					if(center_index.size() >= highest_index.size()){
						size_t common_id = center_index.size();
						for(size_t i=0; i < center_index.size(); i++){
							if(center_index.at(i)==highest_index.at(0)){
								common_id = i;
								break;
							}else continue;
						}
						string common_path;
						if(common_id < center_index.size()){
							for(size_t i=0; i < highest_index.size(); i++){
								if(center_index.at(common_id+i)==highest_index.at(i)){
									common_path += center_index.at(i);
								}else break;
							}
							if(common_path.size() == highest_index.size()){
								map<vector<size_t>, size_t>::iterator it2 = counts.find(it->first);
								assert(it2 != counts.end());
								size_t number = it->second;
								size_t gain = 0;
								gain = number*(center_index.size())-(number+center_index.size()-highest_index.size());
								subtree_gains.insert(make_pair(it->first,gain));
							}
						}
					}			
				}
				for(map<vector<size_t>, size_t>::iterator it = subtree_gains.begin(); it != subtree_gains.end(); it++){
					if(it->second > highest_gain){
						highest_gain = it->second;
						highest_path = it->first;
					}else continue;
				}
			}
	//	}
	/*	for(map<vector<size_t>,size_t>::iterator it = tree_gains.begin(); it != tree_gains.end(); it++){
			cout<< "path"<<endl;
			merged_centers.push_back(it->first);
			for(size_t i = 0 ; i < it->first.size(); i++){
				cout<< it->first.at(i)<< " ";
			}
			cout << " "<<endl;
			cout<< " gain " << it->second << endl;
		}*/
	}

	void merging_centers::merg_alignments(map<string, vector<pw_alignment> > & al_of_a_ccs){
	/*	vector<string> center_index;
		for(map<string, vector<pw_alignment> >::iterator it2=al_of_a_ccs.begin(); it2 != al_of_a_ccs.end();it2++){
			string cent = it2->first;
			center_index.push_back(cent);
		}
		for(size_t i = 0; i < data.numSequences(); i++){
			vector<size_t> successive_centers = centers.get_center(i);
			string center_1 = center_index.at(succ)

		}
		for(size_t i =0; i < merged_centers.size();i++){
			vector<pw_alignment> als;
			vector<size_t> longCenter = merged_centers.at(i);
			size_t left_long_piece =0;
			size_t ref_long_piece ;
			for(size_t j =0; j < longCenter.size(); j++){
				string center = center_index.at(longCenter.at(j));
				vector<string> split;
				strsep(center, ":" , split);
				size_t cent_ref = atoi(split.at(0).c_str());
				size_t cent_left = atoi(split.at(1).c_str());
				map<string, vector<pw_alignment> >::iterator al = al_of_a_ccs.find(center);
				size_t left;
				size_t ref;
				for(size_t k = 0; k < al->second.size();k++){
					pw_alignment & p = al->second.at(k);// find associated ref!
					if(p.getreference1() == cent_ref && cent_left == p.getbegin1()){
						left = p.getbegin2();
						ref = p.getreference2();
					}
					if(p.getreference2() == cent_ref && cent_left == p.getbegin2()){
						left = p.getbegin1();
						ref = p.getreference1();
					}
				}
			}
		}*/
	}



	

#endif
