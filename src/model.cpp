#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP


template<typename T>
void initial_alignment_set<T>::compute(overlap & o) {

	compute_simple_lazy_splits(o);
//	compute_simple(o);
//	compute_vcover_clarkson(o);

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
			std::cout << "avg "<< avg <<std::endl;
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
	std::cout<< "al in lazy split "<< std::endl;
//	al.print();
	if(av_al_gain > 0) {
		splitpoints spl(al, ovrlp, data);
		spl.nonrecursive_splits();
		// sets of alignments that need to be removed and inserted if we want the current alignment 
		std::set<pw_alignment, compare_pw_alignment> remove_als; // remove alignments are pointers to objects contained in the overlap structure
		std::vector<pw_alignment> insert_als;
		spl.split_all(remove_als, insert_als);//TODO

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
		//	std::cout << "INS  gain " <<iav_gain << std::endl;
		//	insert_als.at(i).print();
		//	std::cout << std::endl;
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
	//	al.print();
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
	double random_number;
	for(size_t i=0; i<simmatrix.size(); ++i){//This for loops was added because we observed that AFP has issue in choosing the best center when there are several centers that return the same value.
//TODO replace it with mt 19937 random number generator!
		random_number = std::rand() % 10000;
		random_number = 1 / random_number;
		std::cout<< "test" << random_number << std::endl;
		simmatrix.at(i).at(i) -= random_number;	
	}

}

	template<typename tmodel>
	size_t affpro_clusters<tmodel>::get_sequence_length(size_t ref_idx)const{	
		return 	sequence_lengths.at(ref_idx);
	}
//Finds all the centers the happen on a sequence:
	finding_centers::finding_centers(all_data & d):data(d),AlignmentsFromClustering(data.numSequences()),centersOnSequence(data.numSequences()), centersOfASequence(data.numSequences()){
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (std::map<std::string,std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				if(it->second.size() != 0){
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
	}
	void finding_centers::findMemberOfClusters(map<string,vector<pw_alignment> > & alignmentsOfClusters){
		for(std::map<std::string, std::vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
			if(it->second.size() != 0){
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
	}
	void finding_centers::center_frequency(std::map<std::string,std::vector<pw_alignment> > & alignmentsOfClusters, std::vector<std::map<size_t, std::string> > & centerOnseq){//it basically returns indices of centers on each sequence.
		setOfAlignments(alignmentsOfClusters);
		findMemberOfClusters(alignmentsOfClusters);	
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
							//centersOfASequence.at(i).push_back(cent_index);
							cout<< "cent_index " << cent_index<<endl;
							centersOnSequence.at(i).insert(make_pair(n,center));
							centerOnseq.at(i).insert(make_pair(n,center));
							break;
						}
					}
				}
			}
			size_t last_position;
			size_t first_position;
			std::map<size_t, std::vector<size_t> > AllConnectedOnes;
			vector<size_t> vectorOfcenters;
			for(std::map<size_t , std::string>::iterator it = centersOnSequence.at(i).begin(); it != centersOnSequence.at(i).end();it++){
				if(vectorOfcenters.size()==0){
					last_position = it->first;
					first_position = it->first;
				}
				if(it->first - last_position < 5){
					std::cout << "smaller than 5! "<<std::endl;
					size_t index;
					for(size_t k= 0; k < center_index.size(); k++){
						if( it->second == center_index.at(k)){
							index = k;
							break;
						}
					}
					vectorOfcenters.push_back(index);
					
				}else{
					std::cout<< "bigger than five: " <<std::endl;
					AllConnectedOnes.insert(make_pair(first_position, vectorOfcenters));
					vectorOfcenters.clear();
				}
				std::multimap<size_t, pw_alignment*>::iterator it1=AlignmentsFromClustering.at(i).find(it->first);
				pw_alignment * p = it1->second;
				last_position = it1->first + p->alignment_length()-1;
			}
			for(std::map<size_t, std::vector<size_t> >::iterator it = AllConnectedOnes.begin(); it != AllConnectedOnes.end(); it++){
				if(it->second.size() != 1){
					 centersOfASequence.at(i).insert(make_pair(it->first, it->second));
				}
			}
		}
	}
	std::map< size_t , std::string> finding_centers::get_sequence_centers(size_t& id)const{
		return centersOnSequence.at(id);
	}
	std::string finding_centers::find_center_name(size_t & centerIndex)const{
		return center_index.at(centerIndex);
	}	
	std::map<size_t, std::vector<size_t> >finding_centers::get_center(size_t seq_id)const{
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
		std::cout<< "seq id " << seq_id << std::endl;
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
						std::cout<< "new suffix "<<std::endl;
						for(size_t n =0; n < new_suffix.size();n++){
							std::cout<< new_suffix.at(n)<< " ";
						}
						std::cout << " " << std::endl;
					}
					std::vector<size_t> new_suffix;
					new_suffix.push_back(powerOfTwo.at(31));
					AllSuffixes.at(counter).push_back(new_suffix);
					break;
				}
			}			
			std::cout << "counter "<<counter <<std::endl;
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
		std::cout<< "all the suffixes: "<<std::endl;
		for(size_t i =0; i < suffixes.at(seq_id).size(); i++){
			for(size_t j =0; j < suffixes.at(seq_id).at(i).size(); j++){
				for(size_t  k =0; k < suffixes.at(seq_id).at(i).at(j).size();k++){
					std::cout << suffixes.at(seq_id).at(i).at(j).at(k) << " ";
				}
				std::cout<< " " <<std::endl;
			}
		}
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
			for(size_t i = 0; i < it->second.size(); i++){
				std::cout << it->second.at(i) << std::endl;
			}
		}
	}
	void suffix_tree::find_a_node(size_t& node_number,size_t& parent_node, std::vector<size_t>& node){
		for(size_t i = parent_node; i < nodes.size(); i ++){
		//	std::cout << "size of nodes: " << nodes.size() <<std::endl;
			std::vector<size_t> current_node = nodes.at(i);
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
					std::cout<< "parent " << parent <<std::endl;
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
					std::cout<<"child_node:"<< child_node<< " " << it1->second <<std::endl;
					break;
				}
			}
			//std::cout<<"nodes relation: "<<std::endl;
			//for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
			//	std::cout << it->first <<" "<< it->second << std::endl;
			//}
	}
	void suffix_tree::read_first_parents(std::vector<size_t> & current , std::vector<size_t> & first_parent){
		for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::vector<size_t> parent = it->first;
			if(current.at(0) == parent.at(0)){
				first_parent = it->first;
				std::cout<< "parent index: "<< it->second<<std::endl;
				break;
			}
		}
		std::cout<< "first parent in read function:" <<std::endl;
		for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			std::vector<size_t> f_parent = it->first;
			std::cout << "node index is " << it->second << " ";
			for( size_t j =0; j < f_parent.size() ; j ++){
				std::cout << f_parent.at(j)<< " ";
			}
			std::cout<< " " << std::endl;
		}
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
	void suffix_tree::find_child_nodes(size_t & index, vector<size_t>& childs){
		childs.clear();
		pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > it = nodes_relation.equal_range(index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			childs.push_back(it1->second);
			std::cout << "index "<< index << "children: " << it1->second << std::endl;
		}
	}
	void suffix_tree::first_parent_index(size_t & current_index, size_t & first_parent){
			size_t index = 0;
			for(std::map<std::vector<size_t>, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
				if(it->second <= current_index && index < it->second){
					index = it->second;
					first_parent = it->second;
				}
			}		
	}
	void suffix_tree::create_tree(std::vector<size_t> & center_with_highest_gain, size_t & highest_index ){//makea tree dependent to high_gain. For the first time high_gain == 0 and then it is replaced by new one from merging. For that purpose one may need to difine high_gain vector in main and initialize it with 0 //Is there any better solution than that?
		nodes.clear();
		nodes_relation.clear();
		firstParent.clear();
		for(size_t seq_id =0; seq_id < data.numSequences(); seq_id++){
			if(center_with_highest_gain.size()==0){
				std::cout << "highest is zero! "<<std::endl;
				successive_centers.at(seq_id) = centers.get_center(seq_id);
				create_suffix(seq_id);
			}else{
				std::cout << "when highest is " << highest_index << std::endl;
				update_successive_centers(center_with_highest_gain,highest_index,seq_id);
				create_suffix(seq_id);
			}
			if(suffixes.at(seq_id).size() != 0){
				for(size_t i = 0; i < suffixes.at(seq_id).size(); i++){
				for(size_t q = 0; q < suffixes.at(seq_id).at(i).size();q++){
					std::cout << "suffixes at " << i  << " at " << q <<std::endl;
					std::vector<size_t> first_parent;
					read_first_parents(suffixes.at(seq_id).at(i).at(q),first_parent);
					std::cout<< "f_parent "<<std::endl;
					for(size_t j =0; j < first_parent.size(); j++){
						std::cout << first_parent.at(j)<< " ";
					}
					std::cout << "" << std::endl;
					if(first_parent.size()==0){
						std::cout << "if first parent is empty"<<std::endl;
						nodes.push_back(suffixes.at(seq_id).at(i).at(q));
						firstParent.insert(make_pair(suffixes.at(seq_id).at(i).at(q),nodes.size()-1));	
					}else{
						std::cout << "else "<<std::endl;
						std::vector<size_t> common_part;
						size_t length = suffixes.at(seq_id).at(i).at(q).size();
						if(first_parent.size()<= length){
							length = first_parent.size();
						}
						for(size_t j = 0; j < length ;  j++){
							if(first_parent.at(j)==suffixes.at(seq_id).at(i).at(q).at(j)){
								common_part.push_back(first_parent.at(j));
							}else break;
						}
						std::cout << "common part " << std::endl; 
						for(size_t j =0; j < common_part.size();j++){
							std::cout<< common_part.at(j)<< " ";
						}
						std::cout << " " << std::endl;
						std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(first_parent);
						assert(it!=firstParent.end());
						size_t node_index = it->second;//It s updated later on in a way that always is equal to the parent node index
						std::cout << "node index in make tree: "<<node_index<<std::endl;
						firstParent.erase(it);
						firstParent.insert(make_pair(common_part,node_index));
					//	nodes.at(node_index) = common_part;
						std::vector<size_t> current_parent = first_parent;
						std::vector<size_t> current_string = suffixes.at(seq_id).at(i).at(q);
						bool making_tree = true;
						while(making_tree == true){						
							//The first case:
							if(length == current_parent.size() && common_part.size() == length){//Current parent is shorter than the current suffix, and the current suffix contains all of it
								std::cout << "first case! "<<std::endl;
								std::vector<size_t> other;
								for(size_t j = common_part.size(); j < current_string.size(); j++){
									other.push_back(current_string.at(j));
								}
								std::cout << "other " <<std::endl;
								for(size_t j =0 ; j < other.size(); j++){
									std::cout << other.at(j) << " " ;
								}
								std::cout << " "<<std::endl;
								std::vector<size_t> childs;
								find_child_nodes(node_index,childs);
								std::cout << "child size "<<childs.size()<<std::endl;
								if(childs.size() == 0){
									std::cout << "if it has no child node"<<std::endl;
									if(node_index == nodes.size()-1){//If it is the last node on the tree
										if(other.size() != 0){
											nodes.push_back(other);
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.push_back(new_suffix);
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											if(current_parent.at(current_parent.size()-1) != powerOfTwo.at(31)){
												nodes.push_back(new_suffix);
												nodes_relation.insert(make_pair(node_index,node_index+1));
												std::cout << "# is pushed back!" <<std::endl;
											}else{ // I think i should remove this else!!
												std::cout<< "we are here ! 35" << std::endl;
												making_tree = false;
												break;
											}
										}
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
											std::cout << nodes.size()<< " " << node_index+2 << std::endl;
											if(node_index == nodes.size()-2){
												nodes.push_back(new_suffix);
											}else{
												nodes.at(node_index+2) = new_suffix;
											}
											//The rest of node relation should be updated
											std::multimap<size_t,size_t> intermediate;
											for(size_t shift = nodes.size()-1; shift > node_index;shift--){
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
												}else{//checked!
													//Having same first parents and no child node!
													size_t ItsParent = nodes.size();
													size_t ItsParent1 = nodes.size();
													size_t ItsParent2 = nodes.size();
													find_parent(shift,ItsParent2);
													first_parent_index(shift,ItsParent);
													first_parent_index(node_index,ItsParent1);
													std::cout << "shift "<<shift << " nodeIndex "<<node_index <<" shiftParent: "<< ItsParent << " indexParent " << ItsParent1 << "shift parent "<< ItsParent2<<std::endl;
													if(ItsParent == ItsParent1&& ItsParent != nodes.size()&& ItsParent2 != nodes.size()){
														std::cout<<"Node had no kids but they have same parent" << "shift + 2 " << shift+2 <<std::endl;
														delete_relation(ItsParent2,shift);
														intermediate.insert(make_pair(ItsParent2,shift+2));
													}	
												}
											}
											for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
												nodes_relation.insert(make_pair(it->first, it->second));
											}
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
											std::cout << " node relation for node  " << node_index <<std::endl;
											pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(node_index);
											for(std::multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												std::cout << it1->second << " " ;
											}
											std::cout << " " <<std::endl;
											if(node_index != nodes.size()-2){
												nodes.push_back(second_last_node);
											}
											nodes.push_back(last_node);
											std::cout << "node size " << nodes.size() << std::endl;
											for(size_t shift = node_index+3; shift <nodes.size(); shift++){//was node_index + 1
												std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift));//changed it from shift + 1 to shift+3!
												if(it != firstParent.end() && it->second == shift-2){
													std::cout << "shift - 2 " << shift-2 <<std::endl;
													it->second = shift;
												}
											}
											std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
											if(it != firstParent.end() && it->second == node_index){					
												firstParent.erase(it);
												firstParent.insert(make_pair(common_part,node_index));
												std::cout << "common part is added to the node " << node_index << std::endl;
											}	
											nodes.at(node_index) = common_part;
											std::cout<<"nodes relation1: "<<std::endl;
											for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
												std::cout << it->first <<" "<< it->second << std::endl;
											}
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											if(current_parent.at(current_parent.size()-1) != powerOfTwo.at(31)){
												for(size_t j =nodes.size()-1; j > node_index+1; j--){
													nodes.at(j) = nodes.at(j-1);
												}
												std::cout << "here!"<<std::endl;
												nodes.at(node_index +1) =new_suffix;
												//The rest of node relation should be updated
												std::multimap<size_t,size_t> intermediate;
												for(size_t shift = nodes.size()-1; shift > node_index;shift--){
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
													}else{
												//	if(p1.first == nodes_relation.end())//checked!
														std::cout << "end of the map"<<std::endl;
														size_t ItsParent = nodes.size();// Replaced it with first parent instead
														size_t ItsParent1 = nodes.size();
														size_t ItsParent2 = nodes.size();
														find_parent(shift,ItsParent2);
													//	find_parent(shift,ItsParent);
													//	find_parent(node_index,ItsParent1);
														first_parent_index(shift,ItsParent);
														first_parent_index(node_index,ItsParent1);
														std::cout << "shift "<<shift << " nodeIndex "<<node_index <<" shiftParent: "<< ItsParent << " indexParent " << ItsParent1<< " shift parent "<< ItsParent2 <<std::endl;
														if(ItsParent == ItsParent1&& ItsParent != nodes.size() && ItsParent2 != nodes.size()){
															std::cout<<"Node had no kids but they have same parent"<<std::endl;
															delete_relation(ItsParent2,shift);
															intermediate.insert(make_pair(ItsParent2,shift+1));
														}
													}	
												}
												for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
													nodes_relation.insert(make_pair(it->first, it->second));
												}
												nodes_relation.insert(make_pair(node_index,node_index+1));
												nodes.at(node_index) = common_part;// Seems like an extar thing since we already knew the entire parent is covered. 
												nodes.push_back(last_node);
												for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){//size -1
													std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+1));
													if(it != firstParent.end() && it->second == shift){
														it->second = it->second+1;
													}
												}
												std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
												if(it != firstParent.end() && it->second == node_index){					
													firstParent.erase(it);
													firstParent.insert(make_pair(common_part,node_index));
												}	
												std::cout << "# is pushed back!" <<std::endl;
												for(size_t i = 0; i < nodes.size(); i++){
													std::vector<size_t> node = nodes.at(i);
													std::cout << "node at " << i << " is ";
													for(size_t j =0; j < node.size();j++){
														std::cout << node.at(j) << " " ;
													}
													std::cout << " " <<std::endl;
												}

											}else{// I think i should remove this else!!
												std::cout<< "we are here 35" << std::endl;
												making_tree = false;
												break;
											}
										}

									}
									making_tree = false;
								}else{//if it already has some child nodes, here i need to check all the child nodes for the common context
									size_t biggest_kid = 0;
									bool ItIsNotInAnyOfChildren = false;
									size_t new_parent_index;
									for(size_t j = 0; j < childs.size();j++){
										if(childs.at(j)>biggest_kid){
											biggest_kid = childs.at(j);
										}
									}
									for(size_t j = 0; j < childs.size();j++){
										if(other.size() != 0){
											if(nodes.at(childs.at(j)).at(0)== other.at(0)){
												ItIsNotInAnyOfChildren = true;
												current_parent = nodes.at(childs.at(j));
												new_parent_index = childs.at(j);
												break;
											}
										}
									}
									if(ItIsNotInAnyOfChildren == false){ //In this case current string is added as a new child node.
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
												//nodes.push_back(new_suffix);
												nodes.at(node_index+1)=new_suffix;
											}
										}
										nodes.push_back(last_node);
										std::multimap<size_t,size_t> intermediate;
										for(size_t shift = nodes.size()-2; shift > node_index;shift--){
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
													std::cout <<"shift+1 " << shift + 1<<std::endl;
													intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
												}
											}else{//checked!
												size_t ItsParent = nodes.size();
												size_t ItsParent1 = nodes.size();
												size_t ItsParent2 = nodes.size();
												find_parent(shift,ItsParent2);
												first_parent_index(shift,ItsParent);
												first_parent_index(node_index,ItsParent1);
												if(ItsParent == ItsParent1 && ItsParent2 != nodes.size()){
													delete_relation(ItsParent2,shift);
													nodes_relation.insert(make_pair(ItsParent2,shift+1));
													std::cout <<"shift + 1 " << shift + 1<<std::endl;
												}
											}
											
										}
										for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
										std::cout << "it -> first " << it->first << " it->second "<< it->second << std::endl;
										nodes_relation.insert(make_pair(it->first, it->second));
										}
										std::cout<< "here!"<<std::endl;
										nodes_relation.insert(make_pair(node_index,node_index+1));
										std::cout <<"node index + 1 " << node_index + 1<<std::endl;
										for(size_t shift = node_index+1; shift <nodes.size()-1; shift++){
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
										}
										nodes.at(node_index)= common_part;
										std::cout<<"nodes relation2: "<<std::endl;
										for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
											std::cout << it->first <<" "<< it->second << std::endl;
										}
										making_tree = false;
									}else{//calculate the common_part (This is the only continious condition)
										std::cout<<"continuous one! "<<std::endl;
										nodes.at(node_index) = common_part;//This is just before changing the index to the new one
										std::cout << node_index << std::endl;
										node_index = new_parent_index;
										std::cout << node_index << std::endl;
										current_string = other;
										length = current_string.size();
										if(current_parent.size()<= length){
											length = current_parent.size();
										}
										common_part.clear();
										for(size_t j = 0; j < length ;  j++){
											if(current_parent.at(j)==current_string.at(j)){
												common_part.push_back(current_parent.at(j));
											}else break;
										}
									//	if(node_index < biggest_kid){
									//		size_t ItsParent = nodes.size();
									//		find_parent(node_index,ItsParent);
									//		std::cout<< "node index " << node_index << "its parent: " << ItsParent <<std::endl;					
									//		for(size_t i = node_index + 1 ; i <= biggest_kid; i++){
									//			delete_relation(ItsParent,i);
									//			nodes_relation.insert(make_pair(ItsParent,i+2));
									//			std::cout <<" i + 2 " << i + 2<<std::endl;
									//		}
									//	}
										std::cout<<"common_part:"<<std::endl;
										for(size_t j =0; j < common_part.size(); j++){
											std::cout << common_part.at(j)<< " ";
										}
										std::cout << " " <<std::endl;
									}
								}
							}
							//The second case://The non common part is added as a child node and all old child nodes become children of the non-common part
							else if((length == current_string.size() && common_part.size() == length) || common_part.size() < length ){
								std::vector<size_t> non_common;
								std::vector<size_t> non_common_2;
								for(size_t j =common_part.size(); j < current_parent.size(); j++){
									non_common.push_back(current_parent.at(j));
								}
								if(common_part.size() < length){
									for(size_t j =common_part.size(); j < current_string.size(); j++){
										non_common_2.push_back(current_string.at(j));
									}	
								}else { std::cout <<" there is no non common 2" << std::endl;}
								std::cout << "non_common "<<std::endl;
								assert(non_common.size() != 0);
								for( size_t j =0; j < non_common.size() ; j ++){
									std::cout << non_common.at(j)<< " ";
								}
								std::cout<< " " << std::endl;
								std::cout << "non_common_2"<<std::endl;
								for( size_t j =0; j < non_common_2.size() ; j ++){
									std::cout << non_common_2.at(j)<< " ";
								}
								std::cout<< " " << std::endl;
								vector<size_t> childs;
								find_child_nodes(node_index,childs);
								if(childs.size() == 0){
									if(node_index == nodes.size()-1){//If it is the last node on the tree
									//	nodes.at(node_index) = common_part;
										if(common_part.size() < length){
											nodes.push_back(non_common_2);
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.push_back(new_suffix);
										}
										nodes.push_back(non_common);
									}else if(node_index == nodes.size()-2){
										std::cout<<"node size - 2" <<std::endl;
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
									//	nodes.at(node_index) = common_part;
										if(common_part.size() < length){
											nodes.at(node_index+1)=non_common_2;
										}else{
											std::vector<size_t> new_suffix;
											new_suffix.push_back(powerOfTwo.at(31));
											nodes.at(node_index+1)= new_suffix;
										}
										nodes.push_back(non_common);
										nodes.push_back(last_node);
									}else{//Indices of all the nodes after that should be shifted 
										std::vector<size_t> last_node = nodes.at(nodes.size()-1);
										std::vector<size_t> second_last_node = nodes.at(nodes.size()-2);
										std::cout<< "node index "<< node_index<<std::endl;
										for(size_t j =nodes.size()-1; j > node_index+2; j--){
											nodes.at(j) = nodes.at(j-2);
										}
									//	nodes.at(node_index) = common_part;
										if(common_part.size() < length){
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
										std::multimap<size_t,size_t> intermediate;
										for(size_t shift = nodes.size()-3; shift > node_index;shift--){
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
											}else{//checked!
												size_t ItsParent = nodes.size();
												size_t ItsParent1 = nodes.size();
												size_t ItsParent2 = nodes.size();
												find_parent(shift,ItsParent2); //Replaced it with first parent
											//	find_parent(node_index,ItsParent1);
												first_parent_index(shift,ItsParent);
												first_parent_index(node_index,ItsParent1);
												std::cout << "shift "<< shift << " node_index "<< node_index <<" shift-parent: "<< ItsParent << " index_parent " << ItsParent1<<std::endl;
												if(ItsParent == ItsParent1&& ItsParent != nodes.size() && ItsParent2 != nodes.size()){
													std::cout<<"Node had no kids but they have same parent"<<std::endl;
													delete_relation(ItsParent2,shift);
													intermediate.insert(make_pair(ItsParent2,shift+2));
												}	
											}
										}
										for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
											nodes_relation.insert(make_pair(it->first, it->second));
										}
									}
									nodes.at(node_index) = common_part;
									nodes_relation.insert(make_pair(node_index,node_index+1));
									nodes_relation.insert(make_pair(node_index,node_index+2));
							}else{//It has children
								std::vector<size_t> last_node = nodes.at(nodes.size()-1);// Not sure if it works for the case that node index is equal to the nodes.size()-2
								std::vector<size_t> second_last_node;
								if(nodes.size()-2 != node_index){
									second_last_node = nodes.at(nodes.size()-2);
									std::cout << node_index << " " << nodes.size()-1<< " " << nodes.size()-2 << std::endl;
									for(size_t j =nodes.size()-1; j > node_index+2; j--){
										nodes.at(j) = nodes.at(j-2);
									}
								}else{
									std::cout<< "exception!" <<std::endl;
								}
							//	nodes.at(node_index) = common_part;
								if(common_part.size() < length){
									nodes.at(node_index+1)=non_common_2;
									std::cout<< "HEya!" <<std::endl;
								}else{
									std::vector<size_t> new_suffix;
									new_suffix.push_back(powerOfTwo.at(31));
									nodes.at(node_index+1)= new_suffix;
									std::cout<< "HEYA!" <<std::endl;
								}
								if(nodes.size()-2 != node_index){
									nodes.at(node_index+2)= non_common;
									nodes.push_back(second_last_node);
								}else{
									nodes.push_back(non_common);
								}
								nodes.push_back(last_node);
								for(size_t shift = nodes.size()-2; shift > node_index;shift--){
									pair<std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
									std::vector<size_t> counter;
									std::multimap<size_t,size_t> intermediate;
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
										for(std::multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
											nodes_relation.insert(make_pair(it->first, it->second));
											std::cout << "shift "<< shift<< " it->first " << it->first <<std::endl;
										}
									}else{//checked!
										std::cout << "almost there!" << std::endl;
										size_t ItsParent = nodes.size();
										size_t ItsParent1 = nodes.size();
										size_t ItsParent2 = nodes.size();
										find_parent(shift,ItsParent2);// Replaced it with first parent
									//	find_parent(node_index,ItsParent1);
										first_parent_index(shift,ItsParent);
										first_parent_index(node_index,ItsParent1);
										if(ItsParent == ItsParent1&& ItsParent != nodes.size()&& ItsParent2 != node_index && ItsParent2 != nodes.size()){
											std::cout<<"they have same parent"<<std::endl;
											delete_relation(ItsParent2,shift);
											nodes_relation.insert(make_pair(ItsParent2,shift+2));
										}
									}				
								}
								for(size_t j = 0 ; j < childs.size();j++){
									delete_relation(node_index,childs.at(j));
									nodes_relation.insert(make_pair(node_index+2, childs.at(j)+2));
								}
								std::cout<<"nodes relation4: "<<std::endl;
								for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
									std::cout << it->first <<" "<< it->second << std::endl;
								}
								nodes.at(node_index) = common_part;
								nodes_relation.insert(make_pair(node_index,node_index+1));
								nodes_relation.insert(make_pair(node_index,node_index+2));
							}
							std::map<std::vector<size_t>, size_t> f_par;
							for(size_t shift = node_index+1; shift <nodes.size()-2; shift++){
								std::map<std::vector<size_t>, size_t>::iterator it = firstParent.find(nodes.at(shift+2));
								if(it != firstParent.end() && it->second == shift){
									f_par.insert(make_pair(it->first,it->second+2));
								}
							}
							for(std::map<std::vector<size_t> , size_t>::iterator it = f_par.begin(); it != f_par.end(); it++){
								std::map<std::vector<size_t> , size_t>::iterator it1 = firstParent.find(it->first);
								if(it1 != firstParent.end()){
									it1->second = it->second;
								}
							}
							std::map<std::vector<size_t>,size_t>::iterator it = firstParent.find(nodes.at(node_index));
							if(it != firstParent.end() && it->second == node_index){
								std::cout << "node is first parent!"  << node_index <<std::endl;
								firstParent.erase(it);
								firstParent.insert(make_pair(common_part,node_index));
							}
							std::cout<<"nodes relation: "<<std::endl;
							for(std::multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
								std::cout << it->first <<" "<< it->second << std::endl;
							}
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
					for(size_t j = 0; j < sub_branch.size();j++){
						std::cout << sub_branch.at(j)<< " ";
					}					
					it1->second = it1->second +1;
					std::cout<< "number of happening "<< it1->second << std::endl;
				}
			}			
			}

		}
		cout<< "branch counter: "<<endl;
		for(std::map<std::vector<size_t> , size_t>::iterator it = branch_counter.begin(); it != branch_counter.end(); it++){
			vector<size_t> br = it->first;
			for(size_t i =0; i < br.size();i++){
				cout<< br.at(i)<< " ";
			}
			cout<< " number "<<it->second;
			cout << " " <<endl;
			
		}
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
		std::cout<<"first_parent: "<< std::endl;
		for(size_t i =0 ; i < first_parent.size(); i++){
			std::cout<< first_parent.at(i)<< " " ;
		}
		std::cout << " " <<std::endl;
		return first_parent;
	}
	size_t suffix_tree::get_power_of_two(size_t & power)const{
		size_t power_of_two = powerOfTwo.at(power);
		return power_of_two;
	}
	std::map<size_t, std::vector<size_t> > suffix_tree::get_center_on_a_sequence(size_t & seq_id)const{
		return successive_centers.at(seq_id);
	}
	merging_centers::merging_centers(all_data & d, finding_centers & cent , suffix_tree & t):data(d), centers(cent),tree(t){}
	merging_centers::~merging_centers(){}
	void merging_centers::updating_centers(std::vector<size_t> & center_string, size_t & index){//replace 'center_string' with its new 'index' where ever 'center_sting' occurs on the tree and updates the number of happening.
		//First make a new tree !
		tree.create_tree(center_string,index);
		std::vector<std::vector<size_t> > nodes = tree.get_nodes();
		std::map<std::vector<size_t>, size_t > counts = tree.get_count(); // Updated number of happening for each updated string of centers(Notice that it has node number not the centers indices)
		std::map<std::vector<size_t>, int> intermediate;
		for(std::map<std::vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
			std::cout << "br size is : " << br.size() <<std::endl;
			for(size_t i =0; i < br.size();i++){
				std::cout << br.at(i)<< " ";
			}
			std::cout << " "<< std::endl;
			size_t number = it->second;
			std::vector<size_t> seriesOfCenters;
			for(size_t j = 0 ; j < br.size(); j++){
				for(size_t k =0; k < nodes.at(br.at(j)).size();k++){
					seriesOfCenters.push_back(nodes.at(br.at(j)).at(k));
				}
			}
			size_t power = 31 ;
			if(seriesOfCenters.at(seriesOfCenters.size()-1) == tree.get_power_of_two(power)){
				seriesOfCenters.pop_back();	
			}
			std::map<std::vector<size_t>, int>::iterator it1 = gains.find(seriesOfCenters);
			if(it1 != gains.end()){//It is kept as it is
				intermediate.insert(make_pair(it1->first,it1->second));
			}else{
				size_t number_of_new_center = 0;
				for(size_t i = 0; i < seriesOfCenters.size(); i++){
					if(seriesOfCenters.at(i)== index){
						number_of_new_center = number_of_new_center + 1;
					}
				}
				size_t old_length = seriesOfCenters.size()+ (number_of_new_center*center_string.size());
				size_t gain = number*old_length - (number + seriesOfCenters.size());
				intermediate.insert(make_pair(seriesOfCenters,gain));
			}
		}
		gains.clear();
		for(std::map<std::vector<size_t>,int>::iterator it = intermediate.begin(); it != intermediate.end(); it++){
			gains.insert(make_pair(it->first,it->second));
		}
/*			size_t oldCount = old_count ->second;
			size_t length = centers.size();
			if(center_string.size()<centers.size()){//here center_string is the high path
				length = center_string.size();
			}
			bool StillNeedTobeChecked = true;
			while(StillNeedTobeChecked == true){		
				for(size_t i =0; i < centers.size(); i++){
					size_t first_common_index;
					std::vector<size_t> common_part;
					if(centers.at(i)==center_string.at(0)&& centers.size()-i >= center_string.size()){
						first_common_index = i;
						for(size_t j =0; j < length;j++){
							if(centers.at(j+i)==center_string.at(j)){
								common_part.push_back(centers.at(j));
							}else{
								break;
							}
						}
					}
					if(common_part.size() == center_string.size()){
						std::cout<< "there is a common part "<<std::endl;
					//	std::map<std::vector<size_t> , size_t>::iterator it1 = updated_counts.find(centers);
					//	if(old_count != updated_counts.end() && old_count ->second > 0){//why don't we erase it??
					//		updated_counts.erase(it1);
						//	it1->second = 0;
					//	}
						StillNeedTobeChecked = true;
					//	std::map<std::vector<size_t> , size_t>::iterator oldCount = counts.find(centers);						
						std::vector<size_t> new_center;
						for(size_t j = 0; j < first_common_index ;j++){
							new_center.push_back(centers.at(j));
						}
						new_center.push_back(index);
						for(size_t j = first_common_index + common_part.size(); j < centers.size();j++){
							new_center.push_back(centers.at(j));
						}
				//		centers = new_center;
						std::map<std::vector<size_t> , size_t>::iterator it1 = updated_counts.find(centers);
						if(it1 != updated_counts.end()){
							oldCount = it1->second;
							std::cout<<"old count " << oldCount << std::endl;
							updated_counts.insert(make_pair(new_center,oldCount));
						}else{
							std::cout << "there is something wrong in new counts!" <<std::endl;
						}
						updated_counts.erase(it1);
						centers = new_center;
					}else{
					//	std::cout << "we are here! " <<std::endl;
						StillNeedTobeChecked = false;
					}
				}
			}
			int gain = it->second;
			if(centers.size() != 1 && old_center_size != centers.size()){
				std::cout << oldCount << " " << old_center_size << " " << centers.size()<<std::endl;
				gain = oldCount*(old_center_size)-(oldCount+centers.size());
			}
			intermediate.insert(make_pair(centers,gain));
			std::cout<< "centers " ;
			for(size_t j =0; j < centers.size();j++){
				cout<< centers.at(j) << " " ;
			}
			cout << " " <<endl;
			std::cout << "gain is: " << gain << std::endl;
		}
		for(std::map<std::vector<size_t>, size_t>::iterator new_count = updated_counts.begin(); new_count != updated_counts.end(); new_count ++){
			for(size_t j =0; j < new_count->first.size();j++){
				cout<< new_count->first.at(j) << " " ;
			}
			std::cout << "number of happening " << new_count->second << std::endl;

		}
		size_t gains_size = gains.size();
		std::cout << "updated_count size is " << updated_counts.size() << std::endl;
		std::cout << "gains size " << gains_size <<std::endl;//should be equal to the number of paths
		gains.clear();
		for(std::map<std::vector<size_t>,int>::iterator it = intermediate.begin(); it != intermediate.end(); it++){
			gains.insert(make_pair(it->first,it->second));
		}*/

	}

	void merging_centers::merg_gain_value(){
		std::vector<std::vector<size_t> > nodes = tree.get_nodes();
		cout<< "size: " << nodes.size()<<endl;
		size_t original_center_numbers = centers.get_number_of_centers();
		std::cout << "original center number: "<< original_center_numbers << std::endl;
		std::map<std::vector<size_t> , size_t> counts = tree.get_count();//vector<size_t> shows a path and size_t is its number of happening.
		//calculating the initial gain values:
		for(map<vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
			vector<size_t> br = it->first;//list of nodes of a path
			size_t number = it->second;
			std::vector<size_t> centers;
			int gain = 0;
			for(size_t j = 0 ; j < br.size(); j++){
				for(size_t k =0; k < nodes.at(br.at(j)).size();k++){
					centers.push_back(nodes.at(br.at(j)).at(k));
				}
			}
			size_t power = 31 ;
			if(centers.at(centers.size()-1) == tree.get_power_of_two(power)){
				centers.pop_back();	
			}
			std::cout<< "centers " ;
			for(size_t j =0; j < centers.size();j++){
				cout<< centers.at(j) << " " ;
			}
			cout << " " <<endl;
			std::cout<< "center size: " << centers.size() << "  number of happening: "<< number << std::endl;
			gain = number*(centers.size())-(number+centers.size());
			if(centers.size()!= 0){
				gains.insert(make_pair(centers,gain));//centers-->index of centers
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
		center_numbers = original_center_numbers + 1;//it will be used when ever we are going to insert the next megerd center to the merged_centers map.
		std::cout << "cent number: "<< center_numbers << std::endl;
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
			hi_from_map.push_back(hi->second);//Because we dont want to use the same canter again as the one with the highest gain
			assert(hi != merged_centers.end());
			for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
				if(it->second > highest_gain && it->first != hi_from_map && it->first.size() != 1){
					std::cout<< "gain "<< it->second <<std::endl;
					highest_gain = it->second;
					highest_path = it->first;
				}else continue;
			}
			if(gains.size() != 0){
				std::cout << "check highest path: " <<std::endl;
				for(size_t j = 0; j < highest_path.size(); j++){
					std::cout << highest_path.at(j)<< " ";
				}
				std::cout << " " <<std::endl;
			}
			center_numbers = center_numbers +1;
			merged_centers.insert(make_pair(highest_path, center_numbers));
		}
		std::cout << "final result: " << std::endl;
		for(map<std::vector<size_t>, int>::iterator it = gains.begin(); it != gains.end(); it++){
			for(size_t j =0; j < it->first.size(); j++){
				std::cout<< it->first.at(j)<< " ";
			}
			std::cout<< " gain is " << it->second << std::endl;
		}
		std::cout << "new centers: "<<std::endl;
		for(std::map<std::vector<size_t>, size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){
			if(it != merged_centers.end()){
				for(size_t i = 0; i < it->first.size(); i ++){
					std::cout << it->first.at(i)<< " ";
				}
				std::cout << " its index is " << it->second <<std::endl;
			}else {std::cout << "there is no merged center! " <<std::endl;}
			
		}		
	}
	void merging_centers::adding_new_centers(std::vector<std::vector<std::string> > & long_centers, std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){//Filling in the long centers vector and centersPositionOnASeq
		merg_gain_value();
		size_t biggest_index = 0;
		std::vector<std::string> sequence_of_centers;
		for(std::map<std::vector<size_t>,size_t>::iterator it = merged_centers.begin(); it != merged_centers.end(); it++){
			std::cout << " merged_center: " << it->second << std::endl;
			std::vector<size_t> updated_center = it->first;
			std::vector<size_t> list;
			list = it->first;
			size_t id = list.at(0);	
			bool ThereIsStillABigID = true;
			while (ThereIsStillABigID == true){
				std::cout << "here!" << std::endl;
				std::cout << "number of original centers: "<< centers.get_number_of_centers()<<std::endl;
				for(size_t j =0; j < list.size(); j++){
					id = list.at(j);
					std::cout << "id " << id <<std::endl;
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
			for(size_t j =0; j < list.size(); j++){
				size_t id = list.at(j);
				std::string center = centers.find_center_name(id);
				sequence_of_centers.push_back(center);
			}
			long_centers.push_back(sequence_of_centers);
			for(size_t i = 0; i < data.numSequences(); i++){
				find_new_centers(it->second,sequence_of_centers,i, centersPositionOnASeq);
			}
			sequence_of_centers.clear();
		}
	}
/*	void merging_centers::merg_alignments(vector<vector<std::string> > & long_centers, std::map<std::string, std::vector<pw_alignment> > & al_of_a_ccs, std::map<std::string, std::vector<std::string> > & cluster_result, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers){//new_centers has new centers and their als, it should contains all the centers of all rounds.//TODO has to be fixed!
		size_t counter = 0;
		for(size_t i =0; i < long_centers.size(); i ++){
			std::cout << " i "<< i <<std::endl;
			vector<std::string> long_center = long_centers.at(i);
			std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(long_center);
			if(new_cent == new_centers.end()){
				new_centers.insert(make_pair(long_center,vector<pw_alignment>()));
			}
			std::map<std::vector<std::string>, std::vector<pw_alignment> > temp;
			size_t temp_center_ref;
			size_t temp_center_left;
			for(size_t j = 0; j < long_center.size(); j++){
				std::string center = long_center.at(j);
				std::cout << "center: "<< center << std::endl;
				std::map<std::string, std::vector<std::string> >::iterator it = cluster_result.find(center);
				std::map<std::string, std::vector<pw_alignment> >::iterator all_als = al_of_a_ccs.find(center);
				assert(it != al_of_a_ccs.end());
				std::vector<std::string> center_parts;
				strsep(center, ":" , center_parts);
				unsigned int center_ref = atoi(center_parts.at(0).c_str());
				unsigned int center_left = atoi(center_parts.at(1).c_str());
				std::vector<std::string> current_string_of_centers;
				pw_alignment al;
				if(j == 0){
					std::vector<std::string> first_center;
					first_center.push_back(center);
					temp.insert(make_pair(first_center, all_als->second));
					temp_center_ref = center_ref;
					temp_center_left= center_left;	
					std::cout << "if j is 0 " << std::endl;				
				}else{
					for(size_t m = 0; m < j ; m++){
						current_string_of_centers.push_back(long_center.at(m));
						std::cout << long_center.at(m) <<std::endl;
					}
					std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator temporary = temp.find(current_string_of_centers);
					std::cout << " size of temp->second " << temporary->second.size() <<std::endl;
					if(temporary != temp.end() && temporary->second.size() != 0){
						std::cout << "temporary != temp.end()"<<std::endl;
						std::vector<std::string> members = it->second;
						for(size_t k = 0; k < members.size();k++){
						std::vector<std::string>member_parts;
						strsep(members.at(k),":", member_parts);
						unsigned int member_left = atoi(member_parts.at(0).c_str());
						unsigned int member_ref = atoi(member_parts.at(0).c_str());
					//	for(size_t k = 0; k < als.size(); k++){
					//		pw_alignment & p = als.at(k);
					//		size_t left1,right1,left2,right2;
					//		p.get_lr1(left1, right1);
					//		p.get_lr2(left2, right2);
					//		size_t ref1 = p.getreference1();
					//		size_t ref2 = p.getreference2();
					//		size_t left, right, ref,id;
					//		if(ref1 != center_ref && left1 != center_left){//center is on the ref2 of alignment
								unsigned int closest_one = 0;
								bool closest_left_is_found = false;
								pw_alignment & p1 = temporary->second.at(0);
								for(size_t n =0; n < temporary->second.size(); n++){
									p1 = temporary->second.at(n);
									size_t left_1,right_1,left_2,right_2;
									p1.get_lr1(left_1, right_1);
									p1.get_lr2(left_2, right_2);
									size_t ref_1 = p1.getreference1();
									size_t ref_2 = p1.getreference2();
									if(ref_1 == member_ref && member_left > right_1){//find the closest!


									}
									if(ref_2 == member_ref && member_left > right_2){
									if(ref_1 != temp_center_ref && left_1 != temp_center_left){//center from temp is on ref2 of its al.
										left = left_1;
										right = right_1;
										ref = ref_1;
										id = 0;
									}else{
										left = left_2;
										right = right_2;
										ref = ref_2;
										id = 1;
									}									
									if(ref == ref1 && left < left1){
										if(left >= closest_one){
											closest_one = left;
											closest_left_is_found = true;
										}
									}
								}
								if(closest_left_is_found == true){ //translate and translate back bits
									std::vector<bool> sample_p1 = p1.getsample(id);
									std::vector<bool> sample_temp_center = p1.getsample(1-id);
									std::vector<bool> sample_p = p.getsample(0);
									std::vector<bool> sample_cent = p.getsample(1);
									std::vector<bool> sample1;
									std::vector<bool> sample2;
									for(size_t m = 0 ; m <sample_p1.size(); m++){
										sample1.push_back(sample_p1.at(m));
										sample2.push_back(sample_temp_center.at(m));
									}
									vector<bool> middle_part_of_sample;
									std::cout << "distance between two als " << left1-(right+1) << std::endl;
									for(size_t m = right+1 ; m < left1 ; m++){
										char base = data.getSequence(ref1).at(m);
										vector<bool> bits(3);
										pw_alignment::get_bits(base,bits);
										for(size_t n = 0; n < 3; n++){
											middle_part_of_sample.push_back(bits.at(n));
										}
									}
									std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
									for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
										sample1.push_back(middle_part_of_sample.at(m));
										sample2.push_back(middle_part_of_sample.at(m));
									}
									for(size_t m = 0 ; m < sample_p.size(); m++){
										sample1.push_back(sample_p.at(m));
										sample2.push_back(sample_cent.at(m));
									}
									al.setreference1(ref1);
									al.setbegin1(closest_one);//to solve the problem of reverse or forward just left and right are considered.
									al.setend1(right1);
									if(temp_center_ref == center_ref){
										al.setreference2(ref2);
										al.setbegin2(temp_center_left);
										al.setend2(right2);
									}else{
										if(temp_center_ref > data.numSequences()){
											al.setreference2(temp_center_ref);
											std::cout<< "temp center ref " << temp_center_ref << std::endl;
										}else{
											counter = counter + 1;
											al.setreference2(data.numSequences()+counter);
										}
										al.setbegin2(0);
										al.setend2(p1.alignment_length()+(left1-(right+1))+p.alignment_length());
									}
									al.set_alignment_bits(sample1,sample2);//set bits
									std::cout << "center is on ref 2 " <<std::endl;
									al.print();
									temp_center_ref = al.getreference2();
									size_t Left_2,Right_2;
									al.get_lr2(Left_2, Right_2);
									temp_center_left = Left_2;
								}
							}else{//if center is on ref1 of alignment
								unsigned int closest_one = 0;
								bool closest_left_is_found = false;
								pw_alignment & p1 = temporary->second.at(0);
								for(size_t n =0; n < temporary->second.size(); n++){// if it goes over temp instead it may work
									p1 = temporary->second.at(n);
									size_t left_1,right_1,left_2,right_2;
									p1.get_lr1(left_1, right_1);
									p1.get_lr2(left_2, right_2);
									size_t ref_1 = p1.getreference1();
									size_t ref_2 = p1.getreference2();
									if(ref_1 != temp_center_ref && left_1 != temp_center_left){//center from temp is on ref2 of its al.
										left = left_1;
										right = right_1;
										ref = ref_1;
										id = 0;
									}else{
										left = left_2;
										right = right_2;
										ref = ref_2;
										id = 1;
									}									
									if(ref == ref2 && left < left2){
										if(left >= closest_one){
											closest_one = left;
											closest_left_is_found = true;
										}
									}
								}
								if(closest_left_is_found == true){ //translate and translate back bits
									std::vector<bool> sample_p1 = p1.getsample(id);
									std::vector<bool> sample_temp_center = p1.getsample(1-id);
									std::vector<bool> sample_p = p.getsample(1);
									std::vector<bool> sample_cent = p.getsample(0);
									std::vector<bool> sample1;
									std::vector<bool> sample2;
									for(size_t m = 0 ; m <sample_p1.size(); m++){
										sample1.push_back(sample_temp_center.at(m));
										sample2.push_back(sample_p1.at(m));
									}
									vector<bool> middle_part_of_sample;
									for(size_t m = right+1 ; m < left2 ; m++){
										char base = data.getSequence(ref2).at(m);
										vector<bool> bits(3);
										pw_alignment::get_bits(base,bits);
										for(size_t n = 0; n < 3; n++){
											middle_part_of_sample.push_back(bits.at(n));
										}
									}
									for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
										sample1.push_back(middle_part_of_sample.at(m));
										sample2.push_back(middle_part_of_sample.at(m));
									}
									for(size_t m = 0 ; m < sample_p.size(); m++){
										sample1.push_back(sample_cent.at(m));
										sample2.push_back(sample_p.at(m));
									}
									if(temp_center_ref == center_ref){
										al.setreference1(ref1);
										al.setbegin1(temp_center_left);
										al.setend1(right1);
									}else{
										if(temp_center_ref > data.numSequences()){
											al.setreference1(temp_center_ref);
										}else{
											counter = counter + 1;
											al.setreference1(data.numSequences()+counter);
										}
										al.setbegin1(0);
										al.setend2(p1.alignment_length()+(left2-(right+1))+p.alignment_length());
									}
									al.setreference2(ref2);
									al.setbegin2(closest_one);//to solve the problem of reverse or forward just left and right are considered.
									al.setend2(right2);
									al.set_alignment_bits(sample1,sample2);//set bits
									std::cout << "center is on ref 1 " <<std::endl;
									al.print();
									temp_center_ref = al.getreference1();
									size_t Left_1,Right_1;
									al.get_lr1(Left_1, Right_1);
									temp_center_left = Left_1;
								}
							}
						}
					}
				}//else	(j != 0)
				//Insert the concatenated center and its members to the temp map
				current_string_of_centers.push_back(long_center.at(j));
				std::cout << "size of current string of centers is " << current_string_of_centers.size() <<std::endl;
				for(size_t m = 0; m < current_string_of_centers.size(); m++){
					std::cout << current_string_of_centers.at(m) << " ";
				}
				std::cout << " " <<endl;
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it_temp = temp.find(current_string_of_centers);
				if(it_temp == temp.end()){
					temp.insert(make_pair(current_string_of_centers,std::vector<pw_alignment>()));
					it_temp = temp.find(current_string_of_centers);
				}
				if(al.alignment_length() != 0){
					it_temp->second.push_back(al);
					std::cout << "al: " <<std::endl;
					al.print();
				}
			}
			std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator temporary = temp.find(long_center);
			if(temporary != temp.end()){
				std::map<std::vector<std::string>, std::vector<pw_alignment> >::iterator it1 = new_centers.find(long_center);
				if(it1 != new_centers.end()){
					std::cout<< "add to new center! "<<std::endl;
					it1->second = temporary ->second;
					std::cout << "number of added als: "<< it1->second.size() << std::endl;
				}
			}
		}
	}*/
	void merging_centers::create_alignment(std::vector<std::vector<std::string> > & long_centers, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::vector<std::string> > > & centersPositionOnASeq, std::vector<std::map<size_t , std::string> > & centerOnSequence){
		size_t artificial_ref = data.numSequences()+new_centers.size();
		for(size_t i =0; i < long_centers.size(); i ++){
			vector<std::string> long_center = long_centers.at(i);
			std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(long_center);
			if(new_cent == new_centers.end()){
				new_centers.insert(make_pair(long_center,std::vector<pw_alignment>()));
			}
		}
		for(size_t j = 0 ; j < data.numSequences(); j++){
			for(std::map<size_t , std::vector<std::string> >::iterator it = centersPositionOnASeq.at(j).begin(); it != centersPositionOnASeq.at(j).end(); it++){
				if(it != centersPositionOnASeq.at(j).end()){//If there is a long center on that sequence
					std::vector<bool> sample1;
					std::vector<bool> sample2;
					pw_alignment al;
					al.setreference1(j);
					al.setbegin1(it->first);
					size_t end_of_last_piece=0;
					size_t length;
					size_t position = it->first;
					std::cout << "position " << position << std::endl;
					find_long_center_length(it->second, alignments_in_a_cluster,j,position,length,end_of_last_piece);
					std::cout << "end of last piece "<< end_of_last_piece << " length "<< length << std::endl;
					al.setend1(end_of_last_piece);
					//Up to now first reference of an pw-alignment is set. Center is always considered as the second reference.
					bool centers_are_on_the_same_ref = false;
					std::string first_center = it->second.at(0);
					std::vector<std::string> cent_parts;
					strsep(first_center, ":" , cent_parts);
					unsigned int cent_ref = atoi(cent_parts.at(0).c_str());
					unsigned int cent_left = atoi(cent_parts.at(1).c_str());
					for(size_t i = 0; i < it->second.size();i++){
						std::string center = it->second.at(i);
						std::vector<std::string> center_parts;
						strsep(center, ":" , center_parts);
						unsigned int center_ref = atoi(center_parts.at(0).c_str());
						unsigned int center_left = atoi(center_parts.at(1).c_str());
						if(center_ref == cent_ref){
							centers_are_on_the_same_ref= true;
						}else{
							centers_are_on_the_same_ref=false;
							break;
						}
					}
					if(centers_are_on_the_same_ref==true){
						al.setreference2(cent_ref);
						al.setbegin2(cent_left);
						al.setend2(cent_left + length-1);
					}else{
						std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(it->second);
						if(new_cent->second.size() == 0){
							al.setreference2(artificial_ref);
							artificial_ref = artificial_ref+1;
						}else{
							pw_alignment p1 = new_cent->second.at(0);
							al.setreference2(p1.getreference2());
							assert(p1.getreference2() >=data.numSequences());
						}
						al.setbegin2(0);
						al.setend2(length);
					}
					std::cout<< "begin2 "<< al.getbegin2() << " end2 " << al.getend2() << std::endl;
					//Now second reference is also set.
					//First we should check if the alignment is not there yet
					bool dontAdd = false;
					std::map<vector<std::string>, vector<pw_alignment> >::iterator new_cent = new_centers.find(it->second);
					for(size_t k =0; k < new_cent->second.size();k++){
						pw_alignment p1 = new_cent->second.at(k);
						if(p1.getbegin1()==al.getbegin1()&&p1.getend1()==al.getend1()&&p1.getbegin2()==al.getbegin2()&&p1.getend2()==al.getend2()&&p1.getreference1()==al.getreference1()&&p1.getreference2()==al.getreference2()){
							dontAdd = true;
							break;
						}
					}
					if(dontAdd == false){
						//We are setting samples here:
						size_t right = 0;
						std::cout<< "it second size " << it->second.size() <<std::endl;
						for(size_t i =0; i < it->second.size();i++){
							std::cout<< " center "<< it->second.at(i)<<std::endl;
							std::map<std::string , std::vector<pw_alignment> >::iterator it3 = alignments_in_a_cluster.find(it->second.at(i));
							size_t left_of_a_sample;
							for(std::map<size_t , std::string>::iterator leftOfSample = centerOnSequence.at(j).begin(); leftOfSample !=  centerOnSequence.at(j).end(); leftOfSample++){
							//	std::cout<< "current from the map " << leftOfSample->second << std::endl;
								if (it->second.at(i)== leftOfSample->second && leftOfSample->first >= it->first && leftOfSample->first >= right){
									left_of_a_sample = leftOfSample->first;
									std::cout<< "left is set!"<<std::endl;
									break;
								}
							}
							assert(it3 != alignments_in_a_cluster.end());
							for(size_t k = 0; k < it3->second.size();k++){
								pw_alignment p = it3->second.at(k);
								size_t r1,r2,l1,l2;
								p.get_lr1(l1,r1);
								p.get_lr2(l2,r2);
								std::cout<<  "ref 1 " << p.getreference1() << " ref2 " << p.getreference2()<< " j " << j << " l1 " << l1 << " l2 " << l2 << " left_of_a_sample "<< left_of_a_sample << std::endl;
								if((p.getreference1()== j && l1 == left_of_a_sample)||(p.getreference2() == j && l2 == left_of_a_sample)){
									if(l1 == left_of_a_sample){
										right = r1;
										std::cout << "right1 " << right << std::endl;
									}else{
										right = r2;
										std::cout << "right2 " << right << std::endl;
									}
									std::vector<bool> sample1_p = p.getsample1();
									std::vector<bool> sample2_p = p.getsample2();
									std::cout << " sample1_p.size() " << sample1_p.size()<< " " << sample2_p.size()<< " " << p.alignment_length() <<std::endl;
									std::cout << "pushing back samples " << std::endl;
									for(size_t m =0; m < sample1_p.size(); m++){
										sample1.push_back(sample1_p.at(m));
										sample2.push_back(sample2_p.at(m));
									}
									std::cout << sample1.size() << " " << sample2.size() << std::endl;
									break;							
								}
							}
							if(i != it->second.size()-1){
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
								for(size_t m = right+1 ; m < left_of_next_center ; m++){
									char base = data.getSequence(j).at(m);
									char gap = '-';
									vector<bool> bits(3);
									vector<bool> gap_bit(3);
									pw_alignment::get_bits(base,bits);
									pw_alignment::get_bits(gap,gap_bit);
									for(size_t n = 0; n < 3; n++){
										middle_part_of_sample.push_back(bits.at(n));
										gap_sample.push_back(gap_bit.at(n));
									}
								}
								std::cout<< "size of middle part is " << middle_part_of_sample.size() << std::endl;
								for(size_t m = 0 ; m < middle_part_of_sample.size();m++){
									sample1.push_back(middle_part_of_sample.at(m));
									sample2.push_back(gap_sample.at(m));
								}
							}
						}
						al.set_alignment_bits(sample1,sample2);
						std::cout << "alignment is: " <<std::endl;
						al.print();
						new_cent->second.push_back(al);
					}
				}
			}
		}
	}
//The following function finds long centers of specific sequence:
	void merging_centers::find_new_centers(size_t & center_indices, std::vector<std::string > & current_long_center, size_t & seq_id , std::vector<std::map<size_t, std::vector<std::string> > > & centersPositionOnASeq){
		std::map<size_t,std::vector<size_t> > all_centers = tree.get_center_on_a_sequence(seq_id);
		std::cout<<"sequence is "<< seq_id << std::endl;
		for(std::map<size_t,std::vector<size_t> >::iterator it1 = all_centers.begin(); it1!= all_centers.end();it1++){
				std::cout << center_indices<< std::endl;
				std::cout << " " <<std::endl;
				for(size_t i =0; i < it1->second.size(); i++){
					std::cout << it1->second.at(i)<< " ";
				}
				std::cout << " " <<std::endl;
				if(it1->second.at(0) == center_indices){
					centersPositionOnASeq.at(seq_id).insert(make_pair(it1->first,current_long_center));
					std::cout << "seq "<<seq_id << "position "<<it1->first << std::endl;
				}
		}
	}
	void merging_centers::find_long_center_length(std::vector<std::string> & centers, std::map<std::string,std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & seq_id,size_t & position ,size_t & center_length, size_t & end_of_last_piece){
		end_of_last_piece = 0;
		size_t current_position  = position;
		center_length = 0;
		for(size_t i =0; i < centers.size();i++){
			std::string center = centers.at(i);
			std::vector<std::string> center_parts;
			strsep(center, ":" , center_parts);
			unsigned int center_ref = atoi(center_parts.at(0).c_str());
			unsigned int center_left = atoi(center_parts.at(1).c_str());
			std::map<std::string , std::vector<pw_alignment> >::iterator it = alignments_in_a_cluster.find(centers.at(i));
			assert(it != alignments_in_a_cluster.end());
			std::cout<< "al size "<<it->second.size() << std::endl;
			bool selfAligned = true;
			for(size_t j =0; j < it->second.size();j++){
				pw_alignment p = it->second.at(j);
				size_t r1,r2,l1,l2;
				p.get_lr1(l1,r1);
				p.get_lr2(l2,r2);
				size_t ref1 = p.getreference1();
				size_t ref2 = p.getreference2();
				std::cout << "l1 " << l1 << " r1 "<< r1 << " l2 "<< l2 << " r2 " << r2 << " ref1 "<<ref1 << " ref2 "<< ref2<<std::endl;
				std::cout << "seq id " << seq_id << " cent ref "<< center_ref << " current pos " << current_position << " cent lef "<< center_left << std::endl;
				//If center is not the seq_id
				if((ref1== seq_id && ref2== center_ref && l1 >= current_position && l1 <= current_position + 5 && l2 ==center_left)){ //5 is the number of gaps!
					std::cout<< "on ref 1"<<std::endl;
					center_length += r2-l2+1;
					current_position = r1;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
					}
					break;
				}
				else if(ref2== seq_id && ref1== center_ref && l2 >= current_position && l2<= current_position +5 && l1 ==center_left){
					center_length +=r1-l1+1;
					current_position = r2;
					std::cout << "cur pos " << current_position << " length "<< center_length << std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
					}
					break;
				}
				//If we need to align the same piece against itself
				else if(ref1 == seq_id && seq_id == center_ref && center_left == l1 && center_left >= current_position && center_left <= current_position + 5){
					center_length +=r1-l1+1;
					current_position = r1;
					std::cout << "here!" <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r1;
					}
					break;
				}
				else if(ref2 == seq_id && seq_id == center_ref && center_left == l2 && center_left >= current_position && center_left <= current_position + 5){
					center_length +=r2-l2+1;
					current_position = r2;
					std::cout << "there!" <<std::endl;
					if(i == centers.size()-1){
						end_of_last_piece = r2;
					}
					break;

				}


			}
		}
	}


	

#endif
