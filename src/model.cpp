#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP


// TODO why does ias need an output stream (the gain function should also not need one)
template<typename T>
void initial_alignment_set<T>::compute(overlap & o, ofstream & outs) {

	compute_simple_lazy_splits(o,outs);
//	compute_simple(o,outs);

}

/*

insert al recursively
each insert operation causes the removal of some alignments from overlap
then we split al and the removed alignments (without induced splits)
the gain of all new pieces with positive gain - the gain of all removed pieces is an upper bound of 
the information gain of inserting al. If this upper bound is positive we continue and try to recursively insert all pieces
Afterwards, we have to compute the actual gain. If it is negative we undo the local recursive insert operation
TODO remove outs

*/
template<typename T>
void initial_alignment_set<T>::lazy_split_insert_step(overlap & ovrlp, size_t level, const pw_alignment * al, vector<const pw_alignment*> & inserted_alignments, vector<const pw_alignment *> & removed_alignments, double & local_gain, ofstream & outs) {
// start with computing gain of al
	double gain1;
	double gain2;
	common_model.gain_function(*(al), gain1, gain2,outs);
	double av_al_gain = (gain1 + gain2) / 2 - base_cost;
	// we continue if information gain is possible from the current alignment
	local_gain = 0;
	if(av_al_gain > 0) {
		splitpoints spl(*al, ovrlp, data);
		spl.nonrecursive_splits();
		// sets of alignments that need to be removed and inserted if we want the current alignment (necessary for undo operation after recursion)
		alset remove_als; // remove alignments are pointers to objects contained in the overlap structure
		vector<pw_alignment> insert_als;
		spl.split_all(remove_als, insert_als);

		cout <<"initial: level " << level << " positive gain: " << av_al_gain << " split res rem " << remove_als.size() << " ins " << insert_als.size() << endl;
	
		// we evaluate an upper bound of information gain for the current alignment
		double lost_information_gain = 0;
		for(alset::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			const pw_alignment * ral = *it;
			common_model.gain_function((*ral), gain1, gain2, outs);
			double rav_gain = (gain1 + gain2)/2 - base_cost;
			lost_information_gain += rav_gain;

//			cout << "REM " << endl;
//		       	ral->print();
//			cout << endl;	
		
		}

		multimap<double, pw_alignment> ordered_parts;
		double max_parts_gain = 0;
		for(size_t i=0; i<insert_als.size(); ++i) {
			common_model.gain_function(insert_als.at(i), gain1, gain2, outs);
			double iav_gain = (gain1 + gain2)/2 - base_cost;
			if(iav_gain > 0) {
				max_parts_gain += iav_gain;
				ordered_parts.insert(make_pair(iav_gain, insert_als.at(i)));
			}

//			cout << "INS " << endl;
//			insert_als.at(i).print();
//			cout << endl;
		}

		cout << "splits: level " << level << " on al length " << al->alignment_length() << " gain " << av_al_gain << " causes to remove " << remove_als.size() << " alignments with " << lost_information_gain << 
			" gain. We could insert " << ordered_parts.size() << " small pieces. Upper bound of gain is " << max_parts_gain << endl;

		// here: we actually decide to insert a small part (it does not induce any more splits)
		// TODO rem and ins gets changed, only if it stays good we do insert and return,
		// otherwise go on with recursion
		if(remove_als.empty() && ordered_parts.size() == 1) {
			cout << " PART INSERTED????" << endl;
			pw_alignment alin = insert_als.at(0);
			insert_als.clear();
			splitpoints spl2(alin, ovrlp, data);
			spl2.recursive_splits();
			// TODO check does split all add to remove_als
			spl2.split_all(remove_als, insert_als);
			size_t inserted_counter = 0;
			if(remove_als.empty()) {
				for(size_t i=0; i<insert_als.size(); ++i) {
					double g1;
					double g2;
					common_model.gain_function(insert_als.at(i), g1, g2, outs);
					double avg = (g1 + g2) / 2.0 - base_cost;
					if(avg > 0) {
						pw_alignment * ins =  ovrlp.insert_without_partial_overlap(insert_als.at(i));
						inserted_alignments.push_back(ins);
						local_gain+=avg;
						inserted_counter++;
					}
				}
			



				cout << "insert: " << level << " on al length " << al->alignment_length() << " gain " << av_al_gain << " finally inserted pieces: " << inserted_counter << " real local gain " << local_gain << endl;
				return; 
			// insert pieces for recursion:
			} else {
				ordered_parts.clear();

				for(alset::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
					const pw_alignment * ral = *it;
					common_model.gain_function((*ral), gain1, gain2, outs);
					double rav_gain = (gain1 + gain2)/2 - base_cost;
					lost_information_gain += rav_gain;

		//			cout << "REM2 " << endl;
		//	       		ral->print();
		//			cout << endl;	
		
				}

				max_parts_gain = 0;
				for(size_t i=0; i<insert_als.size(); ++i) {
					common_model.gain_function(insert_als.at(i), gain1, gain2, outs);
					double iav_gain = (gain1 + gain2)/2 - base_cost;
					if(iav_gain > 0) {
						ordered_parts.insert(make_pair(iav_gain, insert_als.at(i)));
						max_parts_gain+=iav_gain;
					}
		//			cout << "1 level " << level << " on al length " << al->alignment_length() << " 
		//			cout << "INS for later" << endl;
		//			insert_als.at(i).print();
		//			cout << endl;


				}

			
				cout << "2nd split: " << level << " on al length " << al->alignment_length() << " gain " << av_al_gain << " new remove " << remove_als.size() << " lost gain " << lost_information_gain <<
				       	" new insert " << insert_als.size() << " gain "<<  max_parts_gain << endl;

			
			}



			

		}
		// if it is possible to increase information gain by inserting this alignment, we 
		// recursively insert all small pieces
		// we start with the high gain pieces, because they have more potential to loose information gain, which can be used to abort the current insert
		double upper_bound_of_gain = max_parts_gain - lost_information_gain;
		if(upper_bound_of_gain > 0) {
			double all_lost_gain = 0; // difference between gain realized by inserting a piece versus the gain of that piece
			double total_recursive_gain = 0; // total gain of all subtrees processed so far

			// all_ins/rem accumulated over all recursive steps below current. If we take current step, those are tranferred to inserted_alignments/removed_alignments
			vector<const pw_alignment*> all_inserted;
			vector<const pw_alignment*> all_removed;

			// remove before recursion
			size_t remove_count = 0;
			for(alset::iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
				const pw_alignment * pwa = *it;
				removed_alignments.push_back(pwa);				
				
	//			cout <<"level " << level << "REMOVE NOW: " << endl;
	//			pwa->print();
	//			cout << " ovrlp size " << ovrlp.size() << endl;
				
				
				ovrlp.remove_alignment_nodelete(pwa);


				remove_count++;
			}

			cout << " before recursion, we have removed " << remove_count << " alignments with gain " << lost_information_gain <<  " upper bound of gain is " << upper_bound_of_gain << endl;
			size_t alnum = 0;
			for(map<double, pw_alignment>::reverse_iterator it = ordered_parts.rbegin(); it!=ordered_parts.rend(); ++it) {
				double ub_of_thisgain = it->first;
				double thisgain;
				vector<const pw_alignment*> this_inserted;
				vector<const pw_alignment*> this_removed;
				pw_alignment thisal = it->second;

			//	thisal.print();
			//	cout << endl;

				cout <<"start recursion: level "<< level << " start rec "<< alnum << " on al length " << thisal.alignment_length() << endl; 
				// this function only changes ovrlp if this was locally good. We may need to undo the changes if they were globally bad
				lazy_split_insert_step(ovrlp, level + 1, &(thisal), this_inserted, this_removed, thisgain, outs);
				double lost_gain = ub_of_thisgain - thisgain; // we have stepped down that much from the previous upper bound of gain
				all_lost_gain += lost_gain;
				total_recursive_gain += thisgain;


				cout <<"end recursion: level "<< level << " end rec "<< alnum << " on al length " << thisal.alignment_length() << " subtree local gain " << thisgain << " ub was " << ub_of_thisgain << " lost gain is " << lost_gain << " total lost gain is " << all_lost_gain << " total gain of this level until now: " << total_recursive_gain <<  " ub is " << upper_bound_of_gain << endl; 
				for(size_t i=0; i<this_inserted.size(); ++i) {
					all_inserted.push_back(this_inserted.at(i));
				}
				for(size_t i=0; i<this_removed.size(); ++i) {
					all_removed.push_back(this_removed.at(i));
				}

				// if we have lost more than we may gain, we undo the entire recursive insertion subtree from here
				if(all_lost_gain > upper_bound_of_gain) {
					// the remove_alignment function requires identical address of the alignment to be removed from overlap
					// it does remove that alignment from overlap and calls delete
					for(size_t i=0; i<all_inserted.size(); ++i) {
						// after recursion we can remove with delete
						ovrlp.remove_alignment(all_inserted.at(i));
					}

					for(size_t i=0; i<all_removed.size(); ++i) {
						// reinsert same pointer
						ovrlp.insert_without_partial_overlap_p((pw_alignment*)all_removed.at(i));
					}

					cout << "undo subtree: level " << level << " on al length " << al->alignment_length() << " gain " << av_al_gain << " ABORT because of low upper bound of gain" << endl;
					// gain of this subtree is 0 because we do not insert
					local_gain = 0;


					return;
				} 
				alnum++;
			} // for ordered_parts

			// we take all alignments underneath of current recursion level if they have more gain than all the one that we had to remove for the current split operation
			if(total_recursive_gain > lost_information_gain) {
				local_gain = total_recursive_gain - lost_information_gain;

				for(size_t i=0; i<all_inserted.size(); ++i) {
					inserted_alignments.push_back(all_inserted.at(i));
				}
				for(size_t i=0; i<all_removed.size(); ++i) {
					removed_alignments.push_back(all_removed.at(i));
				}


	cout << "take subtree: level " << level << " on al length " << al->alignment_length() << " gain " << av_al_gain << " causes to remove " << all_removed.size() << " alignments with " << lost_information_gain << 
			" gain. We want to insert " << all_inserted.size() << " small pieces. TOTAL gain is " << local_gain << endl;

				return;
			
			}
		
		
			// if we have not been able to achieve a global gain, we undo all changes
			for(size_t i=0; i<all_inserted.size(); ++i) {
				// remove with delete
				ovrlp.remove_alignment(all_inserted.at(i));
			}
			for(size_t i=0; i<all_removed.size(); ++i) {
			//	cout << " level "<< level <<" reinsert removed alignment" << endl;
			//	all_removed.at(i).print();
			//	cout << endl;

				ovrlp.insert_without_partial_overlap_p((pw_alignment*)all_removed.at(i));
			}
		
		
		}

	


	}

}



template<typename T>
void initial_alignment_set<T>::compute_simple_lazy_splits(overlap & o,ofstream & outs) {
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;

	for (size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;

		vector<const pw_alignment*> ins_als;
		vector<const pw_alignment*> rem_als;
		lazy_split_insert_step(o,0, al, ins_als, rem_als, gain_of_al, outs);
		
		double gain1;
		double gain2;
		common_model.gain_function(*(al), gain1, gain2,outs);
		double vgain = (gain1 + gain2) / 2 - base_cost;


		double rgain = 0;
		double igain = 0;
		for(size_t j=0; j<ins_als.size(); ++j) {
			common_model.gain_function(*ins_als.at(j), gain1, gain2, outs);
			vgain = (gain1+gain2)/2 - base_cost;
			igain+=vgain;
		}
		for(size_t j=0; j<rem_als.size(); ++j) {
			common_model.gain_function(*rem_als.at(j), gain1, gain2, outs);
			vgain = (gain1+gain2)/2 - base_cost;
			rgain+=vgain;
		}
		

		double checkgain = igain - rgain;

		cout << " al " << i << " gain in " << vgain << " gain out " << gain_of_al << " check gain " << checkgain << endl;
		cout << " out: " << rem_als.size() << " gain " << rgain << " in: " << ins_als.size() << " gain " << igain << endl;

		if(!ins_als.empty()) {
			used++;
			total_gain+=gain_of_al;
			pcs_ins+=ins_als.size();
			pcs_rem+=rem_als.size();
		} else {
			not_used++;
		}
		
		// delete removed alignments
		for(size_t j=0; j<rem_als.size(); ++j) {
			delete rem_als.at(j);
		}

	
	}
	result_gain = total_gain;
	cout << "Input size " << sorted_original_als.size() << " max gain " << max_gain << endl;
	cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << endl;
	cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << endl;



}

// TODO INDEPENDENT SET optimization


/**
	greedy test function 
**/
template<typename T>
void initial_alignment_set<T>::compute_simple(overlap & o,ofstream & outs) {
	size_t used = 0;
	size_t not_used = 0;
	double total_gain = 0;
	size_t pcs_ins = 0;
	size_t pcs_rem = 0;
//	cout<< "sorted alignment size" << sorted_original_als.size()<<endl;	
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;

		// TODO remove
		double gain1, gain2;
		common_model.gain_function(*(al), gain1, gain2,outs);
		gain1-=base_cost;
//		cout << endl<<"at alignment " << i << " length " << al->alignment_length() << " al base gain " << gain1 << endl;
		//al->print();
//		cout << endl;


		alset remove_als;
		vector<pw_alignment> insert_als;

//		// TODO
//		cout << " splitpoints on: " << al << " at " << i <<  " size " <<sorted_original_als.size()<< endl; 
//		al->print();
//		cout << endl;	
//		o.print_all_alignment();


		splitpoints spl(*al, o, data);
		spl.recursive_splits();
		spl.split_all(remove_als, insert_als);
		vector<double> insert_gains(insert_als.size(), 0);

		for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			double g1;
			double g2;
			common_model.gain_function(*(*it), g1, g2,outs);
		//	if(g2<g1) g1 = g2;
			g1-=base_cost;
	//		cout << "r " << (*it)->alignment_length() << " g " << g1 << endl;
			gain_of_al -= g1;

		//	(*it)->print();
		//	cout << endl;

		}	
		for(size_t j=0; j<insert_als.size(); ++j) {
			
			double g1;
			double g2;
			common_model.gain_function(insert_als.at(j), g1, g2,outs);
		//	if(g1<g2) g1 = g2;
			g1-=base_cost;
			insert_gains.at(j) = g1;
			if(g1>0) {
		//		cout << "i " << (insert_als.at(j)).alignment_length() << " g " << g1 << endl;
				gain_of_al += g1;

			//	insert_als.at(j).print();
			//	cout << endl;
			}

		}
//		cout << " al " << i << " rem " << remove_als.size() << " insert " << insert_als.size() << " gain " << gain_of_al << endl;	
		if(gain_of_al>=0) {
			used++;
			for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
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
//	cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << endl;
//	cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << endl;

}



compute_cc::compute_cc(const all_data & dat): als_on_reference(dat.numSequences()), last_pos(dat.numSequences(), 0) {
	max_al_ref_length = 0;
	for(size_t i=0; i<dat.numAlignments(); ++i) {
		const pw_alignment * a = &(dat.getAlignment(i));
		alignments.insert(a);
		add_on_mmaps(a);
	}

}


void compute_cc::add_on_mmaps(const pw_alignment * pwa) {
	size_t ref1 = pwa->getreference1();
	size_t ref2 = pwa->getreference2();
/*	
	cout << "INS " <<endl;
	pwa->print();
	cout << endl;	
*/
	als_on_reference.at(ref1).insert(make_pair(pwa->getbegin1(), pwa));
	als_on_reference.at(ref1).insert(make_pair(pwa->getend1(), pwa));
	als_on_reference.at(ref2).insert(make_pair(pwa->getbegin2(), pwa));
	als_on_reference.at(ref2).insert(make_pair(pwa->getend2(), pwa));

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


void compute_cc::remove_on_mmaps(const pw_alignment * al) {
	size_t ref1 = al->getreference1();
	size_t left, right;
	al->get_lr1(left, right);

	pair<multimap<size_t, const pw_alignment *>::iterator, multimap<size_t, const pw_alignment*>::iterator > l1 = 
		als_on_reference.at(ref1).equal_range(left);
	for(multimap<size_t, const pw_alignment*>::iterator it = l1.first; it!=l1.second; ++it) {
		if(it->second == al) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}
	pair<multimap<size_t, const pw_alignment *>::iterator, multimap<size_t, const pw_alignment*>::iterator  > r1 = 
		als_on_reference.at(ref1).equal_range(left);
	for(multimap<size_t, const pw_alignment*>::iterator it = r1.first; it!=r1.second; ++it) {
		if(it->second == al) {
			als_on_reference.at(ref1).erase(it);
			break;
		}	
	}

	size_t ref2 = al->getreference2();
	al->get_lr2(left, right);
	pair<multimap<size_t, const pw_alignment *>::iterator, multimap<size_t, const pw_alignment*>::iterator  > l2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(multimap<size_t, const pw_alignment*>::iterator it = l2.first; it!=l2.second; ++it) {
		if(it->second == al) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}
	pair<multimap<size_t, const pw_alignment *>::iterator, multimap<size_t, const pw_alignment*>::iterator  > r2 = 
		als_on_reference.at(ref2).equal_range(left);
	for(multimap<size_t, const pw_alignment*>::iterator it = r2.first; it!=r2.second; ++it) {
		if(it->second == al) {
			als_on_reference.at(ref2).erase(it);
			break;
		}	
	}

}


void compute_cc::compute(vector<set< const pw_alignment *, compare_pw_alignment> > & ccs) {
	set <const pw_alignment *, compare_pw_alignment> seen;
	multimap<size_t , set<const pw_alignment *, compare_pw_alignment> > sorter; // sorts ccs to give largest first
	for(set<const pw_alignment *, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = *it;
		set<const pw_alignment *, compare_pw_alignment>::iterator seenal = seen.find(al);
		cout << " seen " << seen.size() << endl;
		if(seenal == seen.end()) {
			cout << " getcc" << endl;
			set<const pw_alignment *, compare_pw_alignment> cc;
			get_cc(al, cc, seen);
			cout << "FOUND CC size " << cc.size() << endl;
			sorter.insert(make_pair(cc.size(), cc));
		}
	
	}
	for(multimap<size_t , set<const pw_alignment *, compare_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); ++it) {
		ccs.push_back(it->second);
	}



}

void compute_cc::get_cc(const pw_alignment * al, set <const pw_alignment *, compare_pw_alignment> & cc, set <const pw_alignment *, compare_pw_alignment> & seen) {
	
	size_t left, right;
	al->get_lr1(left, right);
	cc_step(al->getreference1(), left, right, cc, seen);
	al->get_lr2(left, right);
	cc_step(al->getreference2(), left, right, cc, seen);
}




// TODO further improvements to this function are possible if we store intervals on the references in which all alignments were already processed
void compute_cc::cc_step(size_t ref, size_t left, size_t right, set <const pw_alignment *, compare_pw_alignment> & cc, set <const pw_alignment *, compare_pw_alignment>  & seen ) {
	// search bounds (where could other alignments which overlap with the current one start or end)
	// leftbound: all alignments starting at leftbound or ealier either have the end in the search interval or no overlap with the search interval
	cout << " cc step " << ref << " fr " << left << " to " << right << " seen is " << seen.size() << " we are on " << alignments.size() << " alignments" <<  endl; 
	multimap<size_t, const pw_alignment *>::iterator searchbegin;
	multimap<size_t, const pw_alignment *>::iterator searchend;
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

	set <const pw_alignment *, compare_pw_alignment> seen1;
	set <const pw_alignment *, compare_pw_alignment> seen2;


	// search for overlap first, then do all recursive calls after overlapping alignments have been put to seen 
	// this reduces the maximal recursion level
	size_t numseen = 0;
	for(multimap<size_t, const pw_alignment *>::iterator it = searchbegin; it!=searchend; ++it) {
		const pw_alignment * al = it->second;


//		cout << " at " << it->first << endl;
		set <const pw_alignment *, compare_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) { // if current al not contained in any connected component
//			cout << " not seen" << endl;
			size_t aleft, aright;	
	//		size_t leftmost_point_of_al_on_ref = numeric_limits<size_t>::max(); 
			if(al->getreference1()==ref) {
				al->get_lr1(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					seen1.insert(al);
					cc.insert(al);
			//		cout << "ovlr " << cc.size() << " "  << seen.size() <<  " ref "<< ref << " : " << left << " " << right << " ovrlaps " << endl;
			//		al->print();
				//	al->get_lr2(aleft, aright);
				//	cc_step(al->getreference2(), aleft, aright, cc, seen);
				}
		//		if(aleft < leftmost_point_of_al_on_ref ) {
		//			 leftmost_point_of_al_on_ref = aleft;
		//		}
			}
			if(al->getreference2()==ref) {
				al->get_lr2(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					seen2.insert(al);
					cc.insert(al);
				//	cout << "ovlr " << cc.size() << " "  << seen.size() << " ref "<< ref << " : " << left << " " << right << " ovrlaps " << endl;
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
	cout << " found overlap with " << seen1.size() << " and " << seen2.size() << " alignments, already seen before: " << numseen << endl;


	// now remove all seen alignments to be faster

	for(set<const pw_alignment *>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		remove_on_mmaps(*it);
	}
	for(set<const pw_alignment *>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		remove_on_mmaps(*it);
	}

	size_t debugsum = 0;
	for(size_t i=0; i<als_on_reference.size(); ++i) {
		debugsum+=als_on_reference.at(i).size();
	}
	cout << " mmaps length " <<debugsum << endl;

	for(set<const pw_alignment *>::iterator it = seen1.begin(); it!=seen1.end(); ++it) {
		const pw_alignment * al = *it;
		size_t aleft, aright;	
		al->get_lr2(aleft, aright);
		cc_step(al->getreference2(), aleft, aright, cc, seen);
	} 	
	for(set<const pw_alignment *>::iterator it = seen2.begin(); it!=seen2.end(); ++it) {
		const pw_alignment * al = *it;
		size_t aleft, aright;	
		al->get_lr1(aleft, aright);
		cc_step(al->getreference1(), aleft, aright, cc, seen);
	} 	
	
}
/*
template<typename tmodel>
clustering<tmodel>::clustering(overlap & o, all_data & d,tmodel & m):overl(o),data(d),model(m),als_on_ref(data.numSequences()),gain(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),ava(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),res(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))){
	/*	for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * a = &(data.getAlignment(i));
			alignments.insert(a);
		}*

	}
template<typename tmodel>
clustering<tmodel>::~clustering(){}

template<typename tmodel>
void clustering<tmodel>::als_on_reference(const pw_alignment * p) {
/*	size_t ref1 = p->getreference1();
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
			for(multimap<size_t, pw_alignment*>::iterator it =als_on_ref.at(i).begin(); it!=als_on_ref.at(i).end(); ++it){
				const pw_alignment * p = it-> second;
				//p->print();
				size_t position = it->first;
				cout<<"position: "<<position<<endl;
				model.gain_function(*p,g1,g2);
				gain.at(p->getreference1()).at(p->getreference2()).at(position) += g1;//position has been used to specify the part of reference contains P's sample
				gain.at(p->getreference2()).at(p->getreference1()).at(position) += g2;
				cout<<"gain2: "<< gain.at(p->getreference2()).at(p->getreference1()).at(position)<<endl;
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
			cout<<"ava"<<	ava.at(i).at(j).at(o);
			}
		}
		}
//	}
}

template<typename tmodel>
	void clustering<tmodel>::update_clusters(size_t acc){
	/*	size_t iteration = 100;
		double damp_value = 0.6;
		vector<size_t> examplar(data.numAcc(),0); //examplars of the class of acc, acc will be the center.
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
		}*
		
	}

template<typename tmodel>
	void clustering<tmodel>::update_clusters(){
		vector<vector<vector<double> > >examplar;
		vector<size_t> center;
		for(size_t i = 0; i<data.numSequences(); i++){
			for(size_t n = 0 ; n < gain[i][i].size();n++){
				examplar.at(i).at(i).at(n) = res.at(i).at(i).at(n) + ava.at(i).at(i).at(n);
				if (examplar.at(i).at(i).at(n)>0){
					center.push_back(i);
				}
			}
		//	cout<<"center at: "<< i << endl;

		}
		vector<size_t> idx(data.numSequences(),0);
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
	//		cout << "center of "<< k << " is "<<idx.at(k)<<endl;
		}

	}

*/



template<typename tmodel>
void affpro_clusters<tmodel>::add_alignment(const pw_alignment *al,ofstream& outs) {
	// Get identifiers for both parts of the pairwise alignment
	stringstream sstr1;
//	cout<<"data1 ad in add_al: "<< & dat << endl;	
	size_t left1, right1;
	al->get_lr1(left1, right1);
	sstr1 << al->getreference1()<<":"<<left1;
	stringstream sstr2;
	size_t left2, right2;
	al->get_lr2(left2, right2);
	sstr2 << al->getreference2()<<":"<<left2;
	string ref1name = sstr1.str();
	string ref2name = sstr2.str();

	// Get data indices for both alignment parts/ sequence pieces
	size_t ref1idx = 0;
	size_t ref2idx = 0;
	map<string, size_t>::iterator find1 = sequence_pieces.find(ref1name);
	map<string, size_t>::iterator find2 = sequence_pieces.find(ref2name);
	if(find1 == sequence_pieces.end()) {
		ref1idx = sequence_pieces.size();
		sequence_pieces.insert(make_pair(ref1name, ref1idx));
		sequence_names.push_back(ref1name);
		sequence_lengths.push_back(right1 - left1 + 1);
	} else {
		ref1idx = find1->second;
	}
	if(find2 == sequence_pieces.end()) {
		ref2idx = sequence_pieces.size();
		sequence_pieces.insert(make_pair(ref2name, ref2idx));
		sequence_names.push_back(ref2name);
		sequence_lengths.push_back(right2 - left2 + 1);
	} else {
		ref2idx = find2->second;
	}
//	cout<<"data2 ad in add_al: "<< & dat << endl;	
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
		simmatrix.resize(max, vector<double>(max, -HUGE_VAL));
	}
	double c1;
	double c2;
	double m1; 
	double m2;
//	al->print();
//	cout<<"data3 ad in add_al: "<< & dat << endl;	
//	dat.numAcc();
//	cout << " dat adress " << & dat<< endl;
	model.cost_function(*al, c1, c2, m1, m2,outs);
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

#endif
