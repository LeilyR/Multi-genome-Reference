#include "model.hpp"


#ifndef MODEL_CPP
#define MODEL_CPP

template<typename T>
void initial_alignment_set<T>::compute(overlap & o) {

	compute_simple(o);

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
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;

		// TODO remove
		double gain1, gain2;
		model.gain_function(*(sorted_original_als.at(i)), gain1, gain2);
		gain1-=base_cost;
		cout << endl<< "at alignment " << i << " length " << al->alignment_length() << " al base gain " << gain1 << endl;
		al->print();
		cout << endl;


		alset remove_als;
		vector<pw_alignment> insert_als;
//		cout << " splitpoints on: " << al << " at " << i <<  " size " <<sorted_original_als.size()<< endl; 
//		al->print();
//		cout << endl;	
		splitpoints spl(*al, o, data);
		spl.find_initial_split_point();
		spl.split_all(remove_als, insert_als);

		for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			double g1;
			double g2;
			model.gain_function(*(*it), g1, g2);
		//	if(g2<g1) g1 = g2;
			g1-=base_cost;
			cout << "r " << (*it)->alignment_length() << " r " << g1 << endl;
			gain_of_al -= g1;

			(*it)->print();
			cout << endl;

		}	
		for(size_t j=0; j<insert_als.size(); ++j) {
			double g1;
			double g2;
			model.gain_function(insert_als.at(j), g1, g2);
		//	if(g1<g2) g1 = g2;
			g1-=base_cost;
			cout << "i " << (insert_als.at(j)).alignment_length() << " g " << g1 << endl;
			gain_of_al += g1;

			insert_als.at(j).print();
			cout << endl;

		}
		cout << " al " << i << " rem " << remove_als.size() << " insert " << insert_als.size() << " gain " << gain_of_al << endl;	
		if(gain_of_al>=0) {
			used++;
			for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
				o.remove_alignment(*it);
				pcs_rem++;
			}	
			for(size_t j=0; j<insert_als.size(); ++j) {
				o.insert_without_partial_overlap(insert_als.at(j));
				pcs_ins++;
			}

			total_gain+=gain_of_al;

		} else {
			not_used++;
		}

	
	
	}
	result_gain = total_gain;
	cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << endl;
	cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << endl;

}



compute_cc::compute_cc(const all_data & dat): als_on_reference(dat.numSequences()), last_pos(dat.numSequences(), 0) {
	for(size_t i=0; i<dat.numAlignments(); ++i) {
		const pw_alignment * a = &(dat.getAlignment(i));
		alignments.insert(a);
		add_on_mmaps(a);
	}

}


void compute_cc::add_on_mmaps(const pw_alignment * pwa) {
	size_t ref1 = pwa->getreference1();
	size_t ref2 = pwa->getreference2();

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


void compute_cc::compute(vector<set< const pw_alignment *, compare_pw_alignment> > & ccs) {
	set <const pw_alignment *, compare_pw_alignment> seen;
	multimap<size_t , set<const pw_alignment *, compare_pw_alignment> > sorter;
	for(set<const pw_alignment *, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = *it;
		set<const pw_alignment *, compare_pw_alignment>::iterator seenal = seen.find(al);
	//	cout << " seen " << seen.size() << endl;
		if(seenal == seen.end()) {
	//		cout << " getcc" << endl;
			set<const pw_alignment *, compare_pw_alignment> cc;
			get_cc(al, cc, seen);
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

void compute_cc::cc_step(size_t ref, size_t left, size_t right, set <const pw_alignment *, compare_pw_alignment> & cc, set <const pw_alignment *, compare_pw_alignment>  & seen ) {
	// search bounds (where could other alignments which overlap with the current one start or end)
	size_t leftbound = left;
	if(left > max_al_ref_length) {
		leftbound = left - max_al_ref_length;
	} else {
		leftbound = 0;
	}
	size_t rightbound = right;
	if(right + max_al_ref_length < last_pos.at(ref)) {
		rightbound = right + max_al_ref_length;
	} else {
		rightbound = last_pos.at(ref);
	}

//	cout << " on ref " << ref << " lb " << leftbound << " rb " << rightbound << " l " << left << " r " << right << endl;
	for(multimap<size_t, const pw_alignment *>::iterator it = als_on_reference.at(ref).upper_bound(leftbound); it!=als_on_reference.at(ref).end(); ++it) {
		const pw_alignment * al = it->second;


//		cout << " at " << it->first << endl;
		set <const pw_alignment *, compare_pw_alignment>::iterator seenal = seen.find(al);
		if(seenal == seen.end()) { // if current al not contained in any connected component
//			cout << " not seen" << endl;
			size_t aleft, aright;
			size_t leftmost_point_of_al_on_ref = numeric_limits<size_t>::max(); 
			if(al->getreference1()==ref) {
				al->get_lr1(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					cc.insert(al);
//					cout << "ovlr " << cc.size() << " "  << seen.size() << endl;
					cc_step(al->getreference2(), aleft, aright, cc, seen);
				}
				if(aleft < leftmost_point_of_al_on_ref ) {
					 leftmost_point_of_al_on_ref = aleft;
				}
			}
			if(al->getreference2()==ref) {
				al->get_lr2(aleft, aright);
				if(aright >= left && aleft <= right) {
					seen.insert(al);
					cc.insert(al);
//					cout << "ovlr " << cc.size() << " "  << seen.size() << endl;
					cc_step(al->getreference1(), aleft, aright, cc, seen);
				}
				if(aleft < leftmost_point_of_al_on_ref ) {
					 leftmost_point_of_al_on_ref = aleft;
				}

			}
			if(leftmost_point_of_al_on_ref > rightbound) {
				break;
			}
		}
	}
	
	
}

#endif
