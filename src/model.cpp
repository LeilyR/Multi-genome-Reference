#include "model.hpp"

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

		alset remove_als;
		vector<pw_alignment> insert_als;
		cout << " splitpoints on: " << al << endl; 
		al->print();
	cout << endl;	
		splitpoints spl(*al, o, data);
		spl.find_initial_split_point();
		spl.split_all(remove_als, insert_als);

		for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			double g1;
			double g2;
			model.gain_function(*(*it), g1, g2);
			if(g2<g1) g1 = g2;
			gain_of_al -= g1;
		}	
		for(size_t j=0; j<insert_als.size(); ++j) {
			double g1;
			double g2;
			model.gain_function(insert_als.at(j), g1, g2);
			if(g1<g2) g1 = g2;
			gain_of_al += g1;
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

	cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << endl;
	cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << endl;

}




