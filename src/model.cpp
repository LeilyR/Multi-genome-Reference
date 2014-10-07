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
	cout<< "sorted alignment size" << sorted_original_als.size()<<endl;	
	for(size_t i=0; i<sorted_original_als.size(); ++i) {
		const pw_alignment * al = sorted_original_als.at(i);
		double gain_of_al = 0;

		// TODO remove
		double gain1, gain2;
		common_model.gain_function(*(al), gain1, gain2);
		gain1-=base_cost;
		cout << endl<<"at alignment " << i << " length " << al->alignment_length() << " al base gain " << gain1 << endl;
		al->print();
		cout << endl;


		alset remove_als;
		vector<pw_alignment> insert_als;

		// TODO
//		cout << " splitpoints on: " << al << " at " << i <<  " size " <<sorted_original_als.size()<< endl; 
//		al->print();
//		cout << endl;	
	
//		o.print_all_alignment();


		splitpoints spl(*al, o, data);
		spl.find_initial_split_point();
		spl.split_all(remove_als, insert_als);
		vector<double> insert_gains(insert_als.size(), 0);

		for(alset::const_iterator it = remove_als.begin(); it!=remove_als.end(); ++it) {
			double g1;
			double g2;
			common_model.gain_function(*(*it), g1, g2);
		//	if(g2<g1) g1 = g2;
			g1-=base_cost;
			cout << "r " << (*it)->alignment_length() << " g " << g1 << endl;
			gain_of_al -= g1;

			(*it)->print();
			cout << endl;

		}	
		for(size_t j=0; j<insert_als.size(); ++j) {
			
			double g1;
			double g2;
			common_model.gain_function(insert_als.at(j), g1, g2);
		//	if(g1<g2) g1 = g2;
			g1-=base_cost;
			insert_gains.at(j) = g1;
			if(g1>0) {
				cout << "i " << (insert_als.at(j)).alignment_length() << " g " << g1 << endl;
				gain_of_al += g1;

				insert_als.at(j).print();
				cout << endl;
			}

		}
		cout << " al " << i << " rem " << remove_als.size() << " insert " << insert_als.size() << " gain " << gain_of_al << endl;	
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
	cout << "Used " << used << " alignments with total gain " << total_gain << " not used: " << not_used << endl;
	cout << "pieces removed: " << pcs_rem << " pieces inserted: " << pcs_ins << endl;

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
	multimap<size_t , set<const pw_alignment *, compare_pw_alignment> > sorter; // sorts ccs to give largest first
	for(set<const pw_alignment *, compare_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		const pw_alignment * al = *it;
		set<const pw_alignment *, compare_pw_alignment>::iterator seenal = seen.find(al);
	//	cout << " seen " << seen.size() << endl;
		if(seenal == seen.end()) {
	//		cout << " getcc" << endl;
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
			//		cout << "ovlr " << cc.size() << " "  << seen.size() <<  " ref "<< ref << " : " << left << " " << right << " ovrlaps " << endl;
			//		al->print();
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
				//	cout << "ovlr " << cc.size() << " "  << seen.size() << " ref "<< ref << " : " << left << " " << right << " ovrlaps " << endl;
				//	al->print();
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

template<typename tmodel>
clustering<tmodel>::clustering(overlap & o, all_data & d,tmodel & m):overl(o),data(d),model(m),als_on_ref(data.numSequences()),gain(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),ava(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))),res(data.numSequences(),vector<vector<double> >(data.numSequences(),vector<double>(500000,0))){
	/*	for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * a = &(data.getAlignment(i));
			alignments.insert(a);
		}*/

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
	als_on_ref.at(ref2).insert(make_pair(p->getend2(), p));*/

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
		}*/
		
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
			cout<<"center at: "<< i << endl;

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
			cout << "center of "<< k << " is "<<idx.at(k)<<endl;
		}

	}


template<typename tmodel>
void affpro_clusters<tmodel>::add_alignment(const pw_alignment *al) {

	// Get identifiers for both parts of the pairwise alignment
	stringstream sstr1;
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
	model.cost_function(*al, c1, c2, m1, m2);
	// preferences
	simmatrix.at(ref1idx).at(ref1idx) = -c1 - base_cost;
	simmatrix.at(ref2idx).at(ref2idx) = -c2 - base_cost;
	// similarity (i,j): i has exemplar j, this means modify j to i?
	simmatrix.at(ref2idx).at(ref1idx) = -m2;
	simmatrix.at(ref1idx).at(ref2idx) = -m1;


}


#endif
