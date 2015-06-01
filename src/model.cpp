#include "model.hpp"

#ifndef MODEL_CPP
#define MODEL_CPP

template<typename T>
void initial_alignment_set<T>::compute(overlap & o, ofstream & outs) {

	compute_simple(o,outs);

}




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
		cout << endl;


		alset remove_als;
		vector<pw_alignment> insert_als;

//		// TODO
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
	//	cout << " seen " << seen.size() << endl;
		if(seenal == seen.end()) {
	//		cout << " getcc" << endl;
			set<const pw_alignment *, compare_pw_alignment> cc;
			get_cc(al, cc, seen);
	//		cout << "FOUND CC size " << cc.size() << endl;
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
//	cout << " cc step " << ref << " fr " << left << " to " << right << " seen is " << seen.size() << " we are on " << alignments.size() << " alignments" <<  endl; 
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
//	cout << " found overlap with " << seen1.size() << " and " << seen2.size() << " alignments, already seen before: " << numseen << endl;


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
//	cout << " mmaps length " <<debugsum << endl;

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
	finding_centers::finding_centers(all_data & d):data(d),AlignmentsFromClustering(data.numSequences()),centersOfASequence(data.numSequences(),vector<size_t>()){
	}
	finding_centers::~finding_centers(){}
	void finding_centers::setOfAlignments(map<string,vector<pw_alignment> > & alignmentsOfClusters){//set alignments of each reference
		for(size_t i = 0; i < data.numSequences(); i ++){
			const dnastring sequence = data.getSequence(i);
			for (map<string,vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end();it++){
				assert(it->second.size() != 0);
				string center = it->first;
				vector<string> center_parts;
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
							multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left1);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							}else continue;
						}else if(p->getreference2()== i && p->getreference1()== center_ref&&left1 == center_left){
							multimap<size_t , pw_alignment*>::iterator it2 = AlignmentsFromClustering.at(i).find(left2);
							if(it2 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							}else continue;
						}else continue;
					}
					if(i == center_ref){
						if(p->getreference1()== i && p->getreference2()== i){
							multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(left1);				
							if(it == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left1,p));
							}
							multimap<size_t, pw_alignment*>::iterator it1 = AlignmentsFromClustering.at(i).find(left2);				
							if(it1 == AlignmentsFromClustering.at(i).end()){
								AlignmentsFromClustering.at(i).insert(make_pair(left2,p));
							}
						}else{
							multimap<size_t, pw_alignment*>::iterator it = AlignmentsFromClustering.at(i).find(center_left);
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
		for(map<string, vector<pw_alignment> >::iterator it = alignmentsOfClusters.begin(); it != alignmentsOfClusters.end(); it++){
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
	void finding_centers::center_frequency(map<string,vector<pw_alignment> > & alignmentsOfClusters){//it basically returns indices of centers on each sequence.
		setOfAlignments(alignmentsOfClusters);
		findMemberOfClusters(alignmentsOfClusters);	
		vector<string> center_index;
		for(map<string, vector<pw_alignment> >::iterator it2=alignmentsOfClusters.begin(); it2 != alignmentsOfClusters.end();it2++){
			string cent = it2->first;
			center_index.push_back(cent);
		}
		cout<< "center index size" << center_index.size()<<endl;
		for(size_t i = 0; i< data.numSequences(); i++){
			cout << "sequence: " << i << endl;
			const dnastring & sequence = data.getSequence(i);
			for(size_t n= 0; n < sequence.length(); n++){
				multimap<size_t, pw_alignment*>::iterator it=AlignmentsFromClustering.at(i).find(n);
				if(it != AlignmentsFromClustering.at(i).end()){
					stringstream member;
					member << i << ":" << n;
					map<string,string>::iterator it1 = memberOfCluster.find(member.str());
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
		for(size_t i = 0; i < centersOfASequence.size(); i ++){
			cout<< " centers on sequence " << i << " are " <<endl;
			for(size_t j =0 ; j < centersOfASequence.at(i).size(); j++){
				cout<< centersOfASequence.at(i).at(j) <<endl;
			}
		}
	}
	vector<size_t>  finding_centers::get_center(size_t seq_id)const{
		return centersOfASequence.at(seq_id);
	}
	suffix_tree::suffix_tree(all_data & d, finding_centers & c):data(d),centers(c){
	
	}
	suffix_tree::~suffix_tree(){
	
	}
	void suffix_tree::create_suffix(size_t seq_id){
		suffixes.clear();
		vector<size_t> successive_centers = centers.get_center(seq_id);
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
				cout<< " current " << current_center << " last center "<< last_center<<endl;
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
			cout << " suffixes are " << endl;
			for(size_t i =0; i < suffixes.size();i++){
				string suf = suffixes.at(i);
				for( size_t j =0; j < suf.size() ; j ++){
					cout << int(suf.at(j))<< " ";
				}
					cout<< " " << endl;
			}
		}
	}
	void suffix_tree::find_a_node(size_t& node_number,size_t& parent_node, string& node){
		for(size_t i = parent_node; i < nodes.size(); i ++){
		//	cout << "size of nodes: " << nodes.size() <<endl;
			string current_node = nodes.at(i);
			if(node == current_node){
				cout << "i "<< i << endl;
				node_number = i;
				break;
			}
		}
	}
	void suffix_tree::find_next_sibling(size_t current_node, size_t next_sibling){
	//	

	}
	void suffix_tree::find_next_firstparent(size_t node_index, size_t next_first_parent){
		for(size_t i = node_index+1; i <nodes.size(); i++){
			map<string, size_t>::iterator it = firstParent.find(nodes.at(i));
			if(it != firstParent.end() && it->second == i){
				next_first_parent = i;
			}
		}
	}
/*	void suffix_tree::find_parent_path(size_t node_index,vector<size_t> path){
		size_t parent_position;
		vector<size_t> first_path;
		for(map<vector<size_t>,size_t>::iterator it = branch_counter.begin();it !=branch_counter.end(); it++){
			first_path = it->first;
			for(size_t i =0; i < first_path.size(); i++){
				if(first_path.at(i) == node_index){
					parent_position = i;
					break;
				}
			}
			break;
		}
		for(size_t i =0; i < parent_position; i++){
			path.push_back(first_path.at(i));
		}
		size_t firstparent_index = path.at(0);
		cout<< "first parent: "<< firstparent_index <<endl;
	}*/
	void suffix_tree::make_a_tree(){
		size_t active_length = 0;
		size_t active_node = 0;
		vector<string> first_parent;//To Do: remove it and try to use firstParent map instead of this vector, because they are basically the same thing
		for(size_t seq_id =0; seq_id < data.numSequences(); seq_id++){
			create_suffix(seq_id);
			cout<< "seq id"<<seq_id<<endl;
			for(size_t i = 0; i < suffixes.size(); i++){
				cout<< "suffix i " << i << endl;
				string current = suffixes.at(i);
				size_t temp = active_node;
				size_t last_common_index = 0;
				if(first_parent.size()!= 0){
					for(size_t j= 0; j< first_parent.size(); j++){
						if (current.at(0) == first_parent.at(j).at(0)){
							string current_parent = first_parent.at(j);
							size_t node_index;//current parent node index
						//	vector<size_t> parent_path;
							map<string, size_t>::iterator it = firstParent.find(current_parent);
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
								pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator > p1 = nodes_relation.equal_range(n_index);
								size_t biggest_child =0;
								for(multimap<size_t, size_t>::iterator it = p1.first; it!=p1.second; ++it){
									if(it != nodes_relation.end()){
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
										for(size_t child_node = node_index+1; child_node <= biggest_child;child_node ++){
											if(nodes.at(child_node).at(0)== current.at(first_index_after_parent)){
												adding_a_child ++;
												string updated_current;
												for(size_t update =first_index_after_parent; update< current.size(); update ++){
													updated_current += current.at(update);
												}
												current = updated_current;
												cout << "child node: "<<child_node << endl;
												current_parent = nodes.at(child_node);
												break;
											}else continue;
										}
									}
									if(adding_a_child == 0){//Creating a new child node(*)//inja eshkal dare!
										cout << "here!"<<endl;
										for(size_t shift = node_index+1; shift <nodes.size(); shift++){
											map<string, size_t>::iterator it = firstParent.find(nodes.at(shift));
											if(it != firstParent.end() && it->second == shift){
												it->second = it->second+1;
											}
										}
										string last_node = nodes.at(nodes.size()-1);
										for(size_t shift =nodes.size()-1 ; shift > node_index+1; shift --){
											nodes.at(shift)=nodes.at(shift-1);
										}
										nodes.push_back(last_node);
										string updated_current;
										for(size_t update =first_index_after_parent; update< current.size(); update ++){
											updated_current += current.at(update);
										}
										nodes.at(node_index+1)=updated_current;
										if(biggest_child > 0){
											for(size_t shift = nodes.size()-2; shift > biggest_child;shift--){
												pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												vector<size_t> counter;
												multimap<size_t,size_t> intermediate;
												for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														cout<<"child in map: "<<it1->second<<endl;
													}
												}
												cout<< "count: "<<counter.size() <<endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){																									intermediate.insert(make_pair(shift+1, counter.at(in)+1));					
														cout << "relation in problematic part"<<endl;
	 													cout<< shift+1 << " " << counter.at(in)+1 << endl;
													}
													for(multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
														nodes_relation.insert(make_pair(it->first, it->second));
													}
												}
											}
											nodes_relation.insert(make_pair(node_index,biggest_child+1));
											cout<< "node index: "<< node_index << " biggest_child+1 "<<biggest_child+1 << endl;
											for(size_t m = 0; m < nodes.size(); m++){
												string node = nodes.at(m);
												cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													cout << int(node.at(l)) << " " ;
												}
											}
											cout << " " << endl;	
										}else{//get sure you creat a node with # later!
											for(size_t shift = node_index+1; shift <nodes.size(); shift++){
												map<string, size_t>::iterator it = firstParent.find(nodes.at(shift));
												if(it != firstParent.end() && it->second == shift){
													it->second = it->second+1;
												}
											}
											string last_node = nodes.at(nodes.size()-1);
											for(size_t shift =nodes.size()-1 ; shift > node_index+1; shift --){
												nodes.at(shift)=nodes.at(shift-1);
											}
											nodes.push_back(last_node);
											nodes.at(node_index+1)="#";
											for(size_t shift = nodes.size()-2; shift >=node_index+1;shift--){
												pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												vector<size_t> counter;
												multimap<size_t,size_t> intermediate;
												for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														cout<<"child in map: "<<it1->second<<endl;
													}
												}
												cout<< "count: "<<counter.size() <<endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){																									intermediate.insert(make_pair(shift+2, counter.at(in)+2));					
														cout << "relation in problematic part"<<endl;
	 													cout<< shift+2 << " " << counter.at(in)+2 << endl;
													}
													for(multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
														nodes_relation.insert(make_pair(it->first, it->second));
													}
												}
											}
											nodes_relation.insert(make_pair(node_index,node_index+1));
											nodes_relation.insert(make_pair(node_index,node_index+2));
											cout<< "node index: "<< node_index << " node index+1 "<<node_index+1<< " node index +2  " << node_index +2 << endl;				
										}
										current = updated_current;
										break;
										//need to break while loop here!inam check kon!
									}
								}else{//when we need to break the current parent and make a new branch with extra context, make current kids kids of this new child(**)
									num_kid = 1;
									cout<< "problem in else"<<endl;
									string updated_current;
									cout<< "updated_current: ";
									for(size_t update =first_index_after_parent; update< current.size(); update ++){
										updated_current += current.at(update);
										cout<< int(current.at(update));
									}
									cout<< " " <<endl;
									string extra_part_of_parent;
									cout<< "extra part of parent: ";
									for(size_t update = first_index_after_parent; update<current_parent.size(); update ++){
										extra_part_of_parent += current_parent.at(update);
										cout<< int(current_parent.at(update));
										
									}
									cout << " " << endl;
									cout<< "biggest_child: "<<biggest_child<<endl;
									for(size_t shift = node_index+1; shift <nodes.size(); shift++){
										map<string, size_t>::iterator it = firstParent.find(nodes.at(shift));
										if(it != firstParent.end() && it->second == shift){
											it->second = shift+2;
										}
									}
									if(biggest_child > 0){
									//	for(size_t child_node = node_index+1; child_node <= biggest_child;child_node ++)
									//		nodes.at(child_node) = extra_part_of_parent + nodes.at(child_node);//in ghalate!
									//	
										string last_node = nodes.at(nodes.size()-1);
										string second_last_node = nodes.at(nodes.size()-2);
										for(size_t shift =nodes.size()-1 ; shift >= node_index+1; shift --){
											nodes.at(shift)=nodes.at(shift-2);
										}
										if(current_parent == first_parent.at(j)){
											first_parent.at(j) = common_part;
											map<string, size_t>::iterator it1 = firstParent.find(current_parent);
											if(it1 != firstParent.end()){
												firstParent.erase(it1);
										//	
												map<string, size_t>::iterator it = firstParent.find(common_part);
												if(it == firstParent.end()){
													firstParent.insert(make_pair(common_part,node_index));
												}else{
													it->second = node_index;
												}
											}
										}
										for(size_t shift = nodes.size()-1; shift >= node_index;shift--){
											pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
											vector<size_t> counter;
											multimap<size_t,size_t> intermediate;
											for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												if(it1 != nodes_relation.end()){
													counter.push_back(it1->second);
													cout<<"child in map: "<<it1->second<<endl;
												}
											}
											cout<< "count: "<<counter.size() <<endl;
											if(counter.size() > 0){
												nodes_relation.erase(shift);
												for(size_t in = 0; in < counter.size(); in++){
													intermediate.insert(make_pair(shift+2, counter.at(in)+2));//First insert them to an intermadiate map
													cout << "relation in problematic part1"<<endl;
	 												cout<< shift+2 << " " << counter.at(in)+2 << endl;
												}
												for(multimap<size_t,size_t>::iterator it = intermediate.begin();it != intermediate.end(); it++){
													nodes_relation.insert(make_pair(it->first, it->second));
												}
											}
										}
										nodes_relation.insert(make_pair(node_index,node_index+1));
										nodes_relation.insert(make_pair(node_index,node_index+2));
										cout<< "node index: "<< node_index << " node index+1 "<<node_index+1<< " node index +2  " << node_index +2 << endl;
										nodes.push_back(second_last_node);
										nodes.push_back(last_node);
										nodes.at(node_index+2)= extra_part_of_parent;
										nodes.at(node_index+1)=updated_current;
										nodes.at(node_index)=common_part;
				
									}else{
										string last_node = nodes.at(nodes.size()-1);
										string second_last_node = nodes.at(nodes.size()-2);
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
											map<string, size_t>::iterator it1 = firstParent.find(current_parent);
											if(it1 != firstParent.end()){
												firstParent.erase(it1);
											}
											map<string, size_t>::iterator it = firstParent.find(common_part);
											if(it == firstParent.end()){
												firstParent.insert(make_pair(common_part,node_index));
											}else{
												it->second = node_index;
											}
										}
										for(size_t shift = nodes.size()-3; shift > node_index;shift--){
											pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
											vector<size_t> counter;
											multimap<size_t, size_t> intermediate;
											for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
												if(it1 != nodes_relation.end()){
													counter.push_back(it1->second);
												}
											}
											cout<< "count: "<<counter.size() <<endl;
											if(counter.size() > 0){
												nodes_relation.erase(shift);
												for(size_t in = 0; in < counter.size(); in++){
													intermediate.insert(make_pair(shift+2, counter.at(in)+2));							
													cout << "relation in problematic part2"<<endl;
	 												cout<< shift+2 << " " << counter.at(in)+2 << endl;
												}
												for(multimap<size_t, size_t>::iterator it = intermediate.begin(); it != intermediate.end();it++){
													nodes_relation.insert(make_pair(it->first,it->second));
												}
											}
										}
										nodes_relation.insert(make_pair(node_index,node_index+1));
										nodes_relation.insert(make_pair(node_index,node_index+2));
										for(size_t m = 0; m < nodes.size(); m++){
											string node = nodes.at(m);
											cout << "nodes at " << m << " : ";
											for(size_t l =0; l < node.size();l++){
												cout << int(node.at(l)) << " " ;
											}
										}
										cout << " " << endl;
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
								cout << "node index" << node_index << " , "<< parentInd<<endl;// are they the same? No, node index is from updated parent!
							}//Here is end of while loop!
							cout << "currentsize "<<current.size()<<endl;
							cout<< "current parent size: "<<current_parent.size()<<endl;
							cout<<"num kid: "<<num_kid<<endl;
							last_common_index =0;
							string commonPart;
							if(num_kid ==0){
								for(size_t k = 0 ; k < current.size(); k ++){
									cout<< "parent at " << k << " is " << int(current_parent.at(k))<< " current at " << k << " is " <<int(current.at(k))<<endl;
									if(current.at(k)==current_parent.at(k)){
										last_common_index += 1;
										commonPart += current.at(k);
									}else break;
								}
								cout<<"last common index "<<last_common_index<<endl;
								if(last_common_index != 0){
									cout << "nodes size: " << nodes.size() << "node index" << node_index << endl;
									nodes.at(node_index)= commonPart;
									cout<< "common: " << endl;
									for(size_t c=0; c< commonPart.size(); c++){
										cout<< int(commonPart.at(c))<< " " ;
									}
									cout << " " <<endl;
									if(current_parent == first_parent.at(j)){
										first_parent.at(j) = commonPart;
										map<string, size_t>::iterator it1 = firstParent.find(current_parent);
										if(it1 != firstParent.end()){
											firstParent.erase(it1);
										}
										map<string, size_t>::iterator it = firstParent.find(commonPart);
										if(it == firstParent.end()){
											firstParent.insert(make_pair(commonPart,node_index));
										}else{
											it->second = node_index;
										}
									}
									string other_branch;
									if((commonPart == current&&current_parent.size()>current.size())||commonPart.size() < current.size()){
									//	bool new_relation =false;
										if(commonPart == current&&current_parent.size()>current.size()){
											cout << "all the current is on the parent!"<<endl;
											for(size_t k = last_common_index; k < current_parent.size();k++){
												other_branch += current_parent.at(k);
											}
											cout<< "other1: " << endl;
											for(size_t c=0; c< other_branch.size(); c++){
												cout<< int(other_branch.at(c))<< " " ;
											}
											cout << " " <<endl;	
										//	new_relation = true;
										}
										if(commonPart.size() < current.size()){
											for(size_t k = last_common_index; k < current_parent.size();k++){
												other_branch += current_parent.at(k);
											}
											cout<< "other2: " << endl;
											for(size_t c=0; c< other_branch.size(); c++){
												cout<< int(other_branch.at(c))<< " " ;
											}
											cout << " " <<endl;
										//	new_relation = true;
										}	
										if(node_index == nodes.size()-1){
											cout << "when we are at if"<<endl;
											nodes.push_back(other_branch);
											if(i+last_common_index <= suffixes.size()-1){
												nodes.push_back(suffixes.at(i+last_common_index));
											}else{
												nodes.push_back("#");
											}
											nodes_relation.insert(make_pair(node_index,nodes.size()-1));
											nodes_relation.insert(make_pair(node_index,nodes.size()-2));
											cout<< node_index << " " << nodes.size()-1 << " " << nodes.size()-2 << " " << endl;
											for(size_t m = 0; m < nodes.size(); m++){
												string node = nodes.at(m);
												cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													cout << int(node.at(l)) << " " ;
												}
											}
											cout << " " << endl;	
										}else{// all the nodes in branch counter should be shifted by two, worth it to creat a new function that does that
											cout<< "if we are at else"<<endl;
											string last_node = nodes.at(nodes.size()-1);
											string second_last_node = nodes.at(nodes.size()-2);
										//	find_next_firstparent(node_index,next_first_parent);//If there is no first parent after that it stays 0
											//here is the place to add that shifting function
										//	shift_paths(next_first_parent,2);
											for(size_t shift = node_index+1; shift <nodes.size(); shift++){//inja!
												map<string, size_t>::iterator it = firstParent.find(nodes.at(shift));
												if(it != firstParent.end() && it->second == shift){
													it->second = it->second+2;
													cout<< "it->second + 2 "<<it->second+2 << " " << int(it->first.at(0))<<endl;
												}
											}
											pair < multimap<size_t,size_t>::iterator, multimap<size_t,size_t>::iterator > check_kids = nodes_relation.equal_range(node_index);
											size_t kid_count = 0;
										//	if(new_relation==true){
												for(multimap<size_t,size_t>::iterator it1 = check_kids.first ; it1 != check_kids.second; it1++){
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
												cout<< " i "<< i << " i+ last_common_index "<< i+last_common_index << "commonIndex "<< commonIndex <<endl;
												nodes.at(node_index+1) = suffixes.at(i+last_common_index+commonIndex);
											}else{
												nodes.at(node_index+1) = '#';
												cout<< "node at " << node_index + 1<< " is "  << nodes.at(node_index+1) <<endl;
											}
											for(size_t shift = nodes.size()-2; shift >= node_index+1;shift--){
												cout << "shift is: "<<shift<<endl;
												pair<multimap<size_t,size_t>::iterator , multimap<size_t,size_t>::iterator > p1 = nodes_relation.equal_range(shift);
												vector<size_t> counter;
												multimap<size_t,size_t>intermediate;
												for(multimap<size_t,size_t>::iterator it1 = p1.first ; it1 != p1.second; it1++){
													if(it1 != nodes_relation.end()){
														counter.push_back(it1->second);
														cout<< "parent "<< shift << " , "<< it1->first<<"child in map: "<<it1->second<<endl;
													}
												}
												cout<< "count: "<<counter.size() <<endl;
												if(counter.size() > 0){
													nodes_relation.erase(shift);
													for(size_t in = 0; in < counter.size(); in++){
														intermediate.insert(make_pair(shift+2, counter.at(in)+2));							
														cout << "er"<<endl;
	 													cout<< shift+2 << " " << counter.at(in)+2 << endl;
													}
													for(multimap<size_t, size_t>::iterator it =intermediate.begin();it !=intermediate.end();it++){
														nodes_relation.insert(make_pair(it->first,it->second));
													}
												}
											}
											cout<< "nodeindex: "<< node_index << endl;//always check before inserting
											nodes_relation.insert(make_pair(node_index, node_index+1));
											nodes_relation.insert(make_pair(node_index, node_index+2));
											cout << "node index +1 " << node_index+1 << " node index +2 "<< node_index+2 <<endl;
											if(kid_count> 0){
												cout<< "has sub kids! "<<endl;
												for(size_t i = 0;  i < kid_count; i++){
													nodes_relation.insert(make_pair(node_index+2, node_index+3+i));
													cout << "node index +2 " << node_index+2 << " node index +3 +i"<< node_index+3 + i <<endl;
												}
											}
											for(size_t m = 0; m < nodes.size(); m++){
												string node = nodes.at(m);
												cout << "nodes at " << m << " : ";
												for(size_t l =0; l < node.size();l++){
													cout << int(node.at(l)) << " " ;
													}
											}
											cout << " " << endl;		
										}
									}	
								}
							}
							break;
						}
					}
				}
				if(temp==active_node){//it is used only for adding a new branch to the root.
					cout<< " node size: " << nodes.size() << " current: " << endl;
					for(size_t l =0; l < current.size(); l++){
						cout<< int(current.at(l))<< " ";
					}
					cout << " " << endl;
					first_parent.push_back(current);
					nodes.push_back(current);
				//	vector<size_t> path;
				//	path.push_back(nodes.size()-1);
				//	branch_counter.insert(make_pair(path,1));
					map<string, size_t>::iterator it = firstParent.find(current);
					if(it == firstParent.end()){
						firstParent.insert(make_pair(current,nodes.size()-1));
					}else{
						it->second = nodes.size()-1;
					}
				}
			}
		}
		cout<<"first parents: "<<endl;
		for(size_t k=0; k < first_parent.size(); k++){
			for(size_t i =0; i< first_parent.at(k).size();i++){
				cout<< int(first_parent.at(k).at(i))<< " ";
			}
			cout << " , " ;
		}
		cout<< " "<<endl;
		cout << " first parent map: "<<endl;
		for(map<string, size_t>::iterator it = firstParent.begin();it != firstParent.end(); it++){
			string parent = it->first;
			for(size_t i =0; i< parent.size();i++){
				cout<< int(parent.at(i))<< " ";
			}
			cout << " , " << it->second << endl;
		}
		cout<< "all the nodes:"<<endl;
		for(size_t i = 0; i < nodes.size(); i++){
			string node = nodes.at(i);
			cout << "node at " << i << " is ";
			for(size_t j =0; j < node.size();j++){
				cout << int(node.at(j)) << " " ;
			}
			cout << " " <<endl;
		}
		cout<<"nodes relation: "<<endl;
		for(multimap<size_t , size_t>::iterator it = nodes_relation.begin(); it != nodes_relation.end(); it++){
				cout << it->first <<" "<< it->second << endl;
		}
	}
	void suffix_tree::count_paths(){
		for(size_t seq =0; seq < data.numSequences(); seq++){
			create_suffix(seq);
			for(size_t i = 0; i < suffixes.size(); i++){
				string current = suffixes.at(i);
				vector<size_t> branch;//insert this vector to the branch _counter map
				for(map<string,size_t>::iterator it = firstParent.begin();it!=firstParent.end();it++){
					string first_parent = it->first;
					if(first_parent.at(0)==current.at(0)){
						string updated_current;
						if(first_parent.size()==current.size()){

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
							} else break;
						}
						break;
					}
				}
				for(size_t i = 0; i < branch.size(); i ++){
					vector<size_t>sub_branch;
					cout<< " push back to sub branch: " <<endl;
					for(size_t j = i; j< branch.size();j++){
						sub_branch.push_back(branch.at(j));
						cout<< branch.at(j)<<endl;
						map<vector<size_t>,size_t>::iterator it = branch_counter.find(sub_branch);
						if(it == branch_counter.end()){
							branch_counter.insert(make_pair(sub_branch,0));
							it = branch_counter.find(sub_branch);
						}
						it->second = it->second +1;
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
	vector<string> suffix_tree::get_nodes()const{
		return nodes;
	}
	map<vector<size_t>, size_t> suffix_tree::get_count()const{
		return branch_counter;
	}
	vector<size_t> suffix_tree::get_first_parent()const{
		vector<size_t> first_parent;
		for(map<string,size_t>::const_iterator it = firstParent.begin(); it != firstParent.end(); it++){
			first_parent.push_back(it->second);
		}
		sort(first_parent.begin(),first_parent.end());	
		cout<<"first_parent: "<< endl;
		for(size_t i =0 ; i < first_parent.size(); i++){
			cout<< first_parent.at(i)<< " " ;
		}
		cout << " " <<endl;
		return first_parent;
	}
	merging_centers::merging_centers(finding_center & cent , suffix_tree & t):centers(cent),tree(t){}
	merging_centers::~merging_centers(){}
	void merging_centers::merg_value(){
		vector<string> nodes = tree.get_nodes();
		cout<< "size: " << nodes.size()<<endl;
		map<vector<size_t> , size_t> counts = tree.get_count();
		map<vector<size_t>, size_t> tree_gains;	// just enter the updated gain values
		cout<< "first parent size: "<< tree.get_first_parent().size()<<endl;
		for(size_t i =0; i < tree.get_first_parent().size();i++){
			cout << "i "<< i <<endl;
			size_t current_parent = tree.get_first_parent().at(i);
			map<vector<size_t>, size_t> subtree_gains;//old gain values
			for(map<vector<size_t> , size_t>::iterator it = counts.begin(); it != counts.end(); it++){
				vector<size_t> br = it->first;
				if(br.at(0)== current_parent){
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
					if(gain > 0){
						subtree_gains.insert(make_pair(br,gain));
						cout<< "gain "<<gain <<endl;
					}else{//branches with negative gain values, their gain may change later to a positive one, so I better insert them to the map as well.
						subtree_gains.insert(make_pair(br,0));
						cout<< "gain "<<gain <<endl;	
					}
				}
			}
			size_t highest_gain=0;
			vector<size_t>highest_path;
			for(map<vector<size_t>, size_t>::iterator it = subtree_gains.begin(); it != subtree_gains.end(); it++){
					if(it->second > highest_gain){
						highest_gain = it->second;
						highest_path = it->first;
					}else continue;
			}
			if(highest_gain > 0){
				cout<< "highest path"<<endl;
				for(size_t j = 0; j < highest_path.size(); j++){
					cout<< highest_path.at(j);
				}
				cout<< " " <<endl;
				map<vector<size_t>, size_t>::iterator it = counts.find(highest_path);
				assert(it != counts.end());
				size_t number_of_merged_centers = 0;
				for(size_t j = 0 ; j < it->first.size(); j++){
					number_of_merged_centers += nodes.at(it->first.at(j)).size();
				}
				for(map<vector<size_t>, size_t>::iterator it = subtree_gains.begin(); it != subtree_gains.end(); it++){
					map<vector<size_t>,size_t>::iterator it1 = counts.find(it->first);
					if(it->first.size()>highest_path.size()){
						size_t number = it1->second;
						string centers;
						size_t gain = 0;
						for(size_t j = 0 ; j < it1->first.size(); j++){
							centers += nodes.at(it1->first.at(j));
						}
						gain = number*(centers.size())-(number+centers.size()-number_of_merged_centers);
						tree_gains.insert(make_pair(it1->first,gain));
					}else{
						tree_gains.insert(make_pair(highest_path,highest_gain));
					}
				}
			}			
		}
		for(map<vector<size_t>,size_t>::iterator it = tree_gains.begin(); it != tree_gains.end(); it++){
			cout<< "path"<<endl;
			merged_centers.push_back(it->first);
			for(size_t i = 0 ; i < it->first.size(); i++){
				cout<< it->first.at(i)<< " ";
			}
			cout << " "<<endl;
			cout<< " gain " << it->second << endl;
		}
	}

	void merging_centers::merg_alignments(map<string, vector<pw_alignment> > & al_of_a_ccs){
		vector<string> center_index;
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
		}
	}



	

#endif
