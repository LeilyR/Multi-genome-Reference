#ifndef MODEL_HPP
#define MODEL_HPP

#include "pw_alignment.hpp"
#include"data.hpp"
#include "dynamic_mc.hpp"
#include"overlap.hpp"


//Interval tree:
#include "IntervalTree.hpp"
#include "alignment_index.hpp"

#include <map>
#include <vector>
#include <cassert>
#include <omp.h>

extern "C" {
#include "apcluster.h"
}


typedef  std::set<const pw_alignment*, compare_pw_alignment> alset;


template<typename T, typename overlap_type>
class initial_alignments_from_groups {
	public:
	initial_alignments_from_groups(const all_data & data, const std::vector< std::vector<const pw_alignment *> > & alignment_groups, const T & a_model): data(data), model(a_model){
		std::multimap< double, std::vector<const pw_alignment *> > sorter;
		max_gain = 0;
		for(size_t i=0; i<alignment_groups.size(); ++i) {
			double sumgain = 0;
			for(size_t j=0; j<alignment_groups.at(i).size(); ++j) {
				double gain1, gain2;
				model.gain_function(*(alignment_groups.at(i).at(j)), gain1, gain2);
				double vgain = (gain1+gain2)/2.0;
				sumgain+=vgain;
			
			}
			sorter.insert(std::make_pair(sumgain, alignment_groups.at(i)));
			max_gain+=sumgain;
		}
		for(typename std::multimap< double, std::vector<const pw_alignment *> >::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			sorted_original_als.push_back(rit->second);
		}

	}

	void compute(overlap_type & o, std::string ihead, std::string & info);

	private:
	const all_data & data;
	const T & model;
	std::vector<std::vector< const pw_alignment *> > sorted_original_als; // highest gain first 
	double max_gain;
	double result_gain;
	size_t used_alignments;
	mutable double select_groups_time;

	void lazy_split_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const std::vector<pw_alignment>  & als, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain);
	void lazy_split_full_insert_step(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment  & alin, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain);
	void local_undo(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment > & all_inserted, std::set< pw_alignment , compare_pw_alignment> & all_removed);
	void insert_alignment_sets(std::set<pw_alignment, compare_pw_alignment> & all_ins, std::set<pw_alignment, compare_pw_alignment> & all_rem, std::vector<pw_alignment> & this_ins, std::vector<pw_alignment> & this_rem);
	void all_push_back(std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment > & removed_alignments, std::set<pw_alignment , compare_pw_alignment> & all_inserted, std::set<pw_alignment, compare_pw_alignment> & all_removed );
	void select_from_groups(const std::vector<std::vector<pw_alignment> > & insert_als, const std::vector<std::vector<pw_alignment> > & rem_als_per_ins, std::vector<pw_alignment> & result_ins, std::vector<pw_alignment> & result_rem, double & group_gain) const;
	bool try_rnodes_unselect(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> insert_gain, const std::vector<double> & remove_gain, const std::set<size_t> & rem_to_unselect, vector<bool> & ins_selected, vector<bool> & rem_selected, double & new_gain) const;

	void bb_step(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, const std::vector<bool> & rem_free, const std::vector<bool> & ins_taken, const std::vector<bool> & rem_taken, const double & real_gain, double & best_gain, std::vector<bool> & best_ins, std::vector<bool> & best_rem) const;

	void take_highest_possible(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & gain) const;

	void easy_insert(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, 
		std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & extra_gain) const;

	void ub_gain(const std::vector<std::set<size_t> > & ins_to_rem, const std::vector<std::set<size_t> > & rem_to_ins, const std::vector<double> & insert_gain, const std::vector<double> & remove_gain, 
std::vector<bool> & rem_free, std::vector<bool> & ins_taken, std::vector<bool> & rem_taken, double & ub) const;

};



template<typename T, typename overlap_type>
class initial_alignment_set {
	public:

	initial_alignment_set(const all_data & d, const T & a_model, double base_cost): data(d), common_model(a_model){
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment&> sorter;
		double sumgain = 0;
//		std::cout << "cur "<<std::endl;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment & cur = data.getAlignment(i);
//			cur.print();
			double gain1, gain2;
			common_model.gain_function(cur, gain1, gain2);
			double vgain = (gain1+gain2)/2 - base_cost;
		//	std::cout << " al " << i << " gain1 " << gain1 << std::endl;
			if(vgain > 0.0) {
				sorter.insert(std::make_pair(vgain, cur));
				sumgain+=vgain;
			}
		}

//		sorted_original_als = std::vector<pw_alignment >(sorter.size(), NULL);
//		size_t pos = 0;
		for(std::multimap<double, const pw_alignment &>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment & alit = rit->second;
		//	std::cout << " ral " << alit << std::endl;
		//	alit.print();
		//	std::cout << std::endl;
			sorted_original_als.push_back(&alit);
//			pos++;
		}
		max_gain = sumgain;
	//	assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}
	initial_alignment_set(const all_data & d, const std::set< pw_alignment, compare_pw_alignment> & als, const T & a_model, double base_cost, size_t & num_threads): data(d), common_model(a_model){
		this->base_cost = base_cost;
		this->num_threads = num_threads;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(std::set<pw_alignment , compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); it++) { 
			const pw_alignment * cur = &*it;
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
			double vgain = (gain1+gain2)/2 - base_cost;
			assert(vgain >=0.0);
			sorter.insert(std::make_pair(vgain, cur));

			sumgain+=vgain;
		}

		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); rit++) {
			const pw_alignment * alit = rit->second;
	
			sorted_original_als.push_back(alit);
		}
		max_gain = sumgain;
	}
	initial_alignment_set(const all_data & d, const std::set< const pw_alignment*, compare_pointer_pw_alignment> & als, const T & a_model, double base_cost, size_t & num_threads): data(d), common_model(a_model){
		this->base_cost = base_cost;
		this->num_threads = num_threads;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
	//	std::cout << "calculating vgain! "<<std::endl;
		for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = als.begin(); it!=als.end(); it++) { 
			const pw_alignment * cur = *it;
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
			double vgain = (gain1+gain2)/2 - base_cost;
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al length " << cur->alignment_length() << " gain1 " << gain1 << " gain2 " << gain2 <<  std::endl;
		//	std::cout<< "vgain "<< vgain << std::endl;
			assert(vgain >=0.0);
		//	if(vgain>0.0){
			//	std::cout << " ins " << vgain << " at " << cur << std::endl;
			sorter.insert(std::make_pair(vgain, cur));

				/*
				common_model.gain_function(*cur, gain1, gain2);
				vgain = (gain1+gain2)/2 - base_cost;
				std::cout << " cached gian "<< vgain << std::endl;
				assert(vgain>=0);
				*/

		//	}
			sumgain+=vgain;
		}

	//	sorted_original_als = std::vector<pw_alignment & >(sorter.size()); 
	//	size_t pos = 0;
		for(std::multimap<double, const pw_alignment * >::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); rit++) {
			const pw_alignment * alit = rit->second;
		//	std::cout << " ins2 " << pos << " weight " << rit->first << " at " << alit <<  std::endl;
			
			/*
			double g1, g2;
			common_model.gain_function(*alit, g1, g2);
			double vgain = (g1+g2)/2 - base_cost;
			std::cout << " cached gain " << vgain << std::endl;
			assert(vgain >=0);
*/

			sorted_original_als.push_back(alit);
	//		pos++;
		}
		max_gain = sumgain;
	//	assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}
/*	initial_alignment_set(const all_data & d, std::multimap<double, const pw_alignment*>& sorter, const T & a_model, double base_cost, double& sumgain): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		for(std::multimap<double, const pw_alignment * >::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); rit++) {
			const pw_alignment * alit = rit->second;
			//XXX JUST A TEST! REMOVE IT LATER!
			
			double g1, g2;
			common_model.gain_function(*alit, g1, g2);
			double vgain = (g1+g2)/2 - base_cost;
			std::cout << " cached gain " << vgain << std::endl;
			assert(vgain >=0);
			////// END OF THE TEST
			sorted_original_als.push_back(alit);
		}
		max_gain = sumgain;		
	}*/

	~initial_alignment_set() {}
	void finding_repeats(std::set< const pw_alignment* , compare_pointer_pw_alignment> &, std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> >&, std::set< const pw_alignment* , compare_pointer_pw_alignment> &);

	void compute(overlap_type & o);
	void compute_simple(overlap_type & o);

/*
	This method starts with a weighted INDEPENDENT SET of alignments without partial overlap which is computed
	using a fast and easy to implement VERTEX COVER 2-approximation 
	(Clarkson, KL. A modification of the greedy algorithm for vertex cover. Information Processing Letters. 1983.)
*/
	void compute_vcover_clarkson(overlap_type & o);
	void find_als_weight(std::set<const pw_alignment*>&, std::set<size_t>&, std::vector<const pw_alignment*>&, size_t & ); //It weight the remaining alignments and order them based on their weight
	void compute_simple_lazy_splits(overlap_type & o);
	void lazy_split_insert_step
		(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment & al, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment> & removed_alignments, double & local_gain);
	void lazy_split_full_insert_step
		(overlap_type & ovrlp, size_t level, size_t & rec_calls, const pw_alignment & alin, std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment > & removed_alignments, double & local_gain);
	void all_push_back(std::vector<pw_alignment> & inserted_alignments, vector<pw_alignment > & removed_alignments, std::set<pw_alignment , compare_pw_alignment> & all_inserted, std::set<pw_alignment, compare_pw_alignment> & all_removed );
	void local_undo(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment > & all_inserted, std::set< pw_alignment , compare_pw_alignment> & all_removed);
	void insert_alignment_sets(overlap_type & ovrlp, std::set<pw_alignment, compare_pw_alignment> & all_ins, std::set<pw_alignment, compare_pw_alignment> & all_rem, std::vector<pw_alignment> & this_ins, std::vector<pw_alignment> & this_rem);



	
	double get_max_gain() const {
		return max_gain;
	}
	double get_result_gain() const {
		return result_gain;
	}

	size_t get_used_alignments() const {
		return used_alignments;
	}

	private:
	const all_data & data;
	const T & common_model;
	std::vector<const pw_alignment *> sorted_original_als; // highest gain first 
	double base_cost;
	size_t num_threads;
	double max_gain;
	double result_gain;
	size_t used_alignments;



	void update_remove(size_t index, std::vector<std::set<size_t> > & edges, std::map<size_t, double> & index_to_weight, std::multimap<double, size_t> & weight_to_index);
	void remove_val(size_t index, std::map<size_t, double> & index_to_weight,  std::multimap<double, size_t> & weight_to_index);
	void overlap_graph(std::vector<std::set<size_t> > & edges);
	void search_area(size_t from, size_t to, const std::multimap<size_t, size_t> & data, std::set<size_t> & res);
	
};


/*
template<typename T>
class initial_alignment_set_bb {
	public:
	initial_alignment_set(const all_data & d, const T & a_model, double base_cost,std::ofstream & outs): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * cur = &(data.getAlignment(i));
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2,outs);
			// TODO gain1 and gain2
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al " << i << " gain1 " << gain1 << std::endl;
			if(gain1-base_cost > 0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = std::vector<const pw_alignment*>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::revered from <E2><80><98>void initial_alignment_set<T>::compute(overlap&) [with T = mc_model]<E2><80><99>
main.cpp:rse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
		//	std::cout << " ral " << alit << std::endl;
		//	alit.print();
			std::cout << std::endl;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}
	initial_alignment_set(const all_data & d, const set< const pw_alignment *, compare_pw_alignment> & als, const T & a_model, double base_cost, std::ofstream & outs): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(std::set< const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * cur = *it;
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2,outs);
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al length " << cur->alignment_length() << " gain1 " << gain1 << " gain2 " << gain2 <<  std::endl;
			if(gain1 - base_cost>0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = std::vector<const pw_alignment *>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}

	~initial_alignment_set() {}

	void compute(overlap & o, std::ofstream &);
	void compute_simple(overlap & o,std::ofstream &);

	double get_max_gain() const {
		return max_gain;
	}
	double get_result_gain() const {
		return result_gain;
	}

	private:
	const all_data & data;
	const T & common_model;
	std::vector<const pw_alignment*> sorted_original_als; // highest gain first
	double base_cost;
	double max_gain;
	double result_gain;
};

*/
template<typename overlap_type>
class compute_cc_with_interval_tree{ //XXX this interval tree is designed by Erik Garisson
	public:
	compute_cc_with_interval_tree(const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, size_t num_sequences):intervals(num_sequences),trees(num_sequences){
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = *it;
			add_the_interval(al);
		//	std::cout << "on ref one: "<<std::endl;
			for(size_t j = 0; j < intervals.at(al->getreference1()).size();j++){
				Interval<const pw_alignment*> itv = intervals.at(al->getreference1()).at(j);
			//	itv.PrintInterval();
			}
		//	std::cout << "on ref two: "<<std::endl;
			for(size_t j = 0; j < intervals.at(al->getreference2()).size();j++){
				Interval<const pw_alignment*> itv = intervals.at(al->getreference2()).at(j);
			//	itv.PrintInterval();
			}

			alignments.push_back(al);
		}
		for(size_t i = 0; i < num_sequences;i++){
			IntervalTree<const pw_alignment*> tree;
			tree = IntervalTree<const pw_alignment*>(intervals.at(i));
			trees.at(i)=tree;		
		}
//		for(size_t i =0; i < num_sequences;i++){
//			std::cout<< "on seq: "<< i <<std::endl;
//			for(size_t j = 0; j < intervals.at(i).size();j++){
//				Interval<const pw_alignment*> itv = intervals.at(i).at(j);
//				itv.PrintInterval();
//			}
//		}
	}
	~compute_cc_with_interval_tree(){}
	void add_the_interval(const pw_alignment*);
	void compute(std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> > & ccs);
	void get_cc(const pw_alignment & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & );
	void cc_step(size_t , size_t , size_t , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment* , compare_pointer_pw_alignment>  & );



	private:
	
	std::vector<std::vector<Interval<const pw_alignment*> > >intervals;
	std::vector<IntervalTree<const pw_alignment*> > trees;
	std::vector<const pw_alignment*> alignments;
};
class compute_overlapped_reads_avl {
	public:
	compute_overlapped_reads_avl(const std::vector<pw_alignment> & als_in, size_t num_sequences, size_t num_threads):alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
	//	RIGHT = 0;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(size_t  i =0; i < als_in.size();i++){
			const pw_alignment * al = &als_in.at(i);
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::vector<const pw_alignment*> sorted_als(sorter.size());
		size_t num = 0;
		for(typename std::multimap<size_t, const pw_alignment*>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
		//	al->print();
		//	std::cout << "add one"<<std::endl;
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	}
	~compute_overlapped_reads_avl(){
		delete alind;
	}
	//This two are used for mapping reads against a graph
	void compute_on_second_ref(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs){
		std::cout << "compute CC on " << alignments.size() << std::endl;
		std::set <const pw_alignment*, compare_pointer_pw_alignment > seen;
		std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
		for(typename std::set<const pw_alignment *, compare_pointer_pw_alignment >::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
		//	std::cout << "compute_cc" <<std::endl;
			const pw_alignment * al = *it;
		//	al->print();
			typename std::set< const pw_alignment*, compare_pointer_pw_alignment >::iterator seenal = seen.find(al);
	//		std::cout << " seen " << seen.size() << std::endl;
			if(seenal == seen.end()) {
			//	std::cout << " getcc" << std::endl;
				std::set< const pw_alignment*, compare_pointer_pw_alignment > cc;
				get_cc_on_second_ref(*al, cc, seen);
				std::cout << "FOUND CC size " << cc.size() << std::endl;
			//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::const_iterator it = cc.begin(); it != cc.end();it++){
			//		const pw_alignment * test = *it;
			//		test->print();
			//	}
				sorter.insert(std::make_pair(cc.size(), cc));
			}	
		}
		for(typename std::multimap<size_t , std::set<const pw_alignment *, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
			ccs.push_back(it->second);
		}

	}
	void get_cc_on_second_ref(const pw_alignment & al , std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment *, compare_pointer_pw_alignment> & seen) {
		std::vector<size_t> left(2);
		std::vector<size_t> right(2);
		al.get_lr1(left.at(0), right.at(0));
		al.get_lr2(left.at(1), right.at(1));
		std::vector<size_t>reference(2);
		reference.at(0) = al.getreference1();
		reference.at(1) = al.getreference2();
		cc_step_on_second_ref(&al, reference.at(1), left.at(1), right.at(1), cc, seen);
	}
	void cc_step_on_second_ref(const pw_alignment * al,size_t current, size_t left, size_t right,std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc,std::set <const pw_alignment *,compare_pointer_pw_alignment>  & seen){
		seen.insert(al);
		cc.insert(al);
		std::vector<const pw_alignment *> result;
	
		alind->search_overlap(current, left, right,result);
		for(size_t i =0; i < result.size();i++){
			const pw_alignment * this_result = result.at(i);
			size_t ref = this_result->getreference2();
			size_t l,r;
			this_result->get_lr2(l,r);
			std::set <const pw_alignment* ,compare_pointer_pw_alignment>::iterator it=seen.find(this_result);
			if(it == seen.end()){
				cc_step_on_second_ref(this_result, ref, l, r, cc, seen);
			}
		}
	}

	private:
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments; 
	alignment_index  * alind;
	size_t num_threads;

};

template<typename overlap_type>
class compute_cc_avl {
	public:
	compute_cc_avl(const std::vector<pw_alignment> & als_in, size_t num_sequences, size_t num_threads):alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
	//	RIGHT = 0;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(size_t  i =0; i < als_in.size();i++){
			const pw_alignment * al = &als_in.at(i);
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::vector<const pw_alignment *> sorted_als(sorter.size());
		size_t num = 0;
		for(std::multimap<size_t, const pw_alignment *>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	}

	compute_cc_avl(const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, size_t num_sequences, size_t num_threads): alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
	//	RIGHT = 0;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = *it;
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::vector<const pw_alignment *> sorted_als(sorter.size());
		size_t num = 0;
		for(std::multimap<size_t, const pw_alignment *>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	}
	compute_cc_avl(const std::set<pw_alignment, compare_pw_alignment> & als_in, size_t num_sequences, size_t num_threads): alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
	//	RIGHT = 0;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(std::set< pw_alignment, compare_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = &*it;
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::vector<const pw_alignment *> sorted_als(sorter.size());
		size_t num = 0;
		for(std::multimap<size_t, const pw_alignment *>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	}

	compute_cc_avl(const overlap_type & ovrlp, size_t num_sequences, size_t num_threads){//Is used to seperate the fully overlapped ones from the others
		this->num_threads = num_threads;
		const std::set<pw_alignment, compare_alignment<pw_alignment> > & als = ovrlp.get_all();
		std::cout << " CC AVL build index on " << als.size() << std::endl << std::flush;

		std::vector<const pw_alignment *> vals(als.size());
		size_t num = 0;
		for(std::set<pw_alignment, compare_alignment<pw_alignment> >::const_iterator it = als.begin(); it!=als.end(); it++) {
			const pw_alignment & al = *it;
			vals.at(num) = &al;
			alignments.insert(&al);
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, vals);
		std::cout << " CC AVL created on " << alignments.size() << std::endl;
	}

	~compute_cc_avl(){
		delete alind;
	}
	void compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > &);
	void get_cc(const pw_alignment & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & );
	void cc_step(size_t , size_t , size_t , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment* , compare_pointer_pw_alignment>  & );
	private:
	
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments; 
	alignment_index  * alind;
	size_t num_threads;

};
template<typename overlap_type>
class compute_cc_with_icl{ //XXX Boost icl library
	public:
	compute_cc_with_icl(const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, size_t num_sequences, size_t num_threads):als_on_reference(num_sequences){
		this->num_threads = num_threads;
	//	RIGHT = 0;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = *it;
			alignments.push_back(al);
			add_on_intmap(al); //All the als are added to an interval map.
		}
		std::cout<< "map is filled! "<<std::endl;

	//	for(size_t i =0 ; i< num_sequences;i++){
		//	std::cout << "seq id "<< i << std::endl;
		//	for(boost::icl::interval_map<size_t, std::set<size_t> >::iterator it =als_on_reference.at(i).begin();it!=als_on_reference.at(i).end();it++){
			//	boost::icl::discrete_interval<size_t> itv  = (*it).first;
			//	all_intervals.at(i).insert(itv);
			//	std::set<size_t> overlaps_al = (*it).second;
			//	std::cout << "in interval " << itv << std::endl;
			//	std::cout << "size of overlaps_al "<< overlaps_al.size()<<std::endl;
			//	for(std::set<size_t>::iterator it1 = overlaps_al.begin(); it1 != overlaps_al.end(); it1++){
			//		std::cout << *it1 <<std::endl;
			//		size_t al_id = *it1;
			//		const pw_alignment * pw_al = alignments.at(*it1);
			//		pw_al->print();
			//	}
		//	}
	//	}
	}
	compute_cc_with_icl(const overlap_type & ovrlp, size_t num_sequences, size_t num_threads): als_on_reference(num_sequences){
		this->num_threads = num_threads;
		const std::set<pw_alignment, compare_alignment<pw_alignment> > & als = ovrlp.get_all();
		for(std::set<pw_alignment, compare_alignment<pw_alignment> >::const_iterator it = als.begin(); it!=als.end(); it++) {
			const pw_alignment & al = *it;
			alignments.push_back(&al);
			add_on_intmap(&al);
		}
		std::cout<< "map is filled! "<<std::endl;
	
	}

	compute_cc_with_icl(std::vector<std::pair<size_t,size_t> > & common_int){///XXX it was only for testing the library.
		for(size_t i =0; i < common_int.size();i++){
			boost::icl::discrete_interval<size_t> bounds = boost::icl::construct<boost::icl::discrete_interval<size_t> >(common_int.at(i).first,common_int.at(i).second);
			std::set<size_t> id;
			id.insert(i);
			als_on_ref.add(make_pair(bounds, id));	
		//	reverse_als_on_ref.add(make_pair(id,bounds));
			id_and_bounds.insert(make_pair(i,common_int.at(i)));
		}
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
	~compute_cc_with_icl(){}
	void add_on_intmap(const pw_alignment*);
	void remove_from_intmap(size_t & id, size_t & ref1, size_t & ref2, size_t& left1, size_t & right1, size_t & left2, size_t & right2);
	void compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > &);
	void compute_test(std::vector<std::set<size_t> >&);
	void get_cc(const pw_alignment & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & );
	void cc_step(size_t , size_t , size_t , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment* , compare_pointer_pw_alignment>  & );
	void cc_step_current(size_t & , size_t & , size_t & , std::set<size_t>& , std::set<size_t> & );
	private:
	std::vector<const pw_alignment*> alignments;
	std::vector<boost::icl::interval_map<size_t, std::set<size_t> >  > als_on_reference;
//	boost::icl::interval_map<std::set<size_t>,size_t >reverse_als_on_ref; //XXX It was used only for the simple library test.
	boost::icl::interval_map<size_t, std::set<size_t> >als_on_ref; //XXX It was used only for the simple library test.
	std::map<size_t , std::pair<size_t,size_t> > id_and_bounds;
//	std::vector< boost::icl::interval_set<size_t> >all_intervals; //All the intervals are kept here, while they will be removed gradually from the als_on_reference
	size_t num_threads;


};
template<typename overlap_type>
class compute_cc {
	public:
/*	compute_cc(const std::set<pw_alignment, compare_pw_alignment> & als_in, size_t num_sequences):alignments(als_in), als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		max_al_ref_length = 0;

		for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); ++it) {
			const pw_alignment & al = *it;
			add_on_mmaps(al);
		}
	
	}*/
	compute_cc(const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, size_t num_sequences, size_t num_threads):alignments(als_in), als_on_reference(num_sequences), last_pos(num_sequences, 0) {//XXX This is the one is currently used for creating CCs.
		max_al_ref_length = 0;
		numThreads = num_threads;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++) {
			const pw_alignment * al = *it;
			add_on_mmaps(al);
		}
//		std::cout << "al size " << alignments.size() << std::endl;
	}

	compute_cc(const overlap_type & ovrlp, size_t num_sequences, size_t num_threads): als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		numThreads = num_threads;
		max_al_ref_length = 0;
		const std::set<pw_alignment, compare_alignment<pw_alignment> > & als = ovrlp.get_all();
		for(std::set<pw_alignment, compare_alignment<pw_alignment> >::const_iterator it = als.begin(); it!=als.end(); it++) {
			const pw_alignment & al = *it;
			alignments.insert(&al);
			add_on_mmaps(&al);
		}
	
	}


	compute_cc(const all_data & dat);
	~compute_cc() {}

	void compute(std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> > & ccs); 

	private:
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments; // real objects here, references everywhere else 
	std::vector< std::multimap< size_t, const pw_alignment*> > als_on_reference; // sequence index -> pos on that sequence -> alignment reference
	std::vector<size_t> last_pos;
	size_t max_al_ref_length;
	size_t numThreads;
	void add_on_mmaps(const pw_alignment* pwa);
	void remove_on_mmaps(const pw_alignment* al);
	void get_cc(const pw_alignment & al, std::set <const pw_alignment* , compare_pointer_pw_alignment> & , std::set < const pw_alignment* , compare_pointer_pw_alignment> & );
	void cc_step(size_t ref, size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set < const pw_alignment*, compare_pointer_pw_alignment> & seen );

};

template<typename tmodel, typename overlap_type>
class clustering {
	public:
	clustering(overlap_type &,all_data &, tmodel&);
	~clustering();
	void als_on_reference(const pw_alignment * p);
	void calculate_similarity();//Fill in a matrix with gain values for each two accessions.
	void update_values(); 
	void update_clusters(size_t acc);
	void update_clusters();
	private:
	overlap_type & overl;
	all_data & data;
	tmodel & model;
//	std::set<const pw_alignment *, compare_pw_alignment> alignments;
	std::vector< std::multimap<size_t, pw_alignment *> > als_on_ref;
	std::vector<vector<vector<double> > >gain;//consider it as similarity function, though in your slides you mentioned modification can be the similarity function.
	std::vector<vector<vector<double> >	>ava;//availability		
	std::vector<vector<vector<double> >	>res;//responsibilty
	
};

/*
	TODO 
	can the clustering result become better if we estimate modification costs of those pairs which have no direct pairwise alignment


*/

template<typename tmodel,typename overlap_type>
class affpro_clusters {
	public:

	affpro_clusters(const overlap_type & ovlp, const tmodel & model, double base_cost):model(model), base_cost(base_cost) {

		const std::set<pw_alignment*, compare_pw_alignment> & als = ovlp.get_all();
		for(std::set<pw_alignment*, compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * al = *it;
			all_als.insert(al);
			add_alignment(*al);
		//	al->print();
		}

	}

	affpro_clusters(const std::set<const pw_alignment* , compare_pointer_pw_alignment> & inset, const tmodel & model, double base_cost):model(model), base_cost(base_cost) {
	//	std::cout<<"instd::set size"<<inset.size()<<std::endl;	
		//	std::cout<<"data1 ad in afp: "<< & dat << std::endl;	
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = inset.begin(); it!=inset.end(); ++it) {
		//	std::cout<<"data2 ad in afp: "<< & dat << std::endl;	
			const pw_alignment * al = *it;
		//	std::cout<<"alignment from instd::set: "<<std::endl;
		//	al->print();
		//	dat.numAcc();
		//	std::cout<<"data3 ad in afp: "<< & dat << std::endl;	
			all_als.insert(al);
			add_alignment(*al); 
		//	std::cout<<"data4 ad in afp: "<< & dat << std::endl;
		//	al->print();	
			double g1 ,g2;
			model.gain_function(*al,g1,g2);
			double av_gain = (g1+g2)/2 ;
			assert(av_gain > 0);
		//	assert(g1 > 0);
		//	assert(g2 > 0);


		}
	}


void run(std::map<std::string, std::vector<std::string> > & cluster_result, std::map<std::string, std::vector<pw_alignment> > & local_al_in_a_cluster) {
	// convert to c-style matrix
//	std::cout<< "simmatrix size is "<< simmatrix.size() << std::endl;
	double * data = new double[simmatrix.size() * simmatrix.size()];
	int * result = new int[simmatrix.size()];
	if(simmatrix.size()>2){
		for(size_t i=0; i<simmatrix.size(); ++i) {
			for(size_t j=0; j<simmatrix.size(); ++j) {
				data[i*simmatrix.size() + j] = simmatrix.at(i).at(j);
			//	std::cout<< "simmat at "<<i <<" at "<< j << " " << simmatrix.at(i).at(j)<<std::endl;
			}
		}
		APOPTIONS apoptions;
		apoptions.cbSize = sizeof(APOPTIONS);
		apoptions.lambda = 0.5;
		apoptions.minimum_iterations = 1000;
		apoptions.converge_iterations = 5000;
		apoptions.maximum_iterations = 50000;
		apoptions.nonoise = 1;
		apoptions.progress=NULL; apoptions.progressf=NULL;
		double netsim;	
		int iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
	//	std::cout << "iter " << iter << "result[0] "<< result[0] << " netsim "<< netsim<< std::endl;
		if(iter <= 0 || result[0]==-1) {
			apoptions.minimum_iterations = 10000;
			apoptions.converge_iterations = 15000;
			apoptions.maximum_iterations = 150000;
			apoptions.lambda = 0.6;
			iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
		//	std::cout << "iter " << iter << std::endl;
		}
		if(iter <= 0 || result[0]==-1) {
			apoptions.minimum_iterations = 100000;
			apoptions.converge_iterations = 20000;
			apoptions.maximum_iterations = 250000;
			apoptions.lambda = 0.99;
			iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
		//	std::cout << "iter1 " << iter << "result[0] " <<result[0]<< std::endl;
		}
	} else {
		if(simmatrix.size()==1) {
			result[0] = 0;
		//	std::cout << "result[0] " <<result[0]<< std::endl;

		} else {
		// simpler algorithm for only 2 element
			result[0] = 0;
			double one = simmatrix.at(0).at(0) + simmatrix.at(0).at(1);
			double two = simmatrix.at(1).at(1) + simmatrix.at(1).at(0);
			double separate = simmatrix.at(0).at(0) + simmatrix.at(1).at(1);
		//	std::cout << "one "<< one << " two "<< two << " sep "<< separate << std::endl;
			if(one > two && one > separate) {
				result[0] = 0;
				result[1] = 0;
			} else if(two > one && two > separate) {
				result[0] = 1;
				result[1] = 1;
			} else {
				result[0] = 0;
				result[1] = 1;
			}
		//	std::cout << "2 elements: " << "result[0] " <<result[0] << " result[1] "<< result[1] << std::endl;		
		}
	}
// TODO conv error res -1

	double totalccost = 0;
	for(size_t i=0; i<simmatrix.size(); ++i) {
	//	if(simmatrix.at(i).at(result[i])== -HUGE_VAL){
	//		totalccost = i;
	//	}else{
		totalccost -= simmatrix.at(i).at(i);
	//	}
	}
//	std::cout << "totalcost "<< totalccost<<std::endl;
	double apcost = 0;
//	std::cout<<"size of simmatrix: "<<simmatrix.size()<<std::endl;
	for(size_t i=0; i<simmatrix.size(); ++i) {
		if(result[i]==-1) {
			result[i] = i;
		}
		assert(result[i]>= 0);
		if(i==(size_t)result[i]) {
			apcost-=simmatrix.at(i).at(i);
		//	std::cout << "i is "<< i << std::endl;
		} else {
			if(simmatrix.at(i).at(result[i])== -HUGE_VAL){
				result[i]=i;
				apcost-=simmatrix.at(i).at(i);
			}else{
				apcost-=simmatrix.at(i).at(result[i]);
			}
		}
	}
//	std::cout << "simmatrix size "<< simmatrix.size()<<std::endl;
	for(size_t i=0; i<simmatrix.size(); ++i) {//cluster center may happen whenever i == result[i]
	//	std::cout << sequence_names.at(i) << " res " << i << " is " << result[i] << " ( length " << sequence_lengths.at(i) << ")"<<std::endl;

		if( (size_t)result[i]==i) {
		//	std::cout << " center !!" << simmatrix.at(i).at(i) << std::endl;
			std::map<std::string, std::vector<std::string> >::iterator it=cluster_result.find(sequence_names.at(i));
			if(it==cluster_result.end()){
				cluster_result.insert(make_pair(sequence_names.at(i), std::vector<std::string>()));
			} else {
	//			std::cout << "at else "<<std::endl;
				// it->second.push_back(sequence_names.at(i));
			
			}


		} else {
			/*1. find reverse seq_name if exists compare their gains. For that i might need two containers . One for saving memebers as soon as they happen(map(memebr,center))
			and the other one for getting an al from its sequence names(map(vector(seq1,seq2),al)). 
			 2. Maybe i can also check if the reverse of current center already exists as a member then remove that member. */


		//	std::cout << " " << simmatrix.at(i).at(result[i]) << " : " << simmatrix.at(i).at(i) << std::endl;
			//result[i]is associated one, add them to the map for each center
		/*	std::map<std::string, std::vector<std::string> >::iterator it=cluster_result.find(sequence_names.at(result[i])); //XXX Added it later
			if(it == cluster_result.end()) {
				cluster_result.insert(make_pair(sequence_names.at(result[i]), std::vector<std::string>()));
				it=cluster_result.find(sequence_names.at(result[i]));
			}*/
		//	std::string reverse_seq = make_reverse(sequence_names.at(i));//Reverse of the member
			//Make the al here:
			pw_alignment pal;
		//	std::cout << "center is "<< sequence_names.at(result[i])<<std::endl;
			make_an_alignment(sequence_names.at(result[i]),sequence_names.at(i),pal);//Using result[i] as center and i as member find the corresponding al.
		/*	double gain = model.get_the_gain(pal,sequence_names.at(result[i]));
			std::map<std::string,std::string>::iterator mem = members_of_clusters.find(reverse_seq);
			if(mem != members_of_clusters.end()){
				std::cout << "has reverse!"<<std::endl;
				std::string cent = mem->second;
				std::vector<std::string> temp;
				temp.push_back(cent);
				temp.push_back(reverse_seq);
				std::map<std::vector<std::string>, pw_alignment>::iterator al = alignments.find(temp);
				assert(al != alignments.end());
				pw_alignment p = al->second;
				double rev_gain = model.get_the_gain(p, cent);
				if(rev_gain < gain){
					//Remove it
					alignments.erase(al);
					members_of_clusters.erase(mem);
					//add the new one
					std::vector<std::string> temp1;
					temp1.push_back(sequence_names.at(result[i]));
					temp1.push_back(sequence_names.at(i));
					alignments.insert(std::make_pair(temp1,pal));
					members_of_clusters.insert(std::make_pair(sequence_names.at(i),sequence_names.at(result[i])));
				}else{//Do nothing
			
				}
			}else{
			//	it->second.push_back(sequence_names.at(i)); //XXX Added it later*/
			//add the new one
			std::vector<std::string> temp1;
			temp1.push_back(sequence_names.at(result[i]));
			temp1.push_back(sequence_names.at(i));
			alignments.insert(std::make_pair(temp1,pal));
		//	std::cout << sequence_names.at(result[i]) << " " << sequence_names.at(i) <<std::endl;
			members_of_clusters.insert(std::make_pair(sequence_names.at(i),sequence_names.at(result[i])));
		//	}
			
		}
	}
	//Checks for the redundancy:
//	check_redundancy();
	//Fills in cluster_result and  local_al_in_a_cluster
	write_clusters(cluster_result, local_al_in_a_cluster);
	
	// double apgain = totalccost - apcost;

//	std::cout << "Total sequence cost " << totalccost << " ap clustering cost " << apcost << " gain: " << apgain << std::endl;

	delete [] data;
	delete [] result;
	data = NULL;
	result = NULL;
//	std::cout<< "number of centers: " <<clusterCenter.size()<<std::endl;
}
	size_t get_sequence_length(size_t ref_idx)const; //ref_idx shows the reference that sequence belongs to. It could be either 0 or one.
	
	std::string make_reverse(std::string sequence_name){
		std::vector<std::string> seq_parts;
		strsep(sequence_name, ":" , seq_parts);
		unsigned int dir = atoi(seq_parts.at(0).c_str());
		unsigned int ref = atoi(seq_parts.at(1).c_str());
		unsigned int left = atoi(seq_parts.at(2).c_str());
		if(dir == 0){
			dir = 1;
		}else{
			dir = 0;
		}
		std::stringstream reverse;
		reverse<<dir<<":"<<ref<<":"<<left;
		std::string rev_seq = reverse.str();
		return rev_seq;
	}
	void make_an_alignment(std::string& , std::string & , pw_alignment &);
	void write_clusters(std::map<std::string, std::vector<std::string> > & cluster_result, std::map<std::string, std::vector<pw_alignment> > & local_al_in_a_cluster);
	void check_redundancy();
// simil = neg distance, diagonale pref = neg cost
	private:
//	all_data & dat;
	const tmodel & model;
	double base_cost;
	// TODO do we need distances for all pairs of sequence pieces?
	std::set<const pw_alignment* , compare_pointer_pw_alignment> all_als;
	std::map<std::string, size_t> sequence_pieces; // chr:leftpos -> seq_piece index
	std::vector<std::string> sequence_names;
       	std::vector<size_t> sequence_lengths;	
	std::vector<std::vector<double> > simmatrix;
	std::map<std::string, char> cluster_centers;
	void add_alignment(const pw_alignment  & al);

	std::map<std::string, std::string> members_of_clusters; //first string represents a member and the second string represents its coresponding center
	std::map<std::vector<std::string>,pw_alignment> alignments; //The set represents the name an al's references. (center,other)
};

#include "overlap.cpp"
#include "intervals.cpp"

/*class suffix_tree{
	public:
	suffix_tree(all_data &, finding_centers&);
	~suffix_tree();
	void create_suffix(size_t);
	void update_successive_centers(std::vector<size_t> &,size_t &, size_t & ); // Updating successive center vector by replacing a new merged center 
	void find_parent(size_t&, size_t&);
	void delete_relation(size_t & , size_t& );
	void count_branches();
	std::vector<std::vector<size_t> > get_nodes()const;
	std::map<std::vector<size_t>, size_t> get_count()const;//Returns 'branch counter'
	std::vector<size_t> get_first_parent()const;
	void get_first_parent(std::vector<size_t> & , std::vector<size_t> &);
	void check_first_nodes_after_root(std::vector<size_t> & current, size_t & node_index);
	void find_child_nodes(size_t &, std::vector<size_t> &);
	void first_parent_index(size_t & , size_t &); //current index, its first parent index
	void create_tree(std::vector<size_t> &, size_t &);
	size_t get_power_of_two(size_t &)const;
	void shift_node_relation(size_t & , size_t & , size_t & );
	void shift_first_parent(size_t & , size_t&, std::vector<size_t> &);
	std::map<size_t, std::vector<size_t> > get_center_on_a_sequence(size_t &)const;
	std::vector<std::vector<std::vector<size_t> > > get_suffix_of_sequence(size_t &)const;
	private:
	all_data & data;
	finding_centers & centers;
	std::vector<std::map<size_t,std::vector<size_t> > >successive_centers;//It is updated in each round of making a suffix tree 
	std::vector<std::vector<std::vector<std::vector<size_t> > > >suffixes;//all the suffixes of all the sequences at each around
	std::vector<std::vector<size_t> > nodes;
	std::vector<size_t> powerOfTwo;
	std::multimap<size_t , size_t> nodes_relation; //first size_t shows the parent and second one shows kids
	std::map<std::vector<size_t>, size_t>firstParent; // size_t shows its node index //TODO Can be removed and be replaced by FirstCenterOnFirstNodes and IndexOfFirstNodesAfterRoot
	std::map<std::vector<size_t>,size_t>branch_counter;//vector represents all nodes of a branch , size_t shows the number of that branch is happening
	std::map<size_t,size_t> FirstCenterOnFirstNodes; //size_t -> is the first node on the first row nodes after the root,  size_t -> is the the index of that node
	std::set<size_t> IndexOfFirstNodesAfterRoot; //Index of first row nodes after the root //TODO
};*/


#endif



