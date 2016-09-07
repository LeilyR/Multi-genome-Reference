#ifndef FILTER_ALIGNMENTS_HPP
#define FILTER_ALIGNMENTS_HPP

#include "pw_alignment.hpp"
#include "dynamic_mc.hpp"
#include "alignment_index.hpp"
#include <typeinfo>
#include <iostream>
#include <cmath>

#define FRACTION 0.5
class filter_als{
	public:
	filter_als(dynamic_mc_model & m, std::set<const pw_alignment*,compare_pointer_pw_alignment> & als, size_t num_sequences, size_t num_threads): model(m), alind(NULL){
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = als.begin(); it != als.end(); it++){
			alignments.push_back(*it);
		}
		this ->num_threads = num_threads;
		this->num_sequences = num_sequences;
		alind = new alignment_index(num_sequences, num_threads, alignments);
	}
	~filter_als(){}
	void find_overlapped_references();
	void find_touched_ref(std::vector<const pw_alignment*> & , size_t & , size_t & , size_t & , std::map<const pw_alignment*, size_t> & );
	size_t get_overlap(std::map<const pw_alignment*, size_t> &  , std::map<const pw_alignment*,size_t> & );
	void recalculate_gain_value(const pw_alignment*, size_t &);
	void find_als_with_highest_gain();
	std::set<const pw_alignment*, compare_pointer_pw_alignment> get_filtered_als()const;
	private:
	dynamic_mc_model & model;
	size_t num_threads;
	size_t num_sequences;
	std::vector<const pw_alignment*> alignments;
	std::multimap<double,const pw_alignment*> new_gains;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> filtered_als;
	alignment_index  * alind;




};





class alignment_filter {
	public:



	private:


}


#endif
