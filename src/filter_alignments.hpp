#ifndef FILTER_ALIGNMENTS_HPP
#define FILTER_ALIGNMENTS_HPP

#include "pw_alignment.hpp"
#include "dynamic_mc.hpp"
#include "alignment_index.hpp"
#include "overlap.hpp"
#include "model.hpp"
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


/*

	Idea: Optimist gain
	How much gain could the current alignment generate in the ideal case

	2 concepts:

	1) Homology gain

	seq 1 alalalalal
		olololol


	seq 2 alalalalal
	      ilililil

	seq 3 olololol
	          ilililil

	On seq 1 al overlaps  with ol, on seq 2 il overlaps with al, on seq 3 the overlap of ol and il confirms al
	We intersect ol and il on seq 3, we compute which fractions of al are covered by the intersected parts of ol and il
	This score is added over all the homology confirming alignments and then multiplied by the gain of al


	2) Synteny gain

	aaaa    bbbb     XXXXXX    cccc
	aaaa    bbbb     XXXXXX    cccc

	Alignment xxx gets confirmed by a b and c because they are on the same sequences and confirm the order given by x (take care of alignment directions)
	We assume that x allows for a, b, c, and x to be included in a long center. We estimate how much additional gain per original gain we get on average from using long centers.


	Both optimist gains can quite large on large data sets. To not overdo it, we scale them down. Each alignment gets its original gain doubled without restriction. Furthermore, 
	we increase total original gain at most threefold by linearly sclaling down the remaining optimist gain

*/
class alignment_filter {
	public:
	alignment_filter(dynamic_mc_model & m, const all_data & data, std::set<const pw_alignment*, compare_pointer_pw_alignment> & alignments, const size_t & gap_in_long_centers, const size_t & num_threads) : model(m), data(data), alignments(alignments), alvec(alignments.size()), extra_gain(alignments.size(), 0), gap_in_long_centers(gap_in_long_centers), num_threads(num_threads) {
		size_t ali = 0;
		for( std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
			alvec.at(ali) = *it;
			ali++;
		}


		estimate_optimist_gain_per_gain();

		alind = new alignment_index(data.numSequences(), num_threads, alvec);

		homology_gain();
		synteny_gain();
		filter_als();
	}		
	

	~alignment_filter() {
		delete alind;
	}



	private:
	dynamic_mc_model & model;
	const all_data & data;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> & alignments;
	double optimist_gain_per_gain;
	std::vector<const pw_alignment * > alvec;
	std::vector<double> extra_gain;
	alignment_index * alind;
	const size_t & num_threads;
	size_t gap_in_long_centers;
	double sum_gain, sum_hom, sum_syn;
	
	void overlap_other_side(const pw_alignment * al, size_t ref, const pw_alignment * other_al, size_t & other_ref_of_other_al);
	void estimate_optimist_gain_per_gain();
	double homology_score(const pw_alignment * al, const pw_alignment * al1, const pw_alignment * al2, size_t al1ref, size_t al2ref);
	double best_path_after(const size_t & r1pos, const size_t & r2pos, bool r1forward, bool r2forward, const size_t & rangel, const size_t & ranger, const avl_interval_tree<std::pair<size_t, size_t> > & inrange, std::vector<std::pair< size_t, size_t> > & best_path);
	double synteny_score(size_t alind, const size_t & alref, const size_t & rangel, const size_t & ranger, const avl_interval_tree<std::pair<size_t, size_t> > & inrange);
	void homology_gain();
	void synteny_gain();
	void filter_als();




};


#endif
