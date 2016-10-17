#ifndef OVERLAP_HPP
#define OVERLAP_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include <map>
#include <math.h> 
#include <algorithm>

#include "pw_alignment.hpp"
#include "data.hpp"
#include "alignment_index.hpp"

#include <boost/iostreams/stream.hpp>


class overlap_interval_tree{
public:
	overlap_interval_tree(const all_data&);
	overlap_interval_tree(const overlap_interval_tree & o);
	~overlap_interval_tree();
	void insert_without_partial_overlap(const pw_alignment & p);
	// This function removes an alignment with adress identity to remove from overlap, then deletes remove
	void remove_alignment(const pw_alignment & remove);

	void test_all() const;
	void test_all_part()const;// Change it later to check all those pieces with gap in one sample. Put them all in a std::set and check if the only missed parts of coverage are those parts.
	void test_overlap()const;
	void test_no_overlap_between_ccs(const pw_alignment &, std::set<const pw_alignment*, compare_pointer_pw_alignment> &)const;
	void print_all_alignment() const;
	const pw_alignment * get_al_at_left_end(size_t ref1, size_t ref2, size_t left1, size_t left2) const;
//	std::multimap<size_t, const pw_alignment &>& get_als_on_reference(size_t sequence) ;
//	const std::multimap<size_t, const pw_alignment & >& get_als_on_reference_const(size_t sequence) const ;
	void get_interval_result_on_interval(const size_t& sequence, const size_t& left, const size_t& right, std::vector<const pw_alignment *>&)const;
	void interval_result_on_point(const size_t & sequence, const size_t & split_point, std::vector<const pw_alignment* > & result)const;

	bool checkAlignments(const pw_alignment & p)const;

	const std::set<pw_alignment, compare_pw_alignment> & get_all() const;
	void test_partial_overlap() const;
//	void test_al_on_ref()const;
//	bool check_alignment_address(const pw_alignment & , const pw_alignment * )const;
	static void test_partial_overlap_set(std::set< const pw_alignment *, compare_pw_alignment> & als);
	static void test_partial_overlap_vec(std::vector< const pw_alignment *> & als);
	static bool check_po(size_t l1, size_t r1, size_t l2, size_t r2);

	size_t size() const;
	const all_data & data; // TODO private again?
private:
	std::set<pw_alignment, compare_pw_alignment> alignments;
//	std::vector< std::multimap< size_t, const pw_alignment &> > als_on_reference; // sequence index -> left and right pos on that sequence -> alignment reference
	alignment_index alind;

};

template<typename overlap_type>
class splitpoints_interval_tree {
	public:
	splitpoints_interval_tree(const pw_alignment & , overlap_type &, const all_data &);
	splitpoints_interval_tree(const std::vector<const pw_alignment * > & new_als, overlap_type &, const all_data &);
	~splitpoints_interval_tree();
	void find_initial_split_points(size_t sequence, size_t left, size_t right, bool recursion);
//	void find_initial_split_points_nonrecursive(size_t sequence, size_t left, size_t right);
	void recursive_splits();//initial split points + all induced split points everywhere
	void nonrecursive_splits(); // initial split points only 
	void insert_split_point(size_t sequence, size_t position, bool recursion);//recursively find the other split points
//	void insert_split_point_nonrecursive(size_t sequence, size_t position);//insert without recursion
	void split_all(std::set<pw_alignment, compare_pw_alignment> & remove_alignments, std::vector<pw_alignment> & insert_alignments);
	void splits(const pw_alignment & p,  std::vector<pw_alignment> & insert_alignments);
	bool onlyGapSample(const pw_alignment & p) const;
	void run_initial(bool recursion);

	// std::vector<pw_alignment>  get_insert () const;

	private:
	overlap_type & overl;
	pw_alignment newal;
	std::vector<const pw_alignment * > new_als;
	bool multi_mode;
	
	const all_data & data;
	std::vector<std::set<size_t> > split_points;
	std::vector<pw_alignment> insert_alignments;	

	void initial_splits_on_al(const pw_alignment * al, bool recursion);
	void split_into_alignment_ref(const pw_alignment * al, size_t al_ref, size_t position, bool recursion);
	void split_into_alignment(const pw_alignment * al, size_t sequence, size_t position, bool recursion);

/*
	initial points in std::sets (method with and without sets)
	compute remove and insert alignments
	cost change level 1, abort if no gain
	if there is gain recurse


*/
	
};




#endif
