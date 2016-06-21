#ifndef ALIGNMENT_INDEX_HPP
#define ALIGNMENT_INDEX_HPP

#include "pw_alignment.hpp"
#include "intervals.hpp"
#include <vector>
#include <set>
#include <omp.h>


class alignment_index {
	public:
	alignment_index();
	alignment_index(size_t num_references);
	// parallelized constructor
	alignment_index(size_t num_references, size_t num_threads, std::vector<const pw_alignment *> data);
	~alignment_index();

	void insert(const pw_alignment * al);
	size_t erase(const pw_alignment * al);
	/*
		alignment overlap search function
		* al must point to valid alignment (does not need to be present in index)
		* reference must be 0 or 1
		* we return pointers to all alignments in index that overlap with al on its reference reference
	*/
	void search_overlap(const size_t & reference, const size_t & value, std::vector<const pw_alignment*> & result)const;
	void search_overlap(const pw_alignment & al, const size_t & reference, std::vector<const pw_alignment *> & result) const;
	void search_overlap(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result) const;
	void search_overlap_fraction(const size_t & reference, const size_t & left, const size_t & right, const double & fraction, std::vector<const pw_alignment *> & result) const;

	/*
		same as above plus:
		* all alignments that we just found, are deleted from this reference
		* we return all intervals on all references that were touched by the set of deleted alignments
			touched_intervals structure: reference -> (left, right)
	*/
	void super_search_overlap_and_remove(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result, std::multimap<size_t, std::pair<size_t, size_t> > & touched_intervals);

// this function can maybe be used for finding connected components in parallel
	void super_search_overlap_no_remove(const size_t & reference, const size_t & left, const size_t & right, std::vector<const pw_alignment *> & result, std::multimap<size_t, std::pair<size_t, size_t> > & touched_intervals) const;


	void debug_print() const;	



	private:
	typedef avl_interval_tree<const pw_alignment*> tree_type;
	std::vector<tree_type> trees;

	size_t num_threads;


	void half_insert(const pw_alignment * al, size_t reference);
};

#endif

