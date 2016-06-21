#ifndef INTERVALS_HPP
#define INTERVALS_HPP

#include <stdio.h>
#include <stdlib.h> 
#include <vector>
#include <map>
#include <limits>
#include <iostream>
#include <sstream>
#include <set>

#define DO_REBALANCING 1
#define DO_SLOW_CHECKS 0
#define DEBUG_PRINTS 0

template<typename T>
class avl_interval_tree;

template<typename T>
class it_node {
	typedef T value_type;
	friend class avl_interval_tree<T>;

	public:
	it_node();
	it_node(const size_t & left, const size_t & right, const value_type & v);
	~it_node();

	

	private:
	size_t height; // height of subtree, leaves have height 0. To compute height we only consider paths to the left or right. Center height does not count

	const size_t left; // left, right is both inclusive like we did in pw_alignment
	const size_t right;
	const value_type value;
	

	size_t min_left; // min left position of the subtree below this node
	size_t max_right; // max right of subtree below 


	it_node * str_left;
	it_node * str_right;
	it_node * str_center;
	it_node * parent;



};


template<typename T>
class node_pointer_comp {
	typedef T value_type;
	public:
	bool operator()(const value_type * a, const value_type * b) const {
		return a < b;
	}
	

};


template<typename T>
class avl_interval_tree {
	typedef T value_type;

	public:
	avl_interval_tree();
	~avl_interval_tree();
	
	void insert(const size_t & left, const size_t & right, const value_type  & v);

	void bulk_insert(std::vector< std::pair < std::pair< size_t, size_t > , value_type &>  > & data );
	// erases all intervals with that left and right, returns vector of erased intervals
	void erase_all(const size_t & left, const size_t & right, std::vector<value_type> & erased);
	// erases all intervals with that left and right and value, return number of erased intervals
	size_t erase(const size_t & left, const size_t & right, const value_type & value);
	void overlap_search(const size_t & value, std::vector<value_type> & results)const;
	size_t erase_vector(const std::vector< std::pair< std::pair<size_t, size_t> , value_type > > & todel);
	void overlap_search(const size_t & left, const size_t & right, std::vector<value_type> & results) const;
	void overlap_search_erase(const size_t & left, const size_t & right, std::vector<value_type> & results);
	void join_intervals(std::vector<std::pair<size_t, size_t> > & results) const;
	void overlap_fraction_search(const size_t & left, const size_t & right, const double & fraction, std::vector<value_type> & results) const;

	size_t size() const;
	size_t height() const;
	double list_ratio() const;
	size_t get_num_reba() const;
	size_t get_num_moves() const;
	
	void debug_print() const;


	private:
	typedef it_node<T> node_type;
	typedef typename std::map<node_type *, node_type *, node_pointer_comp<node_type> > node_map_type;
	typedef typename std::set<node_type *, node_pointer_comp<node_type> > node_set_type;
	node_type * root;
	size_t num_nodes;

	size_t num_reba;
	size_t num_moves;
	
	node_set_type unbalanced_nodes; 
#if DO_SLOW_CHECKS
	std::vector< std::pair<std::pair<size_t, size_t>, value_type> > slow_index;
#endif

	// returns a new insertion point. This is either atn or a new parent of atn that was risen by insert switch
	node_type * insert_at_node(node_type * atn, node_type * newn);
	void forest_insert(node_map_type & roots_to_parents, node_set_type & valid_parents); 
	// remove atn from the tree, then delete it, also delete removed node from affected nodes
	void erase_this_node(node_type * atn); 
	void redo_height_up(node_type * n);
	void forest_add(node_type * atn, node_type * insert_node, node_map_type & roots_to_parents, node_set_type & valid_parents); 
	void forest_delete_root(node_type * n, node_set_type & valid_parents); 
	void delete_nodes(const node_set_type & to_Delete); 
	void erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & erased);
	size_t single_erase_at_node(node_type * atn, const size_t & left, const size_t & right, const value_type & value);
	void overlap_search_for_a_value(node_type * atn, const size_t & value, std::vector<value_type> & results)const;

	void overlap_search_erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results, node_set_type & to_Delete );
	void join_intervals_at_node(node_type * atn, std::map<size_t, size_t> & joined_intervals);
	void subtree_max_interval(node_type * atn, size_t & left, size_t & right);
	void overlap_search_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results) const;
	void find_at_node(node_type * atn, const size_t & left, const size_t & right, const value_type & value, node_set_type & result) const;
	void join_intervals_at_node(node_type * atn, std::map<size_t, size_t> & joined_intervals) const;
	void overlap_fraction_search_at_node(node_type * atn, const size_t & left, const size_t & right, const double & fraction, std::vector<value_type> & results) const;
	size_t count_center_at_node(node_type * atn) const;
//	void subtree_max_interval(node_type * atn, size_t & left, size_t & right) const;
	void debug_print_at_node(node_type * atn, const size_t & level, size_t & id, const size_t & parent_id) const;

	// returns new parent node of inserted subtree (because of center insert parent switches)
	node_type * move_subtree(node_type * from, node_type * to);
	void delete_subtree(node_type * n);

	void prune_subtree(node_type * n);
	void left_attach(node_type * parent, node_type * child);
	void right_attach(node_type * parent, node_type * child);
	void center_attach(node_type * parent, node_type * child);
	void delete_node(node_type * n);
	bool redo_height_and_minmax_at_node(node_type * n);
	
	int get_balance(node_type * n) const;
	void balance_up(node_type * n);
	void set_rebalance();
	node_type * balance_node(node_type * n);
	void rotate_lr_to_ll(node_type * n);
	void rotate_rl_to_rr(node_type * n);
	node_type * rotate_left(node_type * n);
	node_type * rotate_right(node_type *n);

	void check_tree() const;
	size_t check_subtree(node_type * n) const;
	void check_tree_at_node(node_type * n,  const size_t & level, size_t & id, const size_t & parent_id, size_t & subtree_size, size_t min_left, size_t max_right, size_t min_in, size_t max_in, bool use_slow_index) const;
	void check_upwards(node_type * n, size_t & minl, size_t & maxr, size_t & intv_l, size_t & intv_r) const;
	bool balance_check_at_node(node_type * n) const;

};

#endif
