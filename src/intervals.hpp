#ifndef INTERVALS_HPP
#define INTERVALS_HPP
#include <stdio.h>
#include <stdlib.h> 
#include <vector>
#include <map>
#include <limits>
#include <sstream>
#include <set>

#define DO_REBALANCING 1
#define DO_SLOW_CHECKS 1

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

	size_t left; // left, right is both inclusive like we did in pw_alignment
	size_t right;
	value_type value;
	

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
	// erases all intervals with that left and right, returns vector of erased intervals
	void erase_all(const size_t & left, const size_t & right, std::vector<value_type> & erased);
	// erases all intervals with that left and right and value, return number of erased intervals
	size_t erase(const size_t & left, const size_t & right, const value_type & value);
	void overlap_search(const size_t & left, const size_t & right, std::vector<value_type> & results);
	void overlap_search_erase(const size_t & left, const size_t & right, std::vector<value_type> & results);
	void join_intervals(std::vector<std::pair<size_t, size_t> > & results);
	
	void debug_print() const;
	

	private:
	typedef it_node<T> node_type;
	typedef typename std::map<node_type *, node_type *, node_pointer_comp<node_type> > node_map_type;
	typedef typename std::set<node_type *, node_pointer_comp<node_type> > node_set_type;
	node_type * root;
	size_t num_nodes;
#if DO_SLOW_CHECKS
	std::vector< std::pair<std::pair<size_t, size_t>, value_type> > slow_index;
#endif

	void insert_at_node(node_type * atn, node_type * newn);
	void forest_insert(node_map_type & roots_to_parents, node_set_type & valid_parents);
	// remove atn from the tree, then delete it
	void erase_this_node(node_type * atn); 
	void redo_height_up(node_type * n);
	void forest_add(node_type * atn, node_type * insert_node, node_map_type & roots_to_parents, node_set_type & valid_parents);
	void forest_delete_root(node_type * n, node_set_type & valid_parents);
	void delete_nodes(const node_set_type & to_Delete);
	void erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & erased);
	size_t single_erase_at_node(node_type * atn, const size_t & left, const size_t & right, const value_type & value);
	void overlap_search_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results);
	void overlap_search_erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results, node_set_type & to_Delete );
	void join_intervals_at_node(node_type * atn, std::map<size_t, size_t> & joined_intervals);
	void subtree_max_interval(node_type * atn, size_t & left, size_t & right);
	void debug_print_at_node(node_type * atn, const size_t & level, size_t & id, const size_t & parent_id) const;


	void move_subtree(node_type * from, node_type * to);
	void delete_subtree(node_type * n);

	void prune_subtree(node_type * n);
	void left_attach(node_type * parent, node_type * child);
	void right_attach(node_type * parent, node_type * child);
	void center_attach(node_type * parent, node_type * child);
	void delete_node(node_type * n);
	void redo_height_and_minmax_at_node(node_type * n);
	void check_tree() const;
	void check_tree_at_node(node_type * n,  const size_t & level, size_t & id, const size_t & parent_id, size_t & subtree_size, size_t min_left, size_t max_right, size_t min_in, size_t max_in) const;


};

 
#endif
