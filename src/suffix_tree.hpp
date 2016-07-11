#ifndef SUFFIX_TREE_HPP
#define SUFFIX_TREE_HPP
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <cmath>
#include <cfloat>

#include "model.hpp"

class suffixTree{
	public:
		suffixTree(size_t & num_seq, finding_centers & cent, std::vector<std::vector<std::vector<size_t> > > &);
		~suffixTree();
		void add_leaves(size_t &, std::vector<size_t>&);//parent,current
		void add_internal_node(size_t &, std::vector<size_t> & );
		void update_node(size_t & , std::vector<size_t> & );
		void split(size_t&,std::vector<size_t> & , std::vector<size_t> &,std::vector<size_t> &, bool &);//split the parent edge(first argument) after the common part(third argument)
		void make_suffix(std::vector<size_t> & ); //Builds suffixes of a sequence for one stirng at the time
		void find_first_edges_after_root(size_t & , std::vector<size_t> &);
		void find_next_parent(size_t & , std::vector<size_t> & , std::vector<size_t> & , size_t & , std::vector<size_t> & );
		void find_children(std::vector<size_t> &,std::vector<std::vector<size_t> >&);
		void find_common_part(std::vector<size_t> & parent , std::vector<size_t> & current, std::vector<size_t> & common_part);
		void make_tree_for_a_seq(size_t & seq);//Builds a suffix tree for one sequence at the time
		void make_tree(std::vector<size_t> & center_with_highest_gain, size_t & highest_index);
		void update_centers(size_t &, std::vector<size_t> & center_with_highest_gain, size_t & highest_index);
		void count_branches();
		void traverse_the_tree(size_t & , std::vector<size_t> & );
		const std::map<size_t, std::vector<size_t> > get_edges()const;
		const std::map<std::vector<size_t>, size_t> get_branches()const;
		const std::vector<std::vector<size_t> > get_current_centers(size_t &)const;
		void print_tree();


	private:
	//	std::multimap<std::pair<std::vector<size_t>, size_t> , std::vector<size_t> > tree; //previous context, current node index, current context
	//	std::multimap<std::vector<size_t>, size_t> edges_relations; //vector<size_t> is parent, size_t is current node index 
		std::multimap<size_t, size_t> edges_relations; //parent, children
		std::map<size_t, std::vector<size_t> > edges; //size_t is current node index , vector<size_t> is current edge content
		std::vector<std::vector<std::vector<size_t> > >words; // successive centers of sequences
		std::vector<std::vector<size_t> >suffixes; // Has suffixes of a single string at the the time
		std::map<size_t, size_t > first_edges; // its first center, its index
		std::map<std::vector<size_t> , size_t> branches;
		size_t last_node_index;
		finding_centers & centers;
		vector<size_t> root;
		size_t num_seq;




};


#endif
