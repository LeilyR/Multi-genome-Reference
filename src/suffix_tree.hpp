#ifndef SUFFIX_TREE_HPP
#define SUFFIX_TREE_HPP
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <cfloat>

#include "model.hpp"

class suffixTree{
	public:
		suffixTree();
		~suffixTree();
		void add(std::vector<size_t> &, std::vector<size_t>&);//parent,current
		void split(std::vector<size_t> & , std::vector<size_t> &,std::vector<size_t> &);//split the parent edge(first argument) after the common part(third argument)
		void make_suffix(std::vector<size_t> & ); //Builds suffixes of a sequence for one stirng at the time
		void find_first_edges_after_root(std::vector<size_t> & , std::vector<size_t> &);
		void find_children(std::vector<size_t> &,std::vector<size_t>&);
		void make_tree(size_t &, bool highest_gain);//Builds a suffix tree for one sequence at the time, size_t ----> id of that sequence

	private:

		std::multimap< std::vector<size_t> , size_t > tree;
		std::multimap<std::pair<std::vector<size_t>, size_t > , std::vector<size_t>& > edges;
		std::vector<std::vector<size_t> > words; // successive centers of a sequence
		std::vector<std::vector<size_t> >suffixes; // Has suffixes of a single string at the the time
		std::map<size_t, std::vector<size_t> > first_edges;
		size_t last_node_index;
		finding_centers & centers;




};


#endif
