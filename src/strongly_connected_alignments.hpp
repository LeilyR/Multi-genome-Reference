#ifndef STRONGLY_CONNECTED_HPP
#define STRONGLY_CONNECTED_HPP

#include "pw_alignment.hpp"

#include "alignment_index.hpp"

//#include "IntervalTree.hpp"

//#include <boost/graph/biconnected_components.hpp>
#include <map>
#include <vector>
#include <cassert>
#include <omp.h>
#include <set>



class detect_overlap_rate{

	public:
	detect_overlap_rate(size_t number , std::set< const pw_alignment*, compare_pointer_pw_alignment> & cc):alind(number){
		this->number = number;
		edge_size =0;
		size_t index=0;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = cc.begin(); it!= cc.end(); it++){
			const pw_alignment * al = *it;
			alind.insert(al);
			alignments.insert(al);
			als_index.insert(std::make_pair(al,index));
			index++;
		}
	}
	~detect_overlap_rate(){ }
	void compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & );
	void get_sc(const pw_alignment & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & );
	void sc_step(size_t , size_t, size_t, size_t, std::set <const pw_alignment* , compare_pointer_pw_alignment>&, size_t& );
	void check_the_strength(std::vector<const pw_alignment*> & , size_t &, size_t &, size_t &, size_t &, std::vector<const pw_alignment*> &, size_t& );
	static void get_adjucent(size_t , std::vector<size_t>);
	size_t get_id (const pw_alignment*)const;
	

	private:
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments;
	std::map<const pw_alignment* , size_t> als_index;
	alignment_index alind;
	size_t number;
	std::multimap<size_t , size_t>edges; //size_t --> node , size_t ----> its adjucent nodes
	size_t edge_size;
	std::set<const pw_alignment* , compare_pointer_pw_alignment> seen;



};

class biconnected_component{
	public:
	biconnected_component(detect_overlap_rate &ov, std::set< const pw_alignment* , compare_pointer_pw_alignment>& als):ov_rate(ov), visited(als.size()), low(als.size()), Depth(als.size()){
		depth = 0;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = als.begin(); it !=als.end();it++){
			alignments.push_back(*it);
		}
	}
	~biconnected_component(){}
	void creat_component();
	void get_articulation_points(size_t node);


	private:
	detect_overlap_rate ov_rate;
	size_t depth;
	std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> > stack;
	std::vector<const pw_alignment*> alignments;//The node indices are their row in the alignments
	std::vector<bool> visited;
	std::vector<size_t> low;
	std::vector<size_t> Depth;
	std::vector<size_t> parent;






};








#endif
