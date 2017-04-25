#ifndef CONNECTED_COMPONENT_HPP
#define CONNECTED_COMPONENT_HPP

#include "pw_alignment.hpp"

#include "alignment_index.hpp"
// #include "overlap.hpp"
#include "dynamic_mc.hpp"
#include "model.hpp"

//#include "IntervalTree.hpp"

//#include <boost/graph/biconnected_components.hpp>
#include <map>
#include <vector>
#include <cassert>
#include <omp.h>
#include <set>

//#define FRACTION 0.5


class compute_cc_avl_fraction {
	public:
	compute_cc_avl_fraction(const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, size_t num_sequences, size_t num_threads): alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
	//	size_t count = 0;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = *it;
		//	std::cout << "al "<< count << " is :"<<std::endl;
		//	al->print();

			als.push_back(al);
			remained_als.insert(al);
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		//	count ++;
		//	if(count == 150) break;
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::cout << "remainder size" << remained_als.size() << std::endl;
		std::vector<const pw_alignment *> sorted_als(sorter.size());
		size_t num = 0;
		for(std::multimap<size_t, const pw_alignment *>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
		//	remained_als.insert(al);
		//	als.push_back(al);
		//	als_index.insert(std::make_pair(al,num));
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	//	alind->debug_print();
	}
	compute_cc_avl_fraction(std::set< pw_alignment, compare_pw_alignment> & als_in, size_t num_sequences, size_t num_threads): alind(NULL) {//Is used to create partially overlapped connected components
		this->num_threads = num_threads;
		std::cout << " CC AVL build index on " << als_in.size() << std::endl << std::flush;
		std::multimap<size_t, const pw_alignment *> sorter;
		for(std::set<pw_alignment, compare_pw_alignment>::iterator it = als_in.begin(); it!=als_in.end(); it++){
			const pw_alignment * al = &*it;
			als.push_back(al);
			remained_als.insert(al);
		//	remainders_als.insert(*it);
			size_t len = al->alignment_length();
			sorter.insert(std::make_pair(len, al));
		}
		std::cout << " Sorting done " << sorter.size()<< std::endl;
		std::cout << "remainder size" << remained_als.size() << std::endl;
		std::vector<const pw_alignment *> sorted_als(sorter.size());
		size_t num = 0;
		for(std::multimap<size_t, const pw_alignment *>::iterator it = sorter.begin(); it!=sorter.end(); ++it) {
			const pw_alignment * al = it->second;
			sorted_als.at(num) = al;
			alignments.insert(al);
		//	remained_als.insert(al);
		//	als.push_back(al);
		//	als_index.insert(std::make_pair(al,num));
			num++;
		}
		alind = new alignment_index(num_sequences, num_threads, sorted_als);
	//	alind->debug_print();
	}

	~compute_cc_avl_fraction(){
		delete alind;
	}
	void non_recursive_compute(std::set< const pw_alignment* , compare_pointer_pw_alignment> &, std::set<const pw_alignment* , compare_pointer_pw_alignment> & , double & );
	void cc_step_non_recursive(const pw_alignment * ,  size_t & , size_t & , size_t &, std::set< const pw_alignment* , compare_pointer_pw_alignment>& , double &  , std::set<const pw_alignment* , compare_pointer_pw_alignment>&);
	void compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > &);
	void get_cc(const pw_alignment & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment*, compare_pointer_pw_alignment> & );
	void cc_step(const pw_alignment &, size_t , size_t , size_t , std::set <const pw_alignment*, compare_pointer_pw_alignment> & , std::set <const pw_alignment* , compare_pointer_pw_alignment>  & );
	void compute_fraction(std::vector<const pw_alignment*> &, const double  ,  const size_t & , const size_t & , const size_t &,std::vector<const pw_alignment*> & );
	void node_adjacents(const pw_alignment*,std::set<const pw_alignment*>&)const;
	void node_adjacents(const pw_alignment*,std::vector<const pw_alignment*>&)const;
	const std::map<const pw_alignment* , std::set<const pw_alignment*> > get_edges()const{
		return edges;
	}
	size_t get_id (const pw_alignment*)const;
	void test_redundancy(const pw_alignment* , std::set< const pw_alignment* , compare_pointer_pw_alignment> & );
	std::set<const pw_alignment*, compare_pointer_pw_alignment> get_remaining_als()const;
//	const std::set<pw_alignment, compare_pw_alignment> get_remainders()const;
	private:
	
	std::set<const pw_alignment*, compare_pointer_pw_alignment> alignments; 
	std::vector<const pw_alignment*> als;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> remained_als; 
//	std::set<pw_alignment, compare_pw_alignment> remainders_als; 

	alignment_index  * alind;
	size_t num_threads;
//	std::map<const pw_alignment* , size_t> als_index;
	std::map<const pw_alignment* , std::set<const pw_alignment*> >edges; 
//	std::set<std::pair<const pw_alignment*, const pw_alignment*> > edges;

	size_t edge_size;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> seen1; 



};
//This class was initially written to find articulation points and make biconnected graph components. IT is able to find the articulations point correctly but fills in the 'stack' with wrong nodes. In case of being used should be fixed.
class biconnected_component{
	public:
	biconnected_component(compute_cc_avl_fraction & cc, std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> >& als):cc_fraction(cc), visited(als.size()), low(als.size()), Depth(als.size()), parent(als.size()),alignments(als.size()),als_index(als.size()){
		depth = 0;
		std::cout << "als size "<< als.size()<<std::endl;
		for(size_t i =0; i < als.size();i++){
			size_t index=0;
			for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = als.at(i).begin(); it !=als.at(i).end();it++){
				alignments.at(i).push_back(*it);
				als_index.at(i).insert(std::make_pair(*it,index));
				index++;
			}
			visited.at(i).resize(index);
			Depth.at(i).resize(index);
			low.at(i).resize(index);
			parent.at(i).resize(index);
		}		
	}
	~biconnected_component(){}
	void creat_biconnected_component(size_t & );
	void get_articulation_points(size_t &, size_t node);


	private:
	compute_cc_avl_fraction & cc_fraction;
	size_t depth;
	std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> > stack;
	std::vector<std::vector<const pw_alignment*> >alignments;//The node indices are their row in the alignments
	std::vector<std::map<const pw_alignment* , size_t> >als_index;
	std::vector<std::vector<bool> >visited;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> temp;
	std::vector<std::vector<size_t> >low;
	std::vector<std::vector<size_t> >Depth;
	std::vector<std::vector<size_t> >parent;
};
//2-edge-connected-component:
//Note that each set of als is a single component graph
class two_edge_cc{
	public:
	two_edge_cc(compute_cc_avl_fraction & cc, std::set< const pw_alignment* , compare_pointer_pw_alignment> & als):cc_fraction(cc){
		count = 0;
		bridge = false;
		std::cout << "als are"<<std::endl;
		size_t count = 0;
		for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = als.begin(); it !=als.end();it++){
			alignments.push_back(*it);
		//	std::cout << "al "<< count << " is :"<<std::endl;
			const pw_alignment* al = *it;
		//	al->print();
			visited_als.insert(std::make_pair(*it,-1));
			als_low.insert(std::make_pair(*it,-1));
			count ++;
	//		if(count == 30) break;
		}
		std::cout<< "visit size "<< visited_als.size() << "low size "<< als_low.size() << " al size "<< alignments.size() <<std::endl;
	}
	two_edge_cc(compute_cc_avl_fraction & cc):cc_fraction(cc){


	}
	void find_bridges(std::vector<std::set<const pw_alignment*, compare_pointer_pw_alignment> >& , std::set<pw_alignment, compare_pw_alignment> &);
	void deep_first_search(const pw_alignment*, const pw_alignment*,std::set<pw_alignment, compare_pw_alignment> &);

//non recurssive:
	void make_two_edge_graph(std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als);
	void find_bridges(const pw_alignment*);
	void remove_bridges(std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als);
	void bfs(const pw_alignment* s, std::map<const pw_alignment* , std::set<const pw_alignment*> > & edges, std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als,std::map<const pw_alignment*, bool> & visited);
	~two_edge_cc(){}
	private:
	compute_cc_avl_fraction & cc_fraction;	
	std::vector<const pw_alignment*> alignments;
	std::map<const pw_alignment*, int> visited_als;
	std::map<const pw_alignment*, int> als_low;
	size_t count;
	bool bridge;

	std::multimap<const pw_alignment*,const pw_alignment*> bridges;
	std::set<const pw_alignment*, compare_pointer_pw_alignment> temp;

	std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> > stack;
	std::set<const pw_alignment*,compare_pointer_pw_alignment> remainder;


	
};

/*
	divide and conquer strategy:

	* get a set of alignments
	* separate into 2-edge connected components
	* components with over 100 alignments become a divide_and_conquer_al_problem (with increased overlap_fraction to further separate them), 
		components with between 3 and 99 alignments get processed into an overlap class
	* aggregate child results with remaining alignments into an overlap object


	STATES for parallelzation
	WAIT_FOR_CC
	* start process alignments into connected components
	WAIT_FOR_CHILDREN
	* connected components are done, waiting till all child processes are done
	COMBINE_RESULTS
	* combining results into single overlap class
	DONE
	

*/
class divide_and_conquer_al_problem;

class divide_and_conquer_alignments {
	public:

	divide_and_conquer_alignments(const all_data & data, const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, const dynamic_mc_model & model, size_t num_sequences, size_t num_threads);

	void run(std::vector<std::vector<pw_alignment> > & overlaps);

	const all_data & data;

	friend class divide_and_conquer_al_problem;
	
	private:
	const dynamic_mc_model & model;
	size_t num_threads;
	size_t num_sequences;
	std::vector<divide_and_conquer_al_problem*> all_problems;
	std::vector<divide_and_conquer_al_problem*> wait_for_cc_problems;
	std::vector<divide_and_conquer_al_problem*> wait_for_overlap_problems;
	std::vector<divide_and_conquer_al_problem*> wait_for_combine_problems;
	size_t waiting_threads;
	size_t running_threads;

	std::vector<const pw_alignment *> all_alignments;
	std::vector<std::string> thread_info;
	
	bool select_next_task(size_t thread);
	bool are_we_done() const;
	void add_cc_problem(divide_and_conquer_al_problem* ap);



};



class divide_and_conquer_al_problem {
	public:
	divide_and_conquer_al_problem(divide_and_conquer_alignments & parent, const std::vector<size_t> & alvec, size_t separation_level, const size_t & num_sequences, const size_t & num_threads): parent(parent), parent_problem(NULL), alvec(alvec), separation_level(separation_level), num_sequences(num_sequences), num_threads(num_threads), results(parent.data) {
		std::cout << " made NEW al prob " << separation_level << " " << num_threads<< std::endl;
		state = WAIT_FOR_CC;
	}
	~divide_and_conquer_al_problem() {
	}
	typedef overlap_interval_tree<pw_alignment> overlap_type;
	typedef dynamic_mc_model use_model;
	typedef initial_alignments_from_groups<dynamic_mc_model,overlap_interval_tree<pw_alignment> > use_ias;
	
	enum state_type {WAIT_FOR_CC, WAIT_FOR_OVERLAP, WAIT_FOR_COMBINE, DONE};

	void run_cc(size_t thread);
	void run_overlap(size_t cl_al_idx, size_t thread);
	void child_done(divide_and_conquer_al_problem * child, const std::vector<std::vector<pw_alignment> > & child_res); 
	void get_results(std::vector<std::vector<pw_alignment> > & child_res);
	void run_combine(size_t thread);
	
	void select_next_alignment_cluster(size_t & idx, bool & found_next, bool & all_done);
	bool are_all_alignment_clusters_done();

	bool is_ready_for_combine();

	void set_parent_problem(divide_and_conquer_al_problem * p);
	bool has_parent() const;
	

	private:
	divide_and_conquer_alignments & parent;
	divide_and_conquer_al_problem * parent_problem;

// input for run_cc
	const std::vector<size_t> alvec;
	size_t separation_level;

// input for run_overlap
	std::vector<std::vector<size_t> > clustered_alignments;
	std::vector<size_t> clustered_al_status; // 0=not yet running, 1=running, 2=completed



// input for run_combine
	std::vector<size_t> leftover_alignments;
	std::vector<divide_and_conquer_al_problem *> child_problems; // if child problems are done, they move their content to clustered_al_results
	std::vector<std::vector<pw_alignment> > clustered_al_results;

// final result
	overlap_type results;



// general
	size_t num_sequences;
	const size_t num_threads;
	state_type state;

	typedef avl_interval_tree<size_t> tree_type;
	// min/max sizes of connected components in order to resolve partial overlap separately
	const size_t MIN_CC_SIZE = 3; 
	const size_t MAX_CC_SIZE = 300; // TODO 100 


	void two_edge_cc(const std::vector<std::pair<size_t, size_t> > & edges, size_t cc_edges, std::vector<std::set<size_t> > & cc_result);
	void bridge_step(size_t & time, const std::map<size_t, std::set<size_t> > & edges, const size_t & node, std::vector<size_t> & parent, std::vector<bool> & visited, std::vector<size_t> & dfsnum, std::vector<size_t> & low, std::vector<std::pair<size_t, size_t> > & bridges);
	void cc_step(const size_t & node, const std::map<size_t, std::set<size_t> > & edges, std::vector<bool> & visited, std::set<size_t> & cur_cc, std::vector<std::pair<size_t, size_t> > & tree_edges);
	void edges_remove(const std::vector<std::pair<size_t, size_t> > & redges, std::map<size_t, std::set<size_t> > & medges);
	void path_remove(std::map<size_t, std::set<size_t> > & edgemap, std::vector<std::set<size_t> > & cc_result);
	void edge_insert(std::map<size_t, std::set<size_t> > & edges, size_t i, size_t j);
	void edge_remove(std::map<size_t, std::set<size_t> > & edges, size_t i, size_t j);
	size_t count_edges(const std::map<size_t, std::set<size_t> > & edges);
	void separation_conditions(double & overlap_fraction, size_t & cc_edges);
	void make_graph(double overlap_fraction, std::vector<std::pair<size_t, size_t> > & edges);
	void make_graph_singlethreaded(double overlap_fraction, std::vector<std::pair<size_t, size_t> > & edges);
	void make_graph(double overlap_fraction, const std::vector<pw_alignment> & als, std::vector<std::pair<size_t, size_t> > & edges);
	void als_to_ccs(std::vector<pw_alignment> & als, std::vector<std::vector<pw_alignment> > & ccs);
};







#endif
