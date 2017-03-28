#ifndef SUFFIX_TREE_HPP
#define SUFFIX_TREE_HPP
#include <vector>
#include <utility>
#include <set>
#include <map>
#include <cmath>
#include <cfloat>
#include "pw_alignment.hpp"
#include"data.hpp"

//#define ALLOWED_GAP 500

class finding_centers{
	public:
	finding_centers(const all_data &);
	~finding_centers();
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > &);
	void findMemberOfClusters(std::map<std::string,std::vector<pw_alignment> > &);
	void center_frequency(std::map<std::string,std::vector<pw_alignment> > &, std::vector<std::map<size_t, std::string> > & );
	void create_long_centers_candidates(std::vector<std::map<size_t, std::string> > & centerOnseq, size_t & gap_in_long_centers);
	void find_seq_centers( size_t & , std::vector<std::map<size_t, std::string> > & centerOnseq);
	void find_long_centers(size_t &, std::map<size_t, std::string> &, size_t & );
	std::map<size_t, std::vector<int> >get_center(size_t &)const;//Returns long centers and their locations
	const std::vector<std::vector<int> > & get_centers(size_t &)const;
	const size_t get_number_of_centers();//Returns the total number of centers
//	std::multimap< size_t, std::string> get_sequence_centers(size_t & )const;//It returns all the centers of a sequence, a center may happen more than once, thus it might be repeated in the vector. size_t shows the position of each center
	std::string find_center_name(int &)const;
	std::vector<size_t> get_long_center_position(size_t & seq_id , std::vector<std::string> & long_center);
	const std::map<std::string, int> get_index()const;
	void add_nonaligned_regions(std::vector<std::map<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::string> > & all_pieces, std::map<std::string, size_t> & non_aligned_right);
	void add_to_alignments(std::vector<std::map<size_t, std::string > > & centerOnSequence, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & ref, size_t & position, std::string & sequence, std::string & center_name , bool direction);
	const size_t find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref)const;
	private://all these containers are local and their containts are replacing after each call of clustering:
	const all_data & data;
	size_t INDEX;
	std::vector< std::map<size_t , pw_alignment> >AlignmentsFromClustering;//It is changing in each call of clustering. contains als on each seq at position size_t
	std::map<std::string,std::string> memberOfCluster; //first string is assocciated member and second one is its center
//	std::vector< std::multimap< size_t, std::string> >centersOnSequence;//all the centers that happen on each sequence and their positions 
	std::vector< std::map< size_t , std::vector<int> > > centersOfASequence;//centers that happen on each sequence and have less than ALLOWED_GAP bases difference on their positions. Fully reveresed are removed. (position, long_cent)
	std::map<size_t ,std::vector<std::vector<int> > > long_centers_of_a_sequence; //Potential long centers. size_t -->seq_id 
	std::vector<std::multimap<std::vector<int>,size_t > > initial_suffixes;
	std::map<std::string,int> center_id;
	std::map<int,std::string> center_index;
	
//	std::map<std::string, int>oriented_index;// centers, their indecies which can be negative if they are used in their reverse direction.
};

class suffixTree{
	public:
		suffixTree(size_t & num_seq, finding_centers & cent, std::vector<std::vector<std::vector<int> > > &);
		~suffixTree();
		void add_leaves(size_t &, std::vector<int>&);//parent,current
		void add_internal_node(size_t &, std::vector<int> & );
		void update_node(size_t & , std::vector<int> & );
		void split(size_t&,std::vector<int> & , std::vector<int> &,std::vector<int> &, bool &);//split the parent edge(first argument) after the common part(third argument)
		void make_suffix(std::vector<int> & ); //Builds suffixes of a sequence for one stirng at the time
		void find_first_edges_after_root(size_t & , std::vector<int> &);
		void find_next_parent(size_t & , std::vector<int> & , std::vector<int> & , size_t & , std::vector<int> & );
		void find_children(std::vector<int> &,std::vector<std::vector<int> >&);
		void find_common_part(std::vector<int> & parent , std::vector<int> & current, std::vector<int> & common_part);
		void make_tree_for_a_seq(size_t & seq);//Builds a suffix tree for one sequence at the time
		void make_tree(std::vector<int> & center_with_highest_gain, int & highest_index);
		void update_centers(size_t &, std::vector<int> & center_with_highest_gain, int & highest_index);
		void count_branches();
		void traverse_the_tree(size_t & , std::vector<int> & );
		const std::map<size_t, std::vector<int> > get_edges()const;
		const std::map<std::vector<size_t>, size_t> get_branches()const;
		const std::vector<std::vector<int> > get_current_centers(size_t &)const;
		void print_tree();


	private:
	//	std::multimap<std::pair<std::vector<size_t>, size_t> , std::vector<size_t> > tree; //previous context, current node index, current context
	//	std::multimap<std::vector<size_t>, size_t> edges_relations; //vector<size_t> is parent, size_t is current node index 
		std::multimap<size_t, size_t> edges_relations; //parent, children
		std::map<size_t, std::vector<int> > edges; //size_t is current node index , vector<int> is current edge content
		std::vector<std::vector<std::vector<int> > >words; // successive centers of sequences
		std::vector<std::vector<int> >suffixes; // Has suffixes of a single string at the the time
		std::map<size_t, size_t > first_edges; // its first center, its index
		std::map<std::vector<size_t> , size_t> branches;
		size_t last_node_index;
		finding_centers & centers;
		vector<int> root;
		size_t num_seq;




};


#endif
