#include "data.hpp"
#include "pw_alignment.hpp"
#include "graph.hpp"
#include "needleman_wunsch.hpp"
#include "dijkstra.hpp"

#define MAXGAP 400
//This class is only used to read the dot file includes the reference graph
class ref_graph{
	public:
	ref_graph( const all_data & d):data(d){
		longname2seqidx = data.getLongname2seqidx();
		std::map<std::string, size_t>::iterator longname = longname2seqidx.begin();
		std::cout << "the begin is "<< longname->first<<std::endl;
	}
	~ref_graph(){}
	void read_dot_file(std::string &, std::string &);
	void add_adjacencies(std::string & , std::string & , std::string & , std::string & );
	void add_adjacencies(std::string & , std::string & );
	void add_predecessor(int & this_node , int & pre_node);
	const std::map<int , std::set<int> > & get_adjacencies()const;
	const std::vector<std::vector<int> > get_predecessor(unsigned int &, bool , size_t &)const;
	const void look_for_predecessor(int & node , size_t & length , size_t & current_length, std::string & acc_name,std::set<int> & visited,std::vector<int> & this_pre_nodes, std::vector<std::vector<int> > & all_pre_nodes)const;
	const std::vector<std::vector<int> > get_successor(unsigned int & ref_id, bool dir , size_t & length)const;
	const void look_for_successor(int & node, size_t & length, size_t & current_length, std::string & acc_name, std::set<int> & visited, std::vector<int> & this_adjacencies, std::vector<std::vector<int> > & all_adjacencies)const;
	void read_gfa_for_adj(std::string &);
	void read_gfa_for_seq(std::string & gfafile , std::ofstream & fasta , std::ofstream & paths);
	const unsigned get_refid(size_t & , int& )const;

	size_t seq_length(std::string &, std::string &);
	std::string seqname(int & );
	void deep_first_search(int &, std::string & , size_t & );
	void look_for_neighbors(int & , std::map<int,bool> &  , std::string &, int & accu_length, std::vector<int> & apath, size_t & );
	void bfs(int & startnode, std::string & refacc, size_t & right_on_ref);
	const vector<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	std::set<int> get_adjacencies(int & node)const{
		std::map<int, std::set<int> >::const_iterator it = adjacencies.find(node);
		assert( it != adjacencies.end());
		return it->second;
	}
	std::set<int> get_predecessor(int& node)const {
		std::set<int> nodes;
		std::cout << "node: "<< node <<std::endl;
		std::map<int, std::set<int> >::const_iterator it = predecessors.find(node);
		if( it != predecessors.end()){
			nodes = it->second;
		}
		return nodes;
	}

	private:
	const all_data & data;
	std::map<std::string, size_t> longname2seqidx;

	std::map<int, std::set<int> > adjacencies;
	std::map<int, std::set<int> > predecessors;

	vector<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;
	std::multimap<int,int> parental_relation; //first one shows the child node ,second one shows the parent
};

/*class deep_first{//It runs deep first algorithm on reference graph to retrun parts of graph that the distance between the first and the last node is less than MAXGAP.
	public:
	deep_first(const all_data & d, std::map<int, std::set<int> > & adjacencies):data(d){
		this -> adjacencies = adjacencies;
	//	for(std::map<int, std::set<int> >::iterator it = adjacencies.begin(); it != adjacencies.end(); it++){
	//		std::cout << it->first << std::endl;
	//	}

	}
	~deep_first(){}
	size_t seq_length(std::string &, std::string &);
	std::string seqname(int & );
	void deep_first_search(int &, std::string & , size_t & );
	void look_for_neighbors(int & , std::map<int,bool> &  , std::string &, int & accu_length, std::vector<int> & apath, size_t & );
	void bfs(int & startnode, std::string & refacc, size_t & right_on_ref);
	const vector<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	private:
	const all_data & data;
	std::map<int, std::set<int> > adjacencies;
	vector<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;
	std::multimap<int,int> parental_relation; //first one shows the child node ,second one shows the parent

};*/

//This class is used to make several components on a read where each component contains aligned regions with less than MAX_GAP base distance
class als_components{

	public:
	als_components(const all_data & d , ref_graph & refg, const dynamic_mc_model & m ,const std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> >  & ccs):data(d), rgraph(refg),model(m){
		index = 2;
		alignments_left.resize(ccs.size());
		alignments_right.resize(ccs.size());

		for(size_t i = 0 ; i < ccs.size(); i++){
			std::cout << "cc at " << i << " : "<<std::endl;
			for(std::set<const pw_alignment* , compare_pointer_pw_alignment >::iterator it = ccs.at(i).begin() ; it != ccs.at(i).end() ; it++){
				const pw_alignment * al = *it;
				alignments_left.at(i).insert(*al);
				alignments_right.at(i).insert(*al);

				al->print();

			}
			std::multiset<pw_alignment , sort_pw_alignment_by_left >::iterator first = alignments_left.at(i).begin();
			pw_alignment p = *first;
			size_t l2,r2;
			p.get_lr2(l2,r2);
			size_t lower = l2;
			std::cout << "lower " << lower <<std::endl;//The least left
			std::multiset<pw_alignment , sort_pw_alignment_by_right>::iterator last = alignments_right.at(i).end(); 
			--last;
			p = *last;
			p.get_lr2(l2,r2);
			size_t upper = r2;
			std::cout << "upper "<< upper<<std::endl;
			boundries.insert(std::make_pair(lower,upper));
			
			
		}
	}
	~als_components(){}
	void merg_components(){//TODO change it with join interval function in avl tree to make it faster
		std::multiset<pw_alignment,sort_pw_alignment_by_left > templ;
		std::multiset<pw_alignment,sort_pw_alignment_by_right > tempr;
		size_t previous_right, current_left;
		size_t counter = 0;
		for(std::map<size_t, size_t>::iterator it = boundries.begin(); it != boundries.end();it++){
			std::cout << it->first <<  " "<< it->second<<std::endl;
			if(it == boundries.begin()){
				previous_right = it->second;
				templ = alignments_left.at(counter);
				tempr = alignments_right.at(counter);
				counter++;
			}else{
				current_left = it->first;
				if(current_left - previous_right < MAXGAP){
					templ.insert(alignments_left.at(counter).begin(), alignments_left.at(counter).end());
					tempr.insert(alignments_right.at(counter).begin(),alignments_right.at(counter).end());
					previous_right = it->second;
					counter++;
				}else{
					std::cout<<"members of a new component"<<std::endl;
					for(std::multiset<pw_alignment,sort_pw_alignment_by_left>::iterator it1 = templ.begin() ; it1 != templ.end();it1++){
						pw_alignment pw = *it1;
						pw.print();
					}
					alignments.push_back(templ);
					ordered_als_by_right.push_back(tempr);
					assert(tempr.size()==templ.size());
					templ = alignments_left.at(counter);
					tempr = alignments_right.at(counter);
					previous_right = it->second;					
					counter++;
				}
			}
		}
		if(templ.size() !=0){//Adding the last one
			assert(tempr.size() != 0);
			assert(tempr.size()==templ.size());
			ordered_als_by_right.push_back(tempr);
			alignments.push_back(templ);
		}
		assert(alignments.size()==ordered_als_by_right.size());
		std::cout << "als size " << alignments.size()<<std::endl;
		adjacencies.resize(alignments.size());
		weight_of_edges.resize(alignments.size());
		indices_nodes.resize(alignments.size());
		node_indices.resize(alignments.size());
	}
	/*Go over the alignments, for each al find all the possible paths on the graph that starts with that alignment and has <=MAXGAP length,
	then find the alignments that their left is closer than MAXGAP bases to the current one.Check them on the graph. If they are on those 
	paths that we already picked then we make new als using needleman wunsch algorithm.
	For those which are not on those paths we check if their left is smaller than beginning of the current alignment component- MAXGAP/2, 
	if so then we try to find the graph paths for it and go further otherwise discard it.*/
	void find_als_on_paths(std::ofstream & output, size_t & refacc, size_t & readacc);
	std::multiset<pw_alignment,sort_pw_alignment_by_left> find_successors(const pw_alignment & p , size_t &);
	void get_subgraph(const pw_alignment &, std::vector<std::vector<int> > &);
	void look_for_successors_on_paths(size_t & comp, const pw_alignment& , std::set<int>& , std::vector<std::vector<int> > & ,std::multiset<pw_alignment,sort_pw_alignment_by_left> &);
	void make_al_on_a_node(size_t & comp, const pw_alignment&, const pw_alignment& ,bool dir, unsigned int & ref1, size_t & left1, size_t & right1, unsigned int & ref2, size_t & left2, size_t & right2, size_t & refacc, size_t & readacc);
	void get_paths(int & node_name , int & current_node_name, std::vector<std::vector<int> > & paths);
	void append_nodes(std::vector<std::vector<int> > & paths, int & name, int & current_node_name, size_t & refacc , std::vector<size_t> & this_refs,std::vector<int> & this_path, std::string & str1);
	void make_alignments(size_t & comp, size_t & begin_on_read, size_t & end_on_read, const pw_alignment & first_al, const pw_alignment & last_al, std::string & read_out, std::string & ref_out, std::vector<size_t> &this_refs, std::vector<int> & this_path, unsigned int & ref2, size_t & first_length, size_t & last_length,size_t &refacc,size_t&readacc);
	void compute_samples(bool tillend, size_t node_length, size_t & current_pos, size_t & read_counter, size_t & ref_counter,std::string & read_in, std::string & read_out, std::string & ref_in, std::string & ref_out);
	void get_reverse_complement(std::string & sequence_in , std::string & sequence_out);
	void add_first_and_last_als(size_t & i, const pw_alignment & p, size_t & , size_t & );
	void looking_for_first_al(size_t & ,const pw_alignment & p,size_t &);
	void looking_for_last_al(size_t & ,const pw_alignment & p,size_t &);
	double get_cost(const pw_alignment & p, size_t & acc1 , size_t & acc2);
	void make_first_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & left_on_current_node,size_t & refacc,size_t & read_id,size_t & current_node, std::vector<pw_alignment> & first_als);
	void make_last_als(std::vector<int> & nodes , std::string & refout , std::string & readout , size_t & first_begin, size_t & last_end, std::string & seq_from_ref , size_t & right_on_current_node,size_t & refacc, size_t & read_id, size_t & current_node, std::vector<pw_alignment> & last_als);
	const std::vector<pw_alignment> find_the_best_als(std::vector<std::vector<pw_alignment> >&, size_t & refacc, size_t & readacc)const;
	void add_adjacencies(size_t & comp, const pw_alignment &, const pw_alignment &,size_t & read_acc, size_t & ref_acc);
	void add_edge_to_begin(size_t & comp, const pw_alignment &,size_t & read_acc, size_t & ref_acc);
	void add_edge_to_end(size_t & comp, const pw_alignment &);
	void add_expensive_edges(size_t & comp,size_t & refacc, size_t & readacc);
	void add_to_maf(const pw_alignment & al, std::ofstream & output, bool & firstal);
	const pw_alignment & get_alignment(size_t & , size_t &)const;
	private:
		const all_data & data;
	//	deep_first & dfirst;
		ref_graph & rgraph;
		const dynamic_mc_model & model;
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_left> >alignments_left;
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_right> >alignments_right;
		std::map<size_t,size_t> boundries;
		std::vector<std::multiset<pw_alignment,sort_pw_alignment_by_left> > alignments;//Includes components of alignments that have more than MAXGAP distance with each other.
		std::vector<std::multiset<pw_alignment , sort_pw_alignment_by_right> >ordered_als_by_right;
	//	std::vector<std::set<const pw_alignment*> >add_to_start;
	//	std::vector<std::set<const pw_alignment*> >add_to_end;
		std::vector<std::vector<size_t> > shortest_path;
		std::vector< std::multimap<size_t , size_t> >adjacencies; //from_al index , to_al index , vector goes over different components
		std::vector<std::multimap<const pw_alignment, size_t, compare_pw_alignment> >node_indices;//al and its index
		std::vector<std::map<std::pair<size_t,size_t>,double> >weight_of_edges; //Edges and their weight. In the case of exisitng adjacencies, weight is always the modification of the end node of an edge 
		std::vector<std::multimap<size_t, const pw_alignment> >indices_nodes;
		size_t index;
		size_t previous_right2;
		size_t previous_right1;
		size_t previous_left1;
		size_t previous_ref;
		bool forward;


};

class map_check{
	public:
	map_check(){}
	~map_check(){}
	void read_graph_maf(std::ifstream & graph_maf);
	void read_alignments(std::ifstream & alignments);//Read all the intial created alignments before mapping any read against a graph
	void check_output(std::ifstream & mapping_maf,const std::string & seqname);
	void check_an_alignment(unsigned int & ref1, const std::string & seqname, size_t & left2 , size_t & right2);
	void read_txt_file_long_center(std::ifstream & , std::string & );
	void read_txt_file(std::ifstream & , std::string &);
	std::multimap<std::string, std::pair<unsigned int, unsigned int> > get_nodes()const{
		return nodes;
	}
	void print_nodes(){
		for(std::multimap<std::string , std::pair<unsigned int ,unsigned int> >::iterator it = nodes.begin() ; it!= nodes.end() ; it++){
			std::vector<std::string> parts;
			strsep(it->first,":",parts);
			std::cout<< "between node " << parts.at(0) << " and "<< parts.at(1) << " : from " << it->second.first << " to " << it->second.second <<std::endl;
		}

	}
	void print_graph(){
		for(std::map<unsigned int , std::vector<std::string> >::iterator it=clusters.begin(); it!= clusters.end(); it++){
			std::cout<< it->first<< std::endl;
			if(it->first==8){
				std::cout<<"if 8 "<<std::endl;
				std::cout<< it->second << std::endl;
			}
		}
	}
	const std::map<unsigned int , std::vector<std::string> > & get_clusters()const;
	const std::multimap<std::string , std::pair<size_t, size_t> > & get_intervals()const;
	const std::map<int,std::pair<int, int> > & get_nonaligned()const;
	const unsigned int get_node_length(unsigned int &)const;
	const std::multimap<size_t , const pw_alignment> & get_als()const{
		return al_from_graph;
	}
	private:
	std::map<unsigned int , std::vector<std::string> > clusters;
	std::multimap<std::string , std::pair<size_t, size_t> > boundries;//Multimap because several part of a sequence can be clustered together in one cluster.
	std::multimap<std::string , std::pair<size_t, size_t> > test_boundries;
	std::multimap<std::string , std::pair<unsigned int ,unsigned int> > nodes;
	std::map<int,std::pair<int, int> > non_aligned; //cluster id, pair(position, length)
	std::map<unsigned int, unsigned int> node_length;
	std::multimap<size_t , const pw_alignment> al_from_graph;

};


class test_sim_reads_mapping{
	public:
		test_sim_reads_mapping(map_check & mp):mp_check(mp){}
		~test_sim_reads_mapping(){}
		void read_sim_maffile(std::ifstream & sim_maffile, std::map<std::string , std::pair<size_t , size_t> > & onreads);
		void this_part_position_on_seq(bool & , size_t & , size_t & , size_t & , size_t &, size_t &, size_t & );
		void check_output(std::ifstream & mapping_maf);
		void check_if_shifted(unsigned int & ref_id, bool & IsMapped);
		void find_position_on_member(unsigned int & center, size_t & begin_on_center , size_t & end_on_center,size_t & begin_on_mem , size_t & end_on_mem,int & position);
	private:
		map_check mp_check;
		std::map<size_t,std::pair<size_t,size_t> > reads_position_on_seq;
		std::map<size_t , const pw_alignment> alignments;

};

class test_reveal{
	public:
		test_reveal(){};
		~test_reveal(){};
		void read_the_result(std::ifstream &);
		void compare_with_path(std::ifstream &);//TOO specific!!
	private:
		std::vector<std::string> nodes;
};
	
#include "needleman_wunsch.cpp"
