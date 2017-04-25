
#ifndef MERGING_CENTERS_HPP
#define MERGING_CENTERS_HPP

#include "pw_alignment.hpp"
#include"data.hpp"
#include "suffix_tree.hpp"

#include <map>
#include <vector>
#include <cassert>
#include <omp.h>


class merging_centers{
	public:
		merging_centers(const all_data & d, finding_centers & cent):data(d), centers(cent), all_current_centers(data.numSequences()){
		
		}
		~merging_centers(){

		}
		void updating_centers(std::vector<int> & , int &);
		void merg_gain_value( suffixTree & );
		void find_highest_gain(std::vector<int> & , int &, std::vector<int> & );
		void adding_new_centers(std::vector<std::vector<std::string> >&, std::vector<std::map<size_t , std::vector<std::string> > >&, std::vector<std::map<size_t, std::string> > &);// Saves new centers in the 'long_centers' vector in main. It gives a unique id to each of them(Their row in the outer vector). vector<string> includes are the centers in the new one. they are saved in the form of ref:left.
		void index_centers(std::map<std::string , std::vector<pw_alignment> > & );
		void add_long_centers_to_map(std::vector<std::vector<std::string> > & , std::map<vector<std::string>, std::vector<pw_alignment> > &); //Initially fills in the new_centers map
		void create_alignment(std::vector<std::vector<std::string> > &, std::map<vector<std::string>, vector<pw_alignment> >&, std::map<std::string , std::vector<pw_alignment> > &,  std::vector<std::map<size_t , std::vector<std::string> > > & , std::vector<std::map<size_t , std::string> > &, size_t &); //Creates all the als using long centers, first ref is always forward!
		void remove_fully_reverse_refs(vector<vector<std::string> > & ,std::map<std::vector<std::string>, std::vector<pw_alignment> > &);
		void find_new_centers(size_t &, std::vector<std::string> & , size_t &, std::vector<std::map<size_t, std::vector<std::string> > >& );//Finds long centers on a certain sequence	
		void find_long_center_length(std::vector<std::string> & , std::map<std::string,std::vector<pw_alignment> > & , size_t & ,size_t & ,size_t &,size_t &, size_t &);
		void set_samples(std::vector<pw_alignment>& ,size_t & , std::vector<bool> & , std::vector<bool> & , size_t & , size_t & , unsigned int & , unsigned int & , size_t &, size_t &);
		void get_merged_als(std::map<std::string, std::vector<pw_alignment> > & merged_als){
			for(std::map<std::string, std::vector<pw_alignment> >::iterator it = merging_als.begin() ; it !=merging_als.end(); it++){
				merged_als.insert(std::make_pair(it->first,it->second));
			}	
		}

	private:
		const all_data & data;
		finding_centers & centers;
		std::map<std::string, int> center_index;
		std::map<std::vector<int> , size_t> merged_centers;//merged centers, vector<int> ----> original indices of megerd centers , size_t ----> its new index
		std::map<std::vector<int>, int >gains; // vector<int> ----> series of centers, int ----> gain of it. It includes the last combination of centers
		std::vector<std::vector<std::vector<int> > >all_current_centers;//All the potential long centers on each sequence (negative values says that their reverse complement is located in that position(in fact it says how the edges should go))
		std::map<std::string, std::vector<pw_alignment> >merging_als;
};

class mixing_centers{
	public:
		mixing_centers(all_data & d):data(d),all_pieces_long_centers(data.numSequences()), all_centers_on_a_seq(data.numSequences()),path_on_a_seq(data.numSequences()){
			
		}
		~mixing_centers(){

		}
		void add_original_centers_to_long_centers(std::map<std::vector<std::string>, std::vector<pw_alignment> > & , std::map<std::string, std::vector<pw_alignment> > & , std::map<std::string, std::vector<pw_alignment> > &, std::vector<std::map<size_t, std::vector<std::string> > > &, std::vector<std::map<size_t, std::string> > & );
		void swap_references(std::map<std::vector<std::string>, std::vector<pw_alignment> > &);
		void add_to_long_center_map(const std::string & center, const pw_alignment &, std::vector<std::map<size_t, std::string> > &, std::map<std::vector<std::string>, std::vector<pw_alignment> > &);
		void get_direction(const pw_alignment & p, size_t & dir);
		void calculate_centers_weight(std::map<std::string, std::vector<pw_alignment> > & , std::map<std::string, unsigned int> &);
		void calculate_long_centers_weight(std::map<std::vector<std::string>, std::vector<pw_alignment> > & , std::map<std::vector<std::string> , unsigned int> &);
		void find_members_of_clusters(std::map<std::string,std::string> &, std::map<std::string, std::vector<pw_alignment> > & );
		void find_adjacencies( std::map<std::string, std::vector<pw_alignment> > & , std::vector<std::multimap<size_t, std::string > > & , std::map<size_t, std::set<size_t> > &);
		void compute_centers_length(std::map<std::string, std::vector<pw_alignment> > & , std::map<std::string , size_t> & );
		const size_t find_right_on_seq(std::string & ,std::map<std::string, std::vector<pw_alignment> > &, size_t & , size_t &, std::map<std::string, size_t> & )const;
		const size_t find_right_on_seq(std::string & center,std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & left , size_t & ref, std::map<std::string, size_t> & non_aligned_right, pw_alignment & al, size_t & gap_on_long_centers)const;
		const size_t find_right_of_long_center_on_seq(std::vector<std::string>& ,  std::map<vector<std::string>, std::vector<pw_alignment> > & , size_t & , size_t &, const std::map<std::string,std::vector<pw_alignment> > &)const;
		void add_nonaligned_regions(std::vector<std::multimap<size_t, std::string > > & , std::map<std::string, std::vector<pw_alignment> > & );
		void remove_inclusive_long_centers(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster , std::vector<std::map<size_t, std::vector<std::string> > > & long_centers_on_seq, std::map<vector<std::string>, std::vector<pw_alignment> > & new_centers, std::map<std::string, size_t> & non_aligned_right, size_t & gap_on_long_centers);
		void make_index(std::vector<std::map<size_t , std::string> > & all_pieces);
		void make_long_cent_index(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::vector<std::string>, std::vector<pw_alignment> > & long_centers, std::map<std::string, size_t> & non_aligned_right);

		void find_adjacencies_with_direction( std::map<int, std::set<int> > &, std::vector<std::map<size_t , std::string> > & all_pieces);
		void find_adjacencies_with_direction_long_centers(std::map<std::string, std::vector<pw_alignment> > & , std::map<int, std::set<int> > &, std::map<std::vector<std::string>, std::vector<pw_alignment> > &, std::map<std::string, size_t> & non_aligned_right);
		void find_direction_of_center(size_t & this_ref, size_t & this_left ,size_t & direction , std::vector<std::string> & center, std::map<std::vector<std::string>, std::vector<pw_alignment> > & long_centers);

		bool find_coordinate(std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, size_t & position , size_t & ref, std::string & center);
		void test(std::map<int, std::set<int> > & adjacencies, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::vector<std::map<size_t , std::string> > & all_pieces, std::map<std::string, size_t>& non_aligned_right);
		const std::map<std::string, size_t> get_center_id()const;
		const std::map<std::vector<std::string>, size_t> get_long_center_id()const;
		const size_t get_right(std::string & , std::map<std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right)const;
                std::vector<std::string> get_reverse(std::vector<std::string>&)const;
		void swap_refs(pw_alignment & al , const std::string & center);
		bool is_forward(const std::string & center , pw_alignment & p);
		const std::vector<std::vector<int> > get_path()const{
			return path_on_a_seq;
		}
	private: 
		const all_data & data;
		std::map<std::string, std::string> non_aligned_context;//The first string is the context and the second string its name(ref:left). It is used for both short or long centers.
		std::map<std::string, std::vector<std::string> > non_aligned_clusters; //The first string is the name of the non aligned piece that is chosen as the "center" of some non laigned regions that have the exact same context as it has, the vector contains the name of all those others with the same context.
	//	std::map<std::string, size_t> non_aligned_right;

		std::vector<std::map<size_t, std::vector<std::string> > >  all_centers_on_a_seq;
		//XXX short centers:
	//	std::vector<std::multimap<size_t , std::string> >all_pieces;//contains short centers and non aligned regions between them.
		std::map<std::string, size_t> center_id;//includes short centers and non aligned regions index

		//XXX long centers:
		std::vector<std::multimap<size_t , std::vector<std::string> > >all_pieces_long_centers;//contains long centers and non aligned regions between them.(pos,center)
		std::map<std::vector<std::string>, size_t> long_center_id;//(center,id)
		//TODO make another container like all_pieces_long_centers and save the content of each piece
		std::vector<std::vector<int> > path_on_a_seq;

};


class write_graph{

	public:
	write_graph(const all_data & d, mixing_centers & m):data(d), mixcenter(m){}
	~write_graph(){}
	void write_graph_dot(std::map<int, std::set<int> > & adjacencies , std::ostream & dotfile);

	void write_graph_fasta(const std::string & graphout,std::ofstream & txtout, std::map< std::string , std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right);
	void write_graph_fasta_with_long_centers(const std::string & graphout, std::ofstream & txtout,std::map<vector<std::string>, std::vector<pw_alignment> > & long_centers, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right);

	void write_maf_record(std::ostream & out, const std::string & src, size_t start, size_t size, char strand, size_t srcSize, const std::string & alignment_part);
	void write_maf_record(std::ostream & out, const all_data & data, const pw_alignment & al, size_t reference);
	void msa_star_alignment(const std::string & center, std::vector<pw_alignment> & alignments);
	void msa_star_al_with_long_centers(const std::vector<std::string> & longCenter, std::vector<pw_alignment> & alignment);
	void write_graph_maf(std::ofstream & gout, const std::map<std::string, std::vector<pw_alignment> > & cluster_result_al, std::map<std::string, size_t> & center_id);
	void write_graph_maf_with_long_centers(std::ofstream & gout ,const std::map<std::vector<std::string>, std::vector<pw_alignment> > & new_centers,std::map<std::vector<std::string>, size_t>& long_center_id);
	void write_graph_gfa(std::map<int, std::set<int> > & adjacencies, std::ostream & graph_gfa, std::map<std::string, std::vector<pw_alignment> > & alignments_in_a_cluster, std::map<std::string, size_t> & non_aligned_right);//TODO
	void make_fasta(std::ostream & graphfasta, std::map<std::string, std::string> & centers) ;
	

	private:
		const all_data & data;
		const mixing_centers & mixcenter;


};
#endif
