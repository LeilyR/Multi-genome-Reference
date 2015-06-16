#ifndef MODEL_HPP
#define MODEL_HPP

#include"data.hpp"
#include "pw_alignment.hpp"
#include <map>
#include <vector>
#include <cassert>

extern "C" {
#include "apcluster.h"
}



typedef  std::set<const pw_alignment*, compare_pw_alignment> alset;

template<typename T>
class initial_alignment_set {
	public:
	initial_alignment_set(const all_data & d, const T & a_model, double base_cost): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * cur = &(data.getAlignment(i));
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
			double vgain = (gain1+gain2)/2 - base_cost;
		//	std::cout << " al " << i << " gain1 " << gain1 << std::endl;
			if(vgain > 0.0) {
				sorter.insert(std::make_pair(vgain, cur));
				sumgain+=vgain;
			}
		}

		sorted_original_als = std::vector<const pw_alignment*>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
		//	std::cout << " ral " << alit << std::endl;
		//	alit.print();
			std::cout << std::endl;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}
	initial_alignment_set(const all_data & d, const std::set< const pw_alignment *, compare_pw_alignment> & als, const T & a_model, double base_cost): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(std::set< const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * cur = *it;
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
			double vgain = (gain1+gain2)/2 - base_cost;
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al length " << cur->alignment_length() << " gain1 " << gain1 << " gain2 " << gain2 <<  std::endl;
			if(vgain>0.0) {
			//	std::cout << " ins " << vgain << " at " << cur << std::endl;
				sorter.insert(std::make_pair(vgain, cur));

				/*
				common_model.gain_function(*cur, gain1, gain2);
				vgain = (gain1+gain2)/2 - base_cost;
				std::cout << " cached gian "<< vgain << std::endl;
				assert(vgain>=0);
				*/

			}
			sumgain+=vgain;
		}

		sorted_original_als = std::vector<const pw_alignment *>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
		//	std::cout << " ins2 " << pos << " weight " << rit->first << " at " << alit <<  std::endl;
			
			/*
			double g1, g2;
			common_model.gain_function(*alit, g1, g2);
			double vgain = (g1+g2)/2 - base_cost;
			std::cout << " cached gain " << vgain << std::endl;
			assert(vgain >=0);
*/

			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}

	~initial_alignment_set() {}

	void compute(overlap & o);
	void compute_simple(overlap & o);

/*
	This method starts with a weighted INDEPENDENT SET of alignments without partial overlap which is computed
	using a fast and easy to implement VERTEX COVER 2-approximation 
	(Clarkson, KL. A modification of the greedy algorithm for vertex cover. Information Processing Letters. 1983.)
*/
	void compute_vcover_clarkson(overlap & o);
	void compute_simple_lazy_splits(overlap & o);
	void lazy_split_insert_step(overlap & ovrlp, size_t level, const pw_alignment * al, std::vector<const pw_alignment*> & inserted_alignments, vector<const pw_alignment*> & removed_alignments, double & local_gain);

	double get_max_gain() const {
		return max_gain;
	}
	double get_result_gain() const {
		return result_gain;
	}

	size_t get_used_alignments() const {
		return used_alignments;
	}

	private:
	const all_data & data;
	const T & common_model;
	std::vector<const pw_alignment*> sorted_original_als; // highest gain first
	double base_cost;
	double max_gain;
	double result_gain;
	size_t used_alignments;


	void update_remove(size_t index, std::vector<std::set<size_t> > & edges, std::map<size_t, double> & index_to_weight, std::multimap<double, size_t> & weight_to_index);
	void remove_val(size_t index, std::map<size_t, double> & index_to_weight,  std::multimap<double, size_t> & weight_to_index);
	void overlap_graph(std::vector<std::set<size_t> > & edges);
	void search_area(size_t from, size_t to, const std::multimap<size_t, size_t> & data, std::set<size_t> & res);
	
};



/*
template<typename T>
class initial_alignment_set_bb {
	public:
	initial_alignment_set(const all_data & d, const T & a_model, double base_cost,std::ofstream & outs): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * cur = &(data.getAlignment(i));
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2,outs);
			// TODO gain1 and gain2
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al " << i << " gain1 " << gain1 << std::endl;
			if(gain1-base_cost > 0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = std::vector<const pw_alignment*>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
		//	std::cout << " ral " << alit << std::endl;
		//	alit.print();
			std::cout << std::endl;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}
	initial_alignment_set(const all_data & d, const set< const pw_alignment *, compare_pw_alignment> & als, const T & a_model, double base_cost, std::ofstream & outs): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		std::multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(std::set< const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * cur = *it;
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2,outs);
		//	if(gain2 > gain1) gain1 = gain2;
		//	std::cout << " al length " << cur->alignment_length() << " gain1 " << gain1 << " gain2 " << gain2 <<  std::endl;
			if(gain1 - base_cost>0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = std::vector<const pw_alignment *>(sorter.size(), NULL);
		size_t pos = 0;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
	//	std::cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << std::endl;
	}

	~initial_alignment_set() {}

	void compute(overlap & o, std::ofstream &);
	void compute_simple(overlap & o,std::ofstream &);

	double get_max_gain() const {
		return max_gain;
	}
	double get_result_gain() const {
		return result_gain;
	}

	private:
	const all_data & data;
	const T & common_model;
	std::vector<const pw_alignment*> sorted_original_als; // highest gain first
	double base_cost;
	double max_gain;
	double result_gain;
};

*/

class compute_cc {
	public:
	compute_cc(const std::set<const pw_alignment *, compare_pw_alignment> & als_in, size_t num_sequences):alignments(als_in), als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		max_al_ref_length = 0;

		for(std::set<const pw_alignment *, compare_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); ++it) {
			const pw_alignment * al = *it;
			add_on_mmaps(al);
		}
	
	}

	compute_cc(const overlap & ovrlp, size_t num_sequences): als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		max_al_ref_length = 0;
		const std::set<pw_alignment*, compare_pw_alignment> & als = ovrlp.get_all();
		for(std::set<pw_alignment*, compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * al = *it;
			alignments.insert(al);
			add_on_mmaps(al);
		}
	
	}


	compute_cc(const all_data & dat);
	~compute_cc() {}

	void compute(std::vector<std::set< const pw_alignment *, compare_pw_alignment> > & ccs); 

	private:
	std::set<const pw_alignment *, compare_pw_alignment> alignments;
	std::vector< std::multimap< size_t, const pw_alignment *> > als_on_reference; // sequence index -> pos on that sequence -> alignment reference
	std::vector<size_t> last_pos;
	size_t max_al_ref_length;

	void add_on_mmaps(const pw_alignment * pwa);
	void remove_on_mmaps(const pw_alignment * al);
	void get_cc(const pw_alignment * al, std::set <const pw_alignment *, compare_pw_alignment> & cc, std::set <const pw_alignment *, compare_pw_alignment> & seen);
	void cc_step(size_t ref, size_t left, size_t right, std::set <const pw_alignment *, compare_pw_alignment> & cc, std::set <const pw_alignment *, compare_pw_alignment> & seen );

};

template<typename tmodel>
class clustering {
	public:
	clustering(overlap &,all_data &, tmodel&);
	~clustering();
	void als_on_reference(const pw_alignment * p);
	void calculate_similarity();//Fill in a matrix with gain values for each two accessions.
	void update_values(); 
	void update_clusters(size_t acc);
	void update_clusters();
	private:
	overlap & overl;
	all_data & data;
	tmodel & model;
//	std::set<const pw_alignment *, compare_pw_alignment> alignments;
	std::vector< std::multimap<size_t, pw_alignment *> > als_on_ref;
	std::vector<vector<vector<double> > >gain;//consider it as similarity function, though in your slides you mentioned modification can be the similarity function.
	std::vector<vector<vector<double> >	>ava;//availability		
	std::vector<vector<vector<double> >	>res;//responsibilty
	
};

/*
	TODO 
	can the clustering result become better if we estimate modification costs of those pairs which have no direct pairwise alignment


*/

template<typename tmodel>
class affpro_clusters {
	public:

	affpro_clusters(const overlap & ovlp, const tmodel & model, double base_cost,std::ofstream& outs):model(model), base_cost(base_cost) {

		const std::set<pw_alignment*, compare_pw_alignment> & als = ovlp.get_all();
		for(std::set<pw_alignment*, compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * al = *it;
			add_alignment(al,outs);
		}

	}

	affpro_clusters(const std::set<const pw_alignment *, compare_pw_alignment> & inset, const tmodel & model, double base_cost,std::ofstream & outs):model(model), base_cost(base_cost) {
	//	std::cout<<"instd::set size"<<inset.size()<<std::endl;	
		//	std::cout<<"data1 ad in afp: "<< & dat << std::endl;	
		for(std::set<const pw_alignment*, compare_pw_alignment>::iterator it = inset.begin(); it!=inset.end(); ++it) {
		//	std::cout<<"data2 ad in afp: "<< & dat << std::endl;	
			const pw_alignment * al = *it;
		//	std::cout<<"alignment from instd::set: "<<std::endl;
		//	al->print();
		//	dat.numAcc();
		//	std::cout<<"data3 ad in afp: "<< & dat << std::endl;	
			add_alignment(al,outs); //Marion : I have an error that I dont understand !
		//	std::cout<<"data4 ad in afp: "<< & dat << std::endl;	

		}
	}


void run(std::map<std::string, std::vector<std::string> > & cluster_result) {
	// convert to c-style matrix
	double * data = new double[simmatrix.size() * simmatrix.size()];
	int * result = new int[simmatrix.size()];
	if(simmatrix.size()>2) {
	for(size_t i=0; i<simmatrix.size(); ++i) {
		for(size_t j=0; j<simmatrix.size(); ++j) {
			data[i*simmatrix.size() + j] = simmatrix.at(i).at(j);
		}
	}

	APOPTIONS apoptions;
	apoptions.cbSize = sizeof(APOPTIONS);
	apoptions.lambda = 0.9;
	apoptions.minimum_iterations = 1;
	apoptions.converge_iterations = 20;
	apoptions.maximum_iterations = 5000;
	apoptions.nonoise = 1;
	apoptions.progress=NULL; apoptions.progressf=NULL;
	double netsim;	
	int iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
//	std::cout << "iter " << iter << std::endl;
	if(iter <= 0 || result[0]==-1) {
		apoptions.converge_iterations = 10;
		apoptions.maximum_iterations = 15000;
		apoptions.lambda = 0.98;
		iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
//		std::cout << "iter " << iter << std::endl;
	}
	if(iter <= 0 || result[0]==-1) {
		apoptions.converge_iterations = 10;
		apoptions.maximum_iterations = 15000;
		apoptions.lambda = 0.995;
		iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
	//	std::cout << "iter " << iter << std::endl;

	}
	} else {
		if(simmatrix.size()==1) {
			result[0] = 0;
		} else {
		// simpler algorithm for only 2 element
			double one = simmatrix.at(0).at(0) + simmatrix.at(1).at(0);
			double two = simmatrix.at(1).at(1) + simmatrix.at(0).at(1);
			double separate = simmatrix.at(0).at(0) + simmatrix.at(1).at(1);

			if(one > two && one > separate) {
				result[0] = 0;
				result[1] = 1;
			} else if(two > one && two > separate) {
				result[0] = 1;
				result[1] = 1;
			} else {
				result[0] = 0;
				result[1] = 1;
			}

		}
	}
// TODO conv error res -1

	double totalccost = 0;
	for(size_t i=0; i<simmatrix.size(); ++i) {
	//	if(simmatrix.at(i).at(result[i])== -HUGE_VAL){
	//		totalccost = i;
	//	}else{
		totalccost -= simmatrix.at(i).at(i);
	//	}
	}

	double apcost = 0;
//	std::cout<<"size of simmatrix: "<<simmatrix.size()<<std::endl;
	for(size_t i=0; i<simmatrix.size(); ++i) {
		if(result[i]==-1) {
			result[i] = i;
		}
		assert(result[i]>= 0);
		if(i==(size_t)result[i]) {
			apcost-=simmatrix.at(i).at(i);
		} else {
			if(simmatrix.at(i).at(result[i])== -HUGE_VAL){
				result[i]=i;
				apcost-=simmatrix.at(i).at(i);
			}else{
				apcost-=simmatrix.at(i).at(result[i]);
			}
		}
	}
	for(size_t i=0; i<simmatrix.size(); ++i) {//cluster center may happen whenever i == result[i]
//		std::cout << sequence_names.at(i) << " res " << i << " is " << result[i] << " ( length " << sequence_lengths.at(i) << ")"<<std::endl;

		if( (size_t)result[i]==i) {
		//	std::cout << " " << simmatrix.at(i).at(i) << std::endl;
			std::map<std::string, std::vector<std::string> >::iterator it=cluster_result.find(sequence_names.at(i));
			if(it==cluster_result.end()) {
				cluster_result.insert(make_pair(sequence_names.at(i), std::vector<std::string>()));
			} else {
				// it->second.push_back(sequence_names.at(i));
			
			}


		} else {
		//	std::cout << " " << simmatrix.at(i).at(result[i]) << " : " << simmatrix.at(i).at(i) << std::endl;
			//result[i]is associated one, add them to the map for each center
			std::map<std::string, std::vector<std::string> >::iterator it=cluster_result.find(sequence_names.at(result[i]));
			if(it == cluster_result.end()) {
				cluster_result.insert(make_pair(sequence_names.at(result[i]), std::vector<std::string>()));
				it=cluster_result.find(sequence_names.at(result[i]));
			}
			it->second.push_back(sequence_names.at(i));

		}
	}

	// double apgain = totalccost - apcost;

//	std::cout << "Total sequence cost " << totalccost << " ap clustering cost " << apcost << " gain: " << apgain << std::endl;

	delete [] data;
	delete [] result;
	data = NULL;
	result = NULL;
//	std::cout<< "number of centers: " <<clusterCenter.size()<<std::endl;
}
	size_t get_sequence_length(size_t ref_idx)const; //ref_idx shows the reference that sequence belongs to. It could be either 0 or one.
	


// simil = neg distance, diagonale pref = neg cost
	private:
//	all_data & dat;
	const tmodel & model;
	double base_cost;
	// TODO do we need distances for all pairs of sequence pieces?
	std::map<std::string, size_t> sequence_pieces; // chr:leftpos -> seq_piece index
	std::vector<std::string> sequence_names;
       	std::vector<size_t> sequence_lengths;	
	std::vector<std::vector<double> > simmatrix;
	std::map<std::string, char> cluster_centers;
	void add_alignment(const pw_alignment *al, std::ofstream &);
};

class finding_centers{
	public:
	finding_centers(all_data &);
	~finding_centers();
	void setOfAlignments(std::map<std::string,std::vector<pw_alignment> > &);
	void findMemberOfClusters(std::map<std::string,std::vector<pw_alignment> > & );
	void center_frequency(std::map<std::string,std::vector<pw_alignment> > &);
	std::vector<size_t> get_center(size_t)const;

	private:
	all_data & data;
	std::vector< std::multimap<size_t , pw_alignment*> >AlignmentsFromClustering;
	std::map<std::string,std::string> memberOfCluster; //first string is assocciated member and second one is its center
	std::vector<std::vector<size_t> >centersOfASequence;//all the centers that happen on each sequence.



};
class suffix_tree{
	public:
	suffix_tree(all_data &, finding_centers&);
	~suffix_tree();
	void create_suffix(size_t);
	void find_a_node(size_t& ,size_t&, std::string&);
	void find_next_sibling(size_t, size_t);
	void find_next_firstparent(size_t, size_t);
	void make_a_tree();
	void insert_node(std::string);
	void count_paths();
	std::vector<std::string> get_nodes()const;
	std::map<std::vector<size_t>, size_t> get_count()const;
	std::vector<size_t> get_first_parent()const;
	
	private:
	all_data & data;
	finding_centers & centers;
	std::vector<std::string>suffixes;//all the suffixes of a sequence at the time
	std::vector<std::string> nodes;
	std::multimap<size_t , size_t> nodes_relation; //first size_t shows the parent and second one shows kids
	std::map<std::string, size_t>firstParent; // size_t shows its index
	std::map<std::vector<size_t>,size_t>branch_counter;//vector represents all nodes of a branch , size_t shows the number of that branch is happening
		
};

class merging_centers{
	public:
		merging_centers(finding_centers &, suffix_tree &);
		~merging_centers();
		void merg_gain_value();
		void merg_alignments(std::map<std::string,std::vector<pw_alignment> > &);
	private:
		finding_centers & centers;
		suffix_tree & tree;
		std::vector<std::vector<size_t> > merged_centers;//merged nodes




};




#endif



