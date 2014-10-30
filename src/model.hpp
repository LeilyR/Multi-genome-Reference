#ifndef MODEL_HPP
#define MODEL_HPP

#include"data.hpp"
#include "pw_alignment.hpp"
#include <map>
#include <vector>

#include "data.hpp"


extern "C" {
#include "apcluster.h"
}

using namespace std;


typedef  set<const pw_alignment*, compare_pw_alignment> alset;

template<typename T>
class initial_alignment_set {
	public:
	initial_alignment_set(const all_data & d, const T & a_model, double base_cost): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * cur = &(data.getAlignment(i));
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
			// TODO gain1 and gain2
		//	if(gain2 > gain1) gain1 = gain2;
			cout << " al " << i << " gain1 " << gain1 << endl;
			if(gain1-base_cost > 0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = vector<const pw_alignment*>(sorter.size(), NULL);
		size_t pos = 0;
		for(multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
		//	cout << " ral " << alit << endl;
		//	alit.print();
			cout << endl;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
		cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << endl;
	}
	initial_alignment_set(const all_data & d, const set< const pw_alignment *, compare_pw_alignment> & als, const T & a_model, double base_cost): data(d), common_model(a_model) {
		this->base_cost = base_cost;
		multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(set< const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * cur = *it;
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
		//	if(gain2 > gain1) gain1 = gain2;
			cout << " al length " << cur->alignment_length() << " gain1 " << gain1 << " gain2 " << gain2 <<  endl;
			if(gain1 - base_cost>0.0) {
				sorter.insert(make_pair(gain1, cur));
			}
			sumgain+=gain1;
		}

		sorted_original_als = vector<const pw_alignment *>(sorter.size(), NULL);
		size_t pos = 0;
		for(multimap<double, const pw_alignment*>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			const pw_alignment * alit = rit->second;
			sorted_original_als.at(pos) = alit;
			pos++;
		}
		max_gain = sumgain - sorter.size() * base_cost;
		assert(pos == sorter.size());
		cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << endl;
	}

	~initial_alignment_set() {}

	void compute(overlap & o);
	void compute_simple(overlap & o);

	double get_max_gain() const {
		return max_gain;
	}
	double get_result_gain() const {
		return result_gain;
	}

	private:
	const all_data & data;
	const T & common_model;
	vector<const pw_alignment*> sorted_original_als; // highest gain first
	double base_cost;
	double max_gain;
	double result_gain;
};




class compute_cc {
	public:
	compute_cc(const set<const pw_alignment *, compare_pw_alignment> & als_in, size_t num_sequences):alignments(als_in), als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		max_al_ref_length = 0;

		for(set<const pw_alignment *, compare_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); ++it) {
			const pw_alignment * al = *it;
			add_on_mmaps(al);
		}
	
	}

	compute_cc(const overlap & ovrlp, size_t num_sequences): als_on_reference(num_sequences), last_pos(num_sequences, 0) {
		max_al_ref_length = 0;
		const set<pw_alignment*, compare_pw_alignment> & als = ovrlp.get_all();
		for(set<pw_alignment*, compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * al = *it;
			alignments.insert(al);
			add_on_mmaps(al);
		}
	
	}


	compute_cc(const all_data & dat);
	~compute_cc() {}

	void compute(vector<set< const pw_alignment *, compare_pw_alignment> > & ccs); 

	private:
	set<const pw_alignment *, compare_pw_alignment> alignments;
	vector< multimap< size_t, const pw_alignment *> > als_on_reference; // sequence index -> pos on that sequence -> alignment reference
	vector<size_t> last_pos;
	size_t max_al_ref_length;

	void add_on_mmaps(const pw_alignment * pwa);
	void get_cc(const pw_alignment * al, set <const pw_alignment *, compare_pw_alignment> & cc, set <const pw_alignment *, compare_pw_alignment> & seen);
	void cc_step(size_t ref, size_t left, size_t right, set <const pw_alignment *, compare_pw_alignment> & cc, set <const pw_alignment *, compare_pw_alignment> & seen );

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
//	set<const pw_alignment *, compare_pw_alignment> alignments;
	vector< multimap<size_t, pw_alignment *> > als_on_ref;
	vector<vector<vector<double> > >gain;//consider it as similarity function, though in your slides you mentioned modification can be the similarity function.
	vector<vector<vector<double> >	>ava;//availability		
	vector<vector<vector<double> >	>res;//responsibilty
	
};



template<typename tmodel>
class affpro_clusters {
	public:

	affpro_clusters(const overlap & ovlp, const tmodel & model, double base_cost): model(model), base_cost(base_cost) {

		const set<pw_alignment*, compare_pw_alignment> & als = ovlp.get_all();
		for(set<pw_alignment*, compare_pw_alignment>::const_iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * al = *it;
			add_alignment(al);
		}

	}

	affpro_clusters(const set<const pw_alignment *, compare_pw_alignment> & inset, const tmodel & model, double base_cost):model(model), base_cost(base_cost) {
		for(set<const pw_alignment*, compare_pw_alignment>::iterator it = inset.begin(); it!=inset.end(); ++it) {
			const pw_alignment * al = *it;
			add_alignment(al);
		}
	}


void run() {
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
	cout << "iter " << iter << endl;
	if(iter <= 0 || result[0]==-1) {
		apoptions.converge_iterations = 10;
		apoptions.maximum_iterations = 15000;
		apoptions.lambda = 0.98;
		iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
		cout << "iter " << iter << endl;
	}
	if(iter <= 0 || result[0]==-1) {
		apoptions.converge_iterations = 10;
		apoptions.maximum_iterations = 15000;
		apoptions.lambda = 0.995;
		iter = apcluster32(data, NULL, NULL, simmatrix.size()*simmatrix.size(), result, &netsim, &apoptions);
		cout << "iter " << iter << endl;

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
		totalccost -= simmatrix.at(i).at(i);
	}

	double apcost = 0;
	for(size_t i=0; i<simmatrix.size(); ++i) {
		if(result[i]==-1) {
			result[i] = i;
		}
		if(i==result[i]) {
			apcost-=simmatrix.at(i).at(i);
		} else {
			apcost-=simmatrix.at(i).at(result[i]);
		}
	}
	for(size_t i=0; i<simmatrix.size(); ++i) {
		cout << sequence_names.at(i) << " res " << i << " is " << result[i] << " ( length " << sequence_lengths.at(i) << ")";

		if(result[i]==i) {
			cout << " " << simmatrix.at(i).at(i) << endl;
		} else {
			cout << " " << simmatrix.at(i).at(result[i]) << " : " << simmatrix.at(i).at(i) << endl;
		}
	}

	double apgain = totalccost - apcost;

	cout << "Total sequence cost " << totalccost << " ap clustering cost " << apcost << " gain: " << apgain << endl;

	delete [] data;
	delete [] result;


}


// simil = neg distance, diagonale pref = neg cost

	private:
	const tmodel & model;
	double base_cost;
	// TODO do we need distances for all pairs of sequence pieces?
	map<string, size_t> sequence_pieces; // chr:leftpos -> seq_piece index
	vector<string> sequence_names;
       	vector<size_t> sequence_lengths;	
	vector<vector<double> > simmatrix;
	map<string, char> cluster_centers;
	void add_alignment(const pw_alignment *al);
};


#endif



