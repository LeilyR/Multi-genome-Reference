#ifndef MODEL_HPP
#define MODEL_HPP

#include"data.hpp"
#include "pw_alignment.hpp"
#include <map>
#include <vector>

#include "data.hpp"

using namespace std;


typedef  set<const pw_alignment*, compare_pw_alignment> alset;

template<typename T>
class initial_alignment_set {
	public:
	initial_alignment_set(const all_data & d, const T & a_model, double base_cost): data(d), common_model(a_model), sorted_original_als(d.numAlignments()) {
		this->base_cost = base_cost;
		multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment * cur = &(data.getAlignment(i));
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
		//	if(gain2 > gain1) gain1 = gain2;

			sorter.insert(make_pair(gain1, cur));

			sumgain+=gain1;
		}

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
	initial_alignment_set(const all_data & d, const set< const pw_alignment *, compare_pw_alignment> & als, const T & a_model, double base_cost): data(d), common_model(a_model), sorted_original_als(als.size()) {
		this->base_cost = base_cost;
		multimap<double, const pw_alignment*> sorter;
		double sumgain = 0;
		for(set< const pw_alignment *, compare_pw_alignment>::iterator it = als.begin(); it!=als.end(); ++it) {
			const pw_alignment * cur = *it;
	
			double gain1, gain2;
			common_model.gain_function(*cur, gain1, gain2);
		//	if(gain2 > gain1) gain1 = gain2;

			sorter.insert(make_pair(gain1, cur));

			sumgain+=gain1;
		}

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
class clustering {
	public:
	clustering(overlap &,all_data &, mc_model&);
	~clustering();
	void als_on_reference(const pw_alignment * p);
	void calculate_similarity();//Fill in a matrix with gain values for each two accessions.
	void update_values(); 
	void update_clusters(size_t acc);
	void update_clusters();
	private:
	overlap & overl;
	all_data & data;
	mc_model & model;
//	set<const pw_alignment *, compare_pw_alignment> alignments;
	vector< multimap<size_t, pw_alignment *> > als_on_ref;
	vector<vector<vector<double> > >gain;//consider it as similarity function, though in your slides you mentioned modification can be the similarity function.
	vector<vector<vector<double> >	>ava;//availability		
	vector<vector<vector<double> >	>res;//responsibilty
	
};



#endif



