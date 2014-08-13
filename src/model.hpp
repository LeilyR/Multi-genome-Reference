#ifndef MODEL_HPP
#define MODEL_HPP


#include "pw_alignment.hpp"
#include "data.hpp"
#include <map>


template<typename T>
class initial_alignment_set {
	public:
	initial_alignment_set(const all_data & d, const T & model): data(d), model(model), sorted_original_als(d.numAlignments()) {
		multimap<double, const pw_alignment&> sorter;
		double sumgain = 0;
		for(size_t i=0; i<data.numAlignments(); ++i) {
			const pw_alignment & cur = data.getAlignment(i);
	
			double gain1, gain2;
			model.gain_function(cur, gain1, gain2);
			if(gain2 > gain1) gain1 = gain2;

			sorter.insert(make_pair(gain1, cur));

			sumgain+=gain1;
		}

		size_t pos = 0;
		for(multimap<double, const pw_alignment&>::reverse_iterator rit = sorter.rbegin(); rit!=sorter.rend(); ++rit) {
			sorted_original_als.at(pos) = &(rit->second);
			pos++;
		}
		cout << " " << sorter.size() << " input alignments, total gain: " << sumgain << " bit " << endl;
	}
	~initial_alignment_set() {}

	void compute(overlap & o);
	void compute_simple(overlap & o);

	private:
	const all_data & data;
	const T & model;
	vector<const pw_alignment*> sorted_original_als; // highest gain first
};








#endif



