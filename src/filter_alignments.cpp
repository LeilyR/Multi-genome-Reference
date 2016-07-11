#ifndef FILTER_ALIGNMENTS_CPP
#define FILTER_ALIGNMENTS_CPP

#include "filter_alignments.hpp"

	void filter_als::find_overlapped_references(){
#pragma omp parallel for num_threads(num_threads)
		for(size_t i =0; i < alignments.size();i++){
			//For each reference of an al check the interval tree and find all the alignments that have more than 50% overlap with that ref. Afterwards check if the obtain result from both ref have overlap with each other, if so then count those pieces. 
				const pw_alignment * p = alignments.at(i);
				size_t r1,r2,l1,l2;
				p->get_lr1(l1,r1);
				p->get_lr2(l2,r2);
				size_t ref1 = p->getreference1();
				size_t ref2 = p->getreference2();
				std::vector<const pw_alignment*> results1;
				std::vector<const pw_alignment*> results2;
				std::map<const pw_alignment*,size_t> not_touched_ref1;
				std::map<const pw_alignment*,size_t> not_touched_ref2;
#pragma omp critical(al_insert)
{
				alind->search_overlap_fraction(ref1, l1, r1, FRACTION, results1);//TODO change th function later (>=50% of p's length should be covered maybe also save the ref that has overlap)
				alind->search_overlap_fraction(ref2, l2, r2, FRACTION, results2);//TODO change the function later

				find_touched_ref(results1, ref1, l1, r1, not_touched_ref1);
				find_touched_ref(results2, ref2, l2, r2, not_touched_ref2);
}
				size_t no_of_overlap = get_overlap(not_touched_ref1,not_touched_ref2);
				recalculate_gain_value(p,no_of_overlap);
		}
		find_als_with_highest_gain();
	}
	void filter_als::find_touched_ref(std::vector<const pw_alignment*> & results, size_t & ref, size_t & left, size_t & right, std::map<const pw_alignment*,size_t> & not_touched_ref){
		for(size_t i =0; i < results.size(); i++){
			const pw_alignment * p = results.at(i);
			size_t r1,r2,l1,l2;
			p->get_lr1(l1,r1);
			p->get_lr2(l2,r2);
			size_t ref1 = p->getreference1();
			size_t ref2 = p->getreference2();
			if(ref1 == ref && l1 <= right && r1 >= left){
				not_touched_ref.insert(std::make_pair(p,1));
			//	if(ref2 != ref || l2 > right || r2< left){ 		
			//		p->print();
			//		std::cout << "ref "<< ref << " left "<< left << " right "<< right << std::endl;
			//	}
			//	assert(ref2 != ref || l2 > right || r2< left);//XXX why is it wrong?!

			}else{
				assert(ref2 == ref && l2 <= right && r2 >= left);
				not_touched_ref.insert(std::make_pair(p,0));

			}

		}
	}

	size_t filter_als::get_overlap(std::map<const pw_alignment*, size_t> & ntouched1 , std::map<const pw_alignment*,size_t> & ntouched2){
		alignment_index alind1(num_sequences);
		for(std::map<const pw_alignment*,size_t>::iterator it=ntouched1.begin();it!= ntouched1.end();it++){//Make a tree for the first one
			alind1.half_insert(it->first, it->second);
		}
		std::vector<const pw_alignment *> result;
		for(std::map<const pw_alignment*,size_t>::iterator it=ntouched2.begin();it!= ntouched2.end();it++){//search on that tree for the second one
			const pw_alignment * al = it->first;
			size_t left,right;
			if(it->second == 0){
				size_t ref1 = al->getreference1();
				al->get_lr1(left,right);
				alind1.search_overlap(ref1, left, right,result);
			}else{
				assert(it->second == 1);
				size_t ref2 = al->getreference2();
				al->get_lr2(left,right);
				alind1.search_overlap(ref2, left, right, result);
			}

		}
		return result.size();
	}
	void filter_als::recalculate_gain_value(const pw_alignment* p,size_t & number){// recalcualte all the als gain value and add them to 'new_gains' 
		double g1 ,g2;
		model.gain_function(*p,g1,g2);
		double av_gain = (g1+g2)/2;
		double new_gain = av_gain*(1+number/2);
		new_gains.insert(std::make_pair(new_gain,p));

	}

	void filter_als::find_als_with_highest_gain(){//if the gain is higher that mean gain i keep it
		size_t counter = 0;
		double median = alignments.size()/2;
		for(std::multimap<double, const pw_alignment*>::reverse_iterator rit = new_gains.rbegin(); rit != new_gains.rend();rit++){
		//	std::cout << rit->first << std::endl;
			if(counter < median){
				filtered_als.insert(rit->second);
				counter++;
			}else break;
		}
	//	std::cout << "filtered size "<< filtered_als.size() << " median " <<median<<std::endl;
		assert(filtered_als.size()==median);

	}
	std::set<const pw_alignment*, compare_pointer_pw_alignment> filter_als::get_filtered_als()const{
		return filtered_als;
	}



#endif
