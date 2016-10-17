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
/*#pragma omp parallel // TODO Try apply it for parallelization of maps! 
{
		size_t cnt =0;
		for(auto element = ntouched2.begin(); element != ntouched2.end(); element++, cnt++){//search on that tree for the second one
			if(cnt%num_threads != 0) continue;
			const pw_alignment * al = element->first;
			size_t left,right;
			if(element->second == 0){
				size_t ref1 = al->getreference1();
				al->get_lr1(left,right);
				alind1.search_overlap(ref1, left, right,result);
			}else{
				assert(element->second == 1);
				size_t ref2 = al->getreference2();
				al->get_lr2(left,right);
				alind1.search_overlap(ref2, left, right, result);
			}
		}
}*/
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



	void alignment_filter::estimate_optimist_gain_per_gain() {

		typedef dynamic_mc_model use_model;
		typedef overlap_interval_tree overlap_type;
		typedef initial_alignment_set<use_model,overlap_type> use_ias;


		std::multimap<double, const pw_alignment *> gain_to_ra;
		const size_t num_als = 30; // TODO this number can be decrease to run the program faster for debugging
		for(size_t i=0; i< 10*num_als; i++) {
			int r = rand();
			size_t ri = (size_t) ( (double)r / (double) RAND_MAX * alvec.size());
			const pw_alignment * al = alvec.at(ri);
			double g1, g2;
			model.gain_function(*al, g1, g2);
			double gain = (g1+g2)/2.0;
			gain_to_ra.insert(std::make_pair(gain, al));
		}

		std::vector<const pw_alignment*> random_als;
		for(std::multimap<double, const pw_alignment *>::reverse_iterator rit = gain_to_ra.rbegin(); rit!=gain_to_ra.rend(); ++rit) {
			random_als.push_back(rit->second);
			if(random_als.size()==num_als) break;
		}

		
		typedef std::multimap<const pw_alignment *, std::vector<const pw_alignment *>, compare_pointer_pw_alignment> maptype;
		maptype rand_to_overl;

		alignment_index small_ind(data.numSequences());
		for(size_t i=0; i<random_als.size(); ++i) {
			small_ind.insert(random_als.at(i));
			std::vector<const pw_alignment*> av;
			rand_to_overl.insert(std::make_pair(random_als.at(i), av));
		}

		
		for(size_t i=0; i<alvec.size(); ++i) {
			std::vector<const pw_alignment *> ires;
			small_ind.search_overlap(*alvec.at(i), 0, ires);
			small_ind.search_overlap(*alvec.at(i), 1, ires);

			for(size_t j=0; j<ires.size(); ++j) {
				typename maptype::iterator fi = rand_to_overl.find(ires.at(j));
				assert(fi!=rand_to_overl.end());
				(fi->second).push_back(alvec.at(i));
			}
		}

		double sum_of_gain = 0;
		double sum_of_optimist_gain = 1; // one bit for free to never divide by 0
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
		for(size_t i=0; i<random_als.size();++i) {
			typename maptype::iterator fi = rand_to_overl.find(random_als.at(i));
			std::vector<const pw_alignment*> ovrlp_with_i = fi->second;
			std::set<const pw_alignment*, compare_pointer_pw_alignment> ovrlp_with_i_set;
			for(size_t j=0; j<ovrlp_with_i.size(); ++j) {
				ovrlp_with_i_set.insert(ovrlp_with_i.at(j));
			}
			size_t local_threads = 1;
			use_ias ias(data, ovrlp_with_i_set, model, 0.0, local_threads); 
			overlap_type ovrlp(data);
			ias.compute(ovrlp);

			const pw_alignment * i_al = fi->first;
			// count number of pieces made out of i_al:
			size_t iref = i_al->getreference1();
			size_t left, right;
			i_al->get_lr1(left, right);
			std::set<size_t> lpoints_in_ial;
			const std::set<pw_alignment, compare_pw_alignment> pieces = ovrlp.get_all();

			size_t al1l, al1r;
			i_al->get_lr1(al1l, al1r);
			for(std::set<pw_alignment, compare_pw_alignment>::const_iterator it = pieces.begin(); it!=pieces.end(); ++it) {
				const pw_alignment & pal = *it;
				if(pal.getreference1()==iref) {
					size_t l, r;
					pal.get_lr1(l,r);
					if( al1l <= l && r <= al1r) {
						lpoints_in_ial.insert(l);
					}
				}
				if(pal.getreference2()==iref) {
					size_t l, r;
					pal.get_lr2(l,r);
					if( al1l <= l && r <= al1r) {
						lpoints_in_ial.insert(l);
					}
				}
			}
			size_t r1 = i_al->getreference1();
			size_t r2 = i_al->getreference2();
			size_t a1 = data.accNumber(r1);
			size_t a2 = data.accNumber(r2);
			double base_cost1 = model.get_alignment_base_cost(a1, a2);
			double base_cost2 = model.get_alignment_base_cost(a2, a1);
			double base_cost = (base_cost1 + base_cost2)/2.0;

			double g1, g2;
			model.gain_function(*i_al, g1, g2);
			double gain = (g1+g2)/2.0;
			double ogpg = (lpoints_in_ial.size() - 1)*base_cost / gain;

#pragma omp critical(sums)
{


			std::cout << " random al " << i << " overlaps with " << fi->second.size() << " other alignments" << std::endl;
			std::cout << " We would cut all these alignments to " << ovrlp.size() << " pieces " << std::endl;
			std::cout << " The first alignment was cut into " << lpoints_in_ial.size() << " pieces " << std::endl;
			std::cout <<  " alignment base cost: " << base_cost << " gain " << gain << " optimist gain per gain " << ogpg << std::endl;




			sum_of_gain += gain;
			if(lpoints_in_ial.size()>1) {
				sum_of_optimist_gain += (lpoints_in_ial.size() - 1) * base_cost;
			}
}
		}
		
		std::cout << " Total gain " << sum_of_gain << " total optimist gain " << sum_of_optimist_gain << ", optimist gain per gain " << sum_of_optimist_gain / sum_of_gain << std::endl; 
		optimist_gain_per_gain = sum_of_optimist_gain / sum_of_gain;

	}


	/*
		requires: other_al overalps with reference ref of al
		output: which reference of other_al does not overlap with al
	*/
	void alignment_filter::overlap_other_side(const pw_alignment * al, size_t ref, const pw_alignment * other_al, size_t & other_ref_of_other_al) {

		size_t oar1 = other_al->getreference1();
		size_t oar2 = other_al->getreference2();
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		if(ref==0) {
			if(r1 == oar1) {
				size_t all, alr;
				size_t othl, othr;
				al->get_lr1(all, alr);
				other_al->get_lr1(othl, othr);
				if(othl<=alr && all<=othr) { // overlap r1 of al with r1 of other_al
					size_t o2l, o2r;
					other_al->get_lr2(o2l, o2r);
					if(o2r<all || o2l > alr) { // no overlap r1 of al with r2 of other_al
						other_ref_of_other_al = 1;
						return;
					}
				}

			}
			if(r1 == oar2) { 
				size_t all, alr;
				size_t othl, othr;
				al->get_lr1(all, alr);
				other_al->get_lr2(othl, othr);
				if(othl<=alr && all<=othr) { // overlap r1 of al with r2 of other_al
					size_t o2l, o2r;
					other_al->get_lr1(o2l, o2r);
					if(o2r<all || o2l > alr) { // no overlap r1 of al with r1 of other al
						other_ref_of_other_al = 0;
						return;
					}
				}

			}
		}
		if(ref==1) {
			if(r1 == oar1) {
				size_t all, alr;
				size_t othl, othr;
				al->get_lr2(all, alr);
				other_al->get_lr1(othl, othr);
				if(othl<=alr && all<=othr) { // overlap r2 of al with r1 of other_al
					size_t o2l, o2r;
					other_al->get_lr2(o2l, o2r);
					if(o2r<all || o2l > alr) { // no overlap r2 of al with r2 of other_al
						other_ref_of_other_al = 1;
						return;
					}
				}

			}
			if(r1 == oar2) {
				size_t all, alr;
				size_t othl, othr;
				al->get_lr2(all, alr);
				other_al->get_lr2(othl, othr);
				if(othl<=alr && all<=othr) { // overlap r2 of al with r2 of other_al
					size_t o2l, o2r;
					other_al->get_lr1(o2l, o2r);
					if(o2r<all || o2l > alr) { // no overlap r2 of al with r1 of other_al
						other_ref_of_other_al = 0;
						return;
					}
				}

			}
		} 
		other_ref_of_other_al = 2;
	}


	double alignment_filter::homology_score(const pw_alignment * al, const pw_alignment * al1, const pw_alignment * al2, size_t al1ref, size_t al2ref) {
		size_t a1l, a1r, a2l, a2r;
		al->get_lr1(a1l, a1r);
		al->get_lr2(a2l, a2r);
		// special case: al within one sequence al1 or al2 overlap with both references, this is no confirmation of al
		if(al->getreference1() == al->getreference2()) {
			// here we look for the part of al1/al2 on the same reference with al
			size_t al1sl, al1sr, al2sl, al2sr; 
			if(1 - al1ref) {
				al1->get_lr2(al1sl, al1sr);
			} else {
				al1->get_lr1(al1sl, al1sr);
			}
			if(1 - al2ref) {
				al2->get_lr2(al2sl, al2sr);
			} else {
				al2->get_lr1(al2sl, al2sr);
			}	
			// al1 overlaps with both references of al 
			if( ((a1l <= al1sr) && (a1r >= al1sl)) && 
			((a2l <= al1sr) && (a2r >= al1sl)) ) return 0;
			// al2 overlaps with both references of al 
			if( ((a1l <= al2sr) && (a1r >= al2sl)) && 
			((a2l <= al2sr) && (a2r >= al2sl)) ) return 0;
			// al1 overlaps with al2 on the same reference as al
			if( ((al1sl <= al2sr) && (al1sr >= al2sl))) return 0; 
		}




		size_t al1reference, al2reference;
		size_t al1l, al1r, al2l, al2r;
		if(al1ref) {
			al1->get_lr2(al1l, al1r);
			al1reference = al1->getreference2();
		} else {
			al1->get_lr1(al1l, al1r);
			al1reference = al1->getreference1();
		}
		if(al2ref) {
			al2->get_lr2(al2l, al2r);
			al2reference = al2->getreference2();
		} else {
			al2->get_lr1(al2l, al2r);
			al2reference = al2->getreference1();
		}
		assert(al1reference == al2reference);
		assert( (al1l<=al2r && al1r>=al2l)   ); // there is some overlap



		pw_alignment intersection1;
		pw_alignment intersection2;
		al1->single_ref_intersect(al1ref, *al2, al2ref, intersection1, intersection2);

		// lr on the same reference as al (opposite to the one given in al1ref/al2ref)
		size_t i1l, i1r;
		size_t i2l, i2r;
		if(al1ref) {
			al1->get_lr1(i1l, i1r);
			assert(al1->getreference1() == al->getreference1());
		} else {
			al1->get_lr2(i1l, i1r);
			assert(al1->getreference2() == al->getreference1());
		}
		if(al2ref) {
			al2->get_lr1(i2l, i2r);
			assert(al2->getreference1() == al->getreference2());
		} else {
			al2->get_lr2(i2l, i2r);
			assert(al2->getreference2() == al->getreference2());
		}

		
		
		double r1_cover = 0; // maximum length fraction of al (ref1) which is covered by intersection1
		double r2_cover = 0; // maximum length fraction of al (ref2) which is covered by intersection1


		double r1_denom = (double) (a1r - a1l + 1);
		double r2_denom = (double) (a2r - a2l + 1);

		if(a1l < i1l) {
			if(a1r < i1r) {
				r1_cover = (double)(a1r - i1l + 1)/r1_denom;
			} else {
				r1_cover = (double)(i1r - i1l + 1)/r1_denom;
			}
		} else if (a1l > i1l) {
			if(a1r < i1r) {
				r1_cover = (double)(a1r - a1l + 1)/r1_denom;
			} else {
				r1_cover = (double)(i1r - a1l + 1)/r1_denom;
			}
		}
		if(a2l < i2l) {
			if(a2r < i2r) {
				r2_cover = (double)(a2r - i2l + 1)/r2_denom;
			} else {
				r2_cover = (double)(i2r - i2l + 1)/r2_denom;
			}
		} else if (a2l > i2l) {
			if(a2r < i2r) {
				r2_cover = (double)(a2r - a2l + 1)/r2_denom;
			} else {
				r2_cover = (double)(i2r - a2l + 1)/r2_denom;
			}
		}



		

		double score = r1_cover * r2_cover;
/*
		if(score>0) {
			std::cout << " Found homology confirmation for " << std::endl;
			al->print();
			std::cout << " a1 " << a1l << " " << a1r << " a2 " << a2l << " " << a2r << std::endl;
			std::cout << " i1 " << i1l << " " << i1r << " i2 " << i2l << " " << i2r << std::endl;
			std::cout << " Could be confirmed by " << std::endl;
			al1->print();
			std::cout << "( on ref "<< al1ref<<") and " << std::endl;
			al2->print();
			std::cout << "( on ref "<< al2ref<< ") Alignment intersection: " << std::endl;
			intersection1.print();
			std::cout << std::endl;
			intersection2.print();
			std::cout << " Intersection covers ref1 of original alignment with " << r1_cover << " and ref2 with " << r2_cover << " homology score is " << score << std::endl;
	


		}
*/

		return score;

	}

	void alignment_filter::homology_gain() {
		typedef avl_interval_tree<size_t> treetype;

#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
		for(size_t i=0; i<alvec.size(); ++i) {
			std::vector<const pw_alignment *> r1res;
			std::vector<const pw_alignment *> r2res;
			const pw_alignment * al = alvec.at(i);
			size_t zero = 0;
			size_t one = 1;
			alind->search_overlap(*al, zero, r1res);
			alind->search_overlap(*al, one, r2res);
			double g1, g2;
			model.gain_function(*al, g1, g2);
			double al_gain = (g1+g2)/2.0;
			
			std::vector< std::pair<const pw_alignment *, size_t> > r1_ovrlp; // alignments which overlap with ref1 of al, second value 0/1 for the reference which does not overlap with al
			for(size_t j=0; j<r1res.size(); ++j) {
				size_t oth1;
				overlap_other_side(al, 0, r1res.at(j), oth1);
				if(oth1<2) {
					r1_ovrlp.push_back(std::make_pair(r1res.at(j), oth1));
				}
			}
			std::vector< std::pair<const pw_alignment *, size_t> > r2_ovrlp; // alignments which overlap with ref2 of al, second value 0/1 for the reference which does not overlap with al
			for(size_t j=0; j<r2res.size(); ++j) {
				size_t oth2;
				overlap_other_side(al, 1, r2res.at(j), oth2);
				if(oth2<2) {
					r2_ovrlp.push_back(std::make_pair(r2res.at(j), oth2));
				}
			}
			// Find references where we have both r1 other sides and r2 other sides
			std::set<size_t> all_references;
			std::multimap<size_t, std::pair<const pw_alignment *, size_t> > ref_to_r1;
			std::multimap<size_t, std::pair<const pw_alignment *, size_t> > ref_to_r2;
			for(size_t j=0; j<r1_ovrlp.size(); ++j) {
				size_t r = r1_ovrlp.at(j).second;
				size_t cref;
				if(r) {
					cref = r1_ovrlp.at(j).first->getreference2();	
				} else {
					cref = r1_ovrlp.at(j).first->getreference1();	
				}
				all_references.insert(cref);
				ref_to_r1.insert(std::make_pair(cref, r1_ovrlp.at(j)));
			}
			for(size_t j=0; j<r2_ovrlp.size(); ++j) {
				size_t r = r2_ovrlp.at(j).second;
				size_t cref;
				if(r) {
					cref = r2_ovrlp.at(j).first->getreference2();	
				} else {
					cref = r2_ovrlp.at(j).first->getreference1();	
				}
				all_references.insert(cref);
				ref_to_r2.insert(std::make_pair(cref, r2_ovrlp.at(j)));
			}
		
			double sum_of_scores = 0;
		//	std::cout << " ovrlps " << r1_ovrlp.size() << " " << r2_ovrlp.size() << " on " << all_references.size() << " references " << std::endl;
			for(std::set<size_t>::iterator it = all_references.begin(); it!=all_references.end(); ++it) {
				size_t cref = *it;
				std::pair< std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator, std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator> eqr1 =
				ref_to_r1.equal_range(cref);
				std::pair< std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator, std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator> eqr2 =
				ref_to_r2.equal_range(cref);

				std::vector< std::pair<const pw_alignment *, size_t> > r2_on_cref;
				for( std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator mit2 = eqr2.first; mit2!=eqr2.second; ++mit2) {
					r2_on_cref.push_back(mit2->second);
				}
				std::vector< std::pair<const pw_alignment *, size_t> > r1_on_cref;
				for( std::multimap<size_t, std::pair<const pw_alignment *, size_t> >::iterator mit1 = eqr1.first; mit1!=eqr1.second; ++mit1) {
					r1_on_cref.push_back(mit1->second);
				}
		//		std::cout << " on ref " << cref << " compare " << r2_on_cref.size() << " against " << r1_on_cref.size() << std::endl;
				treetype tree;
				for(size_t j=0; j<r2_on_cref.size(); ++j) {
					size_t l,r;
					if(r2_on_cref.at(j).second) {
						r2_on_cref.at(j).first->get_lr2(l,r);
					} else {
						r2_on_cref.at(j).first->get_lr1(l,r);
					}	
					tree.insert(l,r,j);
				}

				for(size_t j=0; j<r1_on_cref.size(); ++j) {
					std::vector<size_t> r2res;
					size_t l,r;
					if(r1_on_cref.at(j).second) {
						r1_on_cref.at(j).first->get_lr2(l,r);
					} else {
						r1_on_cref.at(j).first->get_lr1(l,r);
					}	
						
					tree.overlap_search(l, r, r2res);
	
		//			std::cout << " " << r2res.size();

					double maxscore = 0;
					for(size_t k=0; k<r2res.size(); k++) {
						const pw_alignment * r2al = r2_on_cref.at(r2res.at(k)).first;
						size_t r2alref = r2_on_cref.at(r2res.at(k)).second;
						double score = homology_score(al, r1_on_cref.at(j).first, r2al, r1_on_cref.at(j).second, r2alref);
						if(score > maxscore) maxscore = score;
					}
					sum_of_scores += maxscore;
				}

		//		std::cout << std::endl; // TODO 
			}
			
			extra_gain.at(i) = al_gain * sum_of_scores;
//			std::cout << " al " << i << " length " << al->alignment_length() << " gain " << al_gain << " extra " << sum_of_scores << std::endl;


		}



	}

/* DFS algorithm (recursive) to find the best path (highest sum of optimist gains) 
   of alignments after r1pos/r2pos with max distance gap_in_long_centers on r1 and r2
using r1forward/r2forward we can separately control in which direction we go on each reference

*/
double alignment_filter::best_path_after(const size_t & r1pos, const size_t & r2pos, bool r1forward, bool r2forward, const size_t & rangel, const size_t & ranger, const avl_interval_tree<std::pair<size_t, size_t> > & inrange, std::vector<std::pair< size_t, size_t> > & best_path) {

//	std::cout << " search strep from " << r1pos << "  " << r2pos << " dir " << r1forward << " " << r2forward << std::endl;

	std::vector<std::pair<size_t, size_t> > possible_edges;
	if(r1forward) {
		inrange.overlap_search(r1pos+1, ranger, possible_edges);
	} else {
		if(r1pos <= 1) return 0.0;
		inrange.overlap_search(rangel, r1pos - 1, possible_edges);
	}	

	double best_score = 0.0;

	for(size_t i=0; i<possible_edges.size(); ++i) {
		const pw_alignment * nal = alvec.at(possible_edges.at(i).first);
		size_t nalref = possible_edges.at(i).second;

		size_t l1, r1, l2, r2;
		if(nalref) { // reference 2 of nal is ref1
			nal->get_lr2(l1,r1);
			nal->get_lr1(l2,r2);
		} else { 
			nal->get_lr1(l1,r1);
			nal->get_lr2(l2,r2);
		}

		bool ref1ok = false;
		bool ref2ok = false;

//		std::cout << " lr1 " << l1 << " " << r1 << " lr2 " << l2 << " " << r2 << std::endl;

		if(r1forward) {
			if(l1 > r1pos && r1pos + gap_in_long_centers >= l1) ref1ok = true;
		} else {
			if(r1pos <= 1) return 0; // all alignments should be longer than 1, next alignment strictly before r1, r2
			if(r1 < r1pos && r1 + gap_in_long_centers >= r1pos) ref1ok = true;
		}
		if(r2forward) {
			if(l2 > r2pos && r2pos + gap_in_long_centers >= l2) ref2ok = true;
		} else {
			if(r2pos <= 1) return 0; // all alignments should be longer than 1, next alignment strictly before r1, r2
			if(r2 < r2pos && r2 + gap_in_long_centers >= r2pos) ref2ok = true;
		}
		if(ref1ok && ref2ok) {
			double g1, g2;
			model.gain_function(*nal, g1, g2);
			double gain = (g1+g2)/2.0 * extra_gain.at(possible_edges.at(i).first) * optimist_gain_per_gain;

			std::vector<std::pair<size_t, size_t > > local_path;
			size_t nr1pos, nr2pos, nrl, nrr;
			if(r1forward) {
				nr1pos = r1;
				nrl = r1;
				nrr = ranger;
			} else {
				nr1pos = l1;
				nrl = rangel;
				nrr = l1;
			}
			if(r2forward) { 
				nr2pos = r2;
			} else {
				nr2pos = l2;
			}

//			std::cout << " after " << r1pos << " " << r2pos << " dir " << r1forward << " " << r2forward << " we fit " << std::endl;
//			nal->print();
//			std::cout << " gain " << gain << " new pos " << nrl << " " << nrr << std::endl;

			double lscore = best_path_after(nr1pos, nr2pos, r1forward, r2forward, nrl, nrr, inrange, local_path);
			if(lscore + gain > best_score) {
				best_score = lscore + gain;
				local_path.push_back(possible_edges.at(i));
				best_path = local_path;

			}
		}

	}
	return best_score;
}


double alignment_filter::synteny_score(size_t alind, const size_t & alref, const size_t & rangel, const size_t & ranger, const avl_interval_tree<std::pair<size_t, size_t> > & inrange) {



	const pw_alignment * al = alvec.at(alind);


	std::vector<std::pair<size_t, size_t > > path_before; // TODO paths are only needed for debugging
	std::vector<std::pair<size_t, size_t > > path_after;
	
	bool is_same_dir = al->is_same_direction();
	size_t l,r;
	size_t l2, r2;
	if(alref) {
		al->get_lr2(l,r);
		al->get_lr1(l2,r2);
	} else {
		al->get_lr1(l,r);
		al->get_lr2(l2,r2);
	}
	double best_before, best_after;
	if(is_same_dir) {
// before al (on ref alref)
//		std::cout << " search SAME AFTER " << std::endl;
		best_after = best_path_after(r, r2, true, true, r, ranger, inrange, path_after);
// after al (on ref alref)
//		std::cout << " search SAME BEFORE " << std::endl;
		best_before = best_path_after(l, l2, false, false, rangel, l, inrange, path_before);
	} else {
//		std::cout << " search DIFF AFTER " << std::endl;
		best_after = best_path_after(r, l2, true, false, r, ranger, inrange, path_after);
//		std::cout << " search DIFF BEFORE " << std::endl;
		best_before = best_path_after(l, r2, false, true, rangel, l, inrange, path_before);

	}

	double g1, g2;
	model.gain_function(*al, g1, g2);
	double algain = (g1+g2)/2.0 * extra_gain.at(alind) * optimist_gain_per_gain;

	double alscore = algain + best_before + best_after;
	if(best_before > 0 || best_after > 0) {
//		std::cout << " Found synteny confirmation for "<< alref << " algain is " << algain << " total opt gain is " << alscore << std::endl;
//		al->print();

//		std::cout << " path before: " << best_before<< std::endl;
		for(size_t i=0; i<path_before.size(); ++i) {
			const pw_alignment * pal = alvec.at(path_before.at(i).first);
			double pg1, pg2;
			model.gain_function(*pal, pg1, pg2);
			double pgain = (pg2+pg2)/2.0;
			double hgain = pgain * extra_gain.at(path_before.at(i).first);
			double sgain = hgain * optimist_gain_per_gain;
//			std::cout << " ref " << path_before.at(i).second << " orig gain " << pgain << " h gain " << hgain << " s gain " << sgain << std::endl;
//			alvec.at(path_before.at(i).first)->print();

		}
//		std::cout << " path after: " <<best_after << std::endl;
		for(size_t i=0; i<path_after.size(); ++i) {
			const pw_alignment * pal = alvec.at(path_after.at(i).first);
			double pg1, pg2;
			model.gain_function(*pal, pg1, pg2);
			double pgain = (pg2+pg2)/2.0;
			double hgain = pgain * extra_gain.at(path_after.at(i).first);
			double sgain = hgain * optimist_gain_per_gain;
//			std::cout << " ref " << path_after.at(i).second << " orig gain " << pgain << " h gain " << hgain << " s gain " << sgain << std::endl;
//			alvec.at(path_after.at(i).first)->print();
		}

	}

	return alscore;
}

void alignment_filter::synteny_gain() {
	typedef avl_interval_tree<std::pair<size_t, size_t> > treetype;

	// sort first to make multi-threaded index creation easy
	std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > > 
	sorted_als(data.numSequences(), std::vector<std::vector< std::pair<size_t, size_t> > > (data.numSequences())); // ref1->ref2->j-> ( index in alvec, which reference of that al is used)
	for(size_t i=0; i<alvec.size(); ++i) {
		size_t r1 = alvec.at(i)->getreference1();
		size_t r2 = alvec.at(i)->getreference2();
		sorted_als.at(r1).at(r2).push_back(std::make_pair(i, 0));
		sorted_als.at(r2).at(r1).push_back(std::make_pair(i, 1));
	}
	std::vector<std::pair<size_t ,size_t> > comp_refs;
	for(size_t i=0; i<data.numSequences(); ++i) {
		for(size_t j=0; j<data.numSequences(); ++j) {
			comp_refs.push_back(std::make_pair(i, j));
		}
	}

	vector<double> synscores(alvec.size(), 0);
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
	for(size_t i=0; i<comp_refs.size(); ++i) {
		size_t r1 = comp_refs.at(i).first;
		size_t r2 = comp_refs.at(i).second;
/* we build an index (avl_interval_tree) for each pair of reference sequences
   index.at(i).at(j) contains all alignments between i and j (and also j and i), indexed relative to reference i	
*/		
		treetype index;
		for(size_t j=0; j<sorted_als.at(r1).at(r2).size(); ++j) {
			std::pair<size_t, size_t> ald = sorted_als.at(r1).at(r2).at(j);
			size_t l, r;
			const pw_alignment * al = alvec.at(ald.first);
			if(ald.second) {
				al->get_lr2(l,r);
			} else {
				al->get_lr1(l,r);
			}
			index.insert(l, r, ald);
		}
	
		for(size_t j=0; j<sorted_als.at(r1).at(r2).size(); ++j) {
			std::pair<size_t, size_t> ald = sorted_als.at(r1).at(r2).at(j);
			size_t l, r;
			const pw_alignment * al = alvec.at(ald.first);
			if(ald.second) {
				al->get_lr2(l,r);
			} else {
				al->get_lr1(l,r);
			}

			// valid range around al
			size_t len = r - l + 1 + gap_in_long_centers;
			size_t rangel = 0;
			size_t ranger = r + len;
			if(l > len) {
				rangel = l - len;
			}

	
			// find alignments close to al
			std::vector<std::pair<size_t, size_t> > range;
			avl_interval_tree<std::pair<size_t, size_t> > inrange;
			index.overlap_search(rangel, ranger, range);

			// check direction fits with al
			bool samedir = al->is_same_direction();
			for(size_t k=0; k<range.size(); ++k) {
				size_t ind = range.at(k).first;
				if(samedir == alvec.at(ind)->is_same_direction()) {
				if(ald.first != ind) { // dont add al itself. otherwise it could confirm itself if both parts point towards each other on the same sequence
					size_t l, r;
					if(range.at(k).first) {
						alvec.at(ind)->get_lr2(l,r);
					} else {
						alvec.at(ind)->get_lr1(l,r);
					}
					inrange.insert(l,r,range.at(k));
				}
				}
			}

			double score = synteny_score(ald.first, ald.second, rangel, ranger, inrange);
			synscores.at(ald.first)+=score;
		}	
	}

	// compute some statistics, compute result
	sum_gain = 0;
	sum_hom = 0;
	sum_syn = 0;

	double h_gain_below_2 = 1; // 1 instead of 0, to never divide by 0
	double h_gain_above_2 = 1;

	double s_gain_below_2 = 1;
	double s_gain_above_2 = 1;

	for(size_t i=0; i<alvec.size(); ++i) {
		double g1, g2;
		model.gain_function(*(alvec.at(i)), g1, g2);
		double gain = (g1+g2)/2.0;
		double hgain = gain * extra_gain.at(i);
		double sgain = synscores.at(i);
		sum_gain +=gain;

		if(hgain > 2 * gain) {
			h_gain_below_2 += 2*gain;
			h_gain_above_2 += hgain - 2*gain;
		} else {
			h_gain_below_2 += hgain;
		}


		if(sgain > 2*gain) {
			s_gain_below_2 += 2*gain;
			s_gain_above_2 += sgain - 2*gain;	
		} else {
			s_gain_below_2 += sgain;
		}
	}
	double h_extra = 0;
	if(sum_gain > h_gain_below_2) {
		h_extra = sum_gain - h_gain_below_2;
	}
	double s_extra = 0;
	if(sum_gain > s_gain_below_2) {
		s_extra = sum_gain - s_gain_below_2;
	}
	double hfact = (sum_gain + h_extra) / h_gain_above_2;
	double sfact = (sum_gain + s_extra) /s_gain_above_2;
	if(hfact > 1.0) hfact = 1.0;
	if(sfact > 1.0) sfact = 1.0;
	std::cout << " Total sum of orig gain " << sum_gain << " homology gain below factor 2: " << h_gain_below_2 << " above factor 2: " << h_gain_above_2 << " scale by " << hfact << std::endl;
	std::cout <<  "Synteny gain below factor 2: " << s_gain_below_2 << " above factor 2: " << s_gain_above_2 << " scale by " << sfact << std::endl;


	for(size_t i=0; i<alvec.size(); ++i) {
		double g1, g2;
		model.gain_function(*(alvec.at(i)), g1, g2);
		double gain = (g1+g2)/2.0;
		double hgain = gain * extra_gain.at(i);
		double sgain = synscores.at(i);

		if(hgain > 2*gain) {
			hgain = 2*gain + (hgain - 2*gain) * hfact;
		}
		if(sgain > 2*gain) {
			sgain = 2*gain + (sgain - 2*gain) * sfact;
		}

		sum_hom+=hgain;
		sum_syn+=sgain;

		extra_gain.at(i) = gain + hgain + sgain;
	}
	
	double inc_hom = sum_hom / sum_gain;
	double inc_syn = sum_syn / sum_gain;

	std::cout << " increase factor from homology " << inc_hom << " from synteny " << inc_syn << std::endl;

}


void alignment_filter::filter_als() {
	std::multimap<double, size_t> sgtoid;
	std::vector<std::vector<std::multimap<double, size_t> > > sti_on_ref(data.numAcc(), std::vector<std::multimap<double, size_t> >(data.numAcc()));
	std::vector<std::vector<double> > sum_on_acc(data.numAcc(), std::vector<double>(data.numAcc()));
	double sum = 0;
	for(size_t i=0; i<extra_gain.size(); ++i) {
		sgtoid.insert(std::make_pair(extra_gain.at(i), i));
		const pw_alignment * al = alvec.at(i);
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		size_t a1 = data.accNumber(r1);
		size_t a2 = data.accNumber(r2);
		if(a2 < a1) {
			size_t tmp = a1;
			a1 = a2;
			a2 = tmp;
		}
		sti_on_ref.at(a1).at(a2).insert(std::make_pair(extra_gain.at(i), i));
		sum+=extra_gain.at(i);
		sum_on_acc.at(a1).at(a2)+=extra_gain.at(i);
	}
	std::vector<bool> keep(alvec.size(), 0);
	double we_want = 0.9 * sum;
	double we_have = 0;
	for(std::multimap<double, size_t>::reverse_iterator rit = sgtoid.rbegin(); rit!=sgtoid.rend(); ++rit) {
		double g = rit->first;
		size_t id = rit->second;
		we_have+=g;
		keep.at(id) = 1;
		if(we_have>=we_want) {
					const pw_alignment * al = alvec.at(id);
					double g1, g2;
					model.gain_function(*al, g1,  g2);
					double gain = (g1+g2)/2.0;
					std::cout << " last al taken in general list has " << g << " optimist gain " << gain << " original gain "<< al->alignment_length() << " length " << std::endl;
			break;
		}
	} 

	for(size_t a1=0; a1<data.numAcc(); ++a1) {
		for(size_t a2=a1; a2<data.numAcc(); ++a2) {
			double we_want = 0.75 * sum_on_acc.at(a1).at(a2);
			double we_have = 0;
			for(std::multimap<double, size_t>::reverse_iterator rit = sti_on_ref.at(a1).at(a2).rbegin(); rit!=sti_on_ref.at(a1).at(a2).rend(); ++rit) {
				double g = rit->first;
				size_t id = rit->second;
				we_have+=g;
				keep.at(id) = 1;
				if(we_have >= we_want) {
					const pw_alignment * al = alvec.at(id);
					double g1, g2;
					model.gain_function(*al, g1,  g2);
					double gain = (g1+g2)/2.0;
					std::cout << " last al taken between "<< a1 <<" and "<< a2 << " has " << g << " optimist gain " << gain << " original gain "<< al->alignment_length() << " length " << std::endl;



					break;
				}
			}
		}
	}

	double missing_orig_gain = 0;
	std::multimap<double, size_t> missing_map;
	for(size_t i=0; i<alvec.size(); ++i) {
		if(!keep.at(i)) {
			const pw_alignment * al = alvec.at(i);
			double g1, g2;
			model.gain_function(*al, g1,  g2);
			double gain = (g1+g2)/2.0;
			missing_orig_gain += gain;
			missing_map.insert(std::make_pair(gain,  i));
		}
	}
	
	we_want = 0.5 * missing_orig_gain;
	we_have = 0;
	for(std::multimap<double, size_t>::reverse_iterator rit = missing_map.rbegin(); rit!=missing_map.rend(); ++rit) {
		double g = rit->first;
		size_t id = rit->second;
		we_have+=g;
		keep.at(id) = 1;
		if(we_have>=we_want) {
			const pw_alignment * al = alvec.at(id);
			double g1, g2;
			model.gain_function(*al, g1,  g2);
			double gain = (g1+g2)/2.0;
			std::cout << " last al taken in list with original gain has " << g << " optimist gain " << gain << " original gain "<< al->alignment_length() << " length " << std::endl;
			break;
		}
	

	}


	double gain_kept = 0;
	double optgain_kept = 0;
	size_t shortest = std::numeric_limits<size_t>::max();
	double lowest_gain = std::numeric_limits<double>::max();
	size_t longest_erased = 0;
	double best_erased = 0;
	for(size_t i=0; i<alvec.size(); ++i) {
		double g1, g2;
		model.gain_function(*(alvec.at(i)), g1, g2);
		double gain = (g1+g2)/2.0;
		size_t len = alvec.at(i)->alignment_length();
		if(keep.at(i)) {
			gain_kept+=gain;
			optgain_kept+=extra_gain.at(i);
			if(shortest>len) shortest = len;
			if(lowest_gain > gain) lowest_gain = gain;
		} else {
			alignments.erase(alvec.at(i));
			if(gain > best_erased) best_erased = gain;
			if(longest_erased < len) longest_erased = len;
		}
	}
	std::cout << " We have kept " << alignments.size() << " alignments of " << alvec.size() << " original alignments " << std::endl;
	std::cout << " They represent " << gain_kept/sum_gain << " of all gain and " << optgain_kept/sum << " of all optimist gain " << std::endl;
	std::cout << " Shortest kept alignment has length " << shortest << " lowest gain kept is " << lowest_gain << std::endl;
	std::cout << " Longest erased alignment has length " << longest_erased << " highest gain erased is " << best_erased << std::endl;
}

// #include "intervals.cpp"

#endif


