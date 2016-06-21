#include "connected_component.hpp"


	void compute_cc_avl_fraction::compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & ccs){
		std::cout << "compute CC on " << alignments.size() << std::endl;
		std::set <const pw_alignment*, compare_pointer_pw_alignment> seen;
		std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> > sorter; // sorts ccs to give largest first
		for(std::set<const pw_alignment *, compare_pointer_pw_alignment>::iterator it = alignments.begin(); it!=alignments.end(); ++it) {
			const pw_alignment * al = *it;
			al->print();
			std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
			std::cout << " seen size " << seen.size() << std::endl;
			if(seenal == seen.end()) {
				std::cout << " getcc" << std::endl;
				std::set< const pw_alignment*, compare_pointer_pw_alignment> cc;
			//	alind->erase(al);
				seen.insert(al);
				seen1.insert(al);

				cc.insert(al);
				get_cc(*al, cc, seen);
				std::cout << "cc size "<< cc.size()<<std::endl;
				sorter.insert(std::make_pair(cc.size(), cc));
				std::cout << "edges size "<<edges.size()<<std::endl;
				for(std::multimap<const pw_alignment* , const pw_alignment*>::iterator it = edges.begin(); it!=edges.end(); it++){
					std::cout << it->first << " " << it->second << std::endl;
				}
			}
		}
		for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
			if(it->second.size()>1){
				ccs.push_back(it->second);
			}
		}
	//	std::cout << "edges size "<<edges.size()<<std::endl;
	//	for(std::multimap<const pw_alignment* , const pw_alignment*>::iterator it = edges.begin(); it!=edges.end(); it++){
	//		std::cout << it->first << " " << it->second << std::endl;
	//	}
	}
	void compute_cc_avl_fraction::get_cc(const pw_alignment & al , std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment*, compare_pointer_pw_alignment> & seen) {
		std::vector<size_t> left(2);
		std::vector<size_t> right(2);
		al.get_lr1(left.at(0), right.at(0));
		al.get_lr2(left.at(1), right.at(1));
		std::vector<size_t>reference(2);
		reference.at(0) = al.getreference1();
		reference.at(1) = al.getreference2();
#pragma omp parallel for num_threads(num_threads)
		for(size_t i =0; i < 2;i++){
			std::cout << "on ref " << i<<std::endl;
			cc_step(al,reference.at(i), left.at(i), right.at(i), cc, seen);	
		}	
	}

	void compute_cc_avl_fraction::cc_step(const pw_alignment & al, size_t current , size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & cc, std::set <const pw_alignment* , compare_pointer_pw_alignment>  & seen ) {
		std::vector<const pw_alignment*> results;
//#pragma omp critical(seen)
//{
		std::cout << "ref " << current << " l "<< left << " r "<< right << "fraction "<< FRACTION <<std::endl;
//		alind->search_overlap_fraction(current, left, right, FRACTION , results);
		alind->search_overlap(current, left, right, results);

	//	std::cout << "reult size here is "<<results.size()<<std::endl;
		//CHECK FOR 95% OVERLAP AND DELETE THE ALIGNMENT ITSELF FROM RESULTS!
		std::vector<const pw_alignment*> fraction_result;

		compute_fraction(results, FRACTION,current, left, right, fraction_result);
		for(size_t i=0; i<fraction_result.size(); ++i) {
			seen.insert(fraction_result.at(i)); // TODO for now we keep seen, but it should not be necessary as all alignments are deleted from the index
			cc.insert(fraction_result.at(i));
			if(&al != fraction_result.at(i)){
			//	TODO change the container in a way that add it only once
				edges.insert(std::make_pair(&al,fraction_result.at(i)));
			}
		//	alind->erase(results.at(i));//XXX ????? It shouldnt be deleted, because it might itself have overlap with one of the other alignments in results. 
		}
//}
		for(size_t i=0; i<fraction_result.size(); ++i) {
			const pw_alignment* p = fraction_result.at(i);
			std::cout << "al " << &al << " p " << p << std::endl;
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it =seen1.find(p);
			if(it == seen1.end()){
				seen1.insert(p);
		//	if(&al != p){
				size_t ref1 = p->getreference1();
				size_t ref2 = p->getreference2();
				size_t l1,l2,r1,r2;
				p->get_lr1(l1,r1);
				p->get_lr2(l2,r2);
				std::cout << "l1 " << l1 << " r1 " << r1 <<std::endl;
				cc_step(*p,ref1, l1, r1, cc, seen);
				cc_step(*p,ref2, l2, r2, cc, seen);
			}
		//	alind->erase(fraction_result.at(i)); //TODO Then where is the right place for it??
		}
		for(size_t i=0; i<fraction_result.size(); ++i) {
			alind->erase(fraction_result.at(i));
		}
	}
/*	size_t detect_overlap_rate::get_id(const pw_alignment* p)const{
		std::map<const pw_alignment* , size_t>::const_iterator it = als_index.find(p);
		assert(it != als_index.end());
		return it->second;
	}*/
	void compute_cc_avl_fraction::compute_fraction(std::vector<const pw_alignment*> & results, const double  fraction,  const size_t & current , const size_t & left, const size_t & right, std::vector<const pw_alignment*> & fraction_result ){
		for(size_t i =0; i < results.size();i++){
			const pw_alignment * p = results.at(i);
			size_t ref1 = p->getreference1();
			size_t ref2 = p->getreference2();
			size_t l1,l2,r1,r2;
                        p->get_lr1(l1,r1);
                        p->get_lr2(l2,r2);
	 		if( ref1 == current && l1<=right && r1 >= left){
				if(l1<= left){
					double overlap = (r1-left+1)/double(right-l1+1);
					 if(overlap>fraction){
                                      	  fraction_result.push_back(p);
					}
				}else{	
	                                double overlap = (right-l1+1)/double(r1-left+1);
                             		if(overlap>fraction){
                                        	fraction_result.push_back(p);
					}
				}
			}else if(ref2 == current && l2<=right && r2>= left){
				if(l2<= left){
					double overlap = (r2-left+1)/double(right-l2+1);
					 if(overlap>fraction){
                                      	  fraction_result.push_back(p);
					}
				}else{	
	                                double overlap = (right-l2+1)/double(r2-left+1);
                             		if(overlap>fraction){
                                        	fraction_result.push_back(p);
					}
				}
			}
		}
		std::cout << "result size after checking the fraction is "<< results.size() << std::endl;
	}
	void compute_cc_avl_fraction::node_adjucents(const pw_alignment* p ,std::vector<const pw_alignment*> & adj)const{
		std::pair<std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator ,std::multimap<const pw_alignment*,const pw_alignment*>::const_iterator > it = edges.equal_range(p);
		for(std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator it1 = it.first; it1!=it.second; ++it1) {
			adj.push_back(it1->second);
		}
	}
	void biconnected_component::creat_biconnected_component(size_t & at){
		std::cout << "al size at " << at << " is "<< alignments.at(at).size()<<std::endl;
		for(size_t i =0; i < alignments.at(at).size();i++){
			std::cout << "i "<< i << std::endl;
			visited.at(at).at(i) = false;
			std::cout<< "at "<< alignments.at(at).at(i)<<std::endl;
			get_articulation_points(at,i);
		}
		std::cout <<"done! "<<std::endl;
		for(size_t i = 0; i <stack.size();i++){ 
			std::cout<< "stack at " << i <<std::endl;
			for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it= stack.at(i).begin(); it !=stack.at(i).end() ;it++){
				std::cout << *it << std::endl;
			}
		}
		depth = 0;
	}
	void biconnected_component::get_articulation_points(size_t & at, size_t node){
		visited.at(at).at(node) = true;
		Depth.at(at).at(node) = depth;
		low.at(at).at(node) = depth;
		size_t childCount = 0;
		bool isArticulation = false;
		std::vector<const pw_alignment*> adj;
		cc_fraction.node_adjucents(alignments.at(at).at(node), adj);
		std::cout<< alignments.at(at).at(node) <<std::endl;
		for(size_t i =0; i < adj.size();i++){
			std::cout << adj.at(i)<<std::endl;
		}
		std::cout << "adj size "<< adj.size()<<std::endl;
for(std::map<const pw_alignment*,size_t>::iterator it = als_index.at(at).begin() ;it != als_index.at(at).end();it++){
			std::cout<< it->first <<std::endl;

}
		for (size_t ni = 0; ni < adj.size(); ni++){//each ni in adj[node]
			std::map<const pw_alignment*,size_t>::iterator it = als_index.at(at).find(adj.at(ni));
			assert(it != als_index.at(at).end());
			std::cout << it->first << std::endl;
			std::cout << "ni index" << it->second << std::endl;
			std::cout<< "vi size at "<< at << " is "<< visited.at(at).size() <<std::endl;
		        if(visited.at(at).at(it->second) == false){
				temp.insert(alignments.at(at).at(it->second));
				temp.insert(alignments.at(at).at(node));
				std::cout << "inserted " << node << " "<<it->second<<std::endl;
				parent.at(at).at(it->second) = node;
				depth++;
				get_articulation_points(at,it->second);
				childCount = childCount + 1;
				std:: cout << low.at(at).at(it->second) << " "<< Depth.at(at).at(node)<<std::endl;
				if(low.at(at).at(it->second) >= Depth.at(at).at(node)){
					std:: cout << low.at(at).at(it->second) << " "<< Depth.at(at).at(node)<<std::endl;
					if(temp.size() > 2){
						stack.push_back(temp);
					}
					std::cout << "temp size " << temp.size()<<std::endl;
					temp.clear();
          				isArticulation = true;
				        low.at(at).at(node) = std::min(low.at(at).at(node), low.at(at).at(it->second));
				}
	   		}else if( it->second != parent.at(at).at(node)){
					temp.insert(alignments.at(at).at(it->second));
					temp.insert(alignments.at(at).at(node));
					std::cout << "inserted " << node << " "<<it->second<<std::endl;
					std::cout << "here! "<<std::endl;
	      				low.at(at).at(node) = std::min(low.at(at).at(node), Depth.at(at).at(it->second));
			}
		}
    		if( (parent.at(at).at(node) != NULL &&  isArticulation) || (parent.at(at).at(node) == NULL && childCount > 1) ){
			std::cout << node <<" is an articulation point "<< std::endl;
		}
	//	if(temp.size() > 2){
	//		stack.push_back(temp);
	//	}
	//	std::cout << "temp size " << temp.size()<<std::endl;
	//	temp.clear();
		std::cout << "here ! "<<std::endl;
	}
