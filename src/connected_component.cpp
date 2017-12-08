#ifndef CONNECTED_COMPONENT_CPP
#define CONNECTED_COMPONENT_CPP

#include "connected_component.hpp"

	void compute_cc_avl_fraction::non_recursive_compute(std::set< const pw_alignment* , compare_pointer_pw_alignment> & ccs, std::set<const pw_alignment *, compare_pointer_pw_alignment> & remainder, double & fraction){//It creats a graph and keep all the edges in a multimap. Graph might have more than one component
	//	for(size_t i = 0; i < als.size();i++) {
	//		const pw_alignment * al = als.at(i);
		edges.clear();
		std::cout << "edges size1 "<<edges.size()<<std::endl;
	//	size_t count = 0;
//#pragma omp parallel for num_threads(num_threads)
		for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainder.begin(); it!= remainder.end();it++){//TODO time consuming! Need to be parallelized , can simply change the set to a vector (NO you cannot, then you will find the same component many times)
			const pw_alignment * al = *it;
		//	al->print();
		//	std::cout << al << std::endl;
			std::vector<size_t> left(2);
			std::vector<size_t> right(2);
			al->get_lr1(left.at(0), right.at(0));
			al->get_lr2(left.at(1), right.at(1));
			std::vector<size_t>reference(2);
			reference.at(0) = al->getreference1();
			reference.at(1) = al->getreference2();
			for(size_t j =0; j < 2; j++){
			//	std::cout<< al << " on ref "<< j<<std::endl;
				cc_step_non_recursive(al, reference.at(j), left.at(j), right.at(j),ccs,fraction, remainder);	
			}
		}
		
	//	std::cout << "edges size "<<edges.size()<<std::endl;
	//	for(std::map<const pw_alignment* , std::set<const pw_alignment*> >::iterator it = edges.begin(); it!=edges.end(); it++){
	//		std::cout << it->first << " to " ;
	//		for(std::set<const pw_alignment*>::iterator adj = it->second.begin() ; adj != it->second.end(); adj ++){
	//			std::cout << *adj << std::endl;
	//		}
	//	}
		std::cout<< "ccs size "<<ccs.size()<<std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.begin(); it != ccs.end();it ++){
			const pw_alignment * p = *it;
		//	p->print();
	//		std::cout << p << std::endl;
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = remained_als.find(p);
			assert(it2 != remained_als.end());
			remained_als.erase(it2);
		}
		std::cout << "remainded size "<< remained_als.size() <<std::endl;
	}
/*void compute_cc_avl_fraction::non_recursive_compute(std::set< const pw_alignment* , compare_pointer_pw_alignment> & ccs, std::set< pw_alignment , compare_pw_alignment> & remainder, double & fraction){//It creats a graph and keep all the edges in a multimap. Graph might have more than one component
//#pragma omp parallel for num_threads(num_threads)
	//	for(size_t i = 0; i < als.size();i++) {
	//		const pw_alignment * al = als.at(i);
		edges.clear();
		std::cout << "edges size1 "<<edges.size()<<std::endl;
		for(std::set< pw_alignment , compare_pw_alignment>::iterator it = remainder.begin(); it!= remainder.end();it++){
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = remained_als.find(&*it);
			assert(it2 != remained_als.end());

			const pw_alignment * al = &*it;
			std::vector<size_t> left(2);
			std::vector<size_t> right(2);
			al->get_lr1(left.at(0), right.at(0));
			al->get_lr2(left.at(1), right.at(1));
			std::vector<size_t>reference(2);
			reference.at(0) = al->getreference1();
			reference.at(1) = al->getreference2();
			for(size_t j =0; j < 2; j++){
			//	std::cout<< al << " on ref "<< j<<std::endl;
				cc_step_non_recursive(al, reference.at(j), left.at(j), right.at(j),ccs,fraction, remainder);	
			}	
		}
		
		std::cout << "edges size "<<edges.size()<<std::endl;
	//	for(std::multimap<const pw_alignment* , const pw_alignment*>::iterator it = edges.begin(); it!=edges.end(); it++){
	//		std::cout << it->first << " " << it->second << std::endl;
	//	}
		std::cout<< "ccs size "<<ccs.size()<<std::endl;
		std::cout << "remained_als size "<< remained_als.size()<<std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.begin(); it != ccs.end();it ++){
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = remained_als.find(*it);
			assert(it2 != remained_als.end());
			remained_als.erase(*it2);
			remainders_als.erase(**it2);
		}
		std::cout << "remainded size "<< remained_als.size() <<std::endl;
	}*/

	void compute_cc_avl_fraction::cc_step_non_recursive(const pw_alignment  * al,  size_t & ref, size_t & left , size_t & right, std::set< const pw_alignment* , compare_pointer_pw_alignment> & ccs, double & fraction, std::set<const pw_alignment* , compare_pointer_pw_alignment> & remainder){
			std::vector<const pw_alignment*> results;
			std::vector<const pw_alignment*> fraction_result;
		//	std::cout << "left " << left << " right "<< right << std::endl;
			alind->search_overlap_fraction(ref, left, right, fraction , fraction_result);
		//	alind->search_overlap(ref, left, right, results);
		//	compute_fraction(results, fraction,ref, left, right, fraction_result);
#pragma omp parallel for num_threads(num_threads)
			for(size_t i=0; i<fraction_result.size(); ++i) {
				std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainder.find(fraction_result.at(i));
				if(al != fraction_result.at(i) && it != remainder.end()){
#pragma omp critical(edge)
{
				//	std::cout<< "HERE!"<<std::endl;
				//	al->print();
				//	std::cout << "al " << al << " fraction at " << i << " : "<< fraction_result.at(i) <<std::endl;
					std::map<const pw_alignment* , std::set<const pw_alignment*> >::iterator it1 = edges.find(al);
					if(it1 == edges.end()){
						edges.insert(std::make_pair(al,std::set<const pw_alignment*>()));
						it1 = edges.find(al);
					}
					it1->second.insert(fraction_result.at(i));

				//	edges.insert(std::make_pair(al,fraction_result.at(i)));
				//	edges.insert(std::make_pair(fraction_result.at(i),al));
					ccs.insert(al);
				//	ccd.insert(fraction_result.at(i));
}
				}
			}

	}
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
				alind->erase(al);
				seen.insert(al);
				seen1.insert(al);

				cc.insert(al);
				get_cc(*al, cc, seen);
				std::cout << "cc size "<< cc.size()<<std::endl;
				sorter.insert(std::make_pair(cc.size(), cc));
				for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = cc.begin(); it!=cc.end(); it++){
					std::cout << *it << std::endl;
				}
				std::cout << "edges size "<<edges.size()<<std::endl;
			//	for(std::multimap<const pw_alignment* , const pw_alignment*>::iterator it = edges.begin(); it!=edges.end(); it++){
			//		std::cout << it->first << " " << it->second << std::endl;
			//	}
			//	for(std::set<std::pair<const pw_alignment*, const pw_alignment*> >::iterator it = edges.begin();it!=edges.end();it++){
			//		std::pair<const pw_alignment*, const pw_alignment*> pair = *it;
			//		std::cout<< pair.first << " "<<pair.second<<std::endl;
			//	}
			}//TODO walk over seen1 delete them from the tree also clean the seen1 itself.
			for(std::set< const pw_alignment*, compare_pointer_pw_alignment>::iterator it1 = seen1.begin();it1 != seen1.end();it1++){
				alind->erase(*it1);
			}
			seen1.clear();
		}
		for(std::multimap<size_t , std::set<const pw_alignment*, compare_pointer_pw_alignment> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); it++) {
			if(it->second.size()>1){
				ccs.push_back(it->second);
				for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it1 = it->second.begin(); it1 != it->second.end();it1++){
					std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = remained_als.find(*it1);
					assert(it2 != remained_als.end());
					remained_als.erase(it2);
				//	remainders_als.erase(**it2);
				}
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
//		std::cout << "ref " << current << " l "<< left << " r "<< right << "fraction "<< FRACTION <<std::endl;
		std::vector<const pw_alignment*> fraction_result;

//		alind->search_overlap_fraction(current, left, right, FRACTION , fraction_result);
		alind->search_overlap(current, left, right, results);

	//	std::cout << "reult size here is "<<results.size()<<std::endl;
		//CHECK FOR 95% OVERLAP AND DELETE THE ALIGNMENT ITSELF FROM RESULTS!
	//	compute_fraction(results, FRACTION,current, left, right, fraction_result);
		for(size_t i=0; i<fraction_result.size(); ++i) {
			seen.insert(fraction_result.at(i));
			if(&al != fraction_result.at(i)){
			//	TODO change the container in a way that add it only once
			//	std::cout << "al "<< &al << " fraction "<< fraction_result.at(i)<<std::endl;
			/*	std::cout<< edges.size()<<std::endl; //TODO !!!!
				edges.insert(std::make_pair(&al,fraction_result.at(i)));
				std::cout<< edges.size()<<std::endl;
				edges.insert(std::make_pair(fraction_result.at(i),&al));
				std::cout<< edges.size()<<std::endl;*/
			}
		//	alind->erase(results.at(i));//XXX ????? It shouldnt be deleted, because it might itself have overlap with one of the other alignments in results. 
		}
//}
		for(size_t i=0; i<fraction_result.size(); ++i) {
			const pw_alignment* p = fraction_result.at(i);
			std::cout << "al " << &al << " p " << p << std::endl;
			if(&al == p){std::cout << "al is equal to p"<<std::endl;}
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it =seen1.find(p);
			if(it == seen1.end()){
				std::cout<< " p " << p << std::endl;
				cc.insert(fraction_result.at(i));
			//	if(&al != fraction_result.at(i)){
			//		edges.insert(std::make_pair(&al,fraction_result.at(i)));
			//		edges.insert(std::make_pair(fraction_result.at(i),&al));
			//	}

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
	//		alind->erase(fraction_result.at(i)); //TODO Then where is the right place for it??
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
			std::cout << " p " << p << std::endl;
			size_t ref1 = p->getreference1();
			size_t ref2 = p->getreference2();
			size_t l1,l2,r1,r2;
                        p->get_lr1(l1,r1);
                        p->get_lr2(l2,r2);
	 		if( ref1 == current && l1<=right && r1 >= left){
				double overlap = 0.0;
				if(l1<= left){
					if(right > r1){
						overlap = (r1-left+1)/double(right-l1+1);
						std::cout << overlap << std::endl;
					}else{
						overlap = (right-left+1)/double(r1-l1+1);
						std::cout << overlap << std::endl;
					}
					if(overlap>fraction){
						fraction_result.push_back(p);
					}
				}else{//l1 > left	
					if(right<=r1){
	                                	overlap = (right-l1+1)/double(r1-left+1);
					}else{
						overlap= (r1-l1+1)/double(right-left+1);
					}
                             		if(overlap>fraction){
                                        	fraction_result.push_back(p);
					}
				}
			}else if(ref2 == current && l2<=right && r2>= left){
				double overlap = 0.0;
				if(l2<= left){
					if(right > r2){
						overlap = (r2-left+1)/double(right-l2+1);
					}else{
						overlap = (right-left+1)/double(r2-l2+1);
					}
					if(overlap>fraction){
                                      		fraction_result.push_back(p);
					}
				}else{	
					if(right<=r2){
		                                overlap = (right-l2+1)/double(r2-left+1);
					}else{
						overlap= (r2-l2+1)/double(right-left+1);
					}
                             		if(overlap>fraction){
                                        	fraction_result.push_back(p);
					}
				}
			}
		}
		std::cout << "result size after checking the fraction is "<< results.size() << std::endl;
	}
	void compute_cc_avl_fraction::node_adjacents(const pw_alignment* p ,std::set<const pw_alignment*> & adj)const{
	//	std::pair<std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator ,std::multimap<const pw_alignment*,const pw_alignment*>::const_iterator > it = edges.equal_range(p);
	//	for(std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator it1 = it.first; it1!=it.second; ++it1) {
	//		adj.insert(it1->second);
	//	}
		std::map<const pw_alignment* , std::set<const pw_alignment*> >::const_iterator it1 = edges.find(p);
		adj = it1->second;		

	}
	void compute_cc_avl_fraction::node_adjacents(const pw_alignment* p ,std::vector<const pw_alignment*> & adj)const{
//		std::pair<std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator ,std::multimap<const pw_alignment*,const pw_alignment*>::const_iterator > it = edges.equal_range(p);
//		for(std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator it1 = it.first; it1!=it.second; ++it1) {
//			adj.push_back(it1->second);
//		}
		std::map<const pw_alignment* , std::set<const pw_alignment*> >::const_iterator it = edges.find(p);
		for(std::set<const pw_alignment *>::iterator it1 = it->second.begin() ; it1!= it->second.end() ; it1++){
			adj.push_back(*it1);
		}		

	}

	void compute_cc_avl_fraction::test_redundancy(const pw_alignment* p, std::set< const pw_alignment* , compare_pointer_pw_alignment> & cc){
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = cc.begin();it != cc.end();it++){
			const pw_alignment * al = *it;
			if(al == p){
				std::cout << "WARNING! THERE IS A REDUNDANT ALIGNMENT."<<std::endl;
				std::cout << al << " " << p <<std::endl;
				exit(1);
			}
		
		}

	}
	std::set<const pw_alignment*, compare_pointer_pw_alignment> compute_cc_avl_fraction::get_remaining_als()const{	
		return remained_als;
	}
//	const std::set<pw_alignment, compare_pw_alignment> compute_cc_avl_fraction::get_remainders()const{	
//		return remainders_als;
//	}

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
		std::set<const pw_alignment*> adj;
		cc_fraction.node_adjacents(alignments.at(at).at(node), adj);
		std::cout<< alignments.at(at).at(node) <<std::endl;
		for(std::set<const pw_alignment*>::iterator it = adj.begin(); it!= adj.end();it++){
			std::cout << *it <<std::endl;
		}
		std::cout << "adj size "<< adj.size()<<std::endl;
		for(std::map<const pw_alignment*,size_t>::iterator it = als_index.at(at).begin() ;it != als_index.at(at).end();it++){
			std::cout<< it->first <<std::endl;
		}
		for (std::set<const pw_alignment*>::iterator it1 = adj.begin(); it1!= adj.end();it1++){//each ni in adj[node]
			std::map<const pw_alignment*,size_t>::iterator it = als_index.at(at).find(*it1);
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


	void two_edge_cc::find_bridges(std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als){
		for(size_t i =0; i < alignments.size();i++){//TODO add bool , to check if there is no bridge and add temp to stack
		//	if(bridge==false && temp.size()!=0){
			if(temp.size()!=0){
				stack.push_back(temp);
				temp.clear();
			}
			std::map<const pw_alignment*, int>::iterator it = visited_als.find(alignments.at(i));
			assert(it!= visited_als.end());
			if(it->second == -1){
				temp.clear();
				temp.insert(alignments.at(i));
				deep_first_search(alignments.at(i),alignments.at(i), mixed_als);  
			//	non_recursive_deep_first_search(alignments.at(i));  //TODO 1. only find bridges. 2. fill in a map with all the edges in cc class. 3. Then remove these edges
 
			}
		}
		std::cout<< "stack is: "<<std::endl;
		for(size_t i =0; i < stack.size();i++){
			if(stack.at(i).size()>2){
				stacks.push_back(stack.at(i));
			}else{
				for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator it = stack.at(i).begin(); it != stack.at(i).end(); it++){
			//		mixed_als.insert(*it);
				}
			}
		//	std::cout<< "at "<< i << std::endl;
		//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator it = stack.at(i).begin(); it != stack.at(i).end(); it++){
		//		std::cout << *it<<std::endl;
		//	}
		}
	//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator it = mixed_als.begin() ; it!= mixed_als.end();it++){
	//		std::cout << *it << std::endl;
	//	}
	}
	void two_edge_cc::deep_first_search(const pw_alignment* parent, const pw_alignment* child, std::set<pw_alignment, compare_pw_alignment> & mixed_als){
		std::map<const pw_alignment*, int>::iterator vi = visited_als.find(parent);
	//	std::cout << "parent " << parent << std::endl;
		assert(vi!= visited_als.end());
		vi->second = count;
		count++;
		std::map<const pw_alignment*, int>::iterator low = als_low.find(parent);
		assert(low!= als_low.end());
		low->second = vi->second;
		std::vector<const pw_alignment*> adj;
		cc_fraction.node_adjacents(parent, adj);
	//	std::cout<<"parent "<< parent << " adj size "<< adj.size()<<std::endl;
		for(size_t i =0; i < adj.size(); i++){
		//	std::cout<< "adj" <<adj.at(i)<<std::endl;
			std::map<const pw_alignment*, int>::iterator vis = visited_als.find(adj.at(i));
			assert(vis!= visited_als.end());
			if(adj.at(i) == child ) continue;
			if(vis->second == -1){
				deep_first_search(adj.at(i),parent, mixed_als);
				std::map<const pw_alignment*, int>::iterator it = als_low.find(parent);
				assert(it != als_low.end());
				std::map<const pw_alignment*, int>::iterator it1 = als_low.find(adj.at(i));
				assert(it1 != als_low.end());
				it->second = std::min(it->second, it1->second);
				std::map<const pw_alignment*, int>::iterator it2 = visited_als.find(adj.at(i));
				assert(it2 != visited_als.end());
				if(it2->second == it1->second){
					
				//	temp.insert(parent);
				//	stack.push_back(temp);
				//	std::cout << "temp size is "<< temp.size()<<std::endl;
				//	for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator t = temp.begin(); t != temp.end() ;t++){
				//		std::cout << *t <<std::endl;
				//	}
				//	temp.clear();
					bridge = true;
				//	parent->print();
				//	std::cout << " to "<<std::endl;
				//	adj.at(i)->print();
				//	std::cout << " is a bridge."<<std::endl;
				//	std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator t = temp.find(parent);
				//	if(t == temp.end()){
				//		mixed_als.insert(parent);
				//	}
				//	t = temp.find(adj.at(i));
				//	if(t == temp.end()){
				//		mixed_als.insert(adj.at(i));
				//	}
				}else{

				//	std::cout<<"here!"<<std::endl;
					temp.insert(adj.at(i));
				//	std::cout<< "temp size "<<temp.size()<<std::endl;
					bridge = false;

				}
			}else{
				assert(adj.at(i)!= child);
				std::map<const pw_alignment*, int>::iterator it = als_low.find(parent);
				assert(it != als_low.end());
				std::map<const pw_alignment*, int>::iterator it1 = als_low.find(adj.at(i));
				assert(it1 != als_low.end());
				it->second = std::min(it->second, it1->second);
				temp.insert(adj.at(i));
			//	std::cout<<"here!"<<std::endl;
			//	std::cout<< "temp size1 "<<temp.size()<<std::endl;
				bridge = false;

			}			

		}

	}
void two_edge_cc::make_two_edge_graph(std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als){
	size_t count = 0;
	size_t count1 = 0;
	for(size_t i =0; i < alignments.size();i++){
		std::map<const pw_alignment*, int>::iterator it = visited_als.find(alignments.at(i));
		assert(it!= visited_als.end());
		if(it->second == -1){
			count ++;
			find_bridges(alignments.at(i));  
		//	std::cout<< count << std::endl;
		}else{
			count1 ++;
//			std::cout << "seen! "<< count1 <<std::endl;
		}
	}
	std::cout<<"bridges were found! "<<std::endl;
	remove_bridges(stacks,mixed_als);
}
void two_edge_cc::find_bridges(const pw_alignment * parent){
	std::vector<const pw_alignment*> temp_stack;
	temp_stack.push_back(parent);
	std::map<const pw_alignment*, const pw_alignment*> edges;

	while(!temp_stack.empty()){
		const pw_alignment * current_node = temp_stack.back();
		std::map<const pw_alignment*, int>::iterator vi = visited_als.find(current_node);
	//	std::cout << "current node " << current_node << " count "<< count << std::endl;
		assert(vi!= visited_als.end());
		if(vi->second == -1){//Prior to visit the parent node (This condition is correct!)
			vi->second = count;
			count++;
			std::map<const pw_alignment *, int>::iterator low = als_low.find(current_node);
			assert(low!= als_low.end());
			low->second = vi->second;
			std::set<const pw_alignment *> adj;
			cc_fraction.node_adjacents(current_node, adj);
			for(std::set<const pw_alignment *>::iterator this_adj = adj.begin(); this_adj != adj.end(); this_adj++){
				std::map<const pw_alignment *, int>::iterator vis = visited_als.find(*this_adj);
				assert(vis!= visited_als.end());
				if(vis->second == -1){
					temp_stack.push_back(*this_adj);
					edges.insert(std::make_pair(*this_adj,current_node));
				}else{
				//	std::cout << "update "<< std::endl;
					std::map<const pw_alignment *, int>::iterator it = als_low.find(current_node);
					assert(it != als_low.end());
					std::map<const pw_alignment *, int>::iterator it1 = visited_als.find(*this_adj);
					assert(it1 != visited_als.end());
					it->second = std::min(it->second, it1->second);
				}
			}
		}else{ //After visitng the parent node
			temp_stack.pop_back();
			std::set<const pw_alignment *> adj;
			cc_fraction.node_adjacents(current_node, adj);

			for(std::set<const pw_alignment *>::iterator this_adj = adj.begin(); this_adj != adj.end(); this_adj++){
				std::map<const pw_alignment*, const pw_alignment*>::iterator it2 = edges.find(*this_adj);
				if(it2 != edges.end() && it2->second == current_node){
				//	std::cout<<"here! ! !"<<std::endl;
					std::map<const pw_alignment *, int>::iterator it = als_low.find(current_node);
					assert(it != als_low.end());
					std::map<const pw_alignment *, int>::iterator it1 = als_low.find(*this_adj);
					assert(it1 != als_low.end());
					it->second = std::min(it->second, it1->second);
				//	current_node->print();
				//	std::cout<< "to "<<std::endl;
					const pw_alignment * to = *this_adj;
				//	to->print();
					std::map<const pw_alignment *, int>::iterator it3 = visited_als.find(*this_adj);
					assert(it3 != visited_als.end());
					std::map<const pw_alignment *, int>::iterator it4 = visited_als.find(current_node);
					assert(it4 != visited_als.end());
				//	std::cout << "low "<< it1->second << " times " << it3 ->second <<std::endl;
				//	std::cout << "low "<< it->second << " times " << it4 ->second <<std::endl;
				//	if(it3->second == it1->second+1 && it1->second == it4->second)//low adj = vis adj
					if(it1->second > it4->second){
					//	std::cout << "a bridge happened from  "<<std::endl;
					//	current_node->print();
					//	std::cout<< "to "<<std::endl;
						const pw_alignment * to = *this_adj;
					//	to->print();
						bridges.insert(std::make_pair(current_node, to));
					}else{
					}

				}
			}
		}
	}
//	std::cout << "end of while loop"<<std::endl;
}
void two_edge_cc::remove_bridges(std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als){
	std::map<const pw_alignment* , std::set<const pw_alignment*> > edges = cc_fraction.get_edges();
	for(std::multimap<const pw_alignment*,const pw_alignment*>::iterator it = bridges.begin() ; it != bridges.end() ; it++){
		std::map<const pw_alignment* , std::set<const pw_alignment*> >::iterator bridge_edge = edges.find(it->first);
		assert(bridge_edge != edges.end());
		std::set<const pw_alignment*>::iterator adj = bridge_edge->second.find(it->second);
		assert(adj != bridge_edge->second.end());
		bridge_edge->second.erase(adj);
	}
	std::map<const pw_alignment* , std::set<const pw_alignment*> > edges1 = cc_fraction.get_edges();
	if(bridges.size()>0){
		assert(edges1.size() > edges.size());
	}else{
		assert(edges1.size()==edges.size());
	}
	
	std::map<const pw_alignment*, bool> visited;
	//Adding the rest to 'stacks' using breadth first search algorithm to seperate different components of the 2 edge graph
	std::cout << "edges size is "<< edges.size() << std::endl;
	for(std::map<const pw_alignment* , std::set<const pw_alignment*> >:: iterator it = edges.begin() ; it != edges.end() ; it++){
		visited.insert(std::make_pair(it->first, false));
	}
	std::cout << "visited size "<< visited.size() <<std::endl;
	for(std::map<const pw_alignment* , std::set<const pw_alignment*> >:: iterator it = edges.begin() ; it != edges.end() ; it++){
		std::map<const pw_alignment*, bool>::iterator vis = visited.find(it->first);
		assert(vis != visited.end());
		if(vis->second == false){
		//	std::cout << "calling bfs "<<std::endl;
			bfs(it->first, edges,stacks,mixed_als,visited);
		}
	}

	
}
void two_edge_cc::bfs(const pw_alignment* s, std::map<const pw_alignment* , std::set<const pw_alignment*> > & edges, std::vector<std::set<const pw_alignment*,compare_pointer_pw_alignment> >& stacks, std::set<pw_alignment, compare_pw_alignment> & mixed_als,std::map<const pw_alignment*, bool> &visited){
	std::set<const pw_alignment*,compare_pointer_pw_alignment> this_stack;
	std::vector<const pw_alignment*> queue;
	std::map<const pw_alignment*, bool>::iterator vis = visited.find(s);
	assert(vis != visited.end());
	assert(vis->second == false);
	vis->second = true;
	queue.push_back(s);
	this_stack.insert(s);
	while(queue.size() > 0){
		s = queue.back();
		queue.pop_back();
		std::map<const pw_alignment* , std::set<const pw_alignment*> >::iterator it=edges.find(s);
		assert(it != edges.end());
		assert(it->second.size()>0);
		for(std::set<const pw_alignment*>::iterator adj = it->second.begin() ; adj != it->second.end(); adj++){
			std::map<const pw_alignment*, bool>::iterator vis1 = visited.find(*adj);
			assert(vis1 != visited.end());
			if(vis1->second == false){
				vis1->second = true;
				queue.push_back(vis1->first);
				this_stack.insert(vis1->first);
				assert(vis1->first == *adj);
			}
	
		}
	}
	assert(queue.size() == 0);
	assert(this_stack.size()>1);
//	std::cout << "this stack size "<< this_stack.size() << std::endl;
	if(this_stack.size()>2){
		stacks.push_back(this_stack);
	}else{
		assert(this_stack.size() ==2);
		for(std::set<const pw_alignment*,compare_pointer_pw_alignment>::iterator it = this_stack.begin() ; it != this_stack.end() ; it++){
			const pw_alignment* p = *it;
			mixed_als.insert(*p);
		}
	}

}
divide_and_conquer_alignments::divide_and_conquer_alignments(const all_data & data, const std::set<const pw_alignment*, compare_pointer_pw_alignment> & als_in, const dynamic_mc_model & model, size_t num_sequences, size_t num_threads):data(data), model(model), num_threads(num_threads), num_sequences(num_sequences), thread_info(num_threads, "waiting") {
	std::cout << "IN DCA" << std::endl;
	std::vector<size_t> alvec; // input for first problem
	for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::const_iterator it = als_in.begin(); it!=als_in.end(); ++it) {
		alvec.push_back(alvec.size());
		all_alignments.push_back(*it);
	}
//	divide_and_conquer_al_problem * all = new divide_and_conquer_al_problem(*this, alvec, 0.9, num_sequences, num_threads); // first one in parallel
//	all->run_cc();
//	all_problems.push_back(all);
//	wait_for_cc_problems.push_back(all);

}



void divide_and_conquer_alignments::run(std::vector< std::vector< pw_alignment> > & overlaps) {
	// create initial problem containing all alignments (this one runs on multiple threads, all child problems each run on a single thread)
	std::vector<size_t> ali(all_alignments.size());
	for(size_t i=0; i<all_alignments.size(); ++i) ali.at(i) = i;
	divide_and_conquer_al_problem * master_prob = new divide_and_conquer_al_problem(*this, ali, 0, num_sequences, num_threads);
	// create child problems
	master_prob->run_cc(0);
	wait_for_overlap_problems.push_back(master_prob);

	waiting_threads = num_threads;
	running_threads = 0;

// run all child problems and the problems created by them in parallel
#pragma omp parallel for schedule(static) num_threads(num_threads)
	for(size_t tr = 0; tr < num_threads; ++tr) {
		
		while(true) {
			bool got_a_task = select_next_task(tr);
			if(!got_a_task) {
				if(are_we_done()) {
					break;
				} else {
#pragma omp critical(print)
{
					if(num_threads > 1) std::cout << " Thread " << tr << " is waiting for a task " << std::endl;
}
					if(num_threads > 1) sleep(10); // wait for 10 seconds // TODO
				}
			}
		}
	}
	std::cout << " DIVI AND CONQ results " << std::endl;
	master_prob->get_results(overlaps);
	for(size_t i=0; i<overlaps.size(); ++i) {
		std::cout << " set " << i<< " has size " << overlaps.at(i).size() << std::endl;

	}

	delete master_prob;
}
	
bool divide_and_conquer_alignments::select_next_task(size_t thread) {


	divide_and_conquer_al_problem * cc_taken = NULL;
#pragma omp critical(scheduling)
{
	if(!wait_for_cc_problems.empty()) {
		cc_taken = wait_for_cc_problems.back();
		waiting_threads -= 1;
		running_threads += 1;
		wait_for_cc_problems.pop_back();
	}
}
	if(cc_taken!=NULL) {
		cc_taken->run_cc(thread);
#pragma omp critical(scheduling)
{
		waiting_threads += 1;
		running_threads -= 1;
		wait_for_overlap_problems.push_back(cc_taken);
}
		return true;
	}

	divide_and_conquer_al_problem * ol_taken = NULL;
	size_t next_al_cluster = (size_t) -1;
	bool alldone;
#pragma omp critical(scheduling)
{
// ol_taken stays in the queue so that all clusters can be done in parallel
	for(size_t i=0; i<wait_for_overlap_problems.size(); ++i) {
		ol_taken = wait_for_overlap_problems.at(i);
		bool takecl;
		ol_taken->select_next_alignment_cluster(next_al_cluster, takecl, alldone);
		if(!(takecl || alldone)) {
			ol_taken = NULL;

		} else {
			waiting_threads -= 1;
			running_threads += 1;
			break;
		}
	}
} 

	if(ol_taken!=NULL) {
// Here we run on only one set of clustered_alignments
		if(!alldone)
			ol_taken->run_overlap(next_al_cluster, thread);
		

#pragma omp critical(scheduling)
{
	waiting_threads+=1;
	running_threads-=1;
	if(ol_taken->are_all_alignment_clusters_done()) {
		// remove ol_taken from queue
		std::vector<divide_and_conquer_al_problem *> newvec;
		for(size_t i=0; i<wait_for_overlap_problems.size(); ++i) {
			if(wait_for_overlap_problems.at(i)!=ol_taken) newvec.push_back(wait_for_overlap_problems.at(i));
		}
		wait_for_overlap_problems = newvec;
		wait_for_combine_problems.push_back(ol_taken);
	}

}
		return true;
	}

	divide_and_conquer_al_problem * co_taken = NULL;
#pragma omp critical(scheduling)
{
// combine step depends on combine of child problems, therefore we have to search for problems which are ready for combine
	for(size_t i=0; i<wait_for_combine_problems.size(); ++i) {
		if(wait_for_combine_problems.at(i)->is_ready_for_combine()) {
			co_taken = wait_for_combine_problems.at(i);
			std::vector<divide_and_conquer_al_problem *> newvec;
			for(size_t j=0; j<wait_for_combine_problems.size(); ++j) {
				if(co_taken!=wait_for_combine_problems.at(j)) newvec.push_back(wait_for_combine_problems.at(j));
			}
			wait_for_combine_problems = newvec;
			waiting_threads-=1;
			running_threads+=1;
			break;
		}

	}

}
	if(co_taken!=NULL) {

		co_taken->run_combine(thread);

		



#pragma omp critical(scheduling)
{
		waiting_threads+=1;
		running_threads-=1;
// run_combine automatically inserts all results into the parent problem. Therefore the current problem can be deleted
		std::vector<divide_and_conquer_al_problem *> newall;
		for(size_t i=0; i<all_problems.size(); ++i) {
			if(co_taken!=all_problems.at(i)) newall.push_back(all_problems.at(i));
		}
		all_problems = newall;
}		
		return true;
	}


	return false;
}


	
bool divide_and_conquer_alignments::are_we_done() const {

	bool result = false;
#pragma omp critical(scheduling)
{
	if(running_threads==0) {
		if( (( wait_for_cc_problems.empty() &&  
		wait_for_overlap_problems.empty() ) && 
		wait_for_combine_problems.empty() ) ) result = true; 
	}
	std::cout << " CHECK done? Threads running "<< running_threads << " waiting " << waiting_threads  << " Q:  all "<< all_problems.size()<< " cc " << wait_for_cc_problems.size() << " ol " << wait_for_overlap_problems.size() <<  " co " << wait_for_combine_problems.size() << " done " << result <<std::endl;

	for(size_t i=0; i<thread_info.size(); ++i) {
		std::cout << " tr " << i<< " : " << thread_info.at(i) <<std::endl;
	}


}
	if(result) {
		assert(all_problems.size() == 0);
	}

	return result;

}

void divide_and_conquer_alignments::add_cc_problem(divide_and_conquer_al_problem* ap) {

#pragma omp critical(scheduling)
{
	all_problems.push_back(ap);
	wait_for_cc_problems.push_back(ap);
}

}

/*
	graph structure:
	Graph on subsets of all_alignments
	Graph has nodes 0 to alvec.size()-1
	Those ids are used for undirected edges

	node i to alignment: parent.all_alignments.at(alvec.at(i))


*/
// same code with and without multi-threading. Multithreaded code seems not to work when called from another multi-threaded block (even with num_threads=1)
void divide_and_conquer_al_problem::make_graph_singlethreaded(double overlap_fraction, std::vector<std::pair<size_t, size_t> > & edges) {
//	std::cout << " RUN CC single threaded " << std::endl;
	std::map<size_t, size_t> used_refs; // reference sequence number -> index in trees/al_on_seq
	std::vector<std::vector<std::pair<size_t, size_t> > > al_on_seq;
	for(size_t i=0; i<alvec.size(); i++) {
		const pw_alignment * al = parent.all_alignments.at(alvec.at(i));
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		std::map<size_t, size_t>::iterator it1 = used_refs.find(r1);
		if(it1==used_refs.end()) {
			used_refs.insert(std::make_pair(r1, al_on_seq.size()));
			it1=used_refs.find(r1);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}
		std::map<size_t, size_t>::iterator it2 = used_refs.find(r2);
		if(it2==used_refs.end()) {
			used_refs.insert(std::make_pair(r2, al_on_seq.size()));
			it2=used_refs.find(r2);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}


		al_on_seq.at(it1->second).push_back(std::make_pair(0, i));
		al_on_seq.at(it2->second).push_back(std::make_pair(1, i));
	}
	std::vector<tree_type> trees(al_on_seq.size());

//	std::cout << " Using " << trees.size() << " sequences " << std::endl;
		for(size_t s=0; s<al_on_seq.size(); s++) {
			for(size_t i=0; i<al_on_seq.at(s).size(); ++i) {
				std::pair<size_t, size_t> ad = al_on_seq.at(s).at(i);
				const pw_alignment * al = parent.all_alignments.at(alvec.at(ad.second));
				size_t l,r;
				if(ad.first) {
					al->get_lr2(l,r);
				} else {
					al->get_lr1(l,r);
				}
				trees.at(s).insert(l, r, ad.second);
			}	
			std::vector<std::pair<size_t, size_t> > edges_s;
#pragma omp critical(print) 
{
//			std::cout<< " allpairs overlap on " << al_on_seq.at(s).size() << std::endl;
}
			trees.at(s).allpairs_overlap(overlap_fraction, edges_s);
#pragma omp critical(edges) 
{
			for(size_t i=0; i<edges_s.size(); ++i) {
				if(edges_s.at(i).first != edges_s.at(i).second) // alignments can overlap with itself. I that case no edge to keep a simple graph
					edges.push_back(edges_s.at(i));
			}
}
		}
#pragma omp critical(print)
{
//	std::cout << " We found " << edges.size() << " edges " << std::endl;
}



}


void divide_and_conquer_al_problem::make_graph(double overlap_fraction, std::vector<std::pair<size_t, size_t> > & edges) {

 //	std::cout << " RUN CC " << std::endl;
	std::map<size_t, size_t> used_refs; // reference sequence number -> index in trees/al_on_seq
	std::vector<std::vector<std::pair<size_t, size_t> > > al_on_seq;
	for(size_t i=0; i<alvec.size(); i++) {
		const pw_alignment * al = parent.all_alignments.at(alvec.at(i));
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		std::map<size_t, size_t>::iterator it1 = used_refs.find(r1);
		if(it1==used_refs.end()) {
			used_refs.insert(std::make_pair(r1, al_on_seq.size()));
			it1=used_refs.find(r1);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}
		std::map<size_t, size_t>::iterator it2 = used_refs.find(r2);
		if(it2==used_refs.end()) {
			used_refs.insert(std::make_pair(r2, al_on_seq.size()));
			it2=used_refs.find(r2);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}


		al_on_seq.at(it1->second).push_back(std::make_pair(0, i));
		al_on_seq.at(it2->second).push_back(std::make_pair(1, i));
	}
	std::vector<tree_type> trees(al_on_seq.size());

//	std::cout << " Using " << trees.size() << " sequences " << std::endl;
	// same code with and without multi-threading. Multithreaded code seems not to work when called from another multi-threaded block	
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
		for(size_t s=0; s<al_on_seq.size(); s++) {
			for(size_t i=0; i<al_on_seq.at(s).size(); ++i) {
				std::pair<size_t, size_t> ad = al_on_seq.at(s).at(i);
				const pw_alignment * al = parent.all_alignments.at(alvec.at(ad.second));
				size_t l,r;
				if(ad.first) {
					al->get_lr2(l,r);
				} else {
					al->get_lr1(l,r);
				}
				trees.at(s).insert(l, r, ad.second);
			}	
			std::vector<std::pair<size_t, size_t> > edges_s;
#pragma omp critical(print) 
{
//			std::cout<< " allpairs overlap on " << al_on_seq.at(s).size() << std::endl;
}
			trees.at(s).allpairs_overlap(overlap_fraction, edges_s);
#pragma omp critical(edges) 
{
			for(size_t i=0; i<edges_s.size(); ++i) {
				if(edges_s.at(i).first != edges_s.at(i).second) // alignments can overlap with itself. I that case no edge to keep a simple graph
					edges.push_back(edges_s.at(i));
			}
}
		}
#pragma omp critical(print)
{
//	std::cout << " We found " << edges.size() << " edges " << std::endl;
}


}

/*

	Third make graph function:
	make graph of alignments not contained in parent.all_alignments

	single thread implementation only (result graphs are much smaller)

*/
void divide_and_conquer_al_problem::make_graph(double overlap_fraction, const std::vector<pw_alignment> & als, std::vector<std::pair<size_t, size_t> > & edges) {
	std::map<size_t, size_t> used_refs; // reference sequence number -> index in trees/al_on_seq
	std::vector<std::vector<std::pair<size_t, size_t> > > al_on_seq;
	for(size_t i=0; i<als.size(); i++) {
		const pw_alignment * al = &(als.at(i));
		size_t r1 = al->getreference1();
		size_t r2 = al->getreference2();
		std::map<size_t, size_t>::iterator it1 = used_refs.find(r1);
		if(it1==used_refs.end()) {
			used_refs.insert(std::make_pair(r1, al_on_seq.size()));
			it1=used_refs.find(r1);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}
		std::map<size_t, size_t>::iterator it2 = used_refs.find(r2);
		if(it2==used_refs.end()) {
			used_refs.insert(std::make_pair(r2, al_on_seq.size()));
			it2=used_refs.find(r2);
			al_on_seq.push_back(std::vector<std::pair<size_t, size_t> >());
		}


		al_on_seq.at(it1->second).push_back(std::make_pair(0, i));
		al_on_seq.at(it2->second).push_back(std::make_pair(1, i));
	}
	std::vector<tree_type> trees(al_on_seq.size());

//	std::cout << " Using " << trees.size() << " sequences " << std::endl;
	for(size_t s=0; s<al_on_seq.size(); s++) {
		for(size_t i=0; i<al_on_seq.at(s).size(); ++i) {
			std::pair<size_t, size_t> ad = al_on_seq.at(s).at(i);
			const pw_alignment * al = &(als.at(ad.second));
			size_t l,r;
			if(ad.first) {
				al->get_lr2(l,r);
			} else {
				al->get_lr1(l,r);
			}
			trees.at(s).insert(l, r, ad.second);
		}	
		std::vector<std::pair<size_t, size_t> > edges_s;
#pragma omp critical(print) 
{
//		std::cout<< " allpairs overlap on " << al_on_seq.at(s).size() << std::endl;
}
		trees.at(s).allpairs_overlap(overlap_fraction, edges_s);
#pragma omp critical(edges) 
{
		for(size_t i=0; i<edges_s.size(); ++i) {
			if(edges_s.at(i).first != edges_s.at(i).second) // alignments can overlap with itself. I that case no edge to keep a simple graph
				edges.push_back(edges_s.at(i));
		}
}
	}
#pragma omp critical(print)
{
//	std::cout << " We found " << edges.size() << " edges " << std::endl;
}




}


/*
	should only be called from critical(scheduling)

*/
void divide_and_conquer_al_problem::select_next_alignment_cluster(size_t & idx, bool & found_next, bool & all_done) {
	found_next = false;
	all_done = true;

// last one (leftover alignments) first
	if(clustered_al_status.size() > 0) {
		size_t l = clustered_al_status.size() -1 ;
		if(clustered_al_status.at(l) < 2) {
			all_done = false;
		}
		if(clustered_al_status.at(l) == 0) {
			clustered_al_status.at(l) = 1;
			idx = l;
			found_next = true;	
			return;

		}
	}
	for(size_t i=0; i<clustered_al_status.size()-1; ++i) {
		if(clustered_al_status.at(i) < 2) {
			all_done = false;
		}
		if(clustered_al_status.at(i)==0) {
			clustered_al_status.at(i) = 1; // running
			idx = i;
			found_next = true;
			return;
		}
	}
	return;
}

/*
	should only be called from critical(scheduling)

*/
bool divide_and_conquer_al_problem::are_all_alignment_clusters_done() {
	for(size_t i=0; i<clustered_al_status.size(); ++i) {
		if(clustered_al_status.at(i) < 2) {
			return false;
		}
	}

	state = WAIT_FOR_COMBINE;
	return true;
}

/*
	should only be called from critical(scheduling)

*/
bool divide_and_conquer_al_problem::is_ready_for_combine() {
	assert(state == WAIT_FOR_COMBINE);
	return child_problems.empty();
}

void divide_and_conquer_al_problem::set_parent_problem(divide_and_conquer_al_problem * p) {
	parent_problem = p;
}
/*
	
*/
void divide_and_conquer_al_problem::run_cc(size_t thread) {
	assert(state == WAIT_FOR_CC);
	
#pragma omp critical(scheduling)
{
	std::stringstream str;
	str << " run_cc on l " << separation_level << " with " << alvec.size() << " alignments ";
	parent.thread_info.at(thread) = str.str(); 
}
#pragma omp critical(print)
{
//	std::cout << " ON level " << separation_level << " we have " << alvec.size() << " alignments " << std::endl;
}

	

	if(separation_level < 40) {
// how strongly do we need to separate at the current level?
		double overlap_fraction;
		size_t cc_edges;
		separation_conditions(overlap_fraction, cc_edges);

#pragma omp critical(print)
{
	std::cout << " ON level " << separation_level << " overlap " << overlap_fraction << " edges " << cc_edges << " and " << alvec.size() << " alignments, threads: " << num_threads<< std::endl;
}



// create graph of pairwise overlap
		std::vector<std::pair<size_t, size_t> >  edges;
		if(num_threads > 1) 
			make_graph(overlap_fraction, edges);
		else
			make_graph_singlethreaded(overlap_fraction, edges);
// approximate cc_edges-edge connected components
		std::vector<std::set<size_t> > cc_result;	
		two_edge_cc(edges, cc_edges, cc_result);
#pragma omp critical(print)
{
//	std::cout << " ON level " << separation_level << " with " << alvec.size() << " alignments, we made  " << cc_result.size() << " connected components" << std::endl;
}


		std::set<size_t> checkset;
		for(size_t i=0; i<cc_result.size(); ++i) {
			std::set<size_t> & cc = cc_result.at(i);
			if(cc.size() > MAX_CC_SIZE) {
				std::vector<size_t> nals;
				for(std::set<size_t>::iterator it = cc.begin(); it!=cc.end(); ++it) {
					nals.push_back(alvec.at(*it));
					checkset.insert(alvec.at(*it));
				}
#pragma omp critical(print)
{
	std::cout << " ON level " << separation_level << " with " << alvec.size() << " alignments, we made  a child problem with " << nals.size() << " alignments" << std::endl;
}

/*				divide_and_conquer_al_problem * childprob = new divide_and_conquer_al_problem(parent, nals, separation_level + 1, num_sequences, 1);
				childprob->set_parent_problem(this);
				child_problems.push_back(childprob);
				parent.add_cc_problem(childprob);
*/			
			} else if(cc.size() >= MIN_CC_SIZE) {
/*				clustered_alignments.push_back(std::vector<size_t>());
				for(std::set<size_t>::iterator it = cc.begin(); it!=cc.end(); ++it) {
					clustered_alignments.at(clustered_alignments.size()-1).push_back(alvec.at(*it));
					checkset.insert(alvec.at(*it));
				}
*/			} else {
				for(std::set<size_t>::iterator it = cc.begin(); it!=cc.end(); ++it) {
					leftover_alignments.push_back(alvec.at(*it));
					checkset.insert(alvec.at(*it));
				}
			}

		}
//		assert(alvec.size() == checkset.size());
	} else {
		// We somehow cannot separate this alignment cluster further, we have to run slow code on a large data set
		for(size_t i=0; i<alvec.size(); ++i) {
			std::vector<size_t> p(1, i);
			clustered_alignments.push_back(p);
		}
	}

// prepare for next step
//	clustered_al_results.resize(clustered_alignments.size());
	clustered_al_status = std::vector<size_t>(clustered_alignments.size() + 1, 0); // one for leftovers at end

	state = WAIT_FOR_OVERLAP;

	
#pragma omp critical(scheduling)
{
	parent.thread_info.at(thread) = "waiting"; 
}
}



void divide_and_conquer_al_problem::run_overlap(size_t cl_al_idx, size_t thread) {
	assert(state == WAIT_FOR_OVERLAP);


	if(cl_al_idx < clustered_al_status.size()) {
	if(cl_al_idx == clustered_al_status.size() - 1) {
std::stringstream str;
#pragma omp critical(scheduling)
{
	str << separation_level << " run_overlap on group " << cl_al_idx << " on  " << leftover_alignments.size() << " leftover alignments " ;
	parent.thread_info.at(thread) = str.str(); 
}


	// run leftovers first
		std::vector<std::vector<const pw_alignment *> > pvect;
		for(size_t i=0; i<leftover_alignments.size(); ++i) {
			std::vector<const pw_alignment *> nv;
			nv.push_back((parent.all_alignments.at(leftover_alignments.at(i))));
			pvect.push_back(nv);
		}
		use_ias ias(parent.data, pvect, parent.model);
		ias.compute(results, str.str(), parent.thread_info.at(thread));

#pragma omp critical(scheduling)
{
	clustered_al_status.at(cl_al_idx) = 2;
}


	} else {

	std::stringstream str;
#pragma omp critical(scheduling)
{
	str << separation_level << " run_overlap on group " << cl_al_idx << " with " << clustered_alignments.at(cl_al_idx).size() << " alignments ";
	parent.thread_info.at(thread) = str.str(); 
}


//		std::cout << " RUN overlap on " << cl_al_idx << " size " << clustered_alignments.at(cl_al_idx).size() << std::endl;
		assert(clustered_al_status.at(cl_al_idx) == 1);
		// make groups of size 1
		std::vector<std::vector<const pw_alignment *> > algr(clustered_alignments.at(cl_al_idx).size());
		for(size_t i=0; i<clustered_alignments.at(cl_al_idx).size(); ++i) {
			algr.at(i).push_back(parent.all_alignments.at(clustered_alignments.at(cl_al_idx).at(i)));
		}


		use_ias ias(parent.data, algr, parent.model);
		overlap_type ores(parent.data);

		
		ias.compute(ores, str.str(), parent.thread_info.at(thread));

//  slow check
//		std::cout << " Test on " << ores.get_all().size() << std::endl;
//		ores.test_partial_overlap();

		const std::set<pw_alignment, compare_alignment<pw_alignment> > & all = ores.get_all();
		std::vector<pw_alignment> avec;
		for(std::set<pw_alignment, compare_alignment<pw_alignment> >::const_iterator it = all.begin(); it!=all.end(); ++it) {
			avec.push_back(*it);
		}
		std::vector<std::vector<pw_alignment> > ccs;
		als_to_ccs(avec, ccs);
//		std::cout << " in overlap ias made " << avec.size() << " alignments in " << ccs.size() << " connected components"<< std::endl;

#pragma omp critical(scheduling)
{
		for(size_t i=0; i<ccs.size(); ++i) {
			clustered_al_results.push_back(ccs.at(i));
		}
		clustered_al_status.at(cl_al_idx) = 2;
}

	}
	
	}


#pragma omp critical(scheduling)
{
	parent.thread_info.at(thread) = "waiting"; 
}

}



/*
	Tarjan's bridge finding algorithm
Tarjan. Information Processing Letters. 1974.
*/
/*void divide_and_conquer_al_problem::bridge_step(size_t & time, const std::map<size_t, std::set<size_t> > & edges, const size_t & node, std::vector<size_t> & parent, std::vector<bool> & visited, std::vector<size_t> & dfsnum, std::vector<size_t> & low, std::vector<std::pair<size_t, size_t> > & bridges) {
	dfsnum.at(node) = ++time;
	low.at(node) = time;
	visited.at(node) = true;
	std::map<size_t, std::set<size_t> >::const_iterator findn=edges.find(node);
	if(findn!=edges.end()) { 
		const std::set<size_t> & nedges = findn->second;
		for(std::set<size_t>::const_iterator it = nedges.begin(); it!=nedges.end(); ++it) {
			size_t next = *it;
			assert(next != node);
			parent.at(next) = node;
			if(!visited.at(next)) {
				bridge_step(time, edges, next, parent, visited, dfsnum, low, bridges);
				low.at(node) = std::min(low.at(node), low.at(next));
				// found a bridge
				if(low.at(next) > dfsnum.at(node)) {
					bridges.push_back(std::make_pair(node, next));
				}
			} else {
				if(next != parent.at(node)) {
					// found a back-edge:
					low.at(node) = std::min(low.at(node), dfsnum.at(next));
				}
			}
		}
	}	
}*/
/*
	Tarjan's bridge finding algorithm
Tarjan. Information Processing Letters. 1974.

nonrecursive implementation

*/
void divide_and_conquer_al_problem::bridge_step(size_t & time, const std::map<size_t, std::set<size_t> > & edges, const size_t & n, std::vector<size_t> & parent, std::vector<bool> & visited, std::vector<size_t> & dfsnum, std::vector<size_t> & low, std::vector<std::pair<size_t, size_t> > & bridges) {
	std::vector<size_t> nstack;
	nstack.push_back(n);

	while(!nstack.empty()) {
		// take last element of stack, it stays on stack for a second processing after "recursion"
		size_t cur_node = nstack.back();
		if(!visited.at(cur_node)){ // before "recursion" at cur_node
			dfsnum.at(cur_node) = ++time;
			low.at(cur_node) = time;
			visited.at(cur_node) = true;
			std::map<size_t, std::set<size_t> >::const_iterator findn=edges.find(cur_node);
			if(findn!=edges.end()){ 
				const std::set<size_t> & nedges = findn->second;
				for(std::set<size_t>::const_iterator it = nedges.begin(); it!=nedges.end(); ++it){
					size_t next = *it;
					assert(next != cur_node);

					if(!visited.at(next)){
						parent.at(next) = cur_node;
						// put not visited not on the stack (like recursive call)	
						nstack.push_back(next);
					} else{
					//	if(next != parent.at(node)) {
							// found a back-edge:
							low.at(cur_node) = std::min(low.at(cur_node), dfsnum.at(next));
					//	}
					}
				} 
			}


		}else{ // after "recursion" at cur_node. Entire subtree below cur_node was done
			nstack.pop_back(); // remove cur_node from stack, we are done with this subtree
			std::map<size_t, std::set<size_t> >::const_iterator findn=edges.find(cur_node);
			if(findn!=edges.end()) { 
				const std::set<size_t> & nedges = findn->second;
				for(std::set<size_t>::const_iterator it = nedges.begin(); it!=nedges.end(); ++it) {
					size_t next = *it;
					assert(next != cur_node);
					if(parent.at(next) == cur_node){ // after recursion step for cur_node->next
						low.at(cur_node) = std::min(low.at(cur_node), low.at(next));
						// found a bridge
						if(low.at(next) > dfsnum.at(cur_node)) {
							bridges.push_back(std::make_pair(cur_node, next));
						}

					}

				}
			}
		}
	}



}

void divide_and_conquer_al_problem::edge_insert(std::map<size_t, std::set<size_t> > & edges, size_t i, size_t j) {
	std::map<size_t, std::set<size_t> >::iterator findi = edges.find(i);
	if(findi == edges.end()) {
		edges.insert(std::make_pair(i, std::set<size_t>()));
		findi = edges.find(i);
	} 
	std::set<size_t> & jset = findi->second;
	jset.insert(j);
	std::map<size_t, std::set<size_t> >::iterator findj = edges.find(j);
	if(findj == edges.end()) {
		edges.insert(std::make_pair(j, std::set<size_t>()));
		findj = edges.find(j);
	} 
	std::set<size_t> & iset = findj->second;
	iset.insert(i);
}

void divide_and_conquer_al_problem::edge_remove(std::map<size_t, std::set<size_t> > & edges, size_t i, size_t j) {
	std::map<size_t, std::set<size_t> >::iterator findi = edges.find(i);
	if(findi != edges.end()) {
		std::set<size_t> & jset = findi->second;
		jset.erase(j);
	} 
	std::map<size_t, std::set<size_t> >::iterator findj = edges.find(j);
	if(findj != edges.end()) {
		std::set<size_t> & iset = findj->second;
		iset.erase(i);
	}
}

size_t divide_and_conquer_al_problem::count_edges(const std::map<size_t, std::set<size_t> > & edges) {
	size_t res = 0;
	for( std::map<size_t, std::set<size_t> >::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
		res+=(it->second).size();
	}
	return res/2;
}


void divide_and_conquer_al_problem::edges_remove(const std::vector<std::pair<size_t, size_t> > & redges, std::map<size_t, std::set<size_t> > & medges) {
	std::multimap<size_t, size_t> edgemap; // new edgemap

	std::map<size_t, std::set<size_t> > bmap; // map of edges that will be removed
	for(size_t i=0; i<redges.size(); ++i) {
		std::pair<size_t, size_t> b = redges.at(i);
		edge_remove(medges, b.first, b.second);
	}
}




void divide_and_conquer_al_problem::cc_step(const size_t & node, const std::map<size_t, std::set<size_t> > & edges, std::vector<bool> & visited, std::set<size_t> & cur_cc, std::vector<std::pair<size_t, size_t> > & tree_edges) {
	visited.at(node) = true;
	cur_cc.insert(node);
	std::map<size_t, std::set<size_t> >::const_iterator findn = edges.find(node);
	if(findn!=edges.end()) {	
		const std::set<size_t> & nedges = findn->second;

		for(std::set<size_t>::const_iterator it = nedges.begin(); it!=nedges.end(); ++it) {
			size_t next = *it;
			if(!visited.at(next)) {
				tree_edges.push_back(std::make_pair(node, next));
				cc_step(next, edges, visited, cur_cc, tree_edges);
			}
		}
	}
}



void divide_and_conquer_al_problem::path_remove(std::map<size_t, std::set<size_t> > & edgemap, std::vector<std::set<size_t> > & cc_result) {
//	std::cout << " PATH REMOVE on " << count_edges(edgemap) << " edges and " << alvec.size() << " nodes " << std::endl;
	std::vector<bool> visited(alvec.size(), 0);
	std::vector<size_t> parent(alvec.size(), 0);
	std::vector<size_t> dfsnum(alvec.size(), 0);
	std::vector<size_t> low(alvec.size(), 0);
	std::vector<std::pair<size_t, size_t> > bridges;
// find all bridges
	size_t time = 0;
	for(size_t i=0; i<alvec.size(); ++i) {
		if(!visited.at(i)) {
			bridge_step(time, edgemap, i, parent, visited, dfsnum, low, bridges);
		}
	}

//	std::cout << " FOUND " << bridges.size() << " bridges " << std::endl;

// remove bridges from map of edges
	edges_remove(bridges, edgemap);

//	std::cout << " After bridge remove: edges " << count_edges(edgemap) << std::endl;

// find connected components on remaining edges
	std::vector<std::set<size_t> > connected_components(0);
	visited = std::vector<bool>(alvec.size(), 0);
	std::vector<std::pair<size_t, size_t> > tree_edges;
	for(size_t i=0; i<alvec.size(); ++i) {
		if(!visited.at(i)) {
			connected_components.push_back(std::set<size_t>());
			cc_step(i, edgemap, visited, connected_components.at(connected_components.size() -1 ), tree_edges);
		}
	}
	size_t sum = 0;
//	std::cout << " FOUND " << connected_components.size() << " 2-edge connected components ";
	std::multimap<size_t, std::set<size_t> > sorter;
	for(size_t i=0; i<connected_components.size(); ++i) {
		sorter.insert(std::make_pair(connected_components.at(i).size(), connected_components.at(i)));
	}

	cc_result = std::vector<std::set<size_t> >(sorter.size());
	size_t i=0;
	for(std::multimap<size_t, std::set<size_t> >::reverse_iterator it = sorter.rbegin(); it!=sorter.rend(); ++it) {
		std::set<size_t> & cc = it->second;
//		std::cout << " " << cc.size();
		sum+=cc.size();
		cc_result.at(i) = cc;
		i++;
	}
//	std::cout << std::endl; 

	edges_remove(tree_edges, edgemap);

//	std::cout << " EDGES left " << count_edges(edgemap) <<std::endl;
	assert(sum == alvec.size());

}

/*
	compute k-edge connected components
	for k>2 this algorithm computes an approximation where the partition does not have to have the minimal number of sets.
	Each set returned will be k-edge connected

*/
void divide_and_conquer_al_problem::two_edge_cc(const std::vector<std::pair<size_t, size_t> > & edges, size_t cc_edges, std::vector<std::set<size_t> >  & cc_result) {

// put edges in map
	// compute 2-edge cc
	std::vector<size_t> preorder_numbers(alvec.size());
	std::map<size_t, std::set<size_t> > edgemap; // node i -> set of nodes adjacent to e. We save undirected edges in both direction
	for(size_t i=0; i<edges.size(); ++i) {
		std::pair<size_t, size_t> e = edges.at(i);
		edge_insert(edgemap, e.first, e.second);
		edge_insert(edgemap, e.second, e.first);
	}

//	std::cout << " graph has " << edgemap.size() << " nodes and " << edges.size() << " edges " << std::endl;

	assert(cc_edges >= 2);
	for(size_t i=1; i<cc_edges; i++) {
		path_remove(edgemap, cc_result);
	}

}



void divide_and_conquer_al_problem::child_done(divide_and_conquer_al_problem * child, const std::vector<std::vector<pw_alignment> > & child_res) {
	assert(state == WAIT_FOR_COMBINE || state == WAIT_FOR_OVERLAP);
	assert(child->state == DONE);


#pragma omp critical(scheduling) 
{
//	std::cout << " Insert child into vector of ovrlps " << std::endl;
	for(size_t i=0; i<child_res.size(); ++i) {
		clustered_al_results.push_back(child_res.at(i));
	}

	std::vector<divide_and_conquer_al_problem *> newvec;
	for(size_t i=0; i<child_problems.size(); ++i) {
		if(child!=child_problems.at(i)) newvec.push_back(child_problems.at(i));
	}
	child_problems = newvec;

	delete child;
}


}


/*
	Find connected compenents of alignments 

*/
void divide_and_conquer_al_problem::als_to_ccs(std::vector<pw_alignment> & als, std::vector<std::vector<pw_alignment> > & ccs) {

	std::vector<std::pair<size_t, size_t> > edges;	
	make_graph(0.0, als, edges);
	std::map<size_t, std::set<size_t> > ograph;
	for(size_t i=0; i<edges.size(); ++i) {
		edge_insert(ograph, edges.at(i).first, edges.at(i).second);
	}
	
	std::vector<std::set<size_t> > connected_components;
	std::vector<bool> visited(als.size(), 0);
	std::vector<std::pair<size_t, size_t> > tree_edges;
	for(size_t i=0; i<als.size(); ++i) {
		if(!visited.at(i)) {
			connected_components.push_back(std::set<size_t>());
			cc_step(i, ograph, visited, connected_components.at(connected_components.size() -1 ), tree_edges);
		}
	}
	ccs = std::vector<std::vector<pw_alignment> >(connected_components.size());
	for(size_t i=0; i<connected_components.size(); ++i) {
		for(std::set<size_t>::iterator it = connected_components.at(i).begin(); it!=connected_components.at(i).end(); ++it) {
			ccs.at(i).push_back(als.at(*it));

		}
//		std::cout << " FOUND CC size " << ccs.at(i).size() << std::endl;
	}


}
	
void divide_and_conquer_al_problem::run_combine(size_t thread) {
	assert(state==WAIT_FOR_COMBINE);
	size_t in_r_init = results.get_all().size();
#pragma omp critical(scheduling)
{
	std::stringstream str;
	str <<"l "<< separation_level<< " run_combine on " << clustered_al_results.size() << " alignment groups and " << in_r_init << " leftover alignments ";
	parent.thread_info.at(thread) = str.str(); 
}
	std::cout << " IN combine level " << separation_level<< std::endl;
	size_t sum = 0;
	std::vector<std::vector<const pw_alignment *> > pvect(clustered_al_results.size());
	for(size_t i=0; i<clustered_al_results.size(); ++i) {
		std::vector<const pw_alignment *> v(clustered_al_results.at(i).size());
		for(size_t j=0; j<clustered_al_results.at(i).size(); ++j) {
			v.at(j) = &(clustered_al_results.at(i).at(j));
			sum++;
		}
		pvect.at(i) = v;
	}	
//	for(size_t i=0; i<leftover_alignments.size(); ++i) {
//		std::vector<const pw_alignment *> nv;
//		nv.push_back((parent.all_alignments.at(leftover_alignments.at(i))));
//		pvect.push_back(nv);
//	}


std::stringstream str1;
#pragma omp critical(scheduling)
{
	str1 << "l "<< separation_level<< " run_combine on " << clustered_al_results.size() << " alignment groups and " << in_r_init << " leftover alignments, total  "<< sum + in_r_init;
	parent.thread_info.at(thread) = str1.str(); 
}


	std::cout << " co in " << sum << " plus " << leftover_alignments.size() << std::endl;
	use_ias ias(parent.data, pvect, parent.model);
	ias.compute(results, str1.str(), parent.thread_info.at(thread));
#pragma omp critical(scheduling)
{
	std::stringstream str2;
	str2 << "l "<< separation_level<<  " run_combine on " << clustered_al_results.size() << " alignment groups and " << leftover_alignments.size() << " leftover alignments, total  "<< sum + in_r_init + " ias done";
	parent.thread_info.at(thread) = str2.str(); 
}




// slow check for partial overlap TODO 
//	results.test_partial_overlap();

	const std::set<pw_alignment, compare_alignment<pw_alignment> >& rset = results.get_all();
	std::vector<pw_alignment> avec(rset.size());
	size_t i=0;

	double sumg = 0;

	for(std::set<pw_alignment, compare_alignment<pw_alignment> >::const_iterator it = rset.begin(); it!=rset.end(); ++it) {
		avec.at(i) = *it;
		i++;

		double g1,g2;
		parent.model.gain_function(*it, g1,g2);
		sumg += (g1+g2)/2.0;
	}


//	std::vector<std::vector<pw_alignment> > ccs_als;
	als_to_ccs(avec, clustered_al_results);


	std::cout << "l "<< separation_level<<  " combine result has " << avec.size() <<" alignments, gain "<< sumg << " connected components " << clustered_al_results.size()<< std::endl;


//	std::cout << "combine  on " << rset.size() << " alignments, we found " << clustered_al_results.size() << " connected components"<<std::endl;
	
	
	state = DONE;

#pragma omp critical(scheduling)
{
	parent.thread_info.at(thread) = "waiting"; 
}
	if(parent_problem!=NULL) { // if not in master problem
		// attention: deletes current object
		parent_problem->child_done(this, clustered_al_results);
	} else {

	}
}


/* 
In each iteration, we want to separate alignments more strongly (as long as there are clusters that are to large to process them directly)

*/
void divide_and_conquer_al_problem::separation_conditions(double & overlap_fraction, size_t & cc_edges) {
	cc_edges = 2; // 2-edge connected components TODO 
	if(separation_level == 0) {
		overlap_fraction = 0.5; // 0.8;
		return;
	}
	if(separation_level == 1) {
		cc_edges = 3;
		overlap_fraction = 0.7; // 0.95;
		return;
	}
	if(separation_level == 2) {
		cc_edges = 4;
		overlap_fraction = 0.9; // 0.99;
		return;
	}
	if(separation_level == 3) {
		cc_edges = 5;
		overlap_fraction = 0.95; // 0.999;
		return;
	}
	if(separation_level == 4) {
		cc_edges = 6;
		overlap_fraction = 0.99;
		return;
	}
	overlap_fraction = 0.999;
	cc_edges = separation_level + 2; // TODO
	if(separation_level > 10) overlap_fraction = 0.9999;
	if(separation_level > 15) overlap_fraction = 0.99999;
}


bool divide_and_conquer_al_problem::has_parent() const {
	return (parent_problem!=NULL);
}
void divide_and_conquer_al_problem::get_results(std::vector<std::vector<pw_alignment> > & child_res) {
	assert(parent_problem == NULL); // results from master problem
	child_res = clustered_al_results;
}
#include "intervals.cpp"
#include "overlap.cpp"


#endif
