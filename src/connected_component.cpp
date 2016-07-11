#include "connected_component.hpp"

	void compute_cc_avl_fraction::non_recursive_compute(std::set< const pw_alignment* , compare_pointer_pw_alignment> & ccs, std::set<const pw_alignment* , compare_pointer_pw_alignment> & remainder, double & fraction){//It creats a graph and keep all the edges in a multimap. Graph might have more than one component
	//	for(size_t i = 0; i < als.size();i++) {
	//		const pw_alignment * al = als.at(i);
		edges.clear();
		std::cout << "edges size1 "<<edges.size()<<std::endl;
		for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainder.begin(); it!= remainder.end();it++){//TODO time consuming! Need to be parallelized , can simply change the set to a vector
			const pw_alignment * al = *it;
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
	//	std::cout<< "ccs size "<<ccs.size()<<std::endl;
		for(std::set< const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.begin(); it != ccs.end();it ++){
			std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it2 = remained_als.find(*it);
			assert(it2 != remained_als.end());
			remained_als.erase(*it2);
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

	void compute_cc_avl_fraction::cc_step_non_recursive(const pw_alignment*  al,  size_t & ref, size_t & left , size_t & right, std::set< const pw_alignment* , compare_pointer_pw_alignment> & ccs, double & fraction, std::set<const pw_alignment* , compare_pointer_pw_alignment> & remainder){
			std::vector<const pw_alignment*> results;
			std::vector<const pw_alignment*> fraction_result;
			alind->search_overlap_fraction(ref, left, right, fraction , fraction_result);
		//	alind->search_overlap(ref, left, right, results);
		//	compute_fraction(results, fraction,ref, left, right, fraction_result);
			for(size_t i=0; i<fraction_result.size(); ++i) {
				std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = remainder.find(fraction_result.at(i));
				if(al != fraction_result.at(i) && it != remainder.end()){
#pragma omp critical(edge)
{
				//	std::cout << "al " << al << " fraction at " << i << " : "<< fraction_result.at(i) <<std::endl;
					edges.insert(std::make_pair(al,fraction_result.at(i)));
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
				for(std::multimap<const pw_alignment* , const pw_alignment*>::iterator it = edges.begin(); it!=edges.end(); it++){
					std::cout << it->first << " " << it->second << std::endl;
				}
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
					remained_als.erase(*it2);
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
				std::cout<< edges.size()<<std::endl;
				edges.insert(std::make_pair(&al,fraction_result.at(i)));
				std::cout<< edges.size()<<std::endl;
				edges.insert(std::make_pair(fraction_result.at(i),&al));
				std::cout<< edges.size()<<std::endl;
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
	void compute_cc_avl_fraction::node_adjucents(const pw_alignment* p ,std::set<const pw_alignment*> & adj)const{
		std::pair<std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator ,std::multimap<const pw_alignment*,const pw_alignment*>::const_iterator > it = edges.equal_range(p);
		for(std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator it1 = it.first; it1!=it.second; ++it1) {
			adj.insert(it1->second);
		}
		

	}
	void compute_cc_avl_fraction::node_adjucents(const pw_alignment* p ,std::vector<const pw_alignment*> & adj)const{
		std::pair<std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator ,std::multimap<const pw_alignment*,const pw_alignment*>::const_iterator > it = edges.equal_range(p);
		for(std::multimap<const pw_alignment*, const pw_alignment*>::const_iterator it1 = it.first; it1!=it.second; ++it1) {
			adj.push_back(it1->second);
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
		cc_fraction.node_adjucents(alignments.at(at).at(node), adj);
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
		cc_fraction.node_adjucents(parent, adj);
	//	std::cout<<"parent "<< parent << " adj size "<< adj.size()<<std::endl;
		for(size_t i =0; i < adj.size(); i++){
		//	std::cout<< adj.at(i)<<std::endl;
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
				//	std::cout << parent << " to " << adj.at(i) << " is a bridge."<<std::endl;
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
	
