#include "strongly_connected_alignments.hpp"


	void detect_overlap_rate::compute(std::vector<std::set< const pw_alignment* , compare_pointer_pw_alignment> > & strongly_connected){
		for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = alignments.begin(); it!= alignments.end(); it ++){
			const pw_alignment * al = *it;
			std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator seenal = seen.find(al);
			if(seenal == seen.end()){
			//TODO add to adjucents
				std::cout << "get sc " <<std::endl;
				std::set< const pw_alignment*, compare_pointer_pw_alignment> sc;
				sc.insert(al);
				get_sc(*al, sc);
				alind.erase(al);
				std::cout << "sc size " << sc.size() << std::endl;
				if(sc.size()>1){
			//	strongly_connected.push_back(sc);
				}
				for(std::set<const pw_alignment*, compare_pointer_pw_alignment>::iterator it = sc.begin() ; it!= sc.end(); it++){
					std::map<const pw_alignment*, size_t>::iterator it1 = als_index.find(*it);
					assert(it1 != als_index.end());
					std::cout << it1->second<< " ";
				}
				std::cout << " "<<std::endl;
			}
		}
	//	strongly_connected.push_back(sc);//TODO add to adjucents
		for(std::multimap<size_t, size_t>::iterator it=edges.begin(); it != edges.end();it++){
			std::cout << it->first << " " << it->second << std::endl;
		}
	}

	void detect_overlap_rate::get_sc(const pw_alignment & al, std::set <const pw_alignment*, compare_pointer_pw_alignment> & sc ){
		std::vector<size_t> left(2);
		std::vector<size_t> right(2);
		al.get_lr1(left.at(0), right.at(0));
		al.get_lr2(left.at(1), right.at(1));
		std::vector<size_t>reference(2);
		reference.at(0) = al.getreference1();
		reference.at(1) = al.getreference2();
		size_t length = al.alignment_length();
		size_t al_index;
		std::map<const pw_alignment*, size_t>::iterator it = als_index.find(&al);
		assert(it != als_index.end());
		al_index = it->second;
		std::cout << "for index "<< al_index << std::endl;
		for(size_t i =0; i < 2;i++){
			sc_step(length, reference.at(i), left.at(i), right.at(i), sc, al_index);	
		}	
	}
	void detect_overlap_rate::sc_step(size_t length, size_t current_ref , size_t left, size_t right, std::set <const pw_alignment*, compare_pointer_pw_alignment> & sc, size_t & al_index){

		std::vector<const pw_alignment*> results;
		std::cout << "ref "<< current_ref << " left "<< left << " right "<< right << std::endl;
		alind.search_overlap(current_ref,left,right,results);
		std::cout<<"here!"<<std::endl;
		std::vector<const pw_alignment*> strong_als;
		//Check for 95% redundancy
		if(results.size() != 0){
			std::cout<< "here1"<<std::endl;
			check_the_strength(results, length, current_ref,left,right,strong_als, al_index);
			std::cout << "strong al size "<< strong_als.size()<<std::endl;
			for(size_t i=0; i<strong_als.size(); ++i) {
				sc.insert(strong_als.at(i));
			}
			for(size_t i =0; i < strong_als.size();i++){
				const pw_alignment* p = strong_als.at(i);
				std::map<const pw_alignment* , size_t>::iterator it=als_index.find(p);
				assert(it != als_index.end());
				size_t l = p->alignment_length();
				size_t l1,l2,r1,r2;
				p->get_lr1(l1,r1);
				p->get_lr2(l2,r2);
				size_t ref1 = p->getreference1();
				size_t ref2 = p->getreference2();
				sc_step(l, ref1, l1, r1,sc, it->second);
				sc_step(l, ref2, l2, r2,sc, it->second);
				seen.insert(p);
			}//XXX Does it make sense if i remove strong als from tree here? 
			for(size_t i =0; i < strong_als.size();i++){
				const pw_alignment* p = strong_als.at(i);
				alind.erase(p);
			}

			
		}
	}

	void detect_overlap_rate::check_the_strength(std::vector<const pw_alignment*> & results , size_t & length, size_t & ref, size_t & left, size_t & right, std::vector<const pw_alignment*> & strong_als, size_t & node_index){
			for(size_t i =0; i < results.size();i++){
			const pw_alignment * al = results.at(i);
			std::map<const pw_alignment*, size_t>::iterator it = als_index.find(al);
			std::cout << "index " << it->second << " parent index "<< node_index <<std::endl;
			assert(it != als_index.end());
			size_t l1,l2,r1,r2;
			al->get_lr1(l1,r1);
			al->get_lr2(l2,r2);
			if(it->second > node_index && al->getreference1()== ref && l1<=right && right<= r1){
				double overlap = (right-l1)/double((length+al->alignment_length()));
				if(overlap>0.95){
					strong_als.push_back(al);
					edges.insert(std::make_pair(node_index, it->second));
					edges.insert(std::make_pair(it->second,node_index));
				}
			}
			else if(it->second > node_index && al->getreference1()==ref && l1<=left && left<=r1){
				double overlap = (r1-left)/double((length+al->alignment_length()));
				if(overlap>0.95){
					strong_als.push_back(al);
					edges.insert(std::make_pair(node_index, it->second));
					edges.insert(std::make_pair(it->second,node_index));
				}
			}
			else if(it->second > node_index && al->getreference2()==ref && l2<=right && right<= r2){
				double overlap = (right-l2)/double((length+al->alignment_length()));
				if(overlap>0.95){
					strong_als.push_back(al);
					edges.insert(std::make_pair(node_index, it->second));
					edges.insert(std::make_pair(it->second,node_index));
				}
			}
			else if(it->second > node_index && al->getreference2()==ref && left<=r2 && right >=r2){
				double overlap = (r2-left)/double((length+al->alignment_length()));
				if(overlap>0.95){
					strong_als.push_back(al);
					edges.insert(std::make_pair(node_index, it->second));
					edges.insert(std::make_pair(it->second,node_index));

				}
			}	
		}
	}
	size_t detect_overlap_rate::get_id(const pw_alignment* p)const{
		std::map<const pw_alignment* , size_t>::const_iterator it = als_index.find(p);
		assert(it != als_index.end());
		return it->second;
	}

	void biconnected_component::creat_component(){
		for(size_t i =0; i < alignments.size();i++){
			visited.at(i) = false;
			get_articulation_points(i);
		}


	}
	void biconnected_component::get_articulation_points(size_t node){
		visited.at(node) = true;
		Depth.at(node) = depth;
		low.at(node) = depth;
		size_t childCount = 0;
		bool isArticulation = false;
		std::vector<size_t> adj;
	//	std::pair<std::multimap<size_t, size_t>::iterator ,std::multimap<size_t,size_t>::iterator > it = adjucent.equal_range(node);
	//	for(std::multimap<size_t, size_t>::iterator it1 = it.first; it1!=it.second; ++it1) {
	//		adj.push_back(it1->second);
	//	}
		for (size_t ni = 0; ni < adj.size(); ni++){//each ni in adj[node]
		        if(visited.at(ni) == false){
				parent.at(ni) = node;
				depth++;
				get_articulation_points(ni);
				childCount = childCount + 1;
				if(low.at(ni) >= Depth.at(node)){
          				isArticulation = true;
				        low.at(node) = std::min(low.at(node), low.at(ni));
	   			}else if( ni != parent.at(node)){
	      				low.at(node) = std::min(low.at(node), Depth.at(ni));
				}
			}
		}
    		if( (parent.at(node) != NULL &&  isArticulation) || (parent.at(node) == NULL && childCount > 1) ){
			std::cout << "n is an articulation point "<< std::endl;
		}
	
	}
