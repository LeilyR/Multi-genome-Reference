#include "suffix_tree.hpp"


	suffixTree::suffixTree(size_t & num_seq, finding_centers & cent, std::vector<std::vector<std::vector<size_t> > > & words):centers(cent){
		this->num_seq = num_seq;
		this->words =words;
		last_node_index=0;
		edges.insert(make_pair(last_node_index,root));
		last_node_index ++;
	}

	suffixTree::~suffixTree(){}

	void suffixTree::add_leaves(size_t & parent , std::vector<size_t> & new_suffix){
	//	double inf = std::numeric_limits<double>::infinity();
		std::cout<< "last node index "<< last_node_index << " new suffix size "<< new_suffix.size() << std::endl;
		edges_relations.insert(std::make_pair(parent, last_node_index));
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index++;
		std::cout << "relation size "<<edges_relations.size()<< std::endl;
	}
	void suffixTree::add_internal_node(size_t & parent_index , std::vector<size_t> & new_suffix){//equal range for parent
		std::multimap<size_t,size_t>::iterator it2 = edges_relations.find(parent_index);
		assert(it2 != edges_relations.end());
		std::cout << last_node_index << " " << parent_index <<std::endl;
		std::vector<size_t> intermediate;
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			intermediate.push_back(it1->second);
		}
		edges_relations.erase(parent_index);
		it2 = edges_relations.find(parent_index);
		assert(it2 == edges_relations.end());
		std::cout << "child nodes were deleted"<<std::endl;
		edges_relations.insert(std::make_pair(parent_index, last_node_index));
		std::cout<< "intermediate size "<< intermediate.size()<<std::endl;
		for(size_t i =0; i < intermediate.size();i++){
			edges_relations.insert(std::make_pair(last_node_index, intermediate.at(i)));						
			std::cout << "test at "<< i << " is " << intermediate.at(i)<<std::endl;
		}
		edges.insert(std::make_pair(last_node_index,new_suffix));
		last_node_index ++;

	}
	void suffixTree::update_node(size_t & node_index, std::vector<size_t> & new_content){
		std::cout<< "node index "<< node_index << " new content "<< new_content << std::endl;
		std::map<size_t, std::vector<size_t> >::iterator it = edges.find(node_index);
		assert(it != edges.end());
		it->second = new_content;

	}
	void suffixTree::split(size_t & parent_index, std::vector<size_t> & parent , std::vector<size_t> & new_suffix, std::vector<size_t> & common_part, bool & make_a_tree){
		assert(common_part.size() <= parent.size());
		std::vector<size_t> parent_non_common;
		std::cout << "common part size "<< common_part.size() << " parent size "<< parent.size()<<std::endl;
		for(size_t i =common_part.size(); i < parent.size();i++){
			parent_non_common.push_back(parent.at(i));
		}
		if(parent_non_common.size() != 0){//common part shorter than parent's length
			std::cout << "if parent_non_common.size() != 0 "<<std::endl;
			std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
			if(it == edges_relations.end()){
				std::cout<< "add leaves "<<std::endl;
				add_leaves(parent_index,parent_non_common); 
			}else{
				std::cout<< "add internal "<<std::endl;
				add_internal_node(parent_index,parent_non_common);
			}
			std::cout<< "parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() << " " << parent_non_common <<std::endl;
			update_node(parent_index,common_part);//update its content
			std::cout<< "1parent index is "<< parent_index<< " " << new_suffix.size() << " >= " << common_part.size() <<std::endl;
			if(new_suffix.size()>= common_part.size()){
				std::vector<size_t> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
					std::cout << new_suffix.at(i) << " " ;
					current_non_common.push_back(new_suffix.at(i));
				}
				std::cout<< "current non common size is "<< current_non_common <<std::endl;
				if(current_non_common.size() != 0){
					std::cout<< "add leaves1 "<<std::endl;
					add_leaves(parent_index,current_non_common); 
				}
			}
			make_a_tree = false;
		}else{//common part is equal to parent's length - the new suffix is longer than the parent and fully covers it.
			std::cout << "i am here!" <<std::endl;
			assert(parent_non_common.size() == 0);
			if(common_part.size() < new_suffix.size()){ //TRICKY ONE!!!!
				std::vector<size_t> current_non_common;
				for(size_t i =common_part.size(); i < new_suffix.size();i++){
					std::cout<< new_suffix.at(i) << " ";
					current_non_common.push_back(new_suffix.at(i));
				}
				std::cout << " "<<std::endl;
				assert(current_non_common.size()!=0);
				std::multimap<size_t, size_t>::iterator it = edges_relations.find(parent_index);
				if(it != edges_relations.end()){//check on the childrens for a common part
					std::cout<< "has child node"<<std::endl;
					size_t new_parent_index = edges.size();
					find_next_parent(parent_index, current_non_common, common_part, new_parent_index, parent);
					std::cout << new_parent_index <<std::endl;
					if(new_parent_index != edges.size()){
						new_suffix = current_non_common;
						parent_index = new_parent_index;
						std::cout << "here "<< std::endl;
					}else{
						add_leaves(parent_index,current_non_common);
						make_a_tree = false;
					}
				}else{//add it as a leave
					assert(it == edges_relations.end());
					add_leaves(parent_index,current_non_common);
					make_a_tree = false;
				}

			}else{
				assert(common_part.size() == new_suffix.size());
				make_a_tree = false;
			}

		}
	}
	void suffixTree::make_suffix(std::vector<size_t> & word){
		suffixes.clear();
		for(size_t i =0; i < word.size();i++){
			std::cout << word.at(i) << " ";
			std::vector<size_t> suffix;
			for(size_t k =i; k < word.size();k++){
				suffix.push_back(word.at(k));
			}//For convenience add 1<<31 at the end of all the suffixes
			suffix.push_back(1<<31);
			suffixes.push_back(suffix);
		} 		std::cout << " "<<std::endl;	
	
	}
	void suffixTree::find_first_edges_after_root(size_t & deepest_parent, std::vector<size_t> & new_suffix){
		std::map<size_t, size_t >::iterator it = first_edges.find(new_suffix.at(0));
		if(it!= first_edges.end()){
			deepest_parent = it->second;
		}
	
	}
	void suffixTree::find_next_parent(size_t & parent_index, std::vector<size_t> & current_non_common, std::vector<size_t> & common_part, size_t & new_parent_index, std::vector<size_t> & parent){
		common_part.clear();
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t,size_t>::iterator it1 = it.first ; it1 != it.second; it1++){
			std::cout << "child "<< it1->second << std::endl;
			std::map<size_t, std::vector<size_t> >::iterator edge=edges.find(it1->second);
			assert(edge != edges.end());
			std::cout << edge->second.at(0)<< " "<<current_non_common.at(0)<<std::endl;
			if(edge->second.at(0)==current_non_common.at(0)){
				std::cout<< "got it"<< std::endl;
				find_common_part(edge->second, current_non_common,common_part);
				new_parent_index = edge->first;
				parent = edge->second;
				std::cout << "new parent index "<< new_parent_index << std::endl;
				break;
			}
		}


	}
	void suffixTree::find_children(std::vector<size_t> & parent ,std::vector<std::vector<size_t> >& children){


	}
	void suffixTree::find_common_part(std::vector<size_t> & parent , std::vector<size_t> & current, std::vector<size_t> & common_part){
		size_t length = current.size();
		if(length > parent.size()) length = parent.size();
		for(size_t i =0; i < length; i++){
			if(parent.at(i)==current.at(i)){
				common_part.push_back(parent.at(i));
			}else{
				break;
			}
		}
	}
	void suffixTree::make_tree_for_a_seq(size_t & seq){
		std::cout << "at seq "<< seq <<std::endl;
		for(size_t i = 0; i < words.at(seq).size(); i ++){
			std::vector<size_t> this_word = words.at(seq).at(i);//current string of centers on a seq
			make_suffix(this_word);
			std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
			for(size_t j =0; j< suffixes.size(); j++){
				size_t deepest_parent =edges.size();
				find_first_edges_after_root(deepest_parent, suffixes.at(j));
				if(deepest_parent == edges.size()){
					std::cout << "no first parent! "<<std::endl;
					first_edges.insert(std::make_pair(suffixes.at(j).at(0),last_node_index));
					size_t root_index =0;
					add_leaves(root_index,suffixes.at(j));
				}else{
					std::map<size_t, std::vector<size_t> >::iterator it = edges.find(deepest_parent);
					assert(it != edges.end());
					std::cout<< "has a parent"<<std::endl;
					std::vector<size_t> common_part;
					std::vector<size_t> parent = it->second;
					std::vector<size_t> current = suffixes.at(j);
					find_common_part(parent, current,common_part);
					for(size_t n =0; n < common_part.size();n++){
						std::cout<< common_part.at(n) << " ";
					}
					std::cout << " "<<std::endl;
					bool make_a_tree = true;
					size_t parent_index = deepest_parent;
					while(make_a_tree == true){
						split(parent_index,parent,current, common_part, make_a_tree);
					}
				}
			}
		}
//		print_tree();

	}
	void suffixTree::make_tree(std::vector<size_t> & center_with_highest_gain, size_t & highest_index){	
		std::cout << "words size is "<< words.size() <<std::endl;
		for(size_t seq_id = 0; seq_id < num_seq; seq_id++){
			if(center_with_highest_gain.size() == 0){
				std::cout<< "seq id "<<seq_id <<std::endl;
			//	words.at(seq_id) = centers.get_centers(seq_id);
				std::cout << "its words size is "<< words.at(seq_id).size() <<std::endl;
			}else{
				update_centers(seq_id, center_with_highest_gain, highest_index);
			}
			make_tree_for_a_seq(seq_id);
		}
		count_branches();
	}
	void suffixTree::update_centers(size_t & seq_id, std::vector<size_t> & edge_with_highest_gain, size_t & highest_index){
		for(size_t i = 0; i < words.at(seq_id).size(); i++){
			std::vector<size_t> centers = words.at(seq_id).at(i);
			for(size_t k =0; k < centers.size();k++){
				if(centers.at(k) == edge_with_highest_gain.at(0)&& ((centers.size()-k)>= edge_with_highest_gain.size())){
					size_t first_common_index;
					std::vector<size_t> common;
					for(size_t i = 0; i < edge_with_highest_gain.size();i++){
						first_common_index = k;
						if(centers.at(i+k)== edge_with_highest_gain.at(i)){
							common.push_back( edge_with_highest_gain.at(i));
						}else break;
					}
					if(common.size() == edge_with_highest_gain.size()){
						std::vector<size_t> new_center;
						for(size_t m = 0; m < first_common_index ;m++){
							new_center.push_back(centers.at(m));
						}
						new_center.push_back(highest_index);
						for(size_t m = first_common_index + common.size(); m < centers.size();m++){
							new_center.push_back(centers.at(m));
						}
						centers=new_center;
					}
				}
			}
			words.at(seq_id).at(i) = centers;
		}

	}
	void suffixTree::count_branches(){
		for(size_t seq =0; seq < num_seq; seq++){
			for(size_t w = 0; w< words.at(seq).size(); w++){
				std::vector<size_t> this_word = words.at(seq).at(w);//current string of centers on a seq
				make_suffix(this_word);
				std::cout<< "suffixes size "<<suffixes.size()<<std::endl;
				for(size_t j =0; j< suffixes.size(); j++){
					std::cout << "j "<<j <<std::endl;
					std::vector<size_t> branch;
					std::vector<size_t> current = suffixes.at(j);
					size_t deepest_parent =edges.size();
					find_first_edges_after_root(deepest_parent,current);
					assert(deepest_parent != edges.size());
					std::map<size_t, std::vector<size_t> >::iterator it=edges.find(deepest_parent);
					assert(it != edges.end());
					if(it->second.size()==current.size()){
						branch.push_back(deepest_parent);
						std::cout << "all the suffix is on one edge!"<<std::endl;
					}else{//Note that first parent can not be longer than current.
						std::vector<size_t> updated_current;
						std::cout<<"updated current: "<<std::endl;
						for(size_t k =it->second.size() ; k < current.size(); k++){
							std::cout << current.at(k)<< " ";
							updated_current.push_back(current.at(k));
						}
						std::cout << " "<<std::endl;
						size_t parent_index = deepest_parent;
						branch.push_back(parent_index);
						current = updated_current;
						while(current.size()>0){
							traverse_the_tree(parent_index,current);
							std::cout << "current size "<< current.size() <<std::endl;
							branch.push_back(parent_index);							
						}
					}
					std::vector<size_t> sub_branch;
					std::cout<<"branch size "<< branch.size()<<std::endl;
					for(size_t i = 0 ; i < branch.size(); i++){
						std::cout << "i "<< i << " " << branch.at(i) <<std::endl;
						sub_branch.push_back(branch.at(i));
						std::map<std::vector<size_t>,size_t>::iterator it1 = branches.find(sub_branch);
						if(it1 == branches.end()){
							branches.insert(make_pair(sub_branch,0));
							it1 = branches.find(sub_branch);
						}
						it1->second = it1->second +1;
					}

				}
			}
		}
	}
	void suffixTree::traverse_the_tree(size_t & parent_index, std::vector<size_t> & new_suffix){
		std::cout<< "parent index is "<< parent_index <<std::endl;
		std::multimap<size_t, size_t>::iterator check = edges_relations.find(parent_index);
		assert(check != edges_relations.end());
		std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator > it = edges_relations.equal_range(parent_index);
		for(std::multimap<size_t, size_t>::iterator it1 = it.first; it1 !=it.second;it1++){
			std::map<size_t,std::vector<size_t> >::iterator it2 = edges.find(it1->second);
			std::cout<< it1->second <<std::endl;
			assert(it2 != edges.end());
			if(it2->second.at(0)== new_suffix.at(0)){
				parent_index = it1->second;
				std::cout<< "new parent index " << parent_index<<std::endl;
				std::vector<size_t> updated_current;
				std::cout<< "updated:"<<std::endl;
				for(size_t i = it2->second.size(); i < new_suffix.size(); i++){
					std::cout<< new_suffix.at(i)<< " ";
					updated_current.push_back(new_suffix.at(i));
				}		
				std::cout << " "<<std::endl;
				new_suffix = updated_current;
				break;
			}
		}
		


	}
	const std::map<size_t, std::vector<size_t> > suffixTree::get_edges()const{
		return edges;
	}
	const std::map<std::vector<size_t>, size_t > suffixTree::get_branches()const{
		return branches;
	}

	const std::vector<std::vector<size_t> > suffixTree::get_current_centers(size_t & id)const{
		return words.at(id);
	}

	void suffixTree::print_tree(){
		for(std::map<size_t, std::vector<size_t> >::iterator it = edges.begin();it !=edges.end(); it++){
			std::cout<< " node is "<<it->first << " its edge content is ";
			for(size_t i = 0 ; i < it->second.size();i++){
				std::cout << it->second.at(i)<< " ";
			}
			std::cout<< " "<<std::endl;
			std::cout << "its children are "<<std::endl;
			std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> it2 = edges_relations.equal_range(it->first);
			for(std::multimap<size_t,size_t>::iterator it1 = it2.first ; it1 != it2.second; it1++){
				std::cout << it1->second << " ";
			}
			std::cout << " " <<std::endl;
		}
	}



