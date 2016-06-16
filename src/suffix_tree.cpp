#include "suffix_tree.hpp"

//TODO use typename , then one can make a tree using any type of input

	suffixTree::suffixTree(){
		last_node_index=0;
	}

	suffixTree::~suffixTree(){}

	void suffixTree::add(std::vector<size_t> & parent , std::vector<size_t> & new_suffix){


	}
	void suffixTree::split(std::vector<size_t> & parent , std::vector<size_t> & new_suffix, std::vector<size_t> & common_part){
		if(parent.size() >= common_part.size()){
			std::vector<size_t> n_common;
			for(size_t i =commn_part.size(); i < parent.size();i++){
				n_common.push_back(parent.at(i));
			}
			//TODO add it to the map !
		}else{

			//TODO need to check for its child node and update the parent
		}
		if(new_suffix.size()> common_part.size()){
			std::vector<size_t> n_common;
			for(size_t i =commn_part.size(); i < new_suffix.size();i++){
				n_common.push_back(new_suffix.at(i));
			}
			//TODO need to check for its child node and update the parent
			std::vector<std::vector<size_t> > children;
			find_children(parent, children);
			if(children.size()!=0){ //If the current parent has any child node
			//TODO need to check for its child node and update the parent
				

			}else{//If no child node //TODO update 

			}
		}else{
			assert(common_part.size()==new_suffix.size());
		}
	}
	void suffixTree::make_suffix(std::vector<size_t> & word){
		for(size_t i =0; i < word.size();i++){
			std::vector<size_t> suffix;
			for(size_t k =i; k < word.size();k++){
				suffix.push_back()word.at(k);
			}
			suffixs.push_back(suffix);
		} //Add 1<<31 if the last character happens before	
		size_t last_center = word.at(word.size()-1);
		for(size_t j =0; j < word.size()-1; j++){
			if(word.at(j)==last_center){
				for(size_t m =0; m < suffixes.size();m++){
					std::vector<size_t> new_suffix = suffixes.at(m);
					new_suffix.push_back(1<<31);
					suffixes.at(m)= new_suffix;
				}
				break;
			}
		}		
	}
	void suffixTree::find_first_edges_after_root(std::vector<size_t> & deepest_parent, std::vector<size_t> & new_suffix){
		std::map<size_t, std::vector<size_t> >::iterator it = first_edges.find(new_suffix.at(0));
		if(it!= first_edges.end()){
			deepest_parent = it->second;
		}
	
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
	void suffixTree::make_tree(size_t & seq_id, bool highest_gain){
		words = centers.get_centers(seq_id,highest_gain);
		for(size_t i = 0; i < words.size(); i ++){
			std::vecotr<size_t> this_word = words.at(i);//current string of centers on a seq
				make_suffix(this_word);
				for(size_t j =0; j< suffixes.size(); j++){
					vector<size_t> deepest_parent;
					find_first_edges_after_root(deepest_parent, suffixes.at(j));
					if(deepest_parent.size()==0){
						add(deepest_parent,suffixes.at(j));
					}else{
						std::vector<size_t> common_part;
						find_common_part(deepest_parent, suffixes.at(j),common_part);
						std::vector<size_t> parent = deepest_parent;
						while(common_part.size()>0){
							split(parent, common_part);

						}


					}
				}

		
			}
		}
	}



