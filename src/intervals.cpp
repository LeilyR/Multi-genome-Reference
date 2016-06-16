#include "intervals.hpp"

template<typename T>
it_node<T>::it_node(): height(NULL), str_left(NULL), str_right(NULL), str_center(NULL), parent(NULL), min_left(0), max_right(0) {
	assert(0);

}

template<typename T>
it_node<T>::it_node(const size_t & left, const size_t & right, const value_type & v): height(NULL), left(left), right(right), value(v), str_left(NULL), str_right(NULL), str_center(NULL), parent(NULL), min_left(left), max_right(right) { }

template<typename T>
it_node<T>::~it_node() {}








template<typename T>
avl_interval_tree<T>::avl_interval_tree(): root(NULL), num_nodes(0), num_reba(0), num_moves(0)
	 {
}

template<typename T>
avl_interval_tree<T>::~avl_interval_tree() {
	if(root!=NULL) {
		node_type * root_node = root;
		prune_subtree(root);
		delete_subtree(root_node);
	}


}



template<typename T>
void avl_interval_tree<T>::insert(const size_t & left, const size_t & right, const value_type  & v) {
#if DO_SLOW_CHECKS
	assert(slow_index.size() == num_nodes);
#endif
#if DEBUG_PRINTS
	std::cout << " insert l " << left << " r " << right << " d " << v << std::endl;
#endif


	node_type * new_node = new node_type(left, right, v);
	num_nodes++;
	if(root == NULL) {
		root = new_node;
	} else {
		insert_at_node(root, new_node);
#if DO_REBALANCING
		balance_up(new_node);
#endif
	}


#if DO_SLOW_CHECKS
	slow_index.push_back(std::make_pair(std::make_pair(left, right), v));


	assert(slow_index.size() == num_nodes);
#if DEBUG_PRINTS
	std::cout << " After insert check " << std::endl;
#endif
	check_tree();  
#endif
}


template<typename T>
void avl_interval_tree<T>::bulk_insert(std::vector< std::pair < std::pair< size_t, size_t > , value_type & > > & data ) {
/*	multimap< size_t, std::pair < std::pair< size_t, size_t > , value_type & v > > lsorter;
	for(size_t i=0; i<data.size(); ++i) {
		lsorter.insert( std::make_pair( data.at(i).first.first, data.at(i) ) );
	}
	std::vector< std::pair < std::pair< size_t, size_t > , value_type & v > > lsorted_data(lsorter.size());
	size_t num = 0;
	for(multimap< size_t, std::pair < std::pair< size_t, size_t > , value_type & v > >::iterator it = lsorter.begin(); it!=lsorter.end(); ++it) {
		lsorted_data.at(num) = it->second;
		num++;
	}


	std::vector< std::vector< std::pair< size_t, size_t > , value_type & v >  > > level_data; // level , index -> input data element
*/
}

template<typename T>
void avl_interval_tree<T>::erase_all(const size_t & left, const size_t & right, std::vector<value_type> & erased) {

#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " Before erase all check: " << std::endl;
	std::cout << " try to erase: l " << left << " r " << right << std::endl;
#endif
	check_tree();
	assert(num_nodes == slow_index.size());

	
	size_t slow_erased = 0;
	std::set<size_t> to_del;
	for(size_t i=0; i<slow_index.size(); ++i) {
		size_t tl = slow_index.at(i).first.first;
		size_t tr = slow_index.at(i).first.second;
		if((tl == left && tr == right)) {
			to_del.insert(i);
			slow_erased++;
		}
	}
	if(!to_del.empty()) {
		std::vector< std::pair<std::pair<size_t, size_t>, value_type> > new_index(slow_index.size()-to_del.size());
		size_t newat = 0;
		for(size_t i=0; i<slow_index.size(); ++i) {
			std::set<size_t>::iterator findi = to_del.find(i);
			if(findi==to_del.end()) {
				new_index.at(newat) = slow_index.at(i);
				newat++;
			}
		}
		slow_index = new_index;
	}

#endif



	assert(erased.empty());
	if(root!=NULL) {
		erase_at_node(root, left, right, erased);
		num_nodes -= erased.size();
	} 
#if DO_REBALANCING
	set_rebalance();
#endif 

#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " erase all result: " << erased.size() << " nodes were erased" << std::endl;
#endif
	assert(slow_index.size() == num_nodes);
	assert(slow_erased == erased.size());
#if DEBUG_PRINTS
	std::cout << " after erase all check: " << std::endl;
#endif
	check_tree();

#endif

}


template<typename T>
size_t avl_interval_tree<T>::erase(const size_t & left, const size_t & right, const value_type & value) {

#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " Before single single erase check: " << std::endl;
	std::cout << " try to erase: l " << left << " r " << right << " d " << value << std::endl;
#endif
	check_tree();
	assert(num_nodes == slow_index.size());

	
	size_t slow_erased = 0;
	std::set<size_t> to_del;
	for(size_t i=0; i<slow_index.size(); ++i) {
		size_t tl = slow_index.at(i).first.first;
		size_t tr = slow_index.at(i).first.second;
		value_type tv = slow_index.at(i).second;
		if((tl == left && tr == right) && (tv == value)) {
//			std::cout << "tl "<< tl << " tr "<< tr <<std::endl;
			to_del.insert(i);
			slow_erased++;
		}
	}
//	std::cout<< "slow_erased "<< slow_erased << std::endl;
	if(!to_del.empty()) {
		std::vector< std::pair<std::pair<size_t, size_t>, value_type> > new_index(slow_index.size()-to_del.size());
		size_t newat = 0;
		for(size_t i=0; i<slow_index.size(); ++i) {
			std::set<size_t>::iterator findi = to_del.find(i);
			if(findi==to_del.end()) {
				new_index.at(newat) = slow_index.at(i);
				newat++;
			}
		}
		slow_index = new_index;
#if DEBUG_PRINTS
		std::cout << " after deleting " << to_del.size() << " slow index size is " << slow_index.size() << std::endl;
#endif
	}

#endif



	size_t num_erased = 0;
	if(root!=NULL) {
//		std::cout<< "before single erase: "<< left << " "<< right <<std::endl;
		num_erased = single_erase_at_node(root, left, right, value);
		num_nodes -= num_erased; 
	} 
#if DO_REBALANCING
	set_rebalance();
#endif

#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " erase result: " << num_erased << " nodes were erased" << std::endl;
	std::cout << " slow results: " << to_del.size() << " erased, now we have " << slow_index.size() << std::endl;
	std::cout<< "num nodes "<< num_nodes <<std::endl;
#endif
	assert(slow_index.size() == num_nodes);
	assert(slow_erased == num_erased);
#if DEBUG_PRINTS
	std::cout << " after single erase check: " << std::endl;
#endif
	check_tree();

#endif
	return num_erased;
}


template<typename T>
size_t avl_interval_tree<T>::erase_vector(const std::vector< std::pair< std::pair<size_t, size_t> , value_type > > & todel) {
#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " Before erase vector check: " << std::endl;
	for(size_t i=0; i<todel.size(); ++i) {
		size_t left = todel.at(i).first.first;
		size_t right = todel.at(i).first.second;
		value_type value = todel.at(i).second;
		std::cout << " try to erase: l " << left << " r " << right << " d " << value << std::endl;
	}
#endif
	check_tree();
	assert(num_nodes == slow_index.size());

	size_t slow_erased = 0;
	std::set<size_t> to_del;
	for(size_t j=0; j<todel.size(); ++j) {
		size_t left = todel.at(j).first.first;
		size_t right = todel.at(j).first.second;
		value_type value = todel.at(j).second;
		for(size_t i=0; i<slow_index.size(); ++i) {
			size_t tl = slow_index.at(i).first.first;
			size_t tr = slow_index.at(i).first.second;
			value_type tv = slow_index.at(i).second;
			if((tl == left && tr == right) && (tv == value)) {
//			std::cout << "tl "<< tl << " tr "<< tr <<std::endl;
				to_del.insert(i);
				slow_erased++;
			}
		}
	}

//	std::cout<< "slow_erased "<< slow_erased << std::endl;
	if(!to_del.empty()) {
		std::vector< std::pair<std::pair<size_t, size_t>, value_type> > new_index(slow_index.size()-to_del.size());
		size_t newat = 0;
		for(size_t i=0; i<slow_index.size(); ++i) {
			std::set<size_t>::iterator findi = to_del.find(i);
			if(findi==to_del.end()) {
				new_index.at(newat) = slow_index.at(i);
				newat++;
			}
		}
		slow_index = new_index;
#if DEBUG_PRINTS
		std::cout << " after deleting " << to_del.size() << " slow index size is " << slow_index.size() << std::endl;
#endif
	}

#endif


	node_set_type todel_nodes;
	for(size_t i=0; i<todel.size(); ++i) {
		size_t l = todel.at(i).first.first;
		size_t r = todel.at(i).first.second;
		value_type v = todel.at(i).second;
		if(root!=NULL) {
			find_at_node(root, l, r, v, todel_nodes);
		}
	}
	delete_nodes(todel_nodes);
	num_nodes-=todel_nodes.size();
	
#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " erase result: " << todel_nodes.size() << " nodes were erased" << std::endl;
	std::cout << " slow results: " << to_del.size() << " erased, now we have " << slow_index.size() << std::endl;
	std::cout<< "num nodes "<< num_nodes <<std::endl;
#endif
	assert(slow_index.size() == num_nodes);
	assert(slow_erased == todel_nodes.size());
#if DEBUG_PRINTS
	std::cout << " after erase vector check: " << std::endl;
#endif
	check_tree();

#endif
	return todel_nodes.size(); 

}



template<typename T>	
void avl_interval_tree<T>::overlap_search(const size_t & left, const size_t & right, std::vector<value_type> & results) const {
#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " Before overlap search check: " << std::endl;
	std::cout << " check for overlap with l " << left << " r " << right <<  std::endl;
#endif
	check_tree();
#endif


	if(root!=NULL) {
		overlap_search_at_node(root, left, right, results);
	}



#if DO_SLOW_CHECKS
	assert(slow_index.size() == num_nodes);
#if DEBUG_PRINTS
	std::cout << "We have found " << results.size() << " overlapping intervals using the tree "<<std::endl;
#endif
//	for(size_t i=0; i<results.size(); ++i) {
//		std::cout << " v " << results.at(i) << std::endl;
//	}
	std::vector<value_type> slow_results;
	for(size_t i=0; i<slow_index.size(); ++i) {
		size_t tl = slow_index.at(i).first.first;
		size_t tr = slow_index.at(i).first.second;
		value_type tv = slow_index.at(i).second;
		if(tl <= right && tr >= left) {
//			std::cout << " slow found  l " << tl << " r " << tr << " d " << tv << " border " << left<< " r " << right << std::endl;
			slow_results.push_back(tv);
		} else {
//			std::cout << " slow no ovrlp  l " << tl << " r " << tr << " d " << tv << std::endl;

		}
	}
#if DEBUG_PRINTS
	std::cout << "We have found " << slow_results.size() << " overlapping intervals using a simple method " << std::endl;
	std::cout << "result size is "<< results.size() <<std::endl;
#endif
	assert(slow_results.size() == results.size());
	for(size_t i=0; i<slow_results.size(); ++i) {
//		std::cout << " slow res " << slow_results.at(i) << std::endl;
		bool found = false;
		for(size_t j=0; j<results.size(); ++j) {
			if(results.at(j) == slow_results.at(i)) {
				found = true;
				break;
			}
		}
		assert(found);
	}
#endif

}

template<typename T>	
void avl_interval_tree<T>::overlap_search_erase(const size_t & left, const size_t & right, std::vector<value_type> & results) {
#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " Before overlap search erase check: " << std::endl;
	std::cout << " check for overlap with l " << left << " r " << right <<  std::endl;
#endif
	check_tree();
	assert(slow_index.size() == num_nodes);

	std::vector< std::pair<std::pair<size_t, size_t>, value_type> > slow_index_backup = slow_index;
#endif

	if(root!=NULL) {
		node_map_type roots_to_parents;
		node_set_type valid_parents;
		node_set_type to_Delete;
		overlap_search_erase_at_node(root, left, right, results, to_Delete);
		delete_nodes(to_Delete);	
		num_nodes -= to_Delete.size();
		

	}

#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << "We have found " << results.size() << " overlapping intervals using the tree "<<std::endl;
#endif
	for(size_t i=0; i<results.size(); ++i) {
//		std::cout << " v " << results.at(i) << std::endl;
	}
	std::vector<value_type> slow_results;
	slow_index.clear();
	for(size_t i=0; i<slow_index_backup.size(); ++i) {
		size_t tl = slow_index_backup.at(i).first.first;
		size_t tr = slow_index_backup.at(i).first.second;
		value_type tv = slow_index_backup.at(i).second;
		if(tl <= right && tr >= left) {
//			std::cout << " slow found  l " << tl << " r " << tr << " d " << tv << std::endl;
			slow_results.push_back(tv);
		} else {
			slow_index.push_back(slow_index_backup.at(i));
//			std::cout << " slow no ovrlp  l " << tl << " r " << tr << " d " << tv << std::endl;

		}
	}
#if DEBUG_PRINTS
	std::cout << "We have found " << slow_results.size() << " overlapping intervals using a simple method " << std::endl;
#endif
	assert(slow_results.size() == results.size());
	for(size_t i=0; i<slow_results.size(); ++i) {
		bool found = false;
		// check all slow results are in real results
		for(size_t j=0; j<results.size(); ++j) {
			if(results.at(j) == slow_results.at(i)) {
				found = true;
				break;
			}
		}
		assert(found);
	}
	assert(num_nodes == slow_index_backup.size() - slow_results.size());
#if DEBUG_PRINTS
	std::cout << "After overlap search erase check" <<std::endl;
#endif
	check_tree();
#endif




}

template<typename T>
void avl_interval_tree<T>::join_intervals(std::vector<std::pair<size_t, size_t> > & results) const {
	std::map<size_t, size_t> joined_intervals; // left -> right

	if(root!=NULL) {
		join_intervals_at_node(root, joined_intervals);
	}

	for(std::map<size_t, size_t>::iterator it = joined_intervals.begin(); it!=joined_intervals.end(); ++it) {
		results.push_back(*it);
	}

}


template<typename T>
void avl_interval_tree<T>::overlap_fraction_search(const size_t & left, const size_t & right, const double & fraction, std::vector<value_type> & results) const {
	if(root!=NULL) {
		overlap_fraction_search_at_node(root, left, right, fraction, results);
	}
}


template<typename T>
size_t avl_interval_tree<T>::size() const {
	return num_nodes;
}


template<typename T>
size_t avl_interval_tree<T>::height() const {
	if(root==NULL) {
		return 0;
	}
	return root->height;
}


template<typename T>
double avl_interval_tree<T>::list_ratio() const {
	if(root==NULL) {
		return 0;
	}
	size_t num_l = count_center_at_node(root);
	double ratio = (double) num_l / (double) num_nodes;
	return ratio;
}


template<typename T>
size_t avl_interval_tree<T>::get_num_reba() const {
	return num_reba;
}


template<typename T>
size_t avl_interval_tree<T>::get_num_moves() const {
	return num_moves;
}


template <typename T>	
void avl_interval_tree<T>::debug_print() const {
	std::cout << " AVL Interval tree printout: " << std::endl;
	if(root!=NULL) {
		size_t id = 0;
		debug_print_at_node(root, 0, id, 0);
	}

}


template<typename T>
it_node<T> * avl_interval_tree<T>::insert_at_node(node_type * atn, node_type * newn) {
	assert(newn!=NULL);
	assert(newn->str_left==NULL);
	assert(newn->str_right==NULL);
	assert(newn->str_center==NULL);
#if DEBUG_PRINTS
	std::cout << " insert l " << newn->left << " r " << newn->right << " at l " << atn->left << " r " << atn->right << std::endl;;
#endif 

#if DO_SLOW_CHECKS
	size_t minl, maxr, intv_l, intv_r;
	check_upwards(atn, minl, maxr, intv_l, intv_r);

	assert(minl <= newn->min_left);
	assert(maxr >= newn->max_right);
	assert(newn->left <= intv_r);
	assert(newn->right >= intv_l);

#endif


	if(atn==NULL) {
		if(root == NULL) {
			root = newn;
			return root;
		} else {
			insert_at_node(root, newn);
			return root;
		}
	}
	if(newn->right < atn->left) {
// insert into left subtree
		if(atn->str_left==NULL) {
			left_attach(atn, newn);	
			return atn;
		} else {
			insert_at_node(atn->str_left, newn);
			return atn;
		}
	} else if(newn->left > atn->right) {
// insert into right subtree
		if(atn->str_right==NULL) {
			right_attach(atn, newn);
			return atn;
		} else {
			insert_at_node(atn->str_right, newn);
			return atn;
		}
	} else {
// insert into center subtree

	// we swap newn to the position of atn if the new interval is smaller
			// putting smaller intervals at higher positions in the tree makes center subtrees (which cannot be well balanced) shorter


		if( ( (atn->left <= newn->left) && (atn->right >  newn->right) ) ||
		    ( (atn->left <  newn->left) && (atn->right >= newn->right) ) ) {
	
#if DEBUG_PRINTS
			std::cout << "INSERT SWITCH: INS l " << newn->left << " r " << newn->right << " switch with l " << atn->left << " r " << atn->right << std::endl; 
#endif
			node_type * atn_parent = atn->parent;
			if(atn_parent == NULL) {
#if DEBUG_PRINTS
				std::cout << " INSERT SWITCH at root " << std::endl;
#endif
				prune_subtree(atn);
				root = newn;
				move_subtree(atn, root);
				return root;
			} else {
#if DEBUG_PRINTS
				std::cout << " INSERT SWITCH below l " << atn_parent->left << " r " << atn_parent->right << std::endl;
#endif
				if(atn_parent->str_left == atn) {
					prune_subtree(atn);
					left_attach(atn_parent, newn);
					atn_parent = move_subtree(atn, atn_parent);
				} else if(atn_parent->str_center == atn) {
					prune_subtree(atn);
					center_attach(atn_parent, newn);
					atn_parent = move_subtree(atn, atn_parent); 
				} else if (atn_parent->str_right == atn) {
					prune_subtree(atn);
					right_attach(atn_parent, newn);
					atn_parent = move_subtree(atn, atn_parent);
				} 
				return atn_parent;
			}
		} else {
// normal center insert without switch
			if(atn->str_center==NULL) {
				center_attach(atn, newn);
				return atn;
			} else {
				insert_at_node(atn->str_center, newn);
				return atn;
			}
		}
	}	
	assert(0);
}

template<typename T>
void avl_interval_tree<T>::forest_add(node_type * atn, node_type * insert_node, node_map_type & roots_to_parents, node_set_type & valid_parents) {
	roots_to_parents.insert(std::make_pair(atn, insert_node));
	if(insert_node!= NULL){ 
		valid_parents.insert(insert_node);
	}
}

template<typename T>
void avl_interval_tree<T>::redo_height_up(node_type * n) {
	bool n_changed = true;
	while(n != NULL && n_changed){
		n_changed = redo_height_and_minmax_at_node(n);
		n=n->parent;
	}	
}

template<typename T>
void avl_interval_tree<T>::forest_insert(node_map_type & roots_to_parents, node_set_type & valid_parents) {
	for(typename node_map_type::iterator it = roots_to_parents.begin(); it!=roots_to_parents.end(); ++it) {
		node_type * subtree_root = it->first;
//		std::cout << "map size "<< roots_to_parents.size() << std::endl;
		node_type * ins_node = it->second;
		typename node_set_type::iterator find_ian = valid_parents.find(ins_node);
		if(find_ian==valid_parents.end()) {
	//		std::cout<< "ins node null"<<std::endl;
			ins_node = NULL;
		}
		if(ins_node == NULL) ins_node = root;
//		std::cout<< "move"<< subtree_root << " ins_node " << ins_node << std::endl;


		move_subtree(subtree_root, ins_node);
	}
// After all insertions we do some rebalancing (at each inserted subtree root) as all of them propagate upwards, the top part of the tree should be quite good again
#if DO_REBALANCING
	for(typename node_map_type::iterator it = roots_to_parents.begin(); it!=roots_to_parents.end(); ++it) {
		node_type * subtree_root = it->first;
		balance_up(subtree_root);
	}
#endif



}

template<typename T>
void avl_interval_tree<T>::forest_delete_root(node_type * n, node_set_type & valid_parents) {
	valid_parents.erase(n);
}



template<typename T>
void avl_interval_tree<T>::delete_nodes(const node_set_type & to_Delete) {

	node_map_type roots_to_parents;
	node_set_type valid_parents;

	for(typename node_set_type::iterator it = to_Delete.begin(); it!=to_Delete.end(); ++it) {
		node_type * n = *it;
		node_type * np = n->parent;
		node_type * nl = n->str_left;
		node_type * nc = n->str_center;
		node_type * nr = n->str_right;
		if(nl!=NULL) {
			prune_subtree(nl);
			forest_add(nl, np, roots_to_parents, valid_parents);
		}
		if(nc!=NULL) {
			prune_subtree(nc);
			forest_add(nc, np, roots_to_parents, valid_parents);
		}
		if(nr!=NULL) {
			prune_subtree(nr);
			forest_add(nr, np, roots_to_parents, valid_parents);
		}
		if(root == n) {
			prune_subtree(n);
		} else {
			if(np!=NULL) prune_subtree(n); // no double pruning which would trigger assertions
		}
	}
	for(typename node_set_type::iterator it = to_Delete.begin(); it!=to_Delete.end(); ++it) {
		unbalanced_nodes.erase(*it);
		forest_delete_root(*it, valid_parents);
		roots_to_parents.erase(*it);
		delete_node(*it);
	}
	forest_insert(roots_to_parents, valid_parents);
#if DEBUG_PRINTS
	std::cout << " Unbalanced nodes after delete nodes " << unbalanced_nodes.size() << std::endl;
	for(typename node_set_type::iterator it = unbalanced_nodes.begin(); it!=unbalanced_nodes.end(); ++it) {
		node_type * n = *it;
		std::cout << " l " << n->left << " r " << n->right << " v " << n->value << std::endl;
	}

#endif

#if DO_REBALANCING
	set_rebalance();
#endif

}

template<typename T>
void avl_interval_tree<T>::erase_this_node(node_type * atn) {
	node_type * ins_node = atn->parent; // can be NULL if atn is root
	prune_subtree(atn);
	node_type * strl = atn->str_left;
	node_type * strc = atn->str_center;
	node_type * strr = atn->str_right;


// prune all before we start inserting again
	if(strc!=NULL) {
		prune_subtree(strc);
	}
	if(strl!=NULL) {
		prune_subtree(strl);
	}
	if(strr!=NULL) {
		prune_subtree(strr);
	}
//	if(ins_node!=NULL) unbalance_node(ins_node);

	if(strc!=NULL) {
//		std::cout << " center erase this move starts ins_node is " << ins_node << std::endl;


#if DO_SLOW_CHECKS
		size_t minl, maxr, intv_l, intv_r;
		check_upwards(ins_node, minl, maxr, intv_l, intv_r);
//	std::cout << " from above subtree limits l " << minl << " r " << maxr << " all in intv l " << intv_l << " r " << intv_r << std::endl;
		assert(minl <= strc->min_left);
		assert(maxr >= strc->max_right);
		assert(strc->left <= intv_r);
		assert(strc->right >= intv_l);
#endif
		ins_node = move_subtree(strc, ins_node);
//		std::cout << "after center erase check:"<<std::endl;
//		check_subtree(root);
	}
	if(strl!=NULL) {
//		std::cout << " left erase this move starts ins_node is " << ins_node << std::endl;

#if DO_SLOW_CHECKS
	size_t minl, maxr, intv_l, intv_r;
	check_upwards(ins_node, minl, maxr, intv_l, intv_r);
//	std::cout << " from above subtree limits l " << minl << " r " << maxr << " all in intv l " << intv_l << " r " << intv_r << std::endl;
	assert(minl <= strl->min_left);
	assert(maxr >= strl->max_right);
	assert(strl->left <= intv_r);
	assert(strl->right >= intv_l);
#endif

		ins_node = move_subtree(strl, ins_node);
//		std::cout << "after left erase check:"<<std::endl;
//		check_subtree(root);
	}
	if(strr!=NULL) {
//		std::cout << " right erase this move starts ins_node is " << ins_node << std::endl;

#if DO_SLOW_CHECKS
	size_t minl, maxr, intv_l, intv_r;
	check_upwards(ins_node, minl, maxr, intv_l, intv_r);
//	std::cout << " from above subtree limits l " << minl << " r " << maxr << " all in intv l " << intv_l << " r " << intv_r << std::endl;
	assert(minl <= strr->min_left);
	assert(maxr >= strr->max_right);
	assert(strr->left <= intv_r);
	assert(strr->right >= intv_l);
#endif

		ins_node = move_subtree(strr, ins_node);
//		std::cout << "after right erase check:"<<std::endl;
//		check_subtree(root);
	}
//	std::cout << " After erase_this_node check: " << std::endl;
//	check_subtree(root);
#if DO_REBALANCING
	if(ins_node!=NULL) {
		balance_up(ins_node);
	}
#endif
	delete_node(atn);

}


template<typename T>
void avl_interval_tree<T>::erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & erased) {
#if DEBUG_PRINTS
	std::cout << " erasing at l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
#endif
	if(right < atn->left) {
		if(atn->str_left!=NULL) {
			erase_at_node(atn->str_left, left, right, erased);
		} 
	} else if(left > atn->right) {
		if(atn->str_right!=NULL) {
			erase_at_node(atn->str_right, left, right, erased);
		} 
	} else {
		if(left == atn->left && right==atn->right) {
			erased.insert(atn->value);
			erase_this_node(atn);

		}
		if(left <= atn->right && right >= atn->left){
			if(atn->str_center!=NULL) {
				erase_at_node(atn->str_center, left, right, erased);
			} 
		}
	}
}

template<typename T>
size_t avl_interval_tree<T>::single_erase_at_node(node_type * atn, const size_t & left, const size_t & right, const value_type & value) {
#if DEBUG_PRINTS
	std::cout << "atn " << atn << " " <<left << " "<<right <<std::endl;
	std::cout << " erasing at l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
#endif
	size_t num_erased = 0;
	if(right < atn->left) {
		if(atn->str_left!=NULL) {
			num_erased += single_erase_at_node(atn->str_left, left, right, value);
			
		} 
	} else if(left > atn->right) {
		if(atn->str_right!=NULL) {
			num_erased += single_erase_at_node(atn->str_right, left, right, value);
		} 
	} else {
		if(left <= atn->right && right >= atn->left){
			if(atn->str_center!=NULL) {
				num_erased += single_erase_at_node(atn->str_center, left, right, value);
			} 
		}
		if((left == atn->left && right==atn->right) && value == atn->value) {
// Here we have to delete:
#if DEBUG_PRINTS
			std::cout << "ERASE THIS: " << atn->left << " r " << atn->right << " d " << atn->value << " max "<< atn->max_right<< std::endl;
#endif
			erase_this_node(atn);
			num_erased++;
			return num_erased; 
		}
	}
	return num_erased;	
}


template<typename T>	
void avl_interval_tree<T>::overlap_search_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results) const {
#if DEBUG_PRINTS
	std::cout << " ovlrp search for l " << left << " r " << right << " at  l " << atn->left << " r " << atn->right << " d " << atn->value <<std::endl;
#endif
	node_type * strl = atn->str_left;
	node_type * strr = atn->str_right;
	node_type * strc = atn->str_center;

	if(left < atn->left) {
		if(strl!=NULL) {
			if(right >= strl->min_left && left <= strl->max_right) {
				overlap_search_at_node(strl, left, right, results);
			}
		}
	}
	if(right > atn->right) {
		if(atn->str_right!=NULL) {
			if(right >= strr->min_left && left <= strr->max_right) {
				overlap_search_at_node(strr, left, right, results);
			}
		}
	}	
	if(atn->str_center!=NULL) {
			if(right >= strc->min_left && left <= strc->max_right) {
				overlap_search_at_node(strc, left, right, results);
			}
	}
	if(right >= atn->left && left <= atn->right) {
#if DEBUG_PRINTS
		std::cout << " FOUND l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
#endif
		results.push_back(atn->value);
	}
}

template<typename T>
void avl_interval_tree<T>::overlap_search_erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results, node_set_type & to_Delete) {
	
	node_type * strl = atn->str_left;
	node_type * strr = atn->str_right;
	node_type * strc = atn->str_center;
#if DEBUG_PRINTS	
	std::cout << " ovlrp search erase for l " << left << " r " << right << " at  l " << atn->left << " r " << atn->right << " d " << atn->value <<std::endl;
#endif
	if(left < atn->left) {
		if(strl!=NULL) {
			if(right >= strl->min_left && left <= strl->max_right) {
				overlap_search_erase_at_node(strl, left, right, results, to_Delete);
			}
		}
	}
	if(right > atn->right) {
		if(strr!=NULL) {
			if(right >= strr->min_left && left <= strr->max_right) {
				overlap_search_erase_at_node(strr, left, right, results, to_Delete);
			}
		}
	}
	
	if(strc!=NULL) {
		if(right >= strc->min_left && left <= strc->max_right) {
			overlap_search_erase_at_node(strc, left, right, results, to_Delete);
		}
	}
	if(right >= atn->left && left <= atn->right) {
#if DEBUG_PRINTS
		std::cout << " FOUND l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
#endif
		results.push_back(atn->value);
		to_Delete.insert(atn);

	}


}


template<typename T>
void avl_interval_tree<T>::find_at_node(node_type * atn, const size_t & left, const size_t & right, const value_type & value, node_set_type & result) const {
	if( (left == atn->left) && (right == atn->right) && (atn->value == value) ) {
		result.insert(atn);
	}
	if(right < atn->left) {
		if(atn->str_left!=NULL) {
			find_at_node(atn->str_left, left, right, value, result);
		}
	}

	if(left > atn->right) {
		if(atn->str_right!=NULL) {
			find_at_node(atn->str_right, left, right, value, result);
		}
	}

	if( (atn->right >= left) && (atn->left <= right)) {
		if(atn->str_center!=NULL) {
			find_at_node(atn->str_center, left, right, value, result);
		}
	} 
}


template<typename T>
void avl_interval_tree<T>::join_intervals_at_node(node_type * atn, std::map<size_t, size_t> & joined_intervals) const {
	if(atn->str_left!=NULL) {
		join_intervals_at_node(atn->str_left, joined_intervals);
	}	
	if(atn->str_right!=NULL) {
		join_intervals_at_node(atn->str_right, joined_intervals);
	}	
	size_t left = atn->left;
	size_t right = atn->right;
	if(atn->str_center!=NULL) {
		if(left > atn->str_center->min_left) left = atn->str_center->min_left;
		if(right < atn->str_center->max_right) right = atn->str_center->max_right;
//		subtree_max_interval(atn->str_center, left, right);
	}
	joined_intervals.insert(std::make_pair(left, right));
}


template<typename T>
void avl_interval_tree<T>::overlap_fraction_search_at_node(node_type * atn, const size_t & left, const size_t & right, const double & fraction, std::vector<value_type> & results) const {
	size_t minl = atn->min_left;
	size_t maxr = atn->max_right;

	size_t ovl = left;
	size_t ovr = right;
	if(ovl < minl) ovl = minl;
	if(ovr > maxr) ovr = maxr;
	if(ovl > ovr) return;
	double maxfraction = (double)(ovr - ovl + 1) / (double)(right - left + 1);
	if(maxfraction < fraction) return;

// intersection atn and search area
	size_t aleft = atn->left;
	size_t aright = atn->right;
	size_t aovl = aleft;
	size_t aovr = aright;
	if(aovl < left) aovl = left;
	if(aovr > right) aovr = right;
	

// union atn and search area
	size_t aunl = aleft;
	size_t aunr = aright;
	if(left < aunl) aunl = left;
	if(right > aunr) aunr = right;

	if(aovl < aovr) {
		double afraction = (double)(aovr - aovl + 1) / (double)(aunr - aunl + 1);
		if(afraction >= fraction) {
			results.push_back(atn->value);
		}
	}

	node_type * strl = atn->str_left;
	node_type * strr = atn->str_right;
	node_type * strc = atn->str_center;

	if(strl!=NULL) {
		if(right >= strl->min_left && left <= strl->max_right) {
			overlap_fraction_search_at_node(strl, left, right, fraction, results);
		}
	}
	if(strr!=NULL) {
		if(right >= strr->min_left && left <= strr->max_right) {
			overlap_fraction_search_at_node(strl, left, right, fraction, results);
		}
	}
	if(strc!=NULL) {
		if(right >= strc->min_left && left <= strc->max_right) {
			overlap_fraction_search_at_node(strl, left, right, fraction, results);
		}
	}
}


template<typename T>
size_t avl_interval_tree<T>::count_center_at_node(node_type * atn) const {
	size_t res = 0;
	
	node_type * strl = atn->str_left;
	node_type * strr = atn->str_right;
	node_type * strc = atn->str_center;
	if(strl!=NULL) {
		res += count_center_at_node(strl);
	}
	if(strr!=NULL) {
		res += count_center_at_node(strr);
	}
	if(strc!=NULL) {
		res++;
		res += count_center_at_node(strc);
	}

	return res;
}

/* almost the same as min_left  to max_right 


template<typename T>
void avl_interval_tree<T>::subtree_max_interval(node_type * atn, size_t & left, size_t & right) const {
	if(left > atn->left) left = atn->left;
	if(right < atn->right) right = atn->right;

	if(atn->str_left!=NULL) {
		subtree_max_interval(atn->str_left, left, right);
	}
	if(atn->str_center!=NULL) {
		subtree_max_interval(atn->str_center, left, right);
	}
	if(atn->str_right!=NULL) {
		subtree_max_interval(atn->str_right, left, right);
	}
}
*/

template<typename T>
void avl_interval_tree<T>::debug_print_at_node(node_type * atn, const size_t & level, size_t & id, const size_t & parent_id) const {
	std::cout << " node " << id << " on level " << level;
	node_type * p = atn->parent;
	if(p == NULL) {
		std::cout << " is root height " << atn->height<< std::endl;
	} else {
		std::cout << " is a ";
		if(p->str_left==atn) {	
			std::cout << " left ";
		} else if (p->str_center==atn) {
			std::cout << " center ";
		} else if (p->str_right==atn) {
			std::cout << " right ";
		} else {
			assert(0);
		}
		std::cout << "child of " << parent_id << " height " << atn->height << std::endl;
	}

	std::stringstream mrstr;
	if(atn->max_right<std::numeric_limits<size_t>::max()) {
		mrstr << atn->max_right;
	} else {
		mrstr << "MAX";
	}
	std::cout << "l " << atn->left << " r " << atn->right << " data pointer " << atn->value << std::endl;
	std::cout << " saved subtree size l " << atn->min_left << " r " << mrstr.str() << std::endl;

	size_t this_id = id;
	if(atn->str_left!=NULL) {
		debug_print_at_node(atn->str_left, level + 1, ++id, this_id);
	}
	if(atn->str_center!=NULL) {
		debug_print_at_node(atn->str_center, level + 1, ++id, this_id);
	}
	if(atn->str_right!=NULL) {
		debug_print_at_node(atn->str_right, level + 1, ++id, this_id);
	}

}


/* 
   this function takes a subtree (already pruned)
   inserts all its content at to

*/	
template<typename T>
it_node<T> * avl_interval_tree<T>::move_subtree(node_type * from, node_type * to) { 
	num_moves++;
#if DEBUG_PRINTS
	std::cout << "START move subtree function " << std::endl;
	if(from==NULL) {
		std::cout << " FROM is NULL" << std::endl;
	} else {
		std::cout << " FROM l " << from->left << " r " << from->right << " v " << from->value << std::endl;
	}
	if(to==NULL) {
		std::cout << " TO is NULL " << std::endl;
	} else {
		std::cout << " TO l " << to->left << " r " << to->right << " v " << to->value << std::endl;
	}
#endif

#if DO_SLOW_CHECKS
	assert(from!=NULL);
	assert(from->parent == NULL);
#if DEBUG_PRINTS
//	std::cout << "Check move subtree function: "<< from << " to " << to << std::endl;
//	std::cout << "Check of from:" << std::endl;
#endif
//	check_subtree(from);

// prepare to check that from actually fits at to
//	size_t minl, maxr, intv_l, intv_r;
//	check_upwards(to, minl, maxr, intv_l, intv_r);

#endif 


	if(to==NULL) {
		if(root==NULL) {
			root = from;
			return root;;
		}
		to = root;
	}

	
#if DO_SLOW_CHECKS

#if DEBUG_PRINTS
//	std::cout << "Check of to:" << std::endl;
#endif
//	check_subtree(to);

#if DEBUG_PRINTS
//	std::cout << " from above TO: FROM subtree limits l " << minl << " r " << maxr << " all in intv l " << intv_l << " r " << intv_r << std::endl;
#endif

//	assert(minl <= from->min_left);
//	assert(maxr >= from->max_right);
//	assert(from->left <= intv_r);
//	assert(from->right >= intv_l);


#endif




// At first we try to insert the entire subtree at from to to
	if(from->max_right < to->left) { // from fits into left subtree of to
		node_type * tol = to->str_left;
		if(tol == NULL) {
			left_attach(to, from);
			return to;
		} else {
			move_subtree(from, tol);
			return to;
		}
	} 
	if (from->min_left > to->right) { // from fits into right subtree of to
		node_type * tor = to->str_right;
		if(tor==NULL) {
			right_attach(to, from);
			return to;
		} else {
			move_subtree(from, tor);
			return to;
		}
	} 


// If we cannot insert the entire subtree, we insert from and its three subtrees separately
	node_type * strc = from->str_center;
	node_type * strl = from->str_left;
	node_type * strr = from->str_right;

	node_type * newto = to; // this could change the root of the current subtree



	if(strc!=NULL) {
		prune_subtree(strc);
	}
	if(strl!=NULL) {
		prune_subtree(strl);
	}
	if(strr!=NULL) {
		prune_subtree(strr);
	}



	if(strc!=NULL) { // center first, to do center insert switch operations first
/*			std::cout << " MOVE on center"<<std::endl;
	std::cout << "move subtree function " << std::endl;
	if(from==NULL) {
		std::cout << " FROM is NULL" << std::endl;
	} else {
		std::cout << " FROM l " << from->left << " r " << from->right << " v " << from->value << std::endl;
	}
	if(to==NULL) {
		std::cout << " TO is NULL " << std::endl;
	} else {
		std::cout << " TO l " << to->left << " r " << to->right << " v " << to->value << std::endl;
	}
*/
		node_type * tnewto =  move_subtree(strc, newto);
		newto = tnewto;
	}
#if DEBUG_PRINTS		
		std::cout << " INSERT IN MOVE l " << from->left << " r " << from->right << " v " << from->value << std::endl;
#endif
	node_type * tnewto = insert_at_node(newto, from); 
	newto = tnewto;

	if(strl!=NULL) {
/*
			std::cout<< " MOVE on left"<<std::endl;
	std::cout << "move subtree function " << std::endl;
	if(from==NULL) {
		std::cout << " FROM is NULL" << std::endl;
	} else {
		std::cout << " FROM l " << from->left << " r " << from->right << " v " << from->value << std::endl;
	}
	if(to==NULL) {
		std::cout << " TO is NULL " << std::endl;
	} else {
		std::cout << " TO l " << to->left << " r " << to->right << " v " << to->value << std::endl;
	}
*/

		node_type * tnewto = move_subtree(strl, newto);
		newto = tnewto;
	}

	if(strr!=NULL) {
/*
			std::cout << " MOVE on right"<<std::endl;
	std::cout << "move subtree function " << std::endl;
	if(from==NULL) {
		std::cout << " FROM is NULL" << std::endl;
	} else {
		std::cout << " FROM l " << from->left << " r " << from->right << " v " << from->value << std::endl;
	}
	if(to==NULL) {
		std::cout << " TO is NULL " << std::endl;
	} else {
		std::cout << " TO l " << to->left << " r " << to->right << " v " << to->value << std::endl;
		std::cout << " nTO l " << newto->left << " r " << newto->right << " v " << newto->value << std::endl;
	}
*/
		node_type * tnewto = move_subtree(strr, newto);
		newto = tnewto;
	}
	


#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << "Check result of move:" << std::endl;
	std::cout << "move subtree function " << std::endl;
	if(from==NULL) {
		std::cout << " FROM is NULL" << std::endl;
	} else {
		std::cout << " FROM l " << from->left << " r " << from->right << " v " << from->value << std::endl;
	}
	if(to==NULL) {
		std::cout << " TO is NULL " << std::endl;
	} else {
		std::cout << " TO l " << to->left << " r " << to->right << " v " << to->value << std::endl;
	}
#endif

//	if(newto!=NULL) {
//		check_subtree(newto);
//	} else {
//		check_subtree(root);
//	}
#endif

	return newto;
}


template<typename T>
void avl_interval_tree<T>::delete_subtree(node_type * n) {
	assert(n->parent==NULL);
	node_type * strc = n->str_center;
	node_type * strl = n->str_left;
	node_type * strr = n->str_right;

	if(strl!=NULL) {
		prune_subtree(strl);
		delete_subtree(strl);
	}
	if(strr!=NULL) {
		prune_subtree(strr);
		delete_subtree(strr);
	}
	if(strc!=NULL) {
		prune_subtree(strc);
		delete_subtree(strc);
	}
	delete_node(n);
}


template<typename T>
void avl_interval_tree<T>::prune_subtree(node_type * n) {
	assert(n!=NULL);
	if(n->parent==NULL) {
		assert(n==root);
		root = NULL;
	} else {
		node_type * p = n->parent;
		if(p->str_left == n) {
			p->str_left = NULL;
		} else if(p->str_center == n) {
			p->str_center = NULL;
		} else if(p->str_right == n) {
			p->str_right = NULL;
		} else {
			assert(0);
		}
		n->parent = NULL;
		if(p!=NULL) redo_height_up(p); 
	}
}


template<typename T>
void avl_interval_tree<T>::left_attach(node_type * parent, node_type * child) {
	assert(parent->str_left == NULL);
	assert(child->parent == NULL);
	parent->str_left = child;
	child->parent = parent;
	redo_height_up(parent);
}


template<typename T>
void avl_interval_tree<T>::right_attach(node_type * parent, node_type * child) {
	assert(parent->str_right == NULL);
	assert(child->parent == NULL);
	parent->str_right = child;
	child->parent = parent;
	redo_height_up(parent);
}


template<typename T>
void avl_interval_tree<T>::center_attach(node_type * parent, node_type * child) {
	assert(parent->str_center == NULL);
	assert(child->parent == NULL);
	parent->str_center = child;
	child->parent = parent;
	redo_height_up(parent);
}


template<typename T>
void avl_interval_tree<T>::delete_node(node_type * n) {
	assert(n->parent == NULL);
	assert(n->str_left == NULL);
	assert(n->str_right == NULL);
	assert(n->str_center == NULL);
	delete n;
}



/**
	return true if any value at n has changed.
	in that case we have to redo the parent of n
	


**/
template<typename T>	
bool avl_interval_tree<T>::redo_height_and_minmax_at_node(node_type * n) { 
//	std::cout << "n left in redo "<< n->left <<std::endl;
	bool changed = false;
	size_t nh = 0;
	size_t n_min_left = n->left;
	size_t n_max_right = n->right;
	node_type * nl = n->str_left;
	node_type * nr = n->str_right;
	node_type * nc = n->str_center;
	if(nl!=NULL) {
		size_t lh = nl->height + 1;
		if(lh>nh) nh = lh;
		if(nl->min_left < n_min_left) n_min_left = nl->min_left;
		if(nl->max_right > n_max_right) n_max_right = nl->max_right;
	}
	if(nr!=NULL) {
		size_t rh = nr->height + 1;
		if(rh>nh) nh = rh;
		if(nr->min_left < n_min_left) n_min_left = nr->min_left;
		if(nr->max_right > n_max_right) n_max_right = nr->max_right;
	}
	if(nc!=NULL) {

		if(nc->min_left < n_min_left) n_min_left = nc->min_left;
		if(nc->max_right > n_max_right) n_max_right = nc->max_right;

	}
//	std::cout << " redoing node " << n->value << " h " << n->height << " str limits l " << n->min_left << " r " << n->max_right << std::endl; 
	if(nh!=n->height) {
		n->height = nh;
		changed = true;
	}
	if(n_min_left!=n->min_left) {
		n->min_left = n_min_left;
		changed = true;
	}
	if(n_max_right!=n->max_right) {
		n->max_right = n_max_right;
		changed = true;
	}

	int bal = get_balance(n);
	if(bal*bal > 1) {
		unbalanced_nodes.insert(n);
	} else {
		unbalanced_nodes.erase(n);
	}
//	std::cout << " new node values " << n->value << " new h " << n->height << " str limits l " << n->min_left << " r " << n->max_right << std::endl; 
	return changed;
}


template<typename T>
int avl_interval_tree<T>::get_balance(node_type * n) const {
	int bal = 0;
	node_type * strl = n->str_left;
	node_type * strr = n->str_right;
	if(strl!=NULL) {
		bal += strl->height;
	}
	if(strr!=NULL) {
		bal -= strr->height;
	}
	return bal;
}




template<typename T>
void avl_interval_tree<T>::balance_up(node_type * n) {
#if DEBUG_PRINTS
	std::cout << "Balance up at l " << n->left << " r " << n->right << " unb nodes size " << unbalanced_nodes.size() << std::endl;
#endif
	node_type * tn = n;
	while(tn!=NULL) {
		int balance = get_balance(tn);
		if(balance*balance >1) {
			tn = balance_node(tn);
		}
		tn = tn->parent;
	}
}


/** 
	local rebalancing at n
	return new root of the subtree after rebalancing


**/
template<typename T>
it_node<T> * avl_interval_tree<T>::balance_node(node_type * n) {
	assert(DO_REBALANCING);
	num_reba++;
	int balance = get_balance(n);
#if DEBUG_PRINTS
	std::cout << " Balance subtree (bal is " << balance << ") : " << std::endl;
	std::cout << " node is l " << n->left << " r " << n->right << " v " << n->value << std::endl;
#endif
	node_type * res = n;
	if(balance > 1) {
		// left subtree to high
		node_type * strl = n->str_left;
		assert(strl!=NULL);
		int lbal = get_balance(strl);
		if(lbal <0) {
			
#if DEBUG_PRINTS
			std::cout << "LR to LL rotation " << std::endl;
#endif
	
			rotate_lr_to_ll(n);
			balance = get_balance(n);
		}
		if(balance > 1) {
#if DEBUG_PRINTS
			std::cout << "LEFT rotation " << balance << std::endl;
#endif
			res = rotate_left(n);
		}

	}
	if(balance < -1) {
		// right subtree to high
		node_type * strr = n->str_right;
		assert(strr!=NULL);
		int rbal = get_balance(strr);
		if(rbal >0) {
#if DEBUG_PRINTS
			std::cout << "RL to RR rotation " << std::endl;
#endif
			rotate_rl_to_rr(n);
			balance = get_balance(n);
		}
		if(balance < -1) {
#if DEBUG_PRINTS
			std::cout << "RIGHT rotation " << balance << std::endl;
#endif
			res = rotate_right(n);
		}
	}



#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
	std::cout << " After balance check " << std::endl;
	std::cout << " now at node l " << res->left << " r " << res->right << " v " << res->value << std::endl;
	
#endif
	check_subtree(root);
#endif
	return res;

}


#define NUM_FOR_REBALANCE 3

template<typename T>
void avl_interval_tree<T>::set_rebalance() {

	for(size_t nr = 0; nr < NUM_FOR_REBALANCE; ++nr) {
		if(!unbalanced_nodes.empty()) {
			node_type * n = *(unbalanced_nodes.begin());
#if DEBUG_PRINTS
		std::cout << " we will rebalance at l " << n->left << " r " << n->right << " v " << n->value << std::endl;	
		std::cout << " num nodes is " << num_nodes << " num unbalanced nodes is " << unbalanced_nodes.size() << std::endl;
#endif
			balance_up(n);
		}

	}


	
#if DO_SLOW_CHECKS
#if DEBUG_PRINTS
		std::cout << " End of set rebalance check of balance:" << std::endl;
		std::cout << " we have " << num_nodes <<" nodes , " << unbalanced_nodes.size() << " may not be in balance " << std::endl;
#endif
		if(root!=NULL) balance_check_at_node(root); 
#endif


}


template<typename T>
void avl_interval_tree<T>::rotate_lr_to_ll(node_type * n) {
/*
       5
      / d
     3
    a \
       4
      b c
TO 
      5
     / d
    4
   / c
  3
 a b

3 moves down. All intervals in its center subtree are now required to be stricly left of 4. They are guaranteed to be be left of 5. Therefore we do a subtree move to 4
*/



	node_type * node5 = n;
	node_type * node3 = n->str_left;
	assert(node3!=NULL);
	node_type * node4 = node3->str_right;
	assert(node4!=NULL);
	prune_subtree(node3);
	prune_subtree(node4);
	node_type * n3_center = node3->str_center;
	if(n3_center!=NULL) {
		prune_subtree(n3_center);
	}
	node_type * node_b = node4->str_left;
	if(node_b!=NULL) {
		prune_subtree(node_b);
	}
	
	if(node_b!=NULL) {
		right_attach(node3, node_b);
	}
	left_attach(node4, node3);
	left_attach(node5, node4);

//	std::cout << " LR LL rotate before move " << std::endl;
//	check_subtree(node5);


	if(n3_center!=NULL) {
		move_subtree(n3_center, node4); // this can exchange node 4, but we dont care, node5 stays parent of this subtree
	}

	
}


template<typename T>
void avl_interval_tree<T>::rotate_rl_to_rr(node_type * n) {
/*
       3
      a \
         5
        / d
       4
      b c
TO
       3
      a \
         4
        b \
           5
          c d

5 moves down, as above its center subtree has to be moved to 4
*/
	node_type * parent = n->parent;
	node_type * node3 = n;
	node_type * node5 = node3->str_right;
	assert(node5!=NULL);
	node_type * node4 = node5->str_left;
	assert(node4!=NULL);
	prune_subtree(node5);
	prune_subtree(node4);
	node_type * node_c = node4->str_right;
	if(node_c!=NULL) {
		prune_subtree(node_c);
	}
	node_type * n5_center = node5->str_center;
	if(n5_center!=NULL) {
		prune_subtree(n5_center);
	}
	
	if(node_c!=NULL) {
		left_attach(node5, node_c);
	}
	right_attach(node4, node5);
	right_attach(node3, node4);

//	std::cout << " RL RR rotate before move " << std::endl;
//	check_subtree(node3);


	if(n5_center!=NULL) {		
		move_subtree(n5_center, node4);
	}


}


template<typename T>
it_node<T> * avl_interval_tree<T>::rotate_left(node_type * n) {
/*
       5
      / d
     4
    / c
   3
  a b

TO

     4
    / \
   3   5
  a b c d

center of 5 moves to 4

*/
	// remember how we must reattach the new parent of that subtree
	node_type * parent = n->parent;
	bool pleft = false;
	bool pcenter = false;
	bool pright = false;
	if(parent!=NULL) {
		if(parent->str_left==n) pleft = true;
		else if(parent->str_center==n) pcenter = true;
		else if(parent->str_right==n) pright = true;
		else assert(0);
	}

	node_type * node5 = n;
	node_type * node4 = node5->str_left;
	assert(node4!=NULL);
	node_type * node3 = node4->str_left;
	assert(node3!=NULL);
	prune_subtree(node5);
	prune_subtree(node4);
	prune_subtree(node3);
	node_type * node_c = node4->str_right;
	if(node_c!=NULL) {
		prune_subtree(node_c);
	}
	node_type * n5_center = node5->str_center;
	if(n5_center!=NULL) {
		prune_subtree(n5_center);
	}


	if(node_c!=NULL) {
		left_attach(node5, node_c);
	}
	left_attach(node4, node3);
	right_attach(node4, node5);


//	std::cout << " LEFT rotate before move " << std::endl;
//	check_subtree(node4);

	node_type * new_node4 = node4;
	if(n5_center!=NULL) {
		new_node4 = move_subtree(n5_center, node4);
//		if(new_node4 != node4) {
//			unbalance_node(node4); // maybe node4 is not done yet (if balance was more than 2)
//		}
	}
	
	if(pleft) left_attach(parent, new_node4);
	if(pcenter) center_attach(parent, new_node4);
	if(pright) right_attach(parent, new_node4);
	if(parent==NULL) {
		root=new_node4;
	}

	return new_node4;
}


template<typename T>
it_node<T> * avl_interval_tree<T>::rotate_right(node_type *n) {
/*
      3
     a \
        4
       b \
          5
         c d


TO 

       4
      / \
     3   5
    a b c d	

center of 3 has to move to 4

*/
// remember how we must reattach the new parent of that subtree
	node_type * parent = n->parent;
	bool pleft = false;
	bool pcenter = false;
	bool pright = false;
	if(parent!=NULL) {
		if(parent->str_left==n) pleft = true;
		else if(parent->str_center==n) pcenter = true;
		else if(parent->str_right==n) pright = true;
		else assert(0);
	}
	node_type * node3 = n;
	node_type * node4 = node3->str_right;
	assert(node4!=NULL);
	node_type * node5 = node4->str_right;
	assert(node5!=NULL);
	prune_subtree(node3);
	prune_subtree(node4);
	prune_subtree(node5);
	node_type * node_b = node4->str_left;
	if(node_b!=NULL) {
		prune_subtree(node_b);
	}
	node_type * n3_center = node3->str_center;
	if(n3_center!=NULL) {
		prune_subtree(n3_center);
	}

	if(node_b!=NULL) {
		right_attach(node3, node_b);
	}
	left_attach(node4, node3);
	right_attach(node4, node5);


//	std::cout << " RIGHT rotate before move " << std::endl;
//	check_subtree(node4);

	node_type * new_node4 = node4;
	if(n3_center!=NULL) {
		new_node4 = move_subtree(n3_center, node4);
//		if(new_node4!=node4) {
//			unbalance_node(node4);
//		}
	}

	if(pleft) left_attach(parent, new_node4);
	if(pcenter) center_attach(parent, new_node4);
	if(pright) right_attach(parent, new_node4);
	if(parent==NULL) {
		root = new_node4;
	}



	return new_node4;
}


/** Do all kinds of slow checks 

**/
template<typename T>
void avl_interval_tree<T>::check_tree() const {
#if DO_SLOW_CHECKS
	std::cerr << "+++++    Warning: We are running slow checks of AVL interval tree with " << num_nodes << " nodes" << std::endl;

	size_t tsize = 0;
	size_t minl = 0;
	size_t maxr = std::numeric_limits<size_t>::max();
	size_t mini = 0;
	size_t maxi = std::numeric_limits<size_t>::max();
#if DEBUG_PRINTS
	std::cout<< "slow index size " << slow_index.size()<<std::endl;
#endif

	assert(slow_index.size() == num_nodes);
	if(root!=NULL) {
		size_t id = 0;
		check_tree_at_node(root, 0, id, 0, tsize, minl, maxr, mini, maxi, true);
#if DEBUG_PRINTS
		std::cout << " after check num_nodes " << num_nodes << " tree size " << tsize << std::endl;
#endif
		assert(num_nodes == tsize);
	} else {
		assert(num_nodes==0);
	}
	
#endif
}


template<typename T>
size_t avl_interval_tree<T>::check_subtree(node_type * n) const {
#if DEBUG_PRINTS
	std::cout << "######    Warning: We are running slow checks of AVL interval tree on a subtree " <<  std::endl;
#endif
	size_t tsize = 0;
	size_t minl = 0;
	size_t maxr = std::numeric_limits<size_t>::max();
	size_t mini = 0;
	size_t maxi = std::numeric_limits<size_t>::max();
	size_t id = 0;
	check_tree_at_node(n, 0, id, 0, tsize, minl, maxr, mini, maxi, false);
#if DEBUG_PRINTS
	std::cout << "#####   Checked subtree had " << tsize << "nodes " << std::endl;
#endif
	return tsize;
}

// use_slow_index only in method when we choose that option (then we also check for balance)

template<typename T>
void avl_interval_tree<T>::check_tree_at_node(node_type * n,  const size_t & level, size_t & id, const size_t & parent_id, size_t & subtree_size, size_t min_left, size_t max_right, size_t min_in, size_t max_in, bool use_slow_index) const {
#if DO_SLOW_CHECKS
	assert(n!=NULL);
	int bal = get_balance(n);
#if DEBUG_PRINTS
	std::cout << "*** node " << id << " on level " << level;  
#endif
	node_type * p = n->parent;
#if DEBUG_PRINTS
	if(p == NULL) {
		std::cout << " is root height " << n->height<<  " balance " << bal << std::endl;
	} else {
		std::cout << " is a ";
		if(p->str_left==n) {	
			std::cout << " left ";
		} else if (p->str_center==n) {
			std::cout << " center ";
		} else if (p->str_right==n) {
			std::cout << " right ";
		} else {
			assert(0);
		}
		std::cout << "child of " << parent_id << " height " << n->height << " balance " << bal << std::endl;
	}
	std::stringstream mrstr;
	std::stringstream rinstr;
	std::stringstream smrstr;
//	std::cout << "max right "<<max_right <<std::endl;
	if(max_right<std::numeric_limits<size_t>::max()) {
//		std::cout<< "inja"<<std::endl;
		mrstr << max_right;
	} else {
		mrstr << "MAX";
	}
	if(max_in <std::numeric_limits<size_t>::max()) {
		rinstr << max_in;
	} else {
		rinstr << "MAX";
	}
	if(n->max_right < std::numeric_limits<size_t>::max()) {
		smrstr << n->max_right;
	} else {
		smrstr << "MAX";
	}
	std::cout << "   L " << n->left << " R " << n->right << "  data pointer " << n->value << std::endl;
	std::cout << "   saved subtree size l " << n->min_left << " r " << smrstr.str();
	std::cout << "     current borders (computed for check) l " << min_left << " r " << mrstr.str() << " all intv in subtr ovrlp with " << min_in << " " << rinstr.str()<<  std::endl;
#endif
	size_t this_id = id;

	if(use_slow_index) {
	 	bool found_in_slow = false;
//	std::cout << " slow size " << slow_index.size() << std::endl;
		for(size_t i=0; i<slow_index.size(); i++) {
			size_t tl = slow_index.at(i).first.first;
			size_t tr = slow_index.at(i).first.second;
			value_type tv = slow_index.at(i).second;

//		std::cout << " slow ind " << i << " is l " << tl << " r " << tr << " d " << tv << std::endl; 
			if( (tl == n->left && tr == n->right) && tv==n->value) {
				found_in_slow = true;
				break;
			}
		}
		assert(found_in_slow);
	}

	// check saved subtree with bounds from parent nodes
//	std::cout << max_right <<std::endl;
	assert(n->min_left >= min_left);
	assert(n->max_right <= max_right);




	assert(n->left < n->right);
	assert(min_left <= n->left);
	assert(max_right >= n->right);

	assert( min_in <= n->right );
	assert( max_in >= n->left );


	node_type * strc = n->str_center;
	node_type * strl = n->str_left;
	node_type * strr = n->str_right;

	size_t height_check = 0;	
	size_t local_size = 0;
	size_t real_subtree_min = n->left;
	size_t real_subtree_max = n->right;
	
	if(strl!=NULL) {
//		std::cout << "here1"<<std::endl;
		assert(n == strl->parent);
		size_t s = 0;
		size_t minl = min_left;
		size_t maxr = max_right;
		size_t mini = min_in;
		size_t maxi = max_in;
		assert(n->left > 1); // no alignments of length 1 (if there is someting to the left of n, n cannot have left end at 1)
		size_t new_max = n->left - 1;
		if(maxr > new_max) maxr = new_max;
//		std::cout << "maxr "<<maxr<<std::endl;
//		std::cout << "strl "<< strl->max_right << std::endl;
		check_tree_at_node(strl, level+1, ++id, this_id, s, minl, maxr, mini, maxi, use_slow_index);
		local_size+=s;
		size_t new_height = strl->height + 1;
		if(new_height > height_check) height_check = new_height;
		if(real_subtree_min > strl->min_left) real_subtree_min = strl->min_left;
		if(real_subtree_max < strl->max_right) real_subtree_max = strl->max_right;
	}
	if(strr!=NULL) {
//		std::cout << "here2"<<std::endl;

		assert(n == strr->parent);
		size_t s = 0;
		size_t minl = min_left;
		size_t maxr = max_right;
		size_t mini = min_in;
		size_t maxi = max_in;
		size_t new_min = n->right + 1;
		if(minl < new_min) minl = new_min;
		check_tree_at_node(strr, level+1, ++id, this_id,  s, minl, maxr, mini, maxi, use_slow_index);
		local_size+=s;
		size_t new_height = strr->height + 1;
		if(new_height > height_check) height_check = new_height;
		if(real_subtree_min > strr->min_left) real_subtree_min = strr->min_left;
		if(real_subtree_max < strr->max_right) real_subtree_max = strr->max_right;
	}
	if(strc!=NULL) {
//		std::cout << "here3"<<std::endl;

		assert(n == strc->parent);
		size_t s = 0;
		size_t minl = min_left;
		size_t maxr = max_right;
		size_t mini = min_in;
		size_t maxi = max_in;
		size_t new_lo = n->left;
		size_t new_hi = n->right;
		if(new_lo > mini) mini = new_lo;
		if(new_hi < maxi) maxi = new_hi;
		check_tree_at_node(strc, level+1, ++id, this_id,  s, minl, maxr, mini, maxi, use_slow_index);
		local_size+=s;	
		if(real_subtree_min > strc->min_left) real_subtree_min = strc->min_left;
		if(real_subtree_max < strc->max_right) real_subtree_max = strc->max_right;
	}

	subtree_size += 1 + local_size;	
	assert(n->height == height_check);
#if DEBUG_PRINTS
	std::cout << "*   Node " << this_id << " is root of a subtree of size " << subtree_size << " its limits are l " << real_subtree_min << " r " << real_subtree_max << std::endl;
	std::cout << "     saved limits l " << n->min_left << " r " << n->max_right << std::endl;
#endif
	assert(real_subtree_min == n->min_left);
	assert(real_subtree_max == n->max_right);




#endif


	
}


template<typename T>
void avl_interval_tree<T>::check_upwards(node_type * n, size_t & minl, size_t & maxr, size_t & intv_l, size_t & intv_r) const {
	assert(DO_SLOW_CHECKS);
	minl = 0;
	maxr = std::numeric_limits<size_t>::max();
	intv_l = 0;
	intv_r = std::numeric_limits<size_t>::max();

	node_type * parent = NULL;
	node_type * node = n;
	if(node!=NULL) parent = node->parent;
	while(parent!=NULL) {
#if DEBUG_PRINTS
//		std::cout << " parent l " << parent->left << " r " << parent->right << std::endl;
#endif
		if(parent->str_left==node) {
			if(maxr > parent->left -1) maxr = parent->left - 1;
		} else if (parent->str_right == node) {
			if(minl < parent->right+1) minl = parent->right+ 1;
		} else if (parent->str_center == node) {
			if(intv_l < parent->left) intv_l = parent->left;
			if(intv_r > parent->right) intv_r = parent->right;
		} else assert(0);
		


		node = parent;
		parent = node->parent;
	}



}


template<typename T>
bool avl_interval_tree<T>::balance_check_at_node(node_type * n) const {
	assert(DO_SLOW_CHECKS);
	int bal = get_balance(n);
	bool in_balance = true;
	bool subtree_reba = false;

	typename node_set_type::iterator findn = unbalanced_nodes.find(n);
	if(findn==unbalanced_nodes.end()) {
	} else  {
		std::cout << " Node l " << n->left << " r " << n->right << " v " << n->value << " bal " << bal  << " may be not balanced " << std::endl;
		subtree_reba = true;
	}


	if(bal*bal > 1) {
		in_balance = false;
	}
	node_type * nl = n->str_left;
	node_type * nr = n->str_right;
	node_type * nc = n->str_center;
	if(nl!=NULL) {
		bool s = balance_check_at_node(nl);
		subtree_reba = s || subtree_reba;
	}
	if(nr!=NULL) {
		bool s = balance_check_at_node(nr);
		subtree_reba = s || subtree_reba;
	}
	if(nc!=NULL) {
		bool s = balance_check_at_node(nc);
		subtree_reba = s || subtree_reba;
	}
	std::cout << " node l " << n->left << " r " << n->right <<  " height " << n->height <<  " bal " << bal << " in balance: " << in_balance <<  " subtree reba is " << subtree_reba << std::endl;	
	assert(in_balance || subtree_reba);
	return subtree_reba;
}



