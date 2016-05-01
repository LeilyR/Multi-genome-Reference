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
avl_interval_tree<T>::avl_interval_tree(): root(NULL), num_nodes(0) {


}

template<typename T>
avl_interval_tree<T>::~avl_interval_tree() {
	std::cout << " Destroy tree " << std::endl;
	if(root!=NULL) {
		node_type * root_node = root;
		std::cout << " call delete on l " << root->left << " r " << root->right << " d " << root->value << std::endl; 
		prune_subtree(root);
		delete_subtree(root_node);
	}


	std::cout << " tree destroyed " << std::endl;
}

template<typename T>
void avl_interval_tree<T>::insert(const size_t & left, const size_t & right, const value_type  & v) {
#if DO_SLOW_CHECKS
	assert(slow_index.size() == num_nodes);

	std::cout << " insert l " << left << " r " << right << " d " << v << std::endl;
#endif

	node_type * new_node = new node_type(left, right, v);
	num_nodes++;
	if(root == NULL) {
		root = new_node;
	} else {
		insert_at_node(root, new_node);
	}

#if DO_SLOW_CHECKS
	slow_index.push_back(std::make_pair(std::make_pair(left, right), v));


	assert(slow_index.size() == num_nodes);
	std::cout << " After insert check " << std::endl;
	check_tree();
#endif
	
}

template<typename T>
void avl_interval_tree<T>::erase_all(const size_t & left, const size_t & right, std::vector<value_type> & erased) {

#if DO_SLOW_CHECKS
	std::cout << " Before erase all check: " << std::endl;
	std::cout << " try to erase: l " << left << " r " << right << std::endl;
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
#if DO_SLOW_CHECKS
	std::cout << " erase all result: " << erased.size() << " nodes were erased" << std::endl;
	assert(slow_index.size() == num_nodes);
	assert(slow_erased == erased.size());
	std::cout << " after erase all check: " << std::endl;
	check_tree();

#endif

}
template<typename T>
size_t avl_interval_tree<T>::erase(const size_t & left, const size_t & right, const value_type & value) {

#if DO_SLOW_CHECKS
	std::cout << " Before single single erase check: " << std::endl;
	std::cout << " try to erase: l " << left << " r " << right << " d " << value << std::endl;
	check_tree();
	assert(num_nodes == slow_index.size());

	
	size_t slow_erased = 0;
	std::set<size_t> to_del;
	for(size_t i=0; i<slow_index.size(); ++i) {
		size_t tl = slow_index.at(i).first.first;
		size_t tr = slow_index.at(i).first.second;
		value_type tv = slow_index.at(i).second;
		if((tl == left && tr == right) && (tv == value)) {
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

//		std::cout << " after deleting " << to_del.size() << " slow index size is " << slow_index.size() << std::endl;
	}

#endif



	size_t num_erased = 0;
	if(root!=NULL) {
		num_erased = single_erase_at_node(root, left, right, value);
		num_nodes -= num_erased;
	} 
#if DO_SLOW_CHECKS
	std::cout << " erase result: " << num_erased << " nodes were erased" << std::endl;
	std::cout << " slow results: " << to_del.size() << " erased, now we have " << slow_index.size() << std::endl;
	assert(slow_index.size() == num_nodes);
	assert(slow_erased == num_erased);

	std::cout << " after single erase check: " << std::endl;
	check_tree();

#endif
	return num_erased;
}


template<typename T>	
void avl_interval_tree<T>::overlap_search(const size_t & left, const size_t & right, std::vector<value_type> & results) {
#if DO_SLOW_CHECKS
	std::cout << " Before overlap search check: " << std::endl;
	std::cout << " check for overlap with l " << left << " r " << right <<  std::endl;
	check_tree();
#endif


	if(root!=NULL) {
		overlap_search_at_node(root, left, right, results);
	}

#if DO_SLOW_CHECKS
	assert(slow_index.size() == num_nodes);
	std::cout << "We have found " << results.size() << " overlapping intervals using the tree "<<std::endl;
	for(size_t i=0; i<results.size(); ++i) {
		std::cout << " v " << results.at(i) << std::endl;
	}
	std::vector<value_type> slow_results;
	for(size_t i=0; i<slow_index.size(); ++i) {
		size_t tl = slow_index.at(i).first.first;
		size_t tr = slow_index.at(i).first.second;
		value_type tv = slow_index.at(i).second;
		if(tl <= right && tr >= left) {
			std::cout << " slow found  l " << tl << " r " << tr << " d " << tv << " border " << left<< " r " << right << std::endl;
			slow_results.push_back(tv);
		} else {
//			std::cout << " slow no ovrlp  l " << tl << " r " << tr << " d " << tv << std::endl;

		}
	}
	std::cout << "We have found " << slow_results.size() << " overlapping intervals using a simple method " << std::endl;
	assert(slow_results.size() == results.size());
	for(size_t i=0; i<slow_results.size(); ++i) {
		std::cout << " slow res " << slow_results.at(i) << std::endl;
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
	std::cout << " Before overlap search erase check: " << std::endl;
	std::cout << " check for overlap with l " << left << " r " << right <<  std::endl;
	check_tree();

	
	std::vector< std::pair<std::pair<size_t, size_t>, value_type> > slow_index_backup = slow_index;
#endif

	if(root!=NULL) {
		node_map_type roots_to_parents;
		node_set_type valid_parents;
		overlap_search_erase_at_node(root, left, right, results, roots_to_parents, valid_parents);
		forest_insert(roots_to_parents, valid_parents);
	}

#if DO_SLOW_CHECKS
	assert(slow_index.size() == num_nodes);
	std::cout << "We have found " << results.size() << " overlapping intervals using the tree "<<std::endl;
	for(size_t i=0; i<results.size(); ++i) {
		std::cout << " v " << results.at(i) << std::endl;
	}
	std::vector<value_type> slow_results;
	for(size_t i=0; i<slow_index_backup.size(); ++i) {
		size_t tl = slow_index_backup.at(i).first.first;
		size_t tr = slow_index_backup.at(i).first.second;
		value_type tv = slow_index_backup.at(i).second;
		if(tl <= right && tr >= left) {
//			std::cout << " slow found  l " << tl << " r " << tr << " d " << tv << std::endl;
			slow_results.push_back(tv);
		} else {
//			std::cout << " slow no ovrlp  l " << tl << " r " << tr << " d " << tv << std::endl;

		}
	}
	std::cout << "We have found " << slow_results.size() << " overlapping intervals using a simple method " << std::endl;
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

	std::cout << "After overlap search erase check" <<std::endl;
	check_tree();
#endif




}

template<typename T>
void avl_interval_tree<T>::join_intervals(std::vector<std::pair<size_t, size_t> > & results) {
	std::map<size_t, size_t> joined_intervals; // left -> right

	if(root!=NULL) {
		join_intervals_at_node(root, joined_intervals);
	}

	for(std::map<size_t, size_t>::iterator it = joined_intervals.begin(); it!=joined_intervals.end(); ++it) {
		results.push_back(*it);
	}

}

template <typename T>	
void avl_interval_tree<T>::debug_print() const {
	std::cout << " AVL Interval tree printout: " << std::endl;
	if(root!=NULL) {
		size_t id = 0;
		debug_print_at_node(root, 0, id, 0);
	}

}

/*
	insert at assumes that heigt of newn is already correct


*/

template<typename T>
void avl_interval_tree<T>::insert_at_node(node_type * atn, node_type * newn) {
	size_t max_subtr_height = 0;
	if(atn==NULL) {
		if(root == NULL) {
			root = newn;
		} else {
			insert_at_node(root, newn);
			redo_height_and_minmax_at_node(root);
		}
		return;
	}

	if(newn->right < atn->left) {
// insert into left subtree
		if(atn->str_left==NULL) {
			left_attach(atn, newn);
		} else {
			insert_at_node(atn->str_left, newn);
		}
	} else if(newn->left > atn->right) {
// insert into right subtree
		if(atn->str_right==NULL) {
			right_attach(atn, newn);
		} else {
			insert_at_node(atn->str_right, newn);
		}
	} else {
// insert into center subtree
		if(atn->str_center==NULL) {
			center_attach(atn, newn);
		} else {
			insert_at_node(atn->str_center, newn);
		}
	}	
	redo_height_and_minmax_at_node(atn);
}

template<typename T>
void avl_interval_tree<T>::forest_add(node_type * atn, node_type * insert_node, node_map_type & roots_to_parents, node_set_type & valid_parents) {
	roots_to_parents.insert(std::make_pair(atn, insert_node));
	valid_parents.insert(insert_node);
}

template<typename T>
void avl_interval_tree<T>::redo_height_up(node_type * n) {
	while(n != NULL) {
		redo_height_and_minmax_at_node(n);
		n=n->parent;
	}	
}

template<typename T>
void avl_interval_tree<T>::forest_insert(node_map_type & roots_to_parents, node_set_type & valid_parents) {
	for(typename node_map_type::iterator it = roots_to_parents.begin(); it!=roots_to_parents.end(); ++it) {
		node_type * subtree_root = it->first;
		node_type * ins_node = it->second;
		typename node_set_type::iterator find_ian = valid_parents.find(ins_node);
		if(find_ian==valid_parents.end()) {
			ins_node = NULL;
		}
		move_subtree(subtree_root, ins_node);
	}

	// TODO all downwards rotations first, then all upwards rotations

// Correct upwards heights, from all parents where we did insert
	for(typename node_set_type::iterator it = valid_parents.begin(); it!=valid_parents.end(); it++) {

		node_type * n = *it;
		redo_height_up(n); // TODO here or in move_subtree?
			
	}

	for(typename node_map_type::iterator it = roots_to_parents.begin(); it!=roots_to_parents.end(); ++it) {
		node_type * n = it->first;
		redo_height_up(n);
	}

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
		prune_subtree(n);
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
	}
	for(typename node_set_type::iterator it = to_Delete.begin(); it!=to_Delete.end(); ++it) {
		forest_delete_root(*it);
		delete_node(*it);
	}
	forest_insert(roots_to_parents, valid_parents);
}

template<typename T>
void avl_interval_tree<T>::erase_this_node(node_type * atn) {
	node_type * ins_node = atn->parent; // can be NULL if atn is root
	node_map_type roots_to_parents;
	node_set_type valid_parents;
	prune_subtree(atn);
	node_type * strl = atn->str_left;
	node_type * strc = atn->str_center;
	node_type * strr = atn->str_right;
	if(strc!=NULL) {
		prune_subtree(strc);
		forest_add(strc, ins_node, roots_to_parents, valid_parents);
	}
	if(strl!=NULL) {
		prune_subtree(strl);
		forest_add(strl, ins_node, roots_to_parents, valid_parents);
		}
	if(strr!=NULL) {
		prune_subtree(strr);
		forest_add(strr, ins_node, roots_to_parents, valid_parents);
	}
	delete_node(atn);
	forest_insert(roots_to_parents, valid_parents);
}


template<typename T>
void avl_interval_tree<T>::erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & erased) {
	std::cout << " erasing at l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
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
	std::cout << " erasing at l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
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
			std::cout << "ERASE THIS: " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
			erase_this_node(atn);
			num_erased++;
			return num_erased; // no redo height on now invalid atn
		}
	}
	redo_height_and_minmax_at_node(atn);
	return num_erased;	
}


template<typename T>	
void avl_interval_tree<T>::overlap_search_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results) {
	std::cout << " ovlrp search for l " << left << " r " << right << " at  l " << atn->left << " r " << atn->right << " d " << atn->value <<std::endl;
	if(left < atn->left) {
		if(atn->str_left!=NULL) {
			overlap_search_at_node(atn->str_left, left, right, results);
		}
	}
	if(right > atn->right) {
		if(atn->str_right!=NULL) {
			overlap_search_at_node(atn->str_right, left, right, results);
		}
	}	
	if(right >= atn->left && left <= atn->right) {
		if(atn->str_center!=NULL) {
			overlap_search_at_node(atn->str_center, left, right, results);
		}
		std::cout << " FOUND l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
		results.push_back(atn->value);
	}
}

template<typename T>
void avl_interval_tree<T>::overlap_search_erase_at_node(node_type * atn, const size_t & left, const size_t & right, std::vector<value_type> & results, node_set_type & to_Delete) {
	std::cout << " ovlrp search erase for l " << left << " r " << right << " at  l " << atn->left << " r " << atn->right << " d " << atn->value <<std::endl;
	if(left < atn->left) {
		if(atn->str_left!=NULL) {
			overlap_search_erase_at_node(atn->str_left, left, right, results, to_Delete);
		}
	}
	if(right > atn->right) {
		if(atn->str_right!=NULL) {
			overlap_search_erase_at_node(atn->str_right, left, right, results, to_Delete);
		}
	}
	
	if(right >= atn->left && left <= atn->right) {
		if(atn->str_center!=NULL) {
			overlap_search_erase_at_node(atn->str_center, left, right, results, to_Delete);
		}
		std::cout << " FOUND l " << atn->left << " r " << atn->right << " d " << atn->value << std::endl;
		results.push_back(atn->value);
		to_Delete.insert(atn);

	}


}

template<typename T>
void avl_interval_tree<T>::join_intervals_at_node(node_type * atn, std::map<size_t, size_t> & joined_intervals) {
	if(atn->str_left!=NULL) {
		join_intervals_at_node(atn->str_left, joined_intervals);
	}	
	if(atn->str_right!=NULL) {
		join_intervals_at_node(atn->str_right, joined_intervals);
	}	
	size_t left = atn->left;
	size_t right = atn->right;
	if(atn->str_center!=NULL) {
		subtree_max_interval(atn->str_center, left, right);
	}
	joined_intervals.insert(std::make_pair(left, right));
}


template<typename T>
void avl_interval_tree<T>::subtree_max_interval(node_type * atn, size_t & left, size_t & right) {
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

	std::cout << "l " << atn->left << " r " << atn->right << " data pointer " << atn->value << std::endl;
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
		TODO can an invalid insert happen?
		sth below subtree root which does not fit with insert_at_node?


	*/


/* 
   this function takes a subtree (already pruned)
   inserts all its content at to

	TODO we can get much faster here (using bounds of what is contained in the subtree and where it can go)
*/	
template<typename T>
void avl_interval_tree<T>::move_subtree(node_type * from, node_type * to) {
	std::cout << "move subtree on from " << from->value << std::endl;
	assert(from->parent == NULL);
	node_type * strc = from->str_center;
	node_type * strl = from->str_left;
	node_type * strr = from->str_right;
	if(strc!=NULL) { // center first should be more efficient (less rebalancing)
		prune_subtree(strc);
		move_subtree(strc, to);
	}
	if(from->str_left!=NULL) {
		prune_subtree(strl);
		move_subtree(strl, to);
	}
	if(from->str_right!=NULL) {
		prune_subtree(strr);
		move_subtree(strr, to);
	}
	insert_at_node(to, from); 

	redo_height_up(from); // TODO
}


template<typename T>
void avl_interval_tree<T>::delete_subtree(node_type * n) {
	std::cout << " delete subtree l " << n->left << " r " << n->right << " d " << n->value << std::endl;
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
	if(n->parent==NULL) {
		assert(n==root);
		root = NULL;
	} else {
		node_type * p = n->parent;
		if(p->str_left == n) {
			p->height = 0; // if we prune left/right children, height becomes invalid. Calling function should take care of that
			p->str_left = NULL;
		} else if(p->str_center == n) {
			p->str_center = NULL;
		} else if(p->str_right == n) {
			p->height = 0;
			p->str_right = NULL;
		} else {
			assert(0);
		}
		n->parent = NULL;

		redo_height_up(p); // TODO can we get faster here?
	}
}

template<typename T>
void avl_interval_tree<T>::left_attach(node_type * parent, node_type * child) {
	assert(parent->str_left == NULL);
	assert(child->parent == NULL);
	parent->str_left = child;
	child->parent = parent;
}

template<typename T>
void avl_interval_tree<T>::right_attach(node_type * parent, node_type * child) {
	assert(parent->str_right == NULL);
	assert(child->parent == NULL);
	parent->str_right = child;
	child->parent = parent;
}

template<typename T>
void avl_interval_tree<T>::center_attach(node_type * parent, node_type * child) {
	assert(parent->str_center == NULL);
	assert(child->parent == NULL);
	parent->str_center = child;
	child->parent = parent;
}

template<typename T>
void avl_interval_tree<T>::delete_node(node_type * n) {
	assert(n->parent == NULL);
	assert(n->str_left == NULL);
	assert(n->str_right == NULL);
	assert(n->str_center == NULL);
	delete n;
}


template<typename T>	
void avl_interval_tree<T>::redo_height_and_minmax_at_node(node_type * n) {
	size_t nh = 0;
	node_type * nl = n->str_left;
	node_type * nr = n->str_right;
	node_type * nc = n->str_center;
	if(nl!=NULL) {
		size_t lh = nl->height + 1;
		if(lh>nh) nh = lh;
		if(nl->min_left < n->min_left) n->min_left = nl->min_left;
		if(nl->max_right > n->max_right) n->max_right = nl->max_right;
	}
	if(nr!=NULL) {
		size_t rh = nr->height + 1;
		if(rh>nh) nh = rh;
		if(nr->min_left < n->min_left) n->min_left = nr->min_left;
		if(nr->max_right > n->max_right) n->max_right = nr->max_right;
	}
	if(nc!=NULL) {
		if(nc->min_left < n->min_left) n->min_left = nc->min_left;
		if(nc->max_right > n->max_right) n->max_right = nc->max_right;

	}
	n->height = nh;

	std::cout << " redo node " << n->value << " new h " << n->height << " str limits l " << n->min_left << " r " << n->max_right << std::endl; 
}
/** Do all kinds of slow checks 

**/
template<typename T>
void avl_interval_tree<T>::check_tree() const {
	std::cerr << "Warning: We are running slow checks of AVL interval tree with " << num_nodes << " nodes" << std::endl;

	size_t tsize = 0;
	size_t minl = 0;
	size_t maxr = std::numeric_limits<size_t>::max();
	size_t mini = 0;
	size_t maxi = std::numeric_limits<size_t>::max();

	assert(slow_index.size() == num_nodes);

	if(root!=NULL) {
		size_t id = 0;
		check_tree_at_node(root, 0, id, 0, tsize, minl, maxr, mini, maxi);
		std::cout << " after check num_nodes " << num_nodes << " tree size " << tsize << std::endl;
		assert(num_nodes == tsize);
	} else {
		assert(num_nodes==0);
	}
	

}

template<typename T>
void avl_interval_tree<T>::check_tree_at_node(node_type * n,  const size_t & level, size_t & id, const size_t & parent_id, size_t & subtree_size, size_t min_left, size_t max_right, size_t min_in, size_t max_in) const {

	std::cout << " node " << id << " on level " << level;
	node_type * p = n->parent;
	if(p == NULL) {
		std::cout << " is root height " << n->height<< std::endl;
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
		std::cout << "child of " << parent_id << " height " << n->height << std::endl;
	}
	std::stringstream mrstr;
	std::stringstream rinstr;
	std::stringstream smrstr;
	if(max_right<std::numeric_limits<size_t>::max()) {
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
	std::cout << "l " << n->left << " r " << n->right << " data pointer " << n->value << std::endl;
	std::cout << " saved subtree size l " << n->min_left << " r " << smrstr.str();
	std::cout << " current borders l " << min_left << " r " << mrstr.str() << " all intv in subtr ovrlp with " << min_in << " " << rinstr.str()<<  std::endl;
	size_t this_id = id;


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

	// check saved subtree with bounds from parent nodes

// XXX	assert(n->min_left >= min_left);
//	assert(n->max_right <= max_right);




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
		assert(n == strl->parent);
		size_t s = 0;
		size_t minl = min_left;
		size_t maxr = max_right;
		size_t mini = min_in;
		size_t maxi = max_in;
		assert(n->left > 1); // no alignments of length 1 (if there is someting to the left of n, n cannot have left end at 1)
		size_t new_max = n->left - 1;
		if(maxr > new_max) maxr = new_max;
		check_tree_at_node(strl, level+1, ++id, this_id, s, minl, maxr, mini, maxi);
		local_size+=s;
		size_t new_height = strl->height + 1;
		if(new_height > height_check) height_check = new_height;
		if(real_subtree_min > strl->min_left) real_subtree_min = strl->min_left;
		if(real_subtree_max < strl->max_right) real_subtree_max = strl->max_right;
	}
	if(strr!=NULL) {
		assert(n == strr->parent);
		size_t s = 0;
		size_t minl = min_left;
		size_t maxr = max_right;
		size_t mini = min_in;
		size_t maxi = max_in;
		size_t new_min = n->right + 1;
		if(minl < new_min) minl = new_min;
		check_tree_at_node(strr, level+1, ++id, this_id,  s, minl, maxr, mini, maxi);
		local_size+=s;
		size_t new_height = strr->height + 1;
		if(new_height > height_check) height_check = new_height;
		if(real_subtree_min > strr->min_left) real_subtree_min = strr->min_left;
		if(real_subtree_max < strr->max_right) real_subtree_max = strr->max_right;
	}
	if(strc!=NULL) {
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
		check_tree_at_node(strc, level+1, ++id, this_id,  s, minl, maxr, mini, maxi);
		local_size+=s;	
		if(real_subtree_min > strc->min_left) real_subtree_min = strc->min_left;
		if(real_subtree_max < strc->max_right) real_subtree_max = strc->max_right;
	}

	subtree_size += 1 + local_size;	
	assert(n->height == height_check);
	std::cout << "Subtree limits of node " << this_id << " are l " << real_subtree_min << " r " << real_subtree_max << std::endl;
	std::cout << " saved limits l " << n->min_left << " r " << n->max_right << std::endl;
	assert(real_subtree_min == n->min_left);
	assert(real_subtree_max == n->max_right);

}



