#include "data.hpp"
#include "pw_alignment.hpp"
#include "graph.hpp"

#define MAXGAP 400
//This class is used to read the dot file includes reference graph
class ref_graph{//XXX For the moment I use Marion's graph , if I face any further issue I might need to add new functions to this class to make the graph here.
	public:
	ref_graph(Graph & g, const all_data & d):graph(g), data(d){}
	~ref_graph(){}
	void read_dot_file(std::string &, std::string &);
	void add_adjacencies(std::string & , std::string & , std::string & , std::string & );
	void add_adjacencies(std::string & , std::string & );
	const std::map<int , std::set<int> > & get_adjacencies()const;
	private:
	const all_data & data;
	Graph & graph;
	std::map<int, std::set<int> > adjacencies;

};

class deep_first{//It runs deep first algorithm on reference graph to retrun parts of graph that the distance between the first and the last node is less than MAXGAP.
	public:
	deep_first(const all_data & d, std::map<int, std::set<int> > & adjacencies):data(d){
		this -> adjacencies = adjacencies;
	//	for(std::map<int, std::set<int> >::iterator it = adjacencies.begin(); it != adjacencies.end(); it++){
	//		std::cout << it->first << std::endl;
	//	}

	}
	~deep_first(){}
	size_t seq_length(std::string &, std::string &);
	std::string seqname(int & );
	void deep_first_search(int &, std::string & );
	void look_for_neighbors(int & , std::map<int,bool> &  , std::string &, int & accu_length, std::vector<int> & apath, size_t & );
	const vector<vector<int> > get_paths()const;
	const std::set<int> get_nodes()const;
	private:
	const all_data & data;
	std::map<int, std::set<int> > adjacencies;
	vector<vector<int> > paths;
	std::set<int> nodes_on_paths;
	std::map<int, size_t> path_length;
	std::multimap<int,int> parental_relation; //first one shows the child node ,second one shows the parent

};

//This class is used to make several components on a read where each component contains aligned regions with less than MAX_GAP base distance
class als_components{

	public:
	als_components(const all_data & d , deep_first & df,const std::vector<std::set<const pw_alignment* , compare_pointer_pw_alignment> >  & ccs):data(d), dfirst(df){
		alignments_left.resize(ccs.size());
		alignments_right.resize(ccs.size());

		for(size_t i = 0 ; i < ccs.size(); i++){
			std::cout << "cc at " << i << " : "<<std::endl;
			for(std::set<const pw_alignment* , compare_pointer_pw_alignment>::iterator it = ccs.at(i).begin() ; it != ccs.at(i).end() ; it++){
				const pw_alignment * al = *it;
				alignments_left.at(i).insert(*al);
				alignments_right.at(i).insert(*al);

				al->print();

			}
			std::set<pw_alignment , sort_pw_alignment>::iterator first = alignments_left.at(i).begin();
			pw_alignment p = *first;
			size_t l2,r2;
			p.get_lr2(l2,r2);
			size_t lower = l2;
			std::cout << "lower " << lower <<std::endl;//The least left
			std::set<pw_alignment , sort_right_pw_alignment>::iterator last = alignments_right.at(i).end(); 
			--last;
			p = *last;
			p.get_lr2(l2,r2);
			size_t upper = r2;
			std::cout << "upper "<< upper<<std::endl;
			boundries.insert(std::make_pair(lower,upper));
			
			
		}
	}
	~als_components(){}
	void merg_components(){//TODO change it with join interval function in avl tree to make it faster
		std::set<pw_alignment,sort_pw_alignment> temp;
		size_t previous_right, current_left;
		size_t counter = 0;
		for(std::map<size_t, size_t>::iterator it = boundries.begin(); it != boundries.end();it++){
			std::cout << it->first <<  " "<< it->second<<std::endl;
			if(it == boundries.begin()){
				previous_right = it->second;
				temp = alignments_left.at(counter);
				counter++;
			}else{
				current_left = it->first;
				if(current_left - previous_right < MAXGAP){
					temp.insert(alignments_left.at(counter).begin(), alignments_left.at(counter).end());
					previous_right = it->second;
					counter++;
				}else{
					std::cout<<"members of a new component"<<std::endl;
					for(std::set<pw_alignment,sort_pw_alignment>::iterator it1 = temp.begin() ; it1 != temp.end();it1++){
						pw_alignment pw = *it1;
						pw.print();
					}
					alignments.push_back(temp);
					temp = alignments_left.at(counter);
					previous_right = it->second;					
					counter++;
				}
			}
		}
		if(temp.size() !=0){
			alignments.push_back(temp);
		}
		std::cout << "als size " << alignments.size()<<std::endl;
	}
	/*GO over the alignments, for each al find all the possible paths on the graph that starts with that alignment and has <=MAXGAP length,
	then find the alignments that their left is closer than MAXGAP bases to the current one.Check them on the graph. If they are on those 
	paths that we already picked then we make new als using needleman wunsch algorithm.
	For those which are not on those paths we check if their left is smaller than beginning of the current alignment component- MAXGAP/2, 
	if so then we try to find the graph paths for it and go further otherwise discard it.*/
	void find_als_on_paths();
	void get_paths(const pw_alignment &, std::vector<std::vector<int> > &);
	void look_for_neighbors_on_paths(std::set<int>& , std::set<pw_alignment,sort_pw_alignment> &);
	private:
		const all_data & data;
		deep_first & dfirst;
		std::vector<std::set<pw_alignment , sort_pw_alignment> >alignments_left;
		std::vector<std::set<pw_alignment , sort_right_pw_alignment> >alignments_right;
		std::map<size_t,size_t> boundries;
		std::vector<std::set<pw_alignment,sort_pw_alignment> > alignments;//Includes components of alignments that have more than MAXGAP distance with each other.



};

