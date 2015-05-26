#include "graph.hpp"




// ---------------------------
// 		Vertex
// ---------------------------

// -------------------------------------------Functions
//Graph::Vertex::Vertex(){}
Graph::Vertex::Vertex(int id){
	idVertex = id;
	visited = false;
	previous = -1;
	costScore = 0.0; // initialize at none ? or -1 ?
	distance = std::numeric_limits<double>::infinity();
}
//Graph::Vertex::Vertex(Graph::Vertex const& v){}
//Graph::Vertex::~Vertex(){}
// -------------------------------------------Get
const int& Graph::Vertex::getIdVertex()const{return idVertex;}

//-------------
//Split a string
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
 }


//--------------------------------------------------------- Set and Get
void Graph::Vertex::setStartToStart(const int& v){
	startToStart.push_back(v);
}

void Graph::Vertex::setStartToEnd(const int& v){
	startToEnd.push_back(v);
}

void Graph::Vertex::setEndToStart(const int& v){
	endToStart.push_back(v);
}

void Graph::Vertex::setEndToEnd(const int& v){
	endToEnd.push_back(v);
}
void Graph::Vertex::setCostScore(const double& score){
	costScore = score;
}
void Graph::Vertex::setDistance(const double& d){
	distance = d;
}
const std::vector<int>& Graph::Vertex::getStartToStart(){return startToStart;}
const std::vector<int>& Graph::Vertex::getStartToEnd(){return startToEnd;}
const std::vector<int>& Graph::Vertex::getEndToStart(){return endToStart;}
const std::vector<int>& Graph::Vertex::getEndToEnd(){return endToEnd;}
const double& Graph::Vertex::getCostScore(){return costScore;}
const double& Graph::Vertex::getDistance(){return distance;}

// ---------------------------
// 		Graph
// ---------------------------

// -------------------------------------------Functions
//Graph::Graph(){}

//Graph::Graph(Graph const& g){}

//Graph::~Graph(){}
const std::map<int,Graph::Vertex>& Graph::getVertices()const{return vertices;}
void Graph::setScore(const int& name, const double& score){
	//vertices.begin()->second.setCostScore(score);
	std::map<int,Graph::Vertex>::iterator it = vertices.find(name);
	it->second.setCostScore(score);
	//std::pair<std::string,Graph::Vertex> tmp = vertices[i];

}
void Graph::readDotFile(std::string file){
	/*
	Achtung : Need spaces in the dot file
	like that : 1 -> 3 [headport=value] [tailport=value] ;
	For the moment node in .dot correspond to the index in the fasta file.
	So node are interger !
	Really important to respect the format :
	write from left to right to know the sens of the arrow
	and precise the sens of the node with [headport=value] and [tailport=value]
	replace value by e or w for "east" and "west"
	*/
	std::ifstream dotFileIn(file.c_str());
	if(!dotFileIn){
		std::cout << "Error : Cannot open " << file.c_str() << std::endl;
		exit(1);
	}
	else{
		std::string line;
		while(getline(dotFileIn,line)){
			if(line[0] == '/'){
				continue;
			}
			else{
				size_t pos = line.find("->");
				if(pos != std::string::npos){
					int origin,destination;
					std::string arrow, headbox, tailbox;
					std::istringstream iss(line);
					iss >> origin >> arrow >>destination >> headbox >> tailbox;

					std::vector<std::string> tmp1, tmp2;
					tmp1 = split(headbox, '=');
					tmp2 = split(tailbox, '=');

					//case 1 : startToStart
					if((tmp1[1].find("w") != std::string::npos) && (tmp2[1].find("w") != std::string::npos)){
						addEdge(1,origin,destination);
						continue;
					}
					//case 2 : startToEnd
					if((tmp1[1].find("w") != std::string::npos) && (tmp2[1].find("e") != std::string::npos)){
						addEdge(2,origin,destination);
						continue;
					}
					//case 3 : endToStart
					if((tmp1[1].find("e") != std::string::npos) && (tmp2[1].find("w") != std::string::npos)){
						addEdge(3,origin,destination);
						continue;
					}
					//case 4 : endToEnd
					if((tmp1[1].find("e") != std::string::npos) && (tmp2[1].find("e") != std::string::npos)){
						addEdge(4,origin,destination);
						continue;
					}
				}
			}
		}
	}
}

void Graph::addVertex(Graph::Vertex v){
	vertices.insert(std::pair<int,Graph::Vertex>(v.getIdVertex(), v));
}

void Graph::addEdge(int typeOfEdge,int id1, int id2){ // create nodes and add edges

	std::map<int,Graph::Vertex>::iterator it = vertices.find(id1);
	if(it == vertices.end()){ // If vertex does not exist
		Graph::Vertex v1= Graph::Vertex(id1);
		switch(typeOfEdge){
			case 1 :
				v1.setStartToStart(id2);
				break;
			case 2 :
				v1.setStartToEnd(id2);
				break;
			case 3 :
				v1.setEndToStart(id2);
				break;
			case 4 :
				v1.setEndToEnd(id2);
				break;
		}
		addVertex(v1);
	}
	else{ // if exist, update
		switch(typeOfEdge){
			case 1 :
				it->second.setStartToStart(id2);
				break;
			case 2 :
				it->second.setStartToEnd(id2);
				break;
			case 3 :
				it->second.setEndToStart(id2);
				break;
			case 4 :
				it->second.setEndToEnd(id2);
				break;
		}
	}
	if(vertices.find(id2)== vertices.end()){ // if v2 does not exist, creation, else do nothing
		Graph::Vertex v2 = Graph::Vertex(id2);
		addVertex(v2);
	}
}

std::ostream& operator<<(std::ostream &os, Graph & g)
{
	g.printGraph(os);
	return os;
}

void Graph::printGraph(std::ostream &os){
	os << "There is " << vertices.size() << " nodes" << std::endl;

	for (std::map<int,Graph::Vertex>::iterator it=vertices.begin(); it!=vertices.end(); ++it)
		(it->second).printVertex(os) ;
}

std::ostream& operator<<(std::ostream &os, Graph::Vertex const& v)
{
	v.printVertex(os);
	return os;
}

void Graph::Vertex::printVertex(std::ostream &os)const{
	os << " vertex : " << idVertex << " score "<< costScore << "\t";

	if(startToStart.size() != 0)
		os <<"startToStart: ";
	for(unsigned int i=0 ; i < startToStart.size(); ++i){
		 os << startToStart[i] <<"\t";
	}

	if(startToEnd.size() != 0)
		os <<"startToEnd: ";
	for(unsigned int i=0 ; i < startToEnd.size(); ++i){
		os << startToEnd[i] <<"\t";
	}

	if(endToStart.size() != 0)
		os <<"endToStart: ";
	for(unsigned int i=0 ; i < endToStart.size(); ++i){
		os << endToStart[i] <<"\t";
	}

	if(endToEnd.size() != 0)
		os <<"endToEnd: ";
	for(unsigned int i=0 ; i < endToEnd.size(); ++i){
		os << endToEnd[i] <<"\t";
	}
	os << "\n";
}

bool operator== ( const Graph::Vertex &v1, const Graph::Vertex &v2){
        return v1.getIdVertex() == v2.getIdVertex();
}

//void Graph::updateDistance(Graph::Vertex v1, Graph::Vertex v2, std::map<std::string,int>& minScore, std::map<std::string,std::string>& previous){
	/*if(minScore[v2.getIdVertex()] < minScore[v1.getIdVertex()]+ v2.getCostScore()){
					minScore[v2.getIdVertex()] = minScore[v1.getIdVertex()]+ v2.getCostScore();
					previous.insert(std::pair<std::string,std::string>(v2.getIdVertex(), v1.getIdVertex()));

			}
		*/
//}

void Graph::dijkstra(int start,int end){
	//Init
	std::vector<Graph::Vertex> notVisited;
	std::map<int,int> previous;
	std::vector<Graph::Vertex> neighbor;
    for( std::map<int,Graph::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
	    	previous.insert(std::pair<int,int>(it->first,-1)); // pair <node, previous>
	}
	//initialisation depending to the start
	notVisited.push_back(vertices.find(start)->second);
	notVisited.begin()->setDistance(0);

	while (!notVisited.empty()){
		std::sort(notVisited.begin(),notVisited.end());
		if(notVisited.begin()->getIdVertex() == end){
			break;
		}
		else {
			neighbor.clear();
			std::vector<int> tmpneighbor = notVisited.begin()->getEndToStart(); // I need to sort neighbor to find the next one with the smaller cost
			for(std::vector<int>::iterator itt2 =tmpneighbor.begin(); itt2 != tmpneighbor.end(); ++itt2 ){
				neighbor.push_back(vertices.find(*itt2)->second);
				if (notVisited.begin()->getIdVertex()==start){ //first run, so we need to update distance of neighbor
					neighbor.back().setDistance(neighbor.back().getCostScore()); //set the distance to the costScore of the node
					previous[neighbor.back().getIdVertex()] = notVisited.begin()->getIdVertex();
				}
			}
			std::sort(neighbor.begin(), neighbor.end()); // sort neighbor by distance
			for(std::vector<Graph::Vertex>::iterator it2 = neighbor.begin(); it2 != neighbor.end(); ++it2 ){ // for each neighbor, updateDistance
				//update distance
				double weight = it2->getCostScore();
				if(notVisited.begin()->getDistance()+ weight <= it2->getDistance()){
						it2->setDistance(notVisited.begin()->getDistance()+weight);
						previous[it2->getIdVertex()] = notVisited.begin()->getIdVertex();
				}
				// add node to the queue
				if ( std::find(notVisited.begin(), notVisited.end(),*it2)== notVisited.end())
					notVisited.push_back(*it2);
				std::cout <<notVisited.begin()->getIdVertex() <<" --  "<<it2->getIdVertex()<< " distance " << it2->getDistance()<< "( "<<notVisited.begin()->getDistance()<< " "<< weight << std::endl;
			//Before erase neighbor, update vertices
			vertices.find(it2->getIdVertex())->second.setDistance(it2->getDistance());
			}

		}
		notVisited.erase(notVisited.begin());
	}
	if (previous[end] == -1){
		std::cout << "not reachable "<< "start :"<< start << " end "<< end<< std::endl;
		exit(1);

	}
	else{
	// Find Path :
	std::vector<int> path;
	int node = end;
	path.insert(path.begin(),node); // add end at the path
	while(node != start){
		path.insert(path.begin(),previous[node]);
		node = previous[node];
	}
	// print Path
	std::cout << "path ";
	std::copy(path.begin(), path.end(), std::ostream_iterator<int>(std::cout, " "));
	}
}

void Graph::readAlignment(std::string fastaFile, std::string samFile){
	all_data data(fastaFile,samFile);
		// Train the model on all data
		typedef mc_model use_model;
		use_model m(data);

		ofstream outs("encode",std::ofstream::binary);
		m.train(outs);
		std::cout << "nombre alignment "<< data.numAlignments()<<std::endl;
		vector<pw_alignment> pp = data.getAlignments(); // copy of vector of pw_alignments, not good !!
		std::sort(pp.begin(),pp.end());
		for(int itP =0; itP < pp.size(); ++itP){
			std::cout << pp[itP].getbegin1()<<" id "<<pp[itP].getreference1() << " +++ "<< pp[itP].getbegin2()<< " id "<<pp[itP].getreference2() << std::endl;
			double c1;
			double c2;
			double m1;
			double m2;
			//m.cost_function(pp[itP],c1,c2,m1,m2,outs);//TODO It does not work with more than one pw_alignment :$
			//Error : terminate called after throwing an instance of 'std::out_of_range'
			 //			what():  vector<bool>::_M_range_check

			//Verify that the 5' end of the read fit to the graph
			if(pp[0].getbegin2() != 0)
				std::cout << "pb start of the read does not fit" <<std::endl; //call needleman option 2, how to get previous nodes ?
			//TODO look at all possibility. How to do that ? with a Queue ? If overlap, I take one and put the other which overlap in a queue and I look at them later ?
			// Or Segment/interval tree ?
			// verify if there is overlap or gap
			if ((itP+1) < pp.size()){ // compare locus of 2 alignments
				if(pp[itP].getend2() >= pp[itP+1].getbegin2()) // overlap so I should add to a queue one of them
					std::cout << " overlap ! pp[itP].getend2() " << pp[itP].getend2()<< " blop "<<pp[itP].getbegin2()<<" pp[itP+1].getbegin2() "<< pp[itP+1].getbegin2() << " "<<data.get_seq_name(pp[itP].getreference1()) << " " << data.get_seq_name(pp[itP+1].getreference1()) << std::endl;
				else if(pp[itP].getend2()+1 < pp[itP+1].getbegin2()){
					std::cout << " gap ! pp[itP].getend2() "<< pp[itP].getend2()<<" blop "<<pp[itP].getbegin2() <<" pp[itP+1].getbegin2() " <<pp[itP+1].getbegin2()<< data.get_seq_name(pp[itP].getreference1()) << " " << data.get_seq_name(pp[itP+1].getreference1()) << std::endl;
					// take the gap sequence of the read, take the node after pp[itP] that is different from pp[itP+1]
					// call needleman option 1
					// verify that the new node and pp[itP+1] are reachable
				}
				else
					std::cout << "everything fine  ! pp[itP].getend2() "<< pp[itP].getend2()<<" blop "<<pp[itP].getbegin2() <<" pp[itP+1].getbegin2() " <<pp[itP+1].getbegin2()<< data.get_seq_name(pp[itP].getreference1()) << " " << data.get_seq_name(pp[itP+1].getreference1()) << std::endl;
					//dijkstra(pp[itP].getreference1(),pp[itP+1].getreference1());
			}
		}

}
