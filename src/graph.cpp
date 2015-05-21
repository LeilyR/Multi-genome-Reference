#include "graph.hpp"




// ---------------------------
// 		Vertex
// ---------------------------

// -------------------------------------------Functions
//Graph::Vertex::Vertex(){}
Graph::Vertex::Vertex(std::string name){
	idVertex = name;
	visited = false;
	previous = -1;
	costScore = 0; // initialize at none ? or -1 ?
	distance = std::numeric_limits<double>::infinity();
}
//Graph::Vertex::Vertex(Graph::Vertex const& v){}
//Graph::Vertex::~Vertex(){}
// -------------------------------------------Get
const std::string& Graph::Vertex::getIdVertex()const{return idVertex;}

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
void Graph::Vertex::setStartToStart(const std::string& v){
	startToStart.push_back(v);
}

void Graph::Vertex::setStartToEnd(const std::string& v){
	startToEnd.push_back(v);
}

void Graph::Vertex::setEndToStart(const std::string& v){
	endToStart.push_back(v);
}

void Graph::Vertex::setEndToEnd(const std::string& v){
	endToEnd.push_back(v);
}
void Graph::Vertex::setCostScore(const int& score){
	costScore = score;
}
void Graph::Vertex::setDistance(const int& d){
	distance = d;
}
const std::vector<std::string>& Graph::Vertex::getStartToStart(){return startToStart;}
const std::vector<std::string>& Graph::Vertex::getStartToEnd(){return startToEnd;}
const std::vector<std::string>& Graph::Vertex::getEndToStart(){return endToStart;}
const std::vector<std::string>& Graph::Vertex::getEndToEnd(){return endToEnd;}
const int& Graph::Vertex::getCostScore(){return costScore;}
const int& Graph::Vertex::getDistance(){return distance;}

// ---------------------------
// 		Graph
// ---------------------------

// -------------------------------------------Functions
//Graph::Graph(){}

//Graph::Graph(Graph const& g){}

//Graph::~Graph(){}
const std::map<std::string,Graph::Vertex>& Graph::getVertices()const{return vertices;}

void Graph::setScore(const std::string& name, const double& score){
	//vertices.begin()->second.setCostScore(score);
	std::map<std::string,Graph::Vertex>::iterator it = vertices.find(name);
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
					std::string origin, destination, arrow, headbox, tailbox;
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
	vertices.insert(std::pair<std::string,Graph::Vertex>(v.getIdVertex(), v));
}

void Graph::addEdge(int typeOfEdge,std::string name1, std::string name2){ // create nodes and add edges

	std::map<std::string,Graph::Vertex>::iterator it = vertices.find(name1);
	if(it == vertices.end()){ // If vertex does not exist
		Graph::Vertex v1= Graph::Vertex(name1);
		switch(typeOfEdge){
			case 1 :
				v1.setStartToStart(name2);
				break;
			case 2 :
				v1.setStartToEnd(name2);
				break;
			case 3 :
				v1.setEndToStart(name2);
				break;
			case 4 :
				v1.setEndToEnd(name2);
				break;
		}
		addVertex(v1);
	}
	else{ // if exist, update
		switch(typeOfEdge){
			case 1 :
				it->second.setStartToStart(name2);
				break;
			case 2 :
				it->second.setStartToEnd(name2);
				break;
			case 3 :
				it->second.setEndToStart(name2);
				break;
			case 4 :
				it->second.setEndToEnd(name2);
				break;
		}
	}
	if(vertices.find(name2)== vertices.end()){ // if v2 does not exist, creation, else do nothing
		Graph::Vertex v2 = Graph::Vertex(name2);
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

	for (std::map<std::string,Graph::Vertex>::iterator it=vertices.begin(); it!=vertices.end(); ++it)
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
void Graph::updateDistance(Graph::Vertex v1, Graph::Vertex v2, std::map<std::string,int>& minScore, std::map<std::string,std::string>& previous){
	/*if(minScore[v2.getIdVertex()] < minScore[v1.getIdVertex()]+ v2.getCostScore()){
					minScore[v2.getIdVertex()] = minScore[v1.getIdVertex()]+ v2.getCostScore();
					previous.insert(std::pair<std::string,std::string>(v2.getIdVertex(), v1.getIdVertex()));

			}
		*/
}

void Graph::dijkstra(std::string start,std::string end){
	//Init
	std::vector<Graph::Vertex> notVisited;
	std::map<std::string,std::string> previous;
	std::vector<Graph::Vertex> neighbor;
    for( std::map<std::string,Graph::Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it ) {
	    	previous.insert(std::pair<std::string,std::string>(it->first,"na")); // pair <node, previous>
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
			std::vector<std::string> tmpneighbor = notVisited.begin()->getEndToStart(); // I need to sort neighbor to find the next one with the smaller cost
			for(std::vector<std::string>::iterator itt2 =tmpneighbor.begin(); itt2 != tmpneighbor.end(); ++itt2 ){
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
	if (previous[end] == "na"){
		std::cout << "no path "<< std::endl;
		exit(1);
	}
	else{
	// Find Path :
	std::vector<std::string> path;
	std::string node = end;
	path.insert(path.begin(),node); // add end at the path
	while(node != start){
		path.insert(path.begin(),previous[node]);
		node = previous[node];
	}
	// print Path
	std::cout << "path ";
	std::copy(path.begin(), path.end(), std::ostream_iterator<std::string>(std::cout, " "));
	}
}

void Graph::readAlignment(std::string fastaFile, std::string samFile){
	all_data data(fastaFile,samFile);
	//all_data data("/ebio/abt6_projects7/small_projects/mdubarry/Documents/SampleProgram/bin/output/tmpOut.fasta", "/ebio/abt6_projects7/small_projects/mdubarry/Documents/SampleProgram/bin/output/referenceRead.sam");
		// Train the model on all data
		typedef mc_model use_model;
		use_model m(data);

		ofstream outs("encode",std::ofstream::binary);
		m.train(outs);
		std::cout << "nombre alignment "<< data.numAlignments()<<std::endl;
		for(size_t i = 0 ; i < data.numAlignments(); i++){
			const pw_alignment & p = data.getAlignment(i);
			//std::cout << p.getbegin(0)<< " +++ "<< p.getbegin(1) << std::endl;
			//std::cout << p.getend(0) << " ++ " << p.getend(1) << std::endl;
			double c1;
			double c2;
			double m1;
			double m2;

			m.cost_function(p,c1,c2,m1,m2,outs);
			//cout << data.get_seq_name(p.getreference1())<< " c1 : "<< c1 << " m1 : "<< m1 << endl;
			//cout << data.get_seq_name(p.getreference2())<< " c2 : "<< c2 << " m2 : "<< m2 << endl;
			//char buffer[50];

			//sprintf(buffer,"%d",(i+1));
			//std::string name = buffer;
			//TODO pour verifier si il y a gap/overlap trier les alignments suivant leur place sur read
		/*	if ((i+1) < data.numAlignments()){ // compare locus of 2 alignments
				const pw_alignment & p2 = data.getAlignment(i+1);
				std::cout <<"i "<<i <<" data.numAlignments "<< data.numAlignments()<<" " <<p.getend(0) <<" " << p2.getbegin(0) << std::endl;
				if(p.getend(0) >= p2.getbegin(0))
					std::cout << " overlap ! " << data.get_seq_name(p.getreference1()) << " " << data.get_seq_name(p2.getreference1()) << std::endl;
				else if(p.getend(0)+1 < p2.getbegin(0))
					std::cout << " gap ! " << data.get_seq_name(p.getreference1()) << " " << data.get_seq_name(p2.getreference1()) << std::endl;
				else

			}
		 */
			std::cout << "everything fine in locus "<< std::endl;
		}
}
/*
void Graph::buildNewGraph(){

}
*/
