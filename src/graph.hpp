#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <typeinfo>
#include <iterator>

#include "data.hpp"
#include "pw_alignment.hpp"
#include "needleman.hpp"

// Global var
const bool VERBOSE = true;
// ----------------------------
//		Graph
//-----------------------------
class Graph{

private:
	// ----------------------------
	//		Vertex
	//-----------------------------
	class Vertex{

	public:
		//Vertex();
		Vertex(int idVertex);
		//Vertex(Vertex const& v);
		//~Vertex();
		const int& getIdVertex()const;
		friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
		void printVertex(std::ostream &os)const;
		void setStartToStart(const int& id);
		void setStartToEnd(const int& id);
		void setEndToStart(const int& id);
		void setEndToEnd(const int& id);
		void setCostScore(const double& score);
		void setDistance(const double& distance);
		void setPrevious(const int& idPrevious);
		void setStartOnRead(const int& start);
		void setEndOnRead(const int& end);
		void setName(const std::string& name);
		const std::vector<int>& getStartToStart();
		const std::vector<int>& getStartToEnd();
		const std::vector<int>& getEndToStart();
		const std::vector<int>& getEndToEnd();
		const std::vector<int>& getPrevious();
		const double& getCostScore();
		const double& getDistance();
		const int& getStartOnRead();
		bool operator < (const Graph::Vertex& v)const
				    {
				        return (distance < v.distance);
				    }

	private:
		int idVertex;
		std::string name;
		bool visited;
		std::vector<int> previous;
		double costScore;
		double distance;
		std::vector<int> startToStart;
		std::vector<int> startToEnd;
		std::vector<int> endToStart;
		std::vector<int> endToEnd;
		int startOnRead;
		int endOnRead;

	};

public:
	//Graph();
	//Graph(Graph const& g);
	//~Graph();
	void setStart(const int&id, const int&start);
	void setEnd(const int&id, const int&end);
	void readDotFile(std::string dotFile,std::map<std::string, size_t> longname2seqidx); // Graph from a dot file -> Transform to constructor ?
	void addVertex(Vertex v);
	void addEdge(int typeOfEdge,int id1,std::string name1, int id2, std::string name2);
	const std::map<int,Vertex>& getVertices()const;
	void setScore(const int& id, const double& score);
	const double& getScore(const int& id);
	const Vertex& getVertex(const int& id);
	friend std::ostream& operator<<(std::ostream &os, Graph & g);
	friend std::ostream& operator<<(std::ostream &os, Vertex const& v);
	void printGraph(std::ostream &os);
	void updateDistance(Graph::Vertex& v1, Graph::Vertex& v2, std::map<int,int>& previous);
	std::vector<int> dijkstra(int start,int end);
	void initAllCostScore(all_data data, Graph& newGraph);
	void parseData(all_data data, Graph& newGraph);
	friend bool operator== ( const Vertex &v1, const Vertex &v2);
	std::string extractPartOfSeq(dnastring seq, int start, int end);
	void findFinalPath();


private:
	std::map<int,Vertex> vertices;
};

#endif // GRAPH_HPP


