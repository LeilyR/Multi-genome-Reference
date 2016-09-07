#ifndef DALIGNER_HPP
#define DALIGNER_HPP

#include "data.hpp"
#include "pw_alignment.hpp"


class daligner{

	public:
		daligner(){}
		~daligner(){}
		void read_graph(std::string & , std::ostream &);
		void make_dbread(std::string & single_read, size_t & read_id);


	private:






};

#endif
