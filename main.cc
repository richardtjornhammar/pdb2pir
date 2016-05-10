#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "argpa.hh"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

void
fatal(void) {
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" << std::endl;
	std::cout << "PLEASE CONTACT " << author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

int main (int argc, char *argv[] ) {
	int verbose = 0;

//	PREPARING ARGPARSER
	std::pair < int, char ** >	margs;
	margs.first	= argc;
	margs.second	= argv;

	std::vector < std::pair < int, std::pair<std::string, std::string > > >	opts_n_defaults;
	std::pair	< int, std::pair<std::string,std::string > > 		vipss;
	vipss.first=1;	vipss.second.first = "-rfile";
			vipss.second.second = "seq_l.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-pfile";
			vipss.second.second = "seq_s.pir";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-nbins";
			vipss.second.second = "24";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-verbose";
			vipss.second.second = "0";
	opts_n_defaults.push_back(vipss);

//	EXECUTING ARGPARSE
	rich::arg_parser parser;
	std::pair<int, std::vector< std::string > > run_args = parser.parse( margs, opts_n_defaults );
	if( run_args.first ) {
		fatal();
	}

	const char *vs1[]	= {	"---","ala","asn","asp","arg","cys","gln",
					"glu","gly","his","ile","leu","lys","met",
					"pro","phe","ser","thr","trp","tyr","val",
					"unk","ALA","ASN","ASP","ARG","CYS","GLN",
					"GLU","GLY","HIS","ILE","LEU","LYS","MET",
					"PRO","PHE","SER","THR","TRP","TYR","VAL",
					"UNK" };
	int Ncodes		= 43;
	const char *vc1		= "-andrcqeghilkmpfstwyvxANDRCQEGHILKMPFSTWYVX";

	std::ifstream inputf;
	inputf.open(run_args.second[0].c_str());
	if (inputf.is_open()) {
/*
 *		std::cout << "====================" << std::endl;
 *		std::cout << "   FOUND SEQUENCE   " << std::endl;
 *		std::cout << "====================" << std::endl;
 */
		std::string line;
		int rcnt=0;
		int il=0, ncnt=10, ccnt=0;
		std::vector< std::string > vs;
		std::string chn_str;
		while ( getline ( inputf,line ) ) {
			std::string search_for("ATOM");
			std::size_t found = line.find(search_for);
			if ( found != std::string::npos ) {

				std::string search_ca("CA");
				std::size_t found2 = line.find(search_ca);
				if( found2 == std::string::npos)
					continue;

				std::vector< std::string > wrds;
				boost::split( wrds, line, boost::is_any_of(" ") );
				int bHit = 0;
				int I = 0, J = 0;
				char A = line[21];
				int ccct = (int) (A-'A') ;
				for( int i=0; i<wrds.size() ; i++ )  {
					for(int j=0; j<Ncodes ; j++) {
						if(bHit==1)
							continue;
						std::string tmp(vs1[j]);
						std::size_t found3 = wrds[i].find(tmp);
						if(found3 != std::string::npos) {
							if( (wrds[i].size() >3 && wrds[i][0]=='A')
							 || (wrds[i].size()==3) ) {
								bHit=1;
								I=i;
								J=j;
							}
						}
					}
				}

				if( bHit == 1 ) {
					rcnt++;
					if(ccct-ccnt==1){
						// std::cout << std::endl << std::endl;
						ccnt=ccct;
						il=0;
						vs.push_back(chn_str);
						chn_str.clear();
					}
					if( il >= ncnt ) {
						// std::cout << " ";
						chn_str.append(" ");
						il=1;
					} else {
						il++;
					}
					// std::cout << vc1[J];
					std::stringstream ss;
					ss << vc1[J];
					chn_str.append(ss.str());
				}
			}
		}
		vs.push_back(chn_str);
                inputf.close();
//		std::cout << "====================" << std::endl;
//		FINAL IO
        	std::ofstream outputf;
        	outputf.open(run_args.second[1].c_str());
		if(outputf.is_open())
			outputf << "> ; multiplicity=" << 1 << "; \t size=" << rcnt << std::endl;
		std::cout << std::endl << std::endl <<  "> ; multiplicity=" << 1 << "; \t size=" << rcnt << std::endl;
		for(int i=0;i<vs.size();i++)
			std::cout << vs[i] << std::endl << std::endl;
		if(outputf.is_open())
			for(int i=0;i<vs.size();i++)
				outputf << vs[i] << std::endl << std::endl;
	}

	return 0;
}
