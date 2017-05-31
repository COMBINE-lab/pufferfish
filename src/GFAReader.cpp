#include <iostream>
#include <set>
#include "../gfakluge/src/gfakluge.hpp"

bool isdigit2(const string s){
		//const char* s_str = s.c_str();
		string::const_iterator it;
		for (it = s.begin(); it != s.end(); it++){
			if (!std::isdigit(*it)){
				return false;
		   }
	   }
		return true;
	}

int main(int argc, char* argv[]){

    struct config_elem{
        int tname;
        int pos;
        bool ore;
    };

    std::string config_file = argv[1] ;

    std::ofstream configtable("./config_table.text") ;
    std::ofstream tnameTable("./tnames.text");

    gfak::GFAKluge gg;
    gg.parse_gfa_file(config_file) ;

    auto header = gg.get_header() ;
    auto name_to_seq = gg.get_name_to_seq() ;
    auto name_to_path = gg.get_name_to_path() ;

    std::vector<std::string> seqVec ;
    std::vector<int> contiglens ;
    std::map<std::string,size_t> idMap ;
    std::set<std::string> tnames ;
    std::vector<std::string> tnamesVec ;
    std::map<std::string,int> tnameMap ;

    std::map<size_t, std::vector<config_elem> > configMap ;

	//
	int lines = 0;
    size_t seqid = 0;


    std::cout << "\n Number of segments " << name_to_seq.size() ;
    std::cout << "\n Number of paths " << name_to_path.size() ;

    for(auto spair : name_to_seq){
        auto id = spair.first ;
        auto seqelem = spair.second ;


        if(isdigit2(id)){
            //std::cout << id << " " << seqelem.id << "\n" ;\t
            //break ;
            seqVec.push_back(seqelem.sequence) ;
            contiglens.push_back(seqelem.sequence.length()) ;
            idMap[id] = seqid ;
            ++seqid ;
        }
    }

    int tid = 0;

    for(auto ppair : name_to_path){
        auto tname = ppair.first ;
        auto path = ppair.second ;
        auto segnames = path.segment_names ;
        auto orientations = path.orientations ;

        if(tnames.find(tname) == tnames.end()){
            tnames.insert(tname);
            tnamesVec.push_back(tname) ;
            tnameMap[tname] = tid;
            ++tid ;
        }

        if(tname == "ENST00000335137.3|ENSG00000186092.4|OTTHUMG00000001094.2|OTTHUMT00000003223.2|OR4F5-001|OR4F5|918|CDS:1-918|"){
            std::cout << segnames.size() << " " << orientations.size() << "\n";
            for(size_t i = 0; i < segnames.size() ; ++i){
                std::cout << segnames[i] << " \t " << orientations[i] << "\n" ;
            }
            //break ;

        }

        int pos = 0 ;
        for(int i = 0; i < segnames.size() ; ++i){
            auto contig_id = idMap[segnames[i]] ;
            auto curr_len = contiglens[contig_id] ;

            configMap[contig_id].push_back({tnameMap[tname],(pos-i*31),orientations[i]}) ;
            pos += curr_len ;
        }
    }

    for(auto con : configMap){
        //configtable << con.first << "\t" ;
        for(auto celm : con.second){
            configtable << celm.tname << "," << celm.pos << "," << celm.ore << "\t";
        }
        configtable << "\n" ;
    }

    for(auto tname: tnamesVec){
        tnameTable << tname << "\n" ;
    }

    configtable.close() ;
    tnameTable.close() ;

    return 0 ;
}
