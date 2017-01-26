/*
 * created by Kimberly Robasky, Q2 Solutions|EA Genomics
 *
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kraken_headers.hpp"
#include "fqmapper.hpp"
#include "quickfile.hpp"
#include <iostream>
#include <fstream>
#include <unordered_map>

using std::string;
using std::vector;

namespace kraken {

  // Basic constructor
  FQMapper::FQMapper() {
    bzero(call_name_map, sizeof(call_name_map));
    //std::cout << "hi 1" << std::endl; // xxx
  }
  
#define MAX_FASTQ_NAME_LENGTH 256
  typedef char FQ_name[MAX_FASTQ_NAME_LENGTH];
  
  // ptr points to mmap'ed existing file opened in read or read/write mode
  FQMapper::FQMapper(char *ptr, std::string prefix) {

    //std::cout << "hi 2" << std::endl; // xxx
    out_file_prefix = prefix;

    // fill in call_name_map, only open files as needed
    std::stringstream ss(ptr);
    std::string l;
    int i=0;     

    while (ptr != NULL && ptr[i] != '\0') {
      if(ptr[i] == '\n') { i++; continue; }
      vector<uint32_t> nodes;
      vector<std::string> genes;

      // gene name
      // this could be done in one or two lines, but the following code should be faster:
      FQ_name gene; bzero(gene, sizeof(gene));
      int j = 0 ;
      while(ptr != NULL && ptr[i] != '\t' && ptr[i] != '\0') {
	if(ptr[i] == '\n') { i++; continue; }
	gene[j]=ptr[i]; i++; j++;
      }
      gene[j]='\0';

      // Loop over nodes
      while(ptr != NULL && ptr[i] != '\n' && ptr[i] != '\0') {
	j=0; char nodestr[1000]; bzero(nodestr, sizeof(nodestr));

	while(ptr != NULL && ptr[i] != '\n' && ptr[i] != ',' && ptr[i] != '\0') {
	  if(ptr[i] == '\t') { i++; continue; }
	  nodestr[j]=ptr[i]; i++; j++;
	}
	nodestr[j]='\0';
	int node = atoi(nodestr);
	nodes.push_back((uint32_t)node); 
	call_name_map[node].push_back(std::string(gene));

	if(ptr[i] == ',') {i++;}

      }
    }

  }
  
  FQMapper::~FQMapper() {
    for ( auto it = fq_name_ptrs.begin(); it != fq_name_ptrs.end(); ++it ) {
      FP_vec_t *fp_v = &(it->second);

      for(int i = 0; i < (int) fp_v->size(); i++) {
	std::shared_ptr<std::ofstream> fqout = fp_v->at(i);
	fqout->close();
      }
      // delete fp_v; //xxx doesn't work
    }

  }

  void FQMapper::write(DNASequence &dna, uint32_t call, bool Paired_end) {
    int debug = 0;

    // split the sequence, if needed
    DNASequence d1, d2;
    if(Paired_end) { dna.split(&d1, &d2); }

    #pragma omp critical(fwrit_fastq)
      {
        std::vector<std::string>::const_iterator g;
        for(g=call_name_map[call].begin(); g!=call_name_map[call].end(); ++g){
          // open the file if it's not open (open _1,_2 if SPLIT is true)
          FP_vec_t *fp_v(new FP_vec_t);
    
          if(fq_name_ptrs[(*g)].empty()) {
            
            std::shared_ptr<std::ofstream> fqout(new std::ofstream);
            fqout->open(out_file_prefix + (*g) + "_1.fq");
            fp_v->push_back(fqout); 
    
            if(Paired_end) {
              std::shared_ptr<std::ofstream> fqout2(new std::ofstream);
              fqout2->open(out_file_prefix + (*g) + "_2.fq");
              fp_v->push_back(fqout2); 
            }
    
            fq_name_ptrs[(*g)]=*fp_v;
    
          } else {
            fp_v =  &fq_name_ptrs[(*g)];
          }
    
          // ...and write the FASTQ record out to the files
          std::shared_ptr<std::ofstream> fqout = fp_v->at(0);
          if(Paired_end) {
            d1.write(fqout);
//          d1.write(fqout, call, "\\1");
              
            fqout = fp_v->at(1);
            d2.write(fqout);
//          d2.write(fqout, call, "\\2");
              
          } else {
            dna.write(fqout);
//          dna.write(fqout, call, "");
          }
    
          if(debug) {std::cout<<(*g)<< " "; }
        }
        if(debug) {std::cout<< std::endl;}
      } // end of critical section
  }

  void FQMapper::print_node_map() {
    for(int j=0; j< MAX_TREE_NODES; j++) {
      if(call_name_map[j].size() > 0) {
	std::cout << j << ":";
	std::vector<std::string>::const_iterator g;
	for(g=call_name_map[j].begin(); g!=call_name_map[j].end(); ++g){
	  std::cout<<(*g)<< " ";
	}
	std::cout << std::endl;
      }
    }
  }
  
} // namespace

#ifdef FQMAPPER_TEST
  int main (int argc, char **argv) {
    using namespace std;
    using namespace kraken;

    FQMapper FQ_mapper;

    string Map_filename;
    if(argc == 1) {Map_filename = "divtest.txt";}
    else {Map_filename = argv[1];}

    std::cout << "+Reading from: ("<< Map_filename << ")\n";
    bool Populate_memory = false;

    if (Map_filename != "") {
      QuickFile map_file;
      map_file.open_file(Map_filename);
      if (Populate_memory)
	map_file.load_file();
      FQ_mapper = FQMapper(map_file.ptr());
      printf("+file loaded.\n\n+node map:\n");
      FQ_mapper.print_node_map();
    }

    char filename[256]; strncpy(filename,"sim.HLA00001.head.fq",256);
    printf("+Reading sequence from: (%s)\n",filename);
    string file_str(filename);
    DNASequenceReader *reader;
    DNASequence dna;
    reader = new FastqReader(file_str);

    printf("\n+writing a sequence:\n");
    if(reader->is_valid()) {
      dna = reader->next_sequence();
      FQ_mapper.write(dna,5); // write a sequence that's classified to a node that isn't mapped to anything (shouldn't write anything)
      FQ_mapper.write(dna,7);
    }
    printf("+Sequence written.\n");

    printf("\n+done.\n");

    return 0;
  }
#endif
