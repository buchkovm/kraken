/*
 * Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
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

#include <iostream>
#include "kraken_headers.hpp"

#include "seqreader.hpp"


using namespace std;

namespace kraken {


  bool DNASequence::merge(DNASequence d) {
    std::string id1=id;
    std::string id2=d.id;
    id1.resize(id1.length()-2);
    id2.resize(id2.length()-2);

    if(id1 != id2) { 
      return false; 
    }

    seq = seq + "N" + d.seq;
    quals = quals + (char) 30 + d.quals;
    return true;
  }

  void DNASequence::split(DNASequence *d1, DNASequence *d2) {
    uint32_t split_pos = quals.find((char) 30);
    if(split_pos == (uint32_t)std::string::npos) { errx(EX_USAGE, "CODE ERROR: request to split read that was never merged");};
    d1->id = id; d1->id.resize(id.length()-2); // cut off '\1' on the end
    d2->id = id; d2->id.resize(id.length()-2);
    d1->header_line = header_line;
    d2->header_line = header_line;

    d1->quals = quals.substr(0, (size_t)split_pos);
    d2->quals = quals.substr(split_pos + 1, quals.size());
    d1->seq = seq.substr(0, (size_t)split_pos);
    d2->seq = seq.substr(split_pos + 1, seq.size());
  }
  
  void DNASequence::write( std::shared_ptr<std::ofstream> fqout, uint32_t call, std::string rname_suffix) {
    *fqout << "@" << id << "(" << call << ")" << rname_suffix << std::endl; 
    *fqout << seq << std::endl;
    *fqout << "+" << std::endl;
    *fqout << quals << std::endl;
  }


  FastaReader::FastaReader(string filename) {
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }

  DNASequence FastaReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }

    string line;

    if (linebuffer.empty()) {
      getline(file, line); 
    }
    else {
      line = linebuffer;
      linebuffer.clear();
    }

    if (line[0] != '>') {
      warnx("malformed fasta file - expected header char > not found");
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    istringstream seq_id(dna.header_line);
    seq_id >> dna.id;

    ostringstream seq_ss;

    while (file.good()) {
      getline(file, line);
      if (line[0] == '>') {
	linebuffer = line;
	break;
      }
      else {
	seq_ss << line;
      }
    }

    dna.seq = seq_ss.str();

    if (dna.seq.empty()) {
      warnx("malformed fasta file - zero-length record (%s)", dna.id.c_str());
      valid = false;
      return dna;
    }

    return dna;
  }

  bool FastaReader::is_valid() {
    return valid;
  }

  FastqReader::FastqReader(string filename) {

    // is filename a gz file?
    if(filename.substr( filename.length() - 3 ) == ".gz") {
      gzfile.open(filename.c_str());
      std::istream& s = file;
      s.rdbuf(gzfile.rdbuf());
    }
    file.open(filename.c_str());
    if (file.rdstate() & ifstream::failbit) {
      err(EX_NOINPUT, "can't open %s", filename.c_str());
    }
    valid = true;
  }

  DNASequence FastqReader::next_sequence() {
    DNASequence dna;

    if (! valid || ! file.good()) {
      valid = false;
      return dna;
    }

    string line;
    std::getline(file, line); 

    if (line.empty()) {
      valid = false;  // Sometimes FASTQ files have empty last lines
      return dna;
    }
    if (line[0] != '@') {
      if (line[0] != '\r')
        warnx("malformed fastq file - sequence header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    dna.header_line = line.substr(1);
    istringstream line_ss(dna.header_line);
    
    line_ss >> dna.id;

    getline(file, dna.seq);

    getline(file, line);
    if (line.empty() || line[0] != '+') {
      if (line[0] != '\r')
        warnx("malformed fastq file - quality header (%s)", line.c_str());
      valid = false;
      return dna;
    }
    getline(file, dna.quals); 

    return dna;
  }

  bool FastqReader::is_valid() {
    return valid;
  }

  FastqReader::~FastqReader() {
    file.close();
  }

} // namespace
