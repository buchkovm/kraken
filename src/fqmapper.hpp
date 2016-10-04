/*
 * Copyright 2016, Kimberly Robasky <krobasky@gmail.com>
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

#ifndef FQMAPPER_HPP
#define FQMAPPER_HPP

#include "kraken_headers.hpp"
#include "seqreader.hpp"

#include <vector>
#include <memory>
#include <unordered_map>

namespace kraken {

  class FQMapper {

    public:

    
    FQMapper();
    FQMapper(char *ptr, std::string prefix); // ptr points to mmap'ed existing file opened in read or read/write mode

    ~FQMapper(); //closes map file, flushes fastq's

    void write(DNASequence &dna, uint32_t call, bool Paired_end);
    void print_node_map(); // for debugging

    private:
    std::string out_file_prefix;

    // xxx it would be better to malloc these

    //#define MAX_TREE_NODES 349200
#define MAX_TREE_NODES 20000000
    typedef std::vector< std::shared_ptr<std::ofstream> > FP_vec_t; // one vector of fq_name's per call (calls are integers)
    std::unordered_map<std::string, FP_vec_t> fq_name_ptrs;

    std::vector<std::string> call_name_map[MAX_TREE_NODES]; // one vector of fq_name's per call (calls are integers)
    // can a node ever be mapped to more than one fastq? If not, this can be an array of strings instead of an array of vectors, and performance will improve

  };
}

#endif
