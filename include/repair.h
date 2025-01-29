#ifndef REPAIR_H
#define REPAIR_H

#include <string>
#include <fstream>
#include <cstdint>

extern "C" {
    #include "heap.h"
    #include "array.h"
    #include "basics.h"
    #include "records.h"
    #include "hash.h"
}

uint64_t calculate_parse_bytes(std::ifstream& pfile);
void print_all_records();
void print_record(const std::string message, const Trecord* orec);
void print_ref();
void print_phrase_lst();
void print_hash_ranges();
void prepareRef(std::vector<unsigned char>& rtext);
void createMaxHeap(std::ifstream& pfile);
void populatePhrases(std::ifstream& pfile);
void repair(std::ofstream& R, std::ofstream& C);
void phraseBoundaries(int left_elem, int right_elem);
void sourceBoundaries(int left_elem, int right_elem);


#endif  // REPAIR_H