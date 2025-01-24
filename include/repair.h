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
void print_ref(const std::vector<unsigned char>& rvec);
void print_hash_ranges();
void prepareRef(std::vector<unsigned char>& rvec, int* chars, char* map, int size);
void createMaxHeap(const std::vector<unsigned char>& rvec, std::ifstream& pfile);


#endif  // REPAIR_H