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
void print_record(std::string message, Trecord* orec);
void repair(const std::string ref_str, std::ifstream& pfile);


#endif  // REPAIR_H