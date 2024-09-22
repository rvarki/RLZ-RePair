#ifndef REPAIR_H
#define REPAIR_H

#include <string>
#include <fstream>
#include <cstdint>

uint64_t calculate_parse_bytes(std::ifstream& pfile);
void repair(const std::string ref_str, std::ifstream& pfile);


#endif  // REPAIR_H