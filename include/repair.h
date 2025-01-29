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

uint64_t calculateParseBytes(std::ifstream& pfile);
std::string printSymbol(int elem);
void printAllRecords();
void printRecord(const std::string message, const Trecord* orec);
void printMaxPair(int new_symbol, const Trecord* orec);
void printRef();
void printPhraseList();
void printHashRanges();
void prepareRef(std::vector<unsigned char>& rtext);
void createMaxHeap(std::ifstream& pfile);
void populatePhrases(std::ifstream& pfile);
void repair(std::ofstream& R, std::ofstream& C);
void phraseBoundaries(int left_elem, int right_elem);
void sourceBoundaries(int left_elem, int right_elem);
void decreaseFrequency(int left_elem, int right_elem);
void increaseFrequency(int left_elem, int right_elem);


#endif  // REPAIR_H