/*
* RLZ-RePair - RePair compression using RLZ parse
* Copyright (C) 2025-current Rahul Varki
* Licensed under the GNU General Public License v3 or later.
* See the LICENSE file or <https://www.gnu.org/licenses/> for details.
*/

#ifndef REPAIR_H
#define REPAIR_H

#include <string>
#include <fstream>
#include <cstdint>
#include <node.h>

extern "C" {
    #include "heap.h"
    #include "array.h"
    #include "basics.h"
    #include "records.h"
    #include "hash.h"
}

// Debug functiuons
std::string printSymbol(int elem);
void printAllRecords();
void printRecord(const std::string message, const Trecord* orec);
void printMaxPair(int new_symbol, const Trecord* orec);
void printRef();
void printHashRanges();
void printPhrase(PhraseNode* curr_phrase);
void printSupposedPhrase(PhraseNode* curr_phrase, std::list<int>::iterator leftIt);
int invalidSameCharPair(int letter);
bool checkHeap();
int checkPhraseSize();
void printPhraseList();
bool checkExpPairs();
bool checkPhraseBoundaries();
bool checkSourceBoundaries();
double calculateMemoryUsage(size_t curr_mem, size_t prev_mem);

// Necessary functions
uint64_t calculateParseBytes(std::ifstream& pfile);
void prepareRef(std::vector<unsigned char>& rtext);
void createMaxHeap(std::ifstream& pfile);
void populatePhrases(std::ifstream& pfile, int min_threshold);
void buildIntervalTree();
void updateExpPairs(PhraseNode* p, std::list<int>::iterator it, bool leftInsert);
void updateMergeExpPairs(PhraseNode* p, std::list<int>::iterator it);
void reassignExpPairs(PhraseNode* origPhrase, PhraseNode* newPhrase);
void mergeConsecutiveExpPhrases(PhraseNode* curr_phrase, PhraseNode* next_phrase);
void phraseBoundaries(int left_elem, int right_elem);
void sourceBoundaries(int left_elem, int right_elem);
void decreaseFrequency(int left_elem, int right_elem);
void increaseFrequency(int left_elem, int right_elem);
void repair(std::ofstream& R, std::ofstream& C, int num_iter);


#endif  // REPAIR_H