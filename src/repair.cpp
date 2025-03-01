#include <CLI11.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <functional>
#include <utility>
#include <stack>
#include <cassert>
#include <vector> 
#include <deque>
#include "repair.h"
#include "rbintervaltree.h"
#include "doublelinkedlist.h"
#include <chrono>

extern "C" {
    #include "heap.h"
    #include "array.h"
    #include "basics.h"
    #include "records.h"
    #include "hash.h"
}

int verbosity = 0; // The verbosity level set by the user.

extern "C" float factor = 0.75f;

int minsize = 256; // Not too sure of reason but needed for BigRePair data structure
uint64_t psize; // Number of bytes that RLZ seq parse encodes for

Thash Hash; // hash table of pairs
Theap Heap; // special heap of pairs
Trarray Rec; // records
Tpair pair; // pair for freq
int chars[256]; // Indicates what chars are present in both the ref and seq parse.
char map[256]; // How to map back to the original chars.
int alpha; // Number of characters used prior to RePair
int n; // Technically |R| n - alpha gives number of rules
int c = 0; // The number of chars encoded in .C file 
int oid; // Max heap id

// Stores the int version of the reference.
RefLinkedList rlist;

// Array containing pointers to reference nodes for O(1) access to certain positions when needed.
std::vector<RefNode*> rarray;

// Boost hash function
template <typename T>
void hash_combine(std::size_t& seed, const T& value) {
    std::hash<T> hasher;
    seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Hash function for std::pair<int, int>
struct pair_int_hash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        std::size_t seed = 0;
        hash_combine(seed, p.first);  // Combine the hash of the first element
        hash_combine(seed, p.second); // Combine the hash of the second element
        return seed;
    }
};

// Hash table containing the reference ranges of pairs (bi-grams). 
// The vector contains the left endpoint of the range corresponding to the pair in the bi-gram since the range is left endpoint to left endpoint + 1
std::unordered_map<std::pair<int, int>, std::deque<int>, pair_int_hash> hash_ranges;

// Hash table containing info about the pairs within the explicit phrases
struct ExpPair
{
    PhraseNode* exp_phrase;
    std::list<int>::iterator left; 
    std::list<int>::iterator right;
    ExpPair(PhraseNode* phrase, std::list<int>::iterator l, std::list<int>::iterator r) : exp_phrase(phrase), left(l), right(r) {}
};

struct ExpPairHash {
    std::size_t operator()(const ExpPair& exp) const {
        // Combine the hash of the exp_phrase pointer, left, and right
        std::size_t h1 = std::hash<PhraseNode*>{}(exp.exp_phrase);
        std::size_t h2 = std::hash<int*>{}(&(*exp.left));  
        std::size_t h3 = std::hash<int*>{}(&(*exp.right));
        return h1 ^ (h2 << 1) ^ (h3 << 2);  // Combine the hashes
    }
};

struct ExpPairEqual {
    bool operator()(const ExpPair& lhs, const ExpPair& rhs) const {
        return lhs.exp_phrase == rhs.exp_phrase && lhs.left == rhs.left && lhs.right == rhs.right;
    }
};

std::unordered_map<std::pair<int, int>, std::unordered_set<ExpPair, ExpPairHash, ExpPairEqual>, pair_int_hash> exp_pairs; 

// Hash table containing info about phrase boundaries (store the left phrase of the pair). 
std::unordered_map<std::pair<int, int>, std::unordered_set<PhraseNode*>, pair_int_hash> pbound_pairs;

// Hash table containing info about source boundaries.
std::unordered_map<int, std::unordered_set<PhraseNode*>> start_hash;
std::unordered_map<int, std::unordered_set<PhraseNode*>> end_hash;

// List of explicit and non explicit phrases
PhraseLinkedList plist;

// Interval Tree for the non-explicit phrases
RBIntervalTree<PhraseNode*> phrase_tree;

// Timing information
std::chrono::duration<double> prepare_ref_time{0.0};
std::chrono::duration<double> calculate_size_time{0.0};
std::chrono::duration<double> populate_phrase_time{0.0};
std::chrono::duration<double> create_heap_time{0.0};
std::chrono::duration<double> phrase_boundary_time{0.0};
std::chrono::duration<double> source_boundary_time{0.0};
std::chrono::duration<double> build_interval_time{0.0};
std::chrono::duration<double> update_interval_time{0.0};
std::chrono::duration<double> nonexplicit_phrase_time{0.0};
std::chrono::duration<double> explicit_phrase_time{0.0};
std::chrono::duration<double> hash_range_time{0.0};
std::chrono::duration<double> update_bound_hash_time{0.0};
std::chrono::duration<double> total_time{0.0};

/**
 * @brief Calculates the number of bytes encoded in the RLZ parse.
 * 
 * @param[in] pfile [std::ifstream&] the input file stream of the RLZ parse.
 * @return the number of bytes encoded by the RLZ parse.
 */
uint64_t calculateParseBytes(std::ifstream& pfile)
{
    int count = 1;
    uint64_t seq_orig_size = 0;
    uint64_t len;
    // Skip the first uint64_t bytes
    pfile.seekg(sizeof(uint64_t), std::ios::beg);
    while (pfile.read(reinterpret_cast<char*>(&len), sizeof(uint64_t))) {
        // Every second uint64_t bytes encode for the length of the phrase
        if (count % 2 == 0){
            seq_orig_size = seq_orig_size + len;
        }
        count++;
    }

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);

    return seq_orig_size;
}

/**
 * @brief Converts int to char if possible
 * 
 * @param[in] elem [int] an int element.
 * @return string of char or int 
 */
std::string printSymbol(int elem)
{
    if (elem < alpha)
        return std::string(1, static_cast<char>(map[elem]));
    else
        return std::to_string(elem);
}

/**
 * @brief Prints all records in the max heap. Debug purposes only.
 * @return void
 */

void printAllRecords()
{
    int total = 0;
    spdlog::trace("Current records in the heap");
    for (int i = 0; i < Rec.size; i++)
    {
        total += Rec.records[i].freq;
        spdlog::trace("({},{}) {} occs {} heap position {} hash position", printSymbol(Rec.records[i].pair.left), printSymbol(Rec.records[i].pair.right), Rec.records[i].freq, Rec.records[i].hpos, Rec.records[i].kpos);
    }
    spdlog::trace("Total: {}", total);
    spdlog::trace("");
}

/**
 * @brief Prints specific record in the heap. Debug purposes only.
 * 
 * Prints the left and right elements of the record in addition to the frequency.
 * 
 * @param[in] message [std::string] the message to be printed with the record.
 * @param[in] orec [Trecord*] the record content to be printed.
 * @return void
 */

void printRecord(const std::string message, const Trecord* orec)
{
    spdlog::debug("{}",message);
    spdlog::debug("({},{}) {} occs", printSymbol(orec->pair.left), printSymbol(orec->pair.right), orec->freq);
    spdlog::debug("");
}


/**
 * @brief Print max occuring pair. Debug purposes only.
 * 
 * @param[in] new_symbol [int] the new non-terminal symbol to be added.
 * @param[in] orec [Trecord*] the record content to be printed.
 */

void printMaxPair(int new_symbol, const Trecord* orec)
{
    spdlog::debug("Chosen Pair {} = ({},{}) ({} occs)", new_symbol, printSymbol(orec->pair.left), printSymbol(orec->pair.right), orec->freq); 
}

/**
 * @brief Prints the int version of the reference. Debug purposes only.
 * @return void
 */
void printRef()
{
    std::string ref = "";

    RefNode* curr_ref = rlist.getTail();
    while(curr_ref != nullptr)
    {
        curr_ref = rlist.findNearestRef(curr_ref);
        if (curr_ref == nullptr){
            break;
        }

        ref = printSymbol(curr_ref->val) + " " + ref;
        curr_ref = curr_ref->prev;    
    }

    spdlog::trace("Reference string: {}", ref);
    spdlog::trace("");
}

/**
 * @brief Prints the hash table of bi-grams in the reference.
 * @return void
 */
void printHashRanges()
{
    std::string values = "";
    for (const auto& [key, value] : hash_ranges){
        for (int i = 0; i < value.size(); i++){
            values += std::to_string(value[i]);
            values += " ";
        }
        spdlog::trace("Key: ({},{}), Values: {}", printSymbol(key.first), printSymbol(key.second), values);
        values = "";
    }
    spdlog::trace("");
}

/**
 * @brief Prints specific phrase
 */

void printPhrase(PhraseNode* curr_phrase)
{
    std::string content = "";
    if (!(curr_phrase->exp)){
        RefNode* curr_elem = rlist.findNearestRef(curr_phrase->rnode);
        RefNode* first_elem = rlist.findNearestRef(curr_phrase->lnode);
        first_elem = first_elem->prev;
        while(curr_elem != first_elem){
            if (!(curr_elem->deleted))
                content = " " + printSymbol(curr_elem->val) + content;
            curr_elem = curr_elem->prev;
        }

        spdlog::trace("Phrase (Not explicit): {}", content);
    }
    else{
        for (unsigned int i : curr_phrase->content){
            content += printSymbol(i) + " ";
        }
        spdlog::trace("Phrase (explicit): \033[1;31m{}\033[0m", content);
    }
}

/**
 * @brief Print supposed phrase
 */

void printSupposedPhrase(PhraseNode* curr_phrase, std::list<int>::iterator leftIt)
{
    while (leftIt != curr_phrase->content.begin()){
        leftIt = std::prev(leftIt);
    }

    std::string content = "";
    while(leftIt != curr_phrase->content.end())
    {
        content = content + " " + printSymbol(*leftIt);
        leftIt++;
    }
    spdlog::trace("Phrase (explicit): \033[1;31m{}\033[0m", content);
}

/**
 * @brief Calculate number of invalid consecutive same chars in the phrases
 * For example, eeee -> 1 (the 2nd ee)
 * For example, aaa -> 1 (the 2nd aa)
 * For example, iiiii -> 2 (2nd and 4th ii)
 * Assumes run after phrase and source boundary so no crossing boundary
 */

int invalidSameCharPair(int letter)
{
    PhraseNode* curr_phrase = plist.getHead();
    int invalidNexpCount = 0;
    int invalidExpCount = 0;
    int count;
    int countNexp = 0;
    int countExp = 0;
    int left_elem = -1;
    int right_elem = -1;
    while(curr_phrase != nullptr)
    {
        if (!(curr_phrase->exp))
        {
            RefNode* leftNode = rlist.findNearestRef(curr_phrase->lnode);
            RefNode* rightNode = rlist.findNearestRef(curr_phrase->rnode);
            while (leftNode != rightNode->next)
            {
                if (left_elem == -1){
                    left_elem = leftNode->val;
                    leftNode = leftNode->next;
                    continue;
                }
                right_elem = leftNode->val;
                if (left_elem == letter && right_elem == letter){
                    count++;
                    countNexp++;
                }
                else{
                    count = 0;
                }
                if (count > 0 && count % 2 == 0){
                    invalidNexpCount++;
                }
                left_elem = right_elem;
                leftNode = leftNode->next;
            }
            
        }
        else
        {
            auto it = curr_phrase->content.begin();
            while(it != curr_phrase->content.end())
            {
                if (left_elem == -1){
                    left_elem = *it;
                    it++;
                    continue;
                }
                right_elem = *it;
                if (left_elem == letter && right_elem == letter){
                    count++;
                    countExp++;
                }
                else{
                    count = 0;
                }
                if (count > 0 && count % 2 == 0){
                    invalidExpCount++;
                }
                left_elem = right_elem;
                it++;
            }
        } 
        curr_phrase = curr_phrase->next;
    }
    spdlog::trace("Number of Non-explicit Occurences: {}", countNexp);
    spdlog::trace("Invalid Non-explicit Count: {}", invalidNexpCount);
    spdlog::trace("Number of Explicit Occurences: {}", countExp);
    spdlog::trace("Invalid Explicit Count: {}", invalidExpCount);
    return invalidNexpCount + invalidExpCount;
}

/**
 * @brief Checks whether the frequencies in the heap is correct 
 */

bool checkHeap()
{
    std::unordered_map<std::pair<int,int>, int, pair_int_hash> pair_count;
    PhraseNode* curr_phrase = plist.getHead();
    int left_elem = -1;
    int right_elem = -1;
    while(curr_phrase != nullptr)
    {
        if (!(curr_phrase->exp))
        {
            RefNode* leftNode = rlist.findNearestRef(curr_phrase->lnode);
            RefNode* rightNode = rlist.findNearestRef(curr_phrase->rnode);
            while (leftNode != rightNode->next)
            {
                if (left_elem == -1){
                    left_elem = leftNode->val;
                    leftNode = leftNode->next;
                    continue;
                }
                right_elem = leftNode->val;
                pair_count[{left_elem, right_elem}]++;
                left_elem = right_elem;
                leftNode = leftNode->next;
            }
        }
        else
        {
            auto it = curr_phrase->content.begin();
            while(it != curr_phrase->content.end())
            {
                if (left_elem == -1){
                    left_elem = *it;
                    it++;
                    continue;
                }
                right_elem = *it;
                pair_count[{left_elem, right_elem}]++;
                left_elem = right_elem;
                it++;
            }
        } 
        curr_phrase = curr_phrase->next;
    }

    for (auto p : pair_count)
    {
        std::pair<int,int> key = p.first;
        int value = p.second;
        if (value > 1)
        {
            Tpair new_pair;
            new_pair.left = key.first;
            new_pair.right = key.second;
            int nid = searchHash(Hash,new_pair);
            if (nid == -1){
                spdlog::error("Pair ({},{}) does not exist in heap but should!", printSymbol(key.first), printSymbol(key.second));
                return false;
            }
            else{
                Trecord* nrec = &Rec.records[nid];
                if (value != nrec->freq){
                    spdlog::error("Pair ({},{}) | Heap: {} , Actual: {}",  printSymbol(key.first), printSymbol(key.second), nrec->freq, value);
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * @brief Calculates the size of the phrases currently
 */

int checkPhraseSizes()
{
    PhraseNode* curr_phrase = plist.getHead();
    int nexp_length = 0;
    int exp_length = 0;
    int total_length;
    while(curr_phrase != nullptr)
    {
        if (!(curr_phrase->exp)){
            int length = 0;
            RefNode* left_elem = rlist.findNearestRef(curr_phrase->lnode);
            RefNode* right_elem = rlist.findNearestRef(curr_phrase->rnode);
            while (left_elem != right_elem)
            {
                if (left_elem->deleted == false){
                    length++;
                }
                left_elem = left_elem->next;
            }
            length++;
            nexp_length += length;
        }
        else{
            exp_length += curr_phrase->content.size();
        } 
        curr_phrase = curr_phrase->next;
    }
    
    total_length = nexp_length + exp_length;
    spdlog::debug("Non explicit phrase chars: {}", nexp_length);
    spdlog::debug("Explicit phrase chars: {}", exp_length);
    spdlog::debug("Total characters: {}", total_length);
    return total_length;
}


/**
 * @brief Prints the current phrases. Debug purposes only.
 * @return void
 */

void printPhraseList()
{
    std::string content;
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        content = "";
        if (!(curr_phrase->exp)){
            RefNode* curr_elem = rlist.findNearestRef(curr_phrase->rnode);
            RefNode* first_elem = rlist.findNearestRef(curr_phrase->lnode);
            first_elem = first_elem->prev;
            while(curr_elem != first_elem){
                if (!(curr_elem->deleted))
                    content = " " + printSymbol(curr_elem->val) + content;
                curr_elem = curr_elem->prev;
            }

            spdlog::trace("Phrase (Not explicit): {}", content);
        }
        else{
            for (int i : curr_phrase->content){
                content += printSymbol(i) + " ";
            }
            spdlog::trace("Phrase (explicit): \033[1;31m{}\033[0m", content);
        }
        curr_phrase = curr_phrase->next;
    }
    spdlog::trace("");
}

/**
 * @brief Check the exp_pairs stored.
 * 
 * At any point in time the exp_pairs should be up to date
 * 
 */
bool checkExpPairs()
{
    std::unordered_map<std::pair<int, int>, std::unordered_set<ExpPair, ExpPairHash, ExpPairEqual>, pair_int_hash> exp_pairs_tmp(exp_pairs);
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        int sameCharCount = 1;
        if (curr_phrase->exp)
        {
            auto pit = curr_phrase->content.begin();
            auto it = std::next(pit);
            while (it != curr_phrase->content.end())
            {
                if (*it == *pit){
                    sameCharCount++;
                }
                else{
                    sameCharCount = 1;
                }
                if (sameCharCount == 1 || sameCharCount % 2 == 0)
                {
                    if (exp_pairs_tmp[{*pit,*it}].find(ExpPair(curr_phrase, pit, it)) == exp_pairs_tmp[{*pit,*it}].end()){
                        printPhrase(curr_phrase);
                        spdlog::error("Pair missing: ({},{})", printSymbol(*pit), printSymbol(*it));
                        // std::cout << "Address of phrase: " << &(*curr_phrase) << std::endl;
                        // std::cout << "Address of left element pointed to: " << &(*pit) << std::endl;
                        // std::cout << "Address of right element pointed to: " << &(*it) << std::endl;
                        return false;
                    }
                    else{
                        exp_pairs_tmp[{*pit,*it}].erase(ExpPair(curr_phrase, pit, it));
                    }
                }
                else
                {
                    if (exp_pairs_tmp[{*pit,*it}].find(ExpPair(curr_phrase, pit, it)) != exp_pairs_tmp[{*pit,*it}].end()){
                        printPhrase(curr_phrase);
                        spdlog::error("Pair exists that should not exist: ({},{})", printSymbol(*pit), printSymbol(*it));
                        return false;
                    }
                }
                pit = it;
                it = std::next(it);
            }
        }        
        curr_phrase = curr_phrase->next;
    }
    for (const auto& entry : exp_pairs_tmp)
    {
        std::pair<int, int> key = entry.first;
        if (exp_pairs_tmp[key].size() != 0){
            spdlog::error("Key should not exist: ({},{}) : {} pairs", printSymbol(key.first), printSymbol(key.second), exp_pairs_tmp[key].size());
            return false;
        }
    }
    spdlog::debug("Exp pairs is properly maintained");
    return true;
}

/**
 * @brief Check the phrase boundary hash table for correctness
 */
bool checkPhraseBoundaries()
{
    std::unordered_map<std::pair<int, int>, std::unordered_set<PhraseNode*>, pair_int_hash> pbound_pairs_tmp(pbound_pairs);
    PhraseNode* curr_phrase = plist.getHead();
    PhraseNode* next_phrase = curr_phrase->next;
    while(next_phrase != nullptr)
    {
        int left_elem, right_elem;
        if (!curr_phrase->exp){
            left_elem = rlist.findNearestRef(curr_phrase->rnode)->val;
        }
        else{
            left_elem = curr_phrase->content.back();
        }
        if (!next_phrase->exp){
            right_elem = rlist.findNearestRef(next_phrase->lnode)->val;
        }
        else{
            right_elem = next_phrase->content.front();
        }
        pbound_pairs_tmp[{left_elem, right_elem}].erase(curr_phrase);
        curr_phrase = next_phrase;
        next_phrase = next_phrase->next;
    }
    for (const auto& entry : pbound_pairs_tmp)
    {
        std::pair<int, int> key = entry.first;
        if (pbound_pairs_tmp[key].size() != 0){
            spdlog::error("Key should not exist: ({},{}) : {} pairs", printSymbol(key.first), printSymbol(key.second), pbound_pairs_tmp[key].size());
            return false;
        }
    }
    spdlog::debug("The phrase boundaries are correct");
    return true;
}

/**
 * @brief Check the source boundary hash tables for correctness
 */
 bool checkSourceBoundaries()
 {
    std::unordered_map<int, std::unordered_set<PhraseNode*>> start_hash_tmp(start_hash);
    std::unordered_map<int, std::unordered_set<PhraseNode*>> end_hash_tmp(end_hash);
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        if (!curr_phrase->exp){
            int start_elem = rlist.findNearestRef(curr_phrase->lnode)->val;
            int end_elem = rlist.findNearestRef(curr_phrase->rnode)->val;
            start_hash_tmp[start_elem].erase(curr_phrase);
            end_hash_tmp[end_elem].erase(curr_phrase);
        }
        curr_phrase = curr_phrase->next;
    }
    for (const auto& entry : start_hash_tmp)
    {
        int start = entry.first;
        if (start_hash_tmp[start].size() != 0){
            spdlog::error(" {} start count: {} pairs", printSymbol(start), start_hash_tmp[start].size());
            return false;
        }
    }
    for (const auto& entry : end_hash_tmp)
    {
        int end = entry.first;
        if (end_hash_tmp[end].size() != 0){
            spdlog::error(" {} start count: {} pairs", printSymbol(end), end_hash_tmp[end].size());
            return false;
        }
    }
    spdlog::debug("The start and end boundaries are correct");
    return true;
 }

/**
 * @brief Convert the chars of reference to int.
 * This function remaps the characters in the reference sequence to integers. The remapping works as follows:
 * When a new character is encountered in the reference, it is replaced with the corresponding integer from the chars array.
 * Each time a new character is encountered, the alpha variable, which tracks the number of unique characters, is assigned to 
 * the position in the chars array corresponding to the integer version of the character.
 * The alpha variable is then incremented, ensuring that the next new character gets a unique integer.
 * A map array is also maintained to store the reverse mapping, allowing the original character to be retrieved from its integer representation.
 * The logic for this function was taken from RePair/BigRePair.
 * Additionally, currently creates a hash table that stores the positions of unique bi-grams in the reference
 * 
 * @param[in] rtext [std::vector<unsigned char>] The unsigned int representation of the .
 * 
 * @return void
 */

void prepareRef(std::vector<unsigned char>& rtext)
{
    alpha = 0;

    // Initialize the chars array to be -1 to indicate that no characters have been mapped yet.
    for (int i = 0; i < 256; i++){
        chars[i] = -1;
    }

    // Remaps the chars in ref, updates chars array, populates table of ranges of reference bi-grams
    RefNode* prev_elem;
    RefNode* curr_elem;
    std::pair<int, int> ref_pair;
    rarray.reserve(rtext.size());
    for (int i = 0; i < rtext.size(); i++){
        unsigned char x = rtext[i];
        if (chars[x] == -1){
            chars[x] = alpha++;
        }
        curr_elem = rlist.push_back(chars[x]);
        rarray.emplace_back(curr_elem);
        ref_pair.second = curr_elem->val;
        if (i != 0){
            auto hash_range_start = std::chrono::high_resolution_clock::now();
            hash_ranges[ref_pair].push_back(i-1);
            auto hash_range_end = std::chrono::high_resolution_clock::now();
            hash_range_time += hash_range_end - hash_range_start;
        }
        ref_pair.first = ref_pair.second;          
    }

    // Updates the map array to undo remapping.
    for (int i = 0; i < 256; i++){
        if (chars[i] != -1){
            map[chars[i]] = i;
        }
    }

    // Set n to be alpha
    n = alpha;

    // If trace messages are enabled
    if (verbosity == 2){
        spdlog::trace("Reference at the start");
        spdlog::trace("Size: {}", rlist.getSize());
        printRef();
        spdlog::trace("Hash ranges of reference bi-grams at the start");
        printHashRanges();
    }
}

/**
 * @brief Creates max heap of the unique bi-grams in the sequence via the RLZ parse.
 * 
 * We directly are using Heap/Hash/Record data structures in RePair/BigRePair to create the max heap
 * Max heap allows us to quickly find the bi-gram with the highest occurence.
 * 
 * @param[in] pfile [std::infstream&] the RLZ parse file stream
 * 
 * @return void
 */
void createMaxHeap(std::ifstream& pfile)
{
    // Create pair-freq heap (Thanks BigRePair)
    Rec = createRecords(factor, minsize);
    Heap = createHeap(psize, &Rec, factor, minsize);
    Hash = createHash(256 * 256, &Rec);
    assocRecords(&Rec, &Hash, &Heap, NULL);

    uint64_t num_pairs, pos, len;
    pair.left = -1;
    pair.right = -1;
    RefNode* rnode;
    int id;

    // First uint64_t bytes tell how many (pos,len) pairs are stored in the parse file
    pfile.read(reinterpret_cast<char*>(&num_pairs), sizeof(uint64_t));

    // Update the heap with the frequencies of all pairs
    for (uint64_t i = 0; i < 2 * num_pairs; i++)
    {
        if (i % 2 == 0){
            pfile.read(reinterpret_cast<char*>(&pos), sizeof(uint64_t));
        }
        else
        {
            pfile.read(reinterpret_cast<char*>(&len), sizeof(uint64_t));
            for (uint64_t j = 0; j < len; j++)
            {
                if (j == 0){
                    rnode = rarray[pos];
                }
                if (pair.left == -1){
                    pair.left = rnode->val;
                    rnode = rnode->next;
                }
                else{
                    pair.right = rnode->val;
                    id = searchHash(Hash, pair);
                    if (id == -1) // new pair, insert
                    {
                        id = insertRecord(&Rec, pair);
                    }
                    else
                    {
                        incFreq(&Heap, id);
                    }
                    pair.left = pair.right;
                    rnode = rnode->next;
                }
            }
        }
    }
    // Print the records in the heap
    spdlog::trace("Records in the heap at the start (after removing freq 1 records)");
    // Remove frequency 1 records
    purgeHeap(&Heap);

    if (verbosity == 2){
        printAllRecords();
    }

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);
}

/**
 * @brief Populates the phrases in RLZ parse. 
 * 
 * Each phrase created will be a non-explicit phrase. 
 * A non-explicit phrase stores left and right endpoints that references the reference.
 * Range for each phrase is [left,right] 
 * 
 * @param[in] pfile [std::ifstream&] the RLZ parse filestream
 * @return void 
 */

void populatePhrases(std::ifstream& pfile)
{
    uint64_t num_pairs, pos, len;
    PhraseNode* prevPhrase;
    PhraseNode* nextPhrase;

    // First uint64_t bytes tell how many (pos,len) pairs are stored in the parse file
    pfile.read(reinterpret_cast<char*>(&num_pairs), sizeof(uint64_t));

    // Update the heap with the frequencies of all pairs
    for (uint64_t i = 0; i < 2 * num_pairs; i++)
    {
        if (i % 2 == 0){
            pfile.read(reinterpret_cast<char*>(&pos), sizeof(uint64_t));
        }
        else
        {
            pfile.read(reinterpret_cast<char*>(&len), sizeof(uint64_t));
            plist.push_back(rarray[pos], rarray[pos+len-1]);
            nextPhrase = plist.getTail();
            // Only start adding to the phrase boundary hash table after the first phrase is added
            if (i != 1)
            {
                auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                prevPhrase = plist.getTail()->prev;
                pbound_pairs[{prevPhrase->rnode->val, nextPhrase->lnode->val}].insert(prevPhrase);
                auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
            }
            // Add the start and end character of each phrase to the source boundary hash tables
            auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
            start_hash[nextPhrase->lnode->val].insert(nextPhrase);
            end_hash[nextPhrase->rnode->val].insert(nextPhrase);
            auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
            update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
        }
    }

    checkPhraseBoundaries();
    checkSourceBoundaries();

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);

    // Debug
    if (verbosity == 2){
        spdlog::trace("The non-explicit phrases at the start");
        printPhraseList();
    }
}

/**
 * @brief Builds interval tree from the non-explicit phrases.
 */

void buildIntervalTree()
{
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        if (!curr_phrase->exp){
            int lRange = rlist.findNearestRef(curr_phrase->lnode)->pos;
            int rRange = rlist.findNearestRef(curr_phrase->rnode)->pos;
            phrase_tree.insert({lRange, rRange}, curr_phrase);
        }
        curr_phrase = curr_phrase->next;
    }
}

/**
 * @brief Ensures that consecutive same characters >2 do not have overlapping entries in exp_pairs.
 * When adding from left, we might have to update exp_pairs all pairs of exp_pairs. (ie. e + ee -> ee + e)
 * When adding from right, we check if the consecutive same characters are odd or even
 */

void updateExpPairs(PhraseNode* p, std::list<int>::iterator it, bool leftInsert)
{
    int sameCharCount = 1;
    int letter = *it;
    if (leftInsert)
    {
        auto nextIt = std::next(it);
        while(nextIt != p->content.end())
        {
            if (*nextIt != letter){
                break;
            }
            else{
                sameCharCount++;
                if (sameCharCount % 2 == 0){
                    exp_pairs[{letter,letter}].insert(ExpPair(p, it, nextIt));
                }
                else{
                    exp_pairs[{letter,letter}].erase(ExpPair(p, it, nextIt));
                }
            }
            it = nextIt;
            nextIt = std::next(nextIt);
        }
    }
    else
    {
        while (it != p->content.begin())
        {
            if (*it != letter){
                it = std::next(it);
                break;
            }
            it = std::prev(it);
        }
        if (it == p->content.begin() && *it != letter){
            it = std::next(it);
        }
        auto nextIt = std::next(it);
        while(nextIt != p->content.end())
        {
            if (*nextIt != letter){
                break;
            }
            else{
                sameCharCount++;
                if (sameCharCount % 2 == 0){
                    exp_pairs[{letter,letter}].insert(ExpPair(p, it, nextIt));
                }
                else{
                    exp_pairs[{letter,letter}].erase(ExpPair(p, it, nextIt));
                }
            }
            it = nextIt;
            nextIt = std::next(nextIt);
        }       
    }
}

/** 
 * @brief Ensures that consecutive same characters >2 do not have overlapping entries in exp_pairs.
 * Used when merging explicit phrases together
 */

void updateMergeExpPairs(PhraseNode* p, std::list<int>::iterator it)
{
    int letter = *it;
    // Find the first occurence of the letter
    while (it != p->content.begin())
    {
        if (*it != letter){
            it = std::next(it);
            break;
        }
        it = std::prev(it);
    }
    if (it == p->content.begin() && *it != letter){
        it = std::next(it);
    }
    auto nextIt = std::next(it);
    int sameCharCount = 1;
    while(nextIt != p->content.end())
    {
        if (*nextIt != letter){
            break;
        }
        else{
            sameCharCount++;
            if (sameCharCount % 2 == 0){
                exp_pairs[{letter,letter}].insert(ExpPair(p, it, nextIt));
            }
            else{
                exp_pairs[{letter,letter}].erase(ExpPair(p, it, nextIt));
            }
            it = nextIt;
            nextIt = std::next(nextIt);
        }
    }
}

/**
 * @brief Deletes the exp pairs assigned to one exp phrase and assigns to another exp phrase 
 */

void reassignExpPairs(PhraseNode* origPhrase, PhraseNode* newPhrase)
{
    auto it = origPhrase->content.begin();
    auto nextIt = std::next(it);
    int sameCharCount = 1;
    while(nextIt != origPhrase->content.end())
    {
        if (*it == *nextIt){
            sameCharCount++;
        }
        else{
            sameCharCount = 1;
        }

        if (sameCharCount == 1 || sameCharCount % 2 == 0)
        {
            exp_pairs[{*it, *nextIt}].erase(ExpPair(origPhrase, it, nextIt));
            exp_pairs[{*it, *nextIt}].insert(ExpPair(newPhrase, it, nextIt));
        }
        it = nextIt;
        nextIt = std::next(nextIt);
    }
}

/**
 * @brief Process the phrase list for the phrase boundary condition.
 * 
 * If the rightmost elem of a non-explicit phrase and the non-explicit leftmost elem of the adjacent phrase 
 * form the provided bi-gram, both elements are removed from their respective phrases and added together to create 
 * an explicit phrase. If one of the phrases are explicit already, only from the non-explict phrase is the 
 * elem removed (since the other phrase is already explicit).
 * 
 * @param [in] left_elem [int] the left elem of the max occuring bi-gram
 * @param [in] right_elem [int] the right elem of the max occuring bi-gram
 * 
 * @return void
 */

void phraseBoundaries(int left_elem, int right_elem) // TODO: Take care of case where phrase size is 1 with source boundaries
{
    // If the pair does not exist in the phase boundary list then can return immediately
    if (pbound_pairs.find({left_elem, right_elem}) == pbound_pairs.end()){
        return;
    }
    // Iterate through the phrases with the phrase boundary of interest
    std::unordered_set<PhraseNode*> boundaries = pbound_pairs[{left_elem, right_elem}];
    std::unordered_set<PhraseNode*>::iterator it = boundaries.begin();
    while (it != boundaries.end()) 
    {
        PhraseNode* curr_phrase = *it;
        curr_phrase->ltmp = -1; // TODO: Now later on have to fix if not encountered
        curr_phrase->rtmp = -1;
        PhraseNode* next_phrase = curr_phrase->next;
        // If there is a next phrase, check the phrase boundaries.
        if (next_phrase != nullptr)
        {
            // Both phrases not explicit
            if (!(curr_phrase->exp) && !(next_phrase->exp))
            {
                spdlog::debug("Both non-explict");
                if (rlist.findNearestRef(curr_phrase->rnode)->val == left_elem && rlist.findNearestRef(next_phrase->lnode)->val == right_elem)
                {
                    // Indicates whether a non-explicit phrase got deleted.
                    bool deleteCurr = false;
                    bool deleteNext = false;

                    std::list<int> content;
                    content.push_back(left_elem);
                    content.push_back(right_elem);

                    // First remove the offending entries in the tree.
                    //spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
                    //spdlog::trace("Removing ({},{}) from tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);
                    next_phrase->lnode = rlist.findNearestRef(next_phrase->lnode);
                    next_phrase->rnode = rlist.findNearestRef(next_phrase->rnode);

                    // First we delete the phrase boundary entry and the end char of the curr phrase and start char of the next phrase
                    auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    pbound_pairs[{curr_phrase->rnode->val,next_phrase->lnode->val}].erase(curr_phrase);
                    end_hash[curr_phrase->rnode->val].erase(curr_phrase);
                    start_hash[next_phrase->lnode->val].erase(next_phrase);
                    auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;

                    // Update the pointers of the non-explicit phrases
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);
                    next_phrase->lnode = rlist.findForwardRef(next_phrase->lnode);
                    // Insert the explicit phrase to phrase list
                    PhraseNode* expPhrase = plist.insert(next_phrase, content);

                    // Add the new exp pair to exp_pairs
                    auto l = expPhrase->content.begin();
                    auto r = std::prev(expPhrase->content.end());
                    exp_pairs[{*l, *r}].insert(ExpPair(expPhrase, l, r));
                    // Only two characters guaranteed so do not have to check consecutive characters

                    // If the current or next phrases are empty we delete
                    if (curr_phrase->rnode == nullptr || curr_phrase->lnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){
                        deleteCurr = true; 
                        curr_phrase = plist.remove(curr_phrase);
                    }
                    if (next_phrase->rnode == nullptr || next_phrase->lnode == nullptr || next_phrase->lnode->pos > next_phrase->rnode->pos){
                        deleteNext = true;
                        next_phrase = plist.remove(next_phrase);  
                    }

                    // Update the tree with the new entries
                    if (!deleteCurr && !deleteNext){
                        //spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        //spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        end_hash[curr_phrase->rnode->val].insert(curr_phrase);
                        start_hash[next_phrase->lnode->val].insert(next_phrase);
                        // Have to add boundary between the current phrase and new explicit phrase + new explicit phrase and next phrase
                        pbound_pairs[{curr_phrase->rnode->val, expPhrase->content.front()}].insert(curr_phrase);
                        pbound_pairs[{expPhrase->content.back(), next_phrase->lnode->val}].insert(expPhrase);
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    else if (deleteCurr && !deleteNext){
                        //spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        // Update the hash boundary tables.
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        // Since the current phrase is deleted, we update the start hash of the next phrase and not the end hash of the curr phrase
                        start_hash[next_phrase->lnode->val].insert(next_phrase);
                        pbound_pairs[{expPhrase->content.back(), next_phrase->lnode->val}].insert(expPhrase);
                        // If curr phrase is not expPhrase that means that there is a previous phrase and not at start of phrase list 
                        if (curr_phrase != expPhrase){
                            if (!curr_phrase->exp)
                                pbound_pairs[{rlist.findNearestRef(curr_phrase->rnode)->val, expPhrase->content.front()}].insert(curr_phrase);
                            else
                                pbound_pairs[{curr_phrase->content.back(), expPhrase->content.front()}].insert(curr_phrase);
                        }
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    else if (!deleteCurr && deleteNext){
                        //spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        // Since the next phrase is deleted, we update the end hash of the curr phrase and not the start hash of the next phrase
                        end_hash[curr_phrase->rnode->val].insert(curr_phrase);
                        pbound_pairs[{curr_phrase->rnode->val, expPhrase->content.front()}].insert(curr_phrase);
                        // If expPhrase->next is not nullptr there is a next phrase to connect to 
                        if (expPhrase->next != nullptr){
                            if (!expPhrase->next->exp)
                                pbound_pairs[{expPhrase->content.back(), rlist.findNearestRef(expPhrase->next->lnode)->val}].insert(expPhrase);
                            else
                                pbound_pairs[{expPhrase->content.back(), expPhrase->next->content.front()}].insert(expPhrase);
                        }
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    else{ 
                        // If both the current and next_phrase are deleted then have to connect the prev prev phrase with expPhase + expPhrase with next next phrase
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        if (curr_phrase != expPhrase){
                            if (!curr_phrase->exp)
                                pbound_pairs[{rlist.findNearestRef(curr_phrase->rnode)->val, expPhrase->content.front()}].insert(curr_phrase);
                            else
                                pbound_pairs[{curr_phrase->content.back(), expPhrase->content.front()}].insert(curr_phrase);
                        }
                        if (expPhrase->next != nullptr){
                            if (!expPhrase->next->exp)
                                pbound_pairs[{expPhrase->content.back(), rlist.findNearestRef(expPhrase->next->lnode)->val}].insert(expPhrase);
                            else
                                pbound_pairs[{expPhrase->content.back(), expPhrase->next->content.front()}].insert(expPhrase);
                        }
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    continue;
                }
            }
            // Current phrase not explicit and next phrase explicit
            else if (!(curr_phrase->exp) && (next_phrase->exp))
            {
                if ((rlist.findNearestRef(curr_phrase->rnode)->val == left_elem) && (next_phrase->content.front() == right_elem))
                {
                    // Indicates whether a non-explicit phrase got deleted.
                    bool deleteCurr = false;
                    // First remove the offending entry from the tree
                    //spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

                    // Delete from the hash tables
                    auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    pbound_pairs[{curr_phrase->rnode->val,next_phrase->content.front()}].erase(curr_phrase);
                    end_hash[curr_phrase->rnode->val].erase(curr_phrase);
                    auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;

                    // Getting the iterators
                    auto r = next_phrase->content.begin();
                    next_phrase->content.push_front(left_elem);
                    auto l = next_phrase->content.begin();
                    // Add the new exp pair to exp_pairs
                    exp_pairs[{*l, *r}].insert(ExpPair(next_phrase, l, r));
                    // If consecutive characters, then might have to change exp_pairs (i.e e + ee -> ee + e)
                    if (*l == *r){
                        updateExpPairs(next_phrase, l, true);
                    }
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);
                    // If the non-explicit phrase is empty we delete it.
                    if (curr_phrase->rnode == nullptr || curr_phrase->lnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){ 
                        deleteCurr = true;
                        curr_phrase = plist.remove(curr_phrase);
                    }

                    // Update the tree with the new entry
                    if (!deleteCurr){
                        //spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        end_hash[curr_phrase->rnode->val].insert(curr_phrase);
                        pbound_pairs[{curr_phrase->rnode->val, next_phrase->content.front()}].insert(curr_phrase);
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    else{
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        if (curr_phrase != next_phrase){
                            if (!curr_phrase->exp)
                                pbound_pairs[{rlist.findNearestRef(curr_phrase->rnode)->val, next_phrase->content.front()}].insert(curr_phrase);
                            else
                                pbound_pairs[{curr_phrase->content.back(), next_phrase->content.front()}].insert(curr_phrase);
                        }
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    continue;
                }
            }
            // Current phrase explicit and next phrase not explicit
            else if ((curr_phrase->exp) && !(next_phrase->exp))
            {
                if ((curr_phrase->content.back() == left_elem) && (rlist.findNearestRef(next_phrase->lnode)->val == right_elem))
                {
                    // Indicates whether a non-explicit phrase got deleted.
                    bool deleteNext = false;
                    // First remove the offending entry from the tree
                    //spdlog::trace("Removing ({},{}) from tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    next_phrase->lnode = rlist.findNearestRef(next_phrase->lnode);
                    next_phrase->rnode = rlist.findNearestRef(next_phrase->rnode);

                    // Delete from the hash tables
                    auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    pbound_pairs[{curr_phrase->content.back(), next_phrase->lnode->val}].erase(curr_phrase);
                    start_hash[next_phrase->lnode->val].erase(next_phrase);
                    auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;

                    // Getting the iterators
                    auto l = std::prev(curr_phrase->content.end());
                    curr_phrase->content.push_back(right_elem);
                    auto r = std::prev(curr_phrase->content.end());
                    // Add the new exp pair to exp_pairs
                    exp_pairs[{*l, *r}].insert(ExpPair(curr_phrase, l, r));
                    // If consecutive characters, then might have to change exp_pairs (i.e e + ee -> ee + e)
                    if (*l == *r){
                        updateExpPairs(curr_phrase, r, false);
                    }
                    next_phrase->lnode = rlist.findForwardRef(next_phrase->lnode);
                    if (next_phrase->rnode == nullptr || next_phrase->lnode == nullptr || next_phrase->lnode->pos > next_phrase->rnode->pos){
                        deleteNext = true;
                        next_phrase = plist.remove(next_phrase);
                    }

                    // Update the tree with the new entry
                    if (!deleteNext){
                        //spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        start_hash[next_phrase->lnode->val].insert(next_phrase);
                        pbound_pairs[{curr_phrase->content.back(), next_phrase->lnode->val}].insert(curr_phrase);
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    else{
                        // Update the hash boundary tables
                        update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        if (next_phrase->next != nullptr){
                            if (!next_phrase->next->exp)
                                pbound_pairs[{curr_phrase->content.back(), rlist.findNearestRef(next_phrase->next->lnode)->val}].insert(curr_phrase);
                            else
                                pbound_pairs[{curr_phrase->content.back(), next_phrase->next->content.front()}].insert(curr_phrase);
                        }
                        update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                    continue;
                }
            }
            // Both explicit phrases, move content of next phrase into first phrase and delete next phrase
            else{
                // Update the hash boundary tables
                pbound_pairs[{curr_phrase->content.front(), next_phrase->content.back()}].erase(curr_phrase);
                // Do the re-assignment
                auto l = std::prev(curr_phrase->content.end());
                auto r = next_phrase->content.begin();
                reassignExpPairs(next_phrase, curr_phrase); // Reassigns the pairs in next phrase to current phrase
                curr_phrase->content.splice(curr_phrase->content.end(), next_phrase->content);
                next_phrase = plist.remove(next_phrase); // Should be safe to call remove. The pointers within content should still be valid.
                // Add the new exp pair to exp_pairs
                exp_pairs[{*l, *r}].insert(ExpPair(curr_phrase, l, r));
                if (*l == *r){
                    updateMergeExpPairs(curr_phrase, l);
                }
                continue;
            }
            // Update the phrase if it gets to this point
            it++;
        }
    }

    // Debug
    if (verbosity == 2){
        spdlog::trace("Phrase list after phrase boundary condition.");
        printPhraseList();
    }
}

/**
 * @brief Process the phrase list for the source boundary condition.
 * 
 * If the rightmost elem of a non-explicit phrase and the next elem in the reference form the provided bi-gram
 * or the leftmost elem of a non-explicit phrase and the previous elem in the reference form the provided bi-gram
 * then remove the offending elem from the non-explicit phrase and create an explicit phrase.
 * 
 * This function should uphold the commitment that when we create or add to an explicit phrase to the left of the
 * current phrase that we switch the current phrase to the previous phrase (bi-gram formed to the left of the phrase).
 * If we create or add to an explicit phrase to the right (bi-gram formed to the right of the phrase) 
 * then the iteration should maintain the same current phrase.  
 * 
 * @param [in] left_elem [int] the left elem of the max occuring bi-gram
 * @param [in] right_elem [int] the right elem of the max occuring bi-gram
 * 
 * @return void
 * 
 */
void sourceBoundaries(int left_elem, int right_elem)
{
    PhraseNode* curr_phrase = plist.getHead();
    // Iterate through the phrases in the phrase list
    while (curr_phrase != nullptr) 
    {
        // Check if the curr phrase is explicit or not
        if (!(curr_phrase->exp))
        {
            // Check if the current phrase with its leftmost elem can form bi-gram with prev elem on reference
            RefNode* rightElem = rlist.findNearestRef(curr_phrase->lnode);
            RefNode* leftElem = rlist.findNearestRef(rightElem->prev);

            // If the bi-gram can be formed then make the leftmost element of the current phrase an explicit phrase of its own
            if (leftElem != nullptr && rightElem != nullptr && leftElem->val == left_elem && rightElem->val == right_elem){
                PhraseNode* prev_phrase = curr_phrase->prev;
                // If the previous phrase is explicit, we can add directly to it.
                if (prev_phrase != nullptr && prev_phrase->exp){
                    auto l = std::prev(prev_phrase->content.end());
                    prev_phrase->content.push_back(right_elem);
                    auto r = std::prev(prev_phrase->content.end());
                    // Add the new exp pair to exp_pairs
                    exp_pairs[{*l,*r}].insert(ExpPair(prev_phrase, l, r));
                    if (*l == *r){
                        updateExpPairs(prev_phrase, r, false);
                    }
                }
                // We have to create new explicit phrase anyways
                else{
                    std::list<int> content;
                    content.push_back(right_elem);
                    plist.insert(curr_phrase, content);
                    // Update the hash table if only the old previous phrase is not nullptr
                    if (prev_phrase != nullptr){
                        auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                        pbound_pairs[{rlist.findNearestRef(prev_phrase->rnode)->val, prev_phrase->next->content.front()}].insert(prev_phrase);
                        auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                        update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                    }
                }

                // Indicates whether a non-explicit phrase got deleted.
                bool deleteCurr = false;
                // First remove the offending entry from the tree
                //spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                auto update_interval_start = std::chrono::high_resolution_clock::now();
                phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                auto update_interval_end = std::chrono::high_resolution_clock::now();
                update_interval_time += update_interval_end - update_interval_start;

                // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

                // Delete from the hash tables
                auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                if (prev_phrase != nullptr && prev_phrase->exp){
                    pbound_pairs[{prev_phrase->content.back(), curr_phrase->lnode->val}].erase(prev_phrase);
                }
                else if (prev_phrase != nullptr && !prev_phrase->exp){
                    pbound_pairs[{rlist.findNearestRef(prev_phrase->rnode)->val, curr_phrase->lnode->val}].erase(prev_phrase);
                }
                // Always delete from start hash 
                start_hash[rlist.findNearestRef(curr_phrase->lnode)->val].erase(curr_phrase);
                // If the curr phrase is length 1 then have to delete from end hash as well
                if (curr_phrase->lnode == curr_phrase->rnode){
                    end_hash[rlist.findNearestRef(curr_phrase->rnode)->val].erase(curr_phrase);
                }
                auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;

                // Update the information of the current phrase (left node pointer & existence).
                curr_phrase->lnode = rlist.findForwardRef(curr_phrase->lnode);                    
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->lnode == nullptr || curr_phrase->lnode->pos > curr_phrase->rnode->pos){ 
                    deleteCurr = true;
                    curr_phrase = plist.remove(curr_phrase); // Removes the current phrase and sets the curr phrase to be the previous phrase
                }

                // Update the tree with the new entry
                if (!deleteCurr){
                    //spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
                    // Update the hash boundary tables
                    update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    start_hash[curr_phrase->lnode->val].insert(curr_phrase);
                    prev_phrase = curr_phrase->prev;
                    pbound_pairs[{prev_phrase->content.back(), curr_phrase->lnode->val}].insert(prev_phrase);
                    update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                }
                else{
                    // Update the hash boundary tables
                    update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    if (curr_phrase->next != nullptr)
                    {
                        if (!curr_phrase->next->exp)
                            pbound_pairs[{curr_phrase->content.back(), rlist.findNearestRef(curr_phrase->next->lnode)->val}].insert(curr_phrase);
                        else
                            pbound_pairs[{curr_phrase->content.back(), curr_phrase->content.front()}].insert(curr_phrase);
                    }
                    update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                }

                continue;
            }
            
            // Check if the current phrase with its rightmost elem can form bi-gram with right elem on reference
            leftElem = rlist.findNearestRef(curr_phrase->rnode);
            rightElem = rlist.findForwardRef(leftElem);
            
            // If the bi-gram can be formed then make the rightmost element of the current phrase an explicit phrase of its own
            if (leftElem != nullptr && rightElem != nullptr && leftElem->val == left_elem && rightElem->val == right_elem){
                PhraseNode* next_phrase = curr_phrase->next;                    
                // If the next phrase is explicit, we can add directly to it.
                if (next_phrase != nullptr && next_phrase->exp){
                    auto r = next_phrase->content.begin();
                    next_phrase->content.push_front(left_elem);
                    auto l = next_phrase->content.begin();
                    // Add the new exp pair to exp_pairs
                    exp_pairs[{*l,*r}].insert(ExpPair(next_phrase, l, r));
                    if (*l == *r){
                        updateExpPairs(next_phrase, l, true);
                    }
                }
                // We have to create new explicit phrase anyways
                else{
                    std::list<int> content;
                    content.push_back(left_elem);
                    if (next_phrase == nullptr){
                        plist.push_back(content);
                    }
                    else{
                        plist.insert(next_phrase, content); 
                        // Update the hash table if only the old next phrase is not nullptr
                        if (next_phrase != nullptr){
                            auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                            pbound_pairs[{curr_phrase->next->content.back(), rlist.findNearestRef(next_phrase->lnode)->val}].insert(curr_phrase->next);
                            auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                            update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                        }
                    }
                }

                // Indicates whether a non-explicit phrase got deleted.
                bool deleteCurr = false;
                // First remove the offending entry from the tree
                //spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                auto update_interval_start = std::chrono::high_resolution_clock::now();
                phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                auto update_interval_end = std::chrono::high_resolution_clock::now();
                update_interval_time += update_interval_end - update_interval_start;

                // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

                // Delete from the hash tables
                auto update_bound_hash_start = std::chrono::high_resolution_clock::now();
                if (next_phrase != nullptr && next_phrase->exp){
                    pbound_pairs[{curr_phrase->rnode->val, next_phrase->content.front()}].erase(curr_phrase);
                }
                else if (next_phrase != nullptr && !next_phrase->exp){
                    pbound_pairs[{curr_phrase->rnode->val, rlist.findNearestRef(next_phrase->lnode)->val}].erase(curr_phrase);
                }
                // Always delete from end hash
                end_hash[rlist.findNearestRef(curr_phrase->rnode)->val].erase(curr_phrase);
                // If the curr phrase is length 1 then have to delete from start hash as well
                if (curr_phrase->lnode == curr_phrase->rnode){
                    start_hash[rlist.findNearestRef(curr_phrase->lnode)->val].erase(curr_phrase);
                }
                auto update_bound_hash_end = std::chrono::high_resolution_clock::now();
                update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;

                // Update the information of the current phrase (right node pointer & existence).
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);                       
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->rnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){
                    deleteCurr = true; 
                    curr_phrase = plist.remove(curr_phrase); // Removes the current phrase and sets the curr phrase to be the previous phrase
                }

                // Update the tree with the new entry
                if (!deleteCurr){
                    //spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
                    // Update the hash boundary tables
                    update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    end_hash[curr_phrase->rnode->val].insert(curr_phrase);
                    next_phrase = curr_phrase->next;
                    pbound_pairs[{curr_phrase->rnode->val, next_phrase->content.front()}].insert(curr_phrase);
                    update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                }
                else{
                    // Update the hash boundary tables
                    update_bound_hash_start = std::chrono::high_resolution_clock::now();
                    if (curr_phrase != next_phrase->prev)
                    {
                        if (!curr_phrase->exp)
                            pbound_pairs[{rlist.findNearestRef(curr_phrase->rnode)->val, next_phrase->prev->content.front()}].insert(curr_phrase);
                        else
                            pbound_pairs[{curr_phrase->content.back(), next_phrase->prev->content.front()}].insert(curr_phrase);
                    }
                    update_bound_hash_end = std::chrono::high_resolution_clock::now();
                    update_bound_hash_time += update_bound_hash_end - update_bound_hash_start;
                }

                continue; // If it gets to here then the next iteration will have the same curr phrase.
            }
        }
        // If explicit phrase then check surrounding phrases if they are explicit and if so merge.
        else{
            PhraseNode* prev_phrase = curr_phrase->prev;
            // If prev phrase is also explicit merge into the prev phrase
            if (prev_phrase != nullptr && prev_phrase->exp){
                // Update the hash boundary tables
                pbound_pairs[{prev_phrase->content.back(), curr_phrase->content.front()}].erase(prev_phrase);
                auto l = std::prev(prev_phrase->content.end());
                auto r = curr_phrase->content.begin();
                reassignExpPairs(curr_phrase, prev_phrase); // Reassigns the pairs in curr phrase to previous phrase
                prev_phrase->content.splice(prev_phrase->content.end(), curr_phrase->content);
                curr_phrase = plist.remove(curr_phrase); // Deletes the current phrase and sets the current phrase to be the previous phrase
                // Add the new exp pair to exp_pairs
                exp_pairs[{*l, *r}].insert(ExpPair(prev_phrase, l, r));
                if (*l == *r){
                    updateMergeExpPairs(prev_phrase, l);
                }
                continue;
            }
            PhraseNode* next_phrase = curr_phrase->next;
            if (next_phrase != nullptr && next_phrase->exp){
                // Update the hash boundary tables
                pbound_pairs[{curr_phrase->content.back(), next_phrase->content.front()}].erase(curr_phrase);
                auto l = std::prev(curr_phrase->content.end());
                auto r = next_phrase->content.begin();
                reassignExpPairs(next_phrase, curr_phrase); // Reassigns the pairs in next phrase to current phrase
                curr_phrase->content.splice(curr_phrase->content.end(), next_phrase->content);
                plist.remove(next_phrase); // Deletes the next phrase
                // Add the new exp pair to exp_pairs
                exp_pairs[{*l, *r}].insert(ExpPair(curr_phrase, l, r));
                if (*l == *r){
                    updateMergeExpPairs(curr_phrase, l);
                }
                continue; // If it gets to here then the next iteration will have the same curr phrase.
            }
        }
        curr_phrase = curr_phrase->next;
    }
    // Debug
    if (verbosity == 2){
        spdlog::trace("Phrase list after source boundary condition.");
        printPhraseList();
    }
}

/**
 * @brief Decrease the frequency of a pair in the heap.
 * 
 * @param [in] left_elem [int] the left elem of the max occuring bi-gram
 * @param [in] right_elem [int] the right elem of the max occuring bi-gram
 * 
 * @return void
 */
void decreaseFrequency(int left, int right)
{
    // if (verbosity == 2){
    //     spdlog::trace("Decrease frequency: ({},{})", printSymbol(left), printSymbol(right));
    // }
    Tpair new_pair;
    new_pair.left = left;
    new_pair.right = right;
    int id = searchHash(Hash,new_pair);
    if (id != -1 && id != oid){
        decFreq(&Heap,id);
    }
}

/**
 * @brief Increase the frequency of a pair in the heap.
 * 
 * @param [in] left_elem [int] the left elem of the max occuring bi-gram
 * @param [in] right_elem [int] the right elem of the max occuring bi-gram
 * 
 * @return void
 */
void increaseFrequency(int left, int right)
{
    // if (verbosity == 2){
    //     spdlog::trace("Increase frequency: ({},{})", printSymbol(left), printSymbol(right));
    // }
    Tpair new_pair;
    new_pair.left = left;
    new_pair.right = right;
    int id = searchHash(Hash,new_pair);
    if (id == -1){
        id = insertRecord(&Rec,new_pair);
    }
    else{
        incFreq(&Heap,id);
    }
}

/**
 * @brief Runs RePair on the non-explicit and explicit phrases.
 * 
 * RePair replaces the most occuring bi-gram with a new non-terminal symbol (n).
 * Before we do replacement, we have to check and correct the phrase and source boundary conditions.
 * If there is an explicit, we have to do normal pair by pair checking.
 * Non-explicit phrases, we can handle by looking at their ranges.
 * We also have to keep track of the frequency as we do the replacement.
 * 
 * Our goal is to handle the non-explicit phrases in time O(|phrases|). 
 * 
 * @param [in] R [std::ofstream&] The file where we will write the rules.
 * @param [in] C [std::ofstream&] The file where we will write the compressed text.
 * 
 * @return void
 */
void repair(std::ofstream& R, std::ofstream& C)
{
    int start_size = psize;

    // Write alpha to R file
    R.write(reinterpret_cast<const char*>(&alpha), sizeof(int));
    if (!R) {
        spdlog::error("Error writing alpha to R file");
    }
    // Write map to R file
    R.write(map, alpha);
    if (!R) {
        spdlog::error("Error writing map to R file");
    }

    // Build the tree
    auto build_interval_start = std::chrono::high_resolution_clock::now();
    buildIntervalTree(); 
    auto build_interval_end = std::chrono::high_resolution_clock::now();
    build_interval_time += build_interval_end - build_interval_start;

    oid = extractMax(&Heap);
    while (oid != -1)
    {
        Trecord* orec = &Rec.records[oid];
        // When max frequency is 1, RePair ends.
        if (orec->freq == 1){
            break;
        }

        if (verbosity == 1 || verbosity == 2){
            printMaxPair(n, orec);
        }

        // Write pair to R file 
        R.write(reinterpret_cast<const char*>(&(orec->pair)), sizeof(Tpair));
        if (!R) {
            spdlog::error("Error writing max occurrence pair to R file");
        }

        int left_elem = orec->pair.left;
        int right_elem = orec->pair.right;
        
        auto pbound_start = std::chrono::high_resolution_clock::now();
        phraseBoundaries(left_elem, right_elem);
        checkPhraseBoundaries();
        checkSourceBoundaries();
        auto pbound_end = std::chrono::high_resolution_clock::now();
        phrase_boundary_time += pbound_end - pbound_start;

        auto sbound_start = std::chrono::high_resolution_clock::now();
        sourceBoundaries(left_elem, right_elem);
        checkPhraseBoundaries();
        checkSourceBoundaries();
        auto sbound_end = std::chrono::high_resolution_clock::now();
        source_boundary_time += sbound_end - sbound_start;
        break;
        
        // Calculate number of invalid consecutive pairs of chars
        int maxLeft = orec->pair.left;
        int maxRight = orec->pair.right;
        int max_freq = orec->freq;
        int invalidFreq = 0;
        if (verbosity == 2)
        {
            if (maxLeft == maxRight){
                checkPhraseSizes();
                spdlog::trace("Same Char Pair: ({},{})", printSymbol(maxLeft), printSymbol(maxRight));
                invalidFreq = invalidSameCharPair(orec->pair.left);
            }
        }

        std::pair<int,int> max_pair = {left_elem,right_elem};
        // Check if the current max pair is in hash ranges, if is then have to go through non-explicit phrases.
        if (hash_ranges.find(max_pair) != hash_ranges.end())
        {
            spdlog::trace("Going through the non-explicit phrases");
            std::deque<int> ranges = hash_ranges[max_pair];

            auto nexp_start = std::chrono::high_resolution_clock::now();
            for (int curr_range : ranges)
            {
                // Replace in Ref
                RefNode* lref = rarray[curr_range];
                RefNode* rref = rlist.findForwardRef(lref);

                // There is now no guarantee that the pair specified in the hash table actually exists due to lazily handling the hash table.
                if (lref == nullptr || rref == nullptr || lref->deleted || rref->deleted || lref->val != left_elem || rref->val != right_elem){
                    continue;
                }

                if (rarray[curr_range]->deleted){ 
                    continue;
                }  
                
                // Query the tree
                std::vector<PhraseNode*> phrase_results;
                spdlog::trace("Pair to replace: ({},{})", lref->pos, rref->pos);
                phrase_results = phrase_tree.findContained({lref->pos,rref->pos}); // Fully contained within interval
                spdlog::trace("{} non-explicit phrases contain the pair to replace.", phrase_results.size());
                for (int i = 0; i < phrase_results.size(); i++)
                {
                    PhraseNode* nexp_phrase = phrase_results[i];
                    int leftElem =  lref->val;
                    int rightElem = rref->val;
                    int leftleftElem;
                    int rightrightElem;

                    int lRange = rlist.findNearestRef(phrase_results[i]->lnode)->pos;
                    int rRange = rlist.findNearestRef(phrase_results[i]->rnode)->pos;

                    // If range fully contained within range then we can do the decrease and increase frequencies
                    if (lRange != lref->pos && rRange != rref->pos)
                    {
                        leftleftElem = rlist.findNearestRef(lref->prev)->val;
                        rightrightElem = rlist.findForwardRef(rref)->val;
                        decreaseFrequency(leftleftElem, leftElem);
                        increaseFrequency(leftleftElem, n);
                        decreaseFrequency(rightElem, rightrightElem);
                        increaseFrequency(n, rightrightElem);

                    }
                    // If range touches the edges
                    else if (lRange == lref->pos && rRange == rref->pos)
                    {
                        if (nexp_phrase->prev != nullptr && nexp_phrase->prev->exp){
                            leftleftElem = nexp_phrase->prev->content.back();
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        } 
                        else if (nexp_phrase->prev != nullptr && !nexp_phrase->prev->exp) {
                            if (nexp_phrase->prev->rtmp == -1)
                                leftleftElem = rlist.findNearestRef(nexp_phrase->prev->rnode)->val;
                            else
                                leftleftElem = nexp_phrase->prev->rtmp;
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        }
                        nexp_phrase->ltmp = n;
                        if (nexp_phrase->next != nullptr && nexp_phrase->next->exp){
                            rightrightElem = nexp_phrase->next->content.front();
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        } 
                        else if (nexp_phrase->next != nullptr && !nexp_phrase->next->exp) {
                            if (nexp_phrase->next->ltmp == -1)
                                rightrightElem = rlist.findNearestRef(nexp_phrase->next->lnode)->val;
                            else
                                rightrightElem = nexp_phrase->next->ltmp;
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        }
                        nexp_phrase->rtmp = n;
                    }
                    // If range touches left edge, we can only decrease the left pair at the moment.
                    else if (lRange == lref->pos && rRange != rref->pos)
                    {
                        if (nexp_phrase->prev != nullptr && nexp_phrase->prev->exp){
                            leftleftElem = nexp_phrase->prev->content.back();
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        } 
                        else if (nexp_phrase->prev != nullptr && !nexp_phrase->prev->exp) {
                            if (nexp_phrase->prev->rtmp == -1)
                                leftleftElem = rlist.findNearestRef(nexp_phrase->prev->rnode)->val;
                            else
                                leftleftElem = nexp_phrase->prev->rtmp;
                            nexp_phrase->ltmp = n;
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        }

                        rightrightElem = rlist.findForwardRef(rref)->val;
                        decreaseFrequency(rightElem, rightrightElem);
                        increaseFrequency(n, rightrightElem);
                    }
                    // If range touches right edge
                    else if (lRange != lref->pos && rRange == rref->pos)
                    {
                        leftleftElem = rlist.findNearestRef(lref->prev)->val;
                        decreaseFrequency(leftleftElem, leftElem);
                        increaseFrequency(leftleftElem, n);

                        if (nexp_phrase->next != nullptr && nexp_phrase->next->exp){
                            rightrightElem = nexp_phrase->next->content.front();
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        } 
                        else if (nexp_phrase->next != nullptr && !nexp_phrase->next->exp) {
                            if (nexp_phrase->next->ltmp == -1)
                                rightrightElem = rlist.findNearestRef(nexp_phrase->next->lnode)->val;
                            else
                                rightrightElem = nexp_phrase->next->ltmp;
                            nexp_phrase->rtmp = n;
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        }
                    }
                    else{
                        spdlog::error("Should not logically happen!");
                    }
                }
                // Lazy delete the adjacent pairs effected by the merge
                std::pair<int,int> hash_pair;
                RefNode* llref = rlist.findNearestRef(lref->prev);
                RefNode* rrref = rlist.findForwardRef(rref);
                auto hash_range_start = std::chrono::high_resolution_clock::now(); 
                if (llref != nullptr){
                    hash_pair.first = llref->val;
                    hash_pair.second = n;
                    hash_ranges[hash_pair].push_back(llref->pos);
                }
                if (rrref != nullptr){
                    hash_pair.first = n;
                    hash_pair.second = rrref->val;
                    hash_ranges[hash_pair].push_back(lref->pos);
                }
                auto hash_range_end = std::chrono::high_resolution_clock::now();
                hash_range_time += hash_range_end - hash_range_start;
                // Replace the pair in the reference.
                rlist.replacePair(n, lref, rref);
            }
            auto hash_range_start = std::chrono::high_resolution_clock::now();
            hash_ranges.erase(max_pair); // Delete the max pair in the hash table since we have done all the replacements in the nexp phrases.
            auto hash_range_end = std::chrono::high_resolution_clock::now();
            hash_range_time += hash_range_end - hash_range_start;
            auto nexp_end = std::chrono::high_resolution_clock::now();
            nonexplicit_phrase_time += nexp_end - nexp_start;
        }
        
        // Have to process the explicit phrases
        auto exp_start = std::chrono::high_resolution_clock::now();
        if (exp_pairs.find(max_pair) != exp_pairs.end())
        {
            spdlog::trace("Going through the explicit phrases");
            std::unordered_set<ExpPair,ExpPairHash,ExpPairEqual> max_pairs = exp_pairs[max_pair];
            for (auto it = max_pairs.begin(); it != max_pairs.end(); it++)
            {
                PhraseNode* curr_phrase = it->exp_phrase;
                auto leftIt = it->left;
                auto rightIt = it->right;
                std::list<int>::iterator leftleftIt;
                std::list<int>::iterator rightrightIt;

                // Some assertions
                assert(leftIt <= rightIt);
                assert(leftIt >= curr_phrase->content.begin() && rightIt < curr_phrase->content.end());
                assert(leftIt >= curr_phrase->content.begin() && rightIt < curr_phrase->content.end());

                bool rightAddConsecutive = false; // Adding from the right causes consecutive same chars 
                bool leftAddConsecutive = false; // Adding from the left causes consecutive same chars
                bool rightDeleteConsecutive = false;  // Deleting from the right causes shift of consecutive same chars 
                bool leftDeleteConsecutive = false; // Deleting from the left causes shift of consecutive same chars

                if (*leftIt == left_elem && *rightIt == right_elem)
                {
                    // Have to change the frequencies in max heap.
                    // If the two elements are not at the beginning or end of the phrase then replace happens fully within the phrase
                    if (leftIt != curr_phrase->content.begin() && rightIt != std::prev(curr_phrase->content.end()))
                    {
                        leftleftIt = std::prev(leftIt);
                        rightrightIt = std::next(rightIt);
                        // Decrease frequency of left pair effected by merge.
                        decreaseFrequency(*leftleftIt, *leftIt);
                        exp_pairs[{*leftleftIt, *leftIt}].erase(ExpPair(curr_phrase, leftleftIt, leftIt));
                        // Increase frequency of new pair.
                        increaseFrequency(*leftleftIt, n);
                        exp_pairs[{*leftleftIt,n}].insert(ExpPair(curr_phrase, leftleftIt, leftIt));
                        if (*leftleftIt == n){
                            rightAddConsecutive = true;
                        }
                        else if (*leftleftIt == *leftIt)
                        {
                            rightDeleteConsecutive = true;
                        }
                        // Decrease frequency of right pair effected by merge.
                        decreaseFrequency(*rightIt, *rightrightIt);
                        exp_pairs[{*rightIt, *rightrightIt}].erase(ExpPair(curr_phrase, rightIt, rightrightIt));
                        // Increase frequency of new pair.
                        increaseFrequency(n, *rightrightIt);
                        exp_pairs[{n, *rightrightIt}].insert(ExpPair(curr_phrase, leftIt, rightrightIt));
                        if (*rightrightIt == n){
                            leftAddConsecutive = true;
                        }
                        else if (*rightrightIt == *rightIt){
                            leftDeleteConsecutive = true;
                        }
                    }
                    // If both elments are at the beginning and end of the phrases, then have to look at the end of the prev and start of the next
                    else if (leftIt == curr_phrase->content.begin() && rightIt == std::prev(curr_phrase->content.end()))
                    {
                        // For the left pair effected have to look at previous phrase
                        if (curr_phrase != plist.getHead())
                        {
                            PhraseNode* prev_phrase = curr_phrase->prev;
                            if (prev_phrase->exp)
                            {
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(prev_phrase->content.back(), *leftIt);
                                // Increase frequency of new pair.
                                increaseFrequency(prev_phrase->content.back(), n);
                            }
                            else
                            {
                                int leftleftElem = rlist.findNearestRef(prev_phrase->rnode)->val;
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(leftleftElem, *leftIt);
                                // Increase frequency of new pair.
                                increaseFrequency(leftleftElem, n);
                            }
                        }
                        // For the right pair effected have to look at next phrase
                        if (curr_phrase != plist.getTail())
                        {
                            PhraseNode* next_phrase = curr_phrase->next;
                            if (next_phrase->exp)
                            {
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*rightIt, next_phrase->content.front());
                                // Increase frequency of new pair
                                increaseFrequency(n, next_phrase->content.front());
                            }
                            else{
                                int rightrightElem = rlist.findNearestRef(next_phrase->lnode)->val;
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*rightIt, rightrightElem);
                                // Increase frequency of new pair
                                increaseFrequency(n, rightrightElem);
                            }
                        }
                    }
                    // If the left elem is at the beginning of phrase and the right elem is in middle of phrase
                    else if (leftIt == curr_phrase->content.begin() && rightIt != std::prev(curr_phrase->content.end()))
                    {
                        // For the left pair effected have to look at previous phrase
                        if (curr_phrase != plist.getHead())
                        {
                            PhraseNode* prev_phrase = curr_phrase->prev;
                            if (prev_phrase->exp)
                            {
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(prev_phrase->content.back(), *leftIt);
                                // Increase frequency of new pair.
                                increaseFrequency(prev_phrase->content.back(), n);
                            }
                            else
                            {
                                int leftleftElem = rlist.findNearestRef(prev_phrase->rnode)->val;
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(leftleftElem, *(leftIt));
                                // Increase frequency of new pair.
                                increaseFrequency(leftleftElem, n);
                            }
                        }
                        // Decrease frequency of right pair effected by merge.
                        rightrightIt = std::next(rightIt);
                        decreaseFrequency(*rightIt, *rightrightIt);
                        exp_pairs[{*rightIt, *rightrightIt}].erase(ExpPair(curr_phrase, rightIt, rightrightIt));
                        // Increase frequency of new pair.
                        increaseFrequency(n, *rightrightIt);
                        exp_pairs[{n, *rightrightIt}].insert(ExpPair(curr_phrase, leftIt, rightrightIt));
                        if (*rightrightIt == n){
                            leftAddConsecutive = true;
                        }
                        else if (*rightrightIt == *rightIt){
                            leftDeleteConsecutive = true;
                        }
                    }
                    // If the left elem is in the middle of the phrase and the right elem is at the end of the phrase
                    else if (leftIt != curr_phrase->content.begin() && rightIt == std::prev(curr_phrase->content.end()))
                    {
                        // Decrease frequency of left pair effected by merge.
                        leftleftIt = std::prev(leftIt);
                        decreaseFrequency(*leftleftIt, *leftIt);
                        exp_pairs[{*leftleftIt, *leftIt}].erase(ExpPair(curr_phrase, leftleftIt, leftIt));
                        // Increase frequency of new pair.
                        increaseFrequency(*leftleftIt, n);
                        exp_pairs[{*leftleftIt, n}].insert(ExpPair(curr_phrase, leftleftIt, leftIt));
                        if (*leftleftIt == n){
                            rightAddConsecutive = true;
                        }
                        else if (*leftleftIt == *leftIt){
                            rightDeleteConsecutive = true;
                        }

                        // For the right pair effected have to look at next phrase
                        if (curr_phrase != plist.getTail())
                        {
                            PhraseNode* next_phrase = curr_phrase->next;
                            if (next_phrase->exp)
                            {
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*rightIt, next_phrase->content.front());
                                // Increase frequency of new pair
                                increaseFrequency(n, next_phrase->content.front());
                            }
                            else
                            {
                                int rightrightElem = rlist.findNearestRef(next_phrase->lnode)->val;
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*rightIt, rightrightElem);
                                // Increase frequency of new pair
                                increaseFrequency(n, rightrightElem);
                            }
                        }
                    } 
                    // Replace current elem val with n and remove next elem 
                    *leftIt = n;
                    curr_phrase->content.erase(rightIt);

                    // Potentially fix consecutive same char pairs in exp_pairs
                    if (leftAddConsecutive){
                        updateExpPairs(curr_phrase, leftIt, true);
                    }
                    if (leftDeleteConsecutive){
                        updateExpPairs(curr_phrase, std::next(leftIt), true);
                    }
                    if (rightAddConsecutive){
                        updateExpPairs(curr_phrase, leftIt, false);
                    }
                    if (rightDeleteConsecutive){
                        updateExpPairs(curr_phrase, std::prev(leftIt), false);
                    }
                }
            }
            exp_pairs.erase(max_pair);
        }
        auto exp_end = std::chrono::high_resolution_clock::now();
        explicit_phrase_time += exp_end - exp_start;

        // Check the number of chars replaced is correct
        if (verbosity == 2)
        {
            spdlog::trace("----------------------------");
            spdlog::trace("Max Pair: ({},{})", printSymbol(maxLeft), printSymbol(maxRight));
            int phrase_length = checkPhraseSizes();
            int num_pairs_replaced = start_size - phrase_length;
            max_freq = max_freq - invalidFreq;
            spdlog::trace("Number of occurences: {}", max_freq);
            spdlog::trace("Number replaced: {}", num_pairs_replaced);
            if (num_pairs_replaced != max_freq){
                spdlog::error("Something is wrong");
            } 
            start_size = phrase_length;
            spdlog::trace("----------------------------");
        }

        // Remove old record
        removeRecord(&Rec,oid);

        // Remove frequency 1 records
        purgeHeap(&Heap);

        // Get next value in max heap
        oid = extractMax(&Heap);

        // Update n
        n++;

        // Debug
        if (verbosity == 2){
            spdlog::trace("");
            spdlog::trace("*** Information after bi-gram replacement ***");
            printRef();
            printPhraseList();
            //printAllRecords();
            checkExpPairs();
            checkHeap();
            //phrase_tree.printTree();
            spdlog::trace("*********************************************");
        }
    }

    // Write the final integers to the C file
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        if (!(curr_phrase->exp)){
            // Go backwards via the prev pointers but write in forward direction
            std::stack<int> nexp_stack;
            RefNode* end_elem = rlist.findNearestRef(curr_phrase->rnode);
            RefNode* first_elem = rlist.findNearestRef(rlist.findNearestRef(curr_phrase->lnode)->prev);
            while(end_elem != first_elem){
                nexp_stack.push(end_elem->val);
                end_elem = rlist.findNearestRef(end_elem->prev);
            }
            while (!nexp_stack.empty()) {
                c++;
                unsigned int i = nexp_stack.top();
                C.write(reinterpret_cast<const char*>(&i), sizeof(unsigned int));
                nexp_stack.pop();
            }
        }
        else{
            for (unsigned int i : curr_phrase->content){
                c++;
                C.write(reinterpret_cast<const char*>(&i), sizeof(unsigned int));
            }
        }
        curr_phrase = curr_phrase->next;
    }
}


int main(int argc, char *argv[])
{
    CLI::App app("rlz - Run RePair with the RLZ parse.\n\nImplemented by Rahul Varki");
    std::string ref_file;
    std::string rlz_parse;
    std::string version = "Version: 1.0.0";

    app.add_option("-r,--ref", ref_file, "The reference file used to create the RLZ parse")->configurable()->required();
    app.add_option("-p,--parse", rlz_parse, "The RLZ parse of the sequence file (.rlz)")->configurable()->required();
    app.add_option("-v,--verbosity", verbosity, "Set verbosity level (0 = none, 1 = basic, 2 = detailed)")->default_val(0);
    app.set_version_flag("--version", version);
    app.footer("Example usage:\n"
           "  Compress: ./repair --ref reference.fasta --parse sequence.rlz\n");
    app.description("Run RePair on RLZ parse");
    CLI11_PARSE(app, argc, argv);

    if (verbosity == 2) {
        spdlog::set_level(spdlog::level::trace);
    }
    else if (verbosity == 1){
        spdlog::set_level(spdlog::level::debug);
    }
    else if (verbosity == 0){ 
        spdlog::set_level(spdlog::level::info);
    }
    else{
        spdlog::error("Set verbosity level to be 0,1,2");
    }

    spdlog::debug("The reference file provided: {}", ref_file);
    spdlog::debug("The RLZ parse file provided: {}", rlz_parse);
    auto total_time_start = std::chrono::high_resolution_clock::now();

    // Opening Ref file 
    std::ifstream rfile(ref_file, std::ios::binary | std::ios::in);
    if (!rfile) {
        spdlog::error("Error opening {}", ref_file);
        std::exit(EXIT_FAILURE);
    }

    // Get Ref file size
    rfile.seekg(0, std::ios::end);
    size_t rsize = rfile.tellg();
    rfile.seekg(0, std::ios::beg);

    // Read the Ref file into a vector of unsigned char
    std::vector<unsigned char> rtext(rsize);
    rfile.read(reinterpret_cast<char*>(rtext.data()), rsize);

    // Opening RLZ parse
    std::ifstream pfile(rlz_parse, std::ios::binary | std::ios::in);
    if (!pfile) {
        spdlog::error("Error opening {}", rlz_parse);
        std::exit(EXIT_FAILURE);
    }

    // Create R file
    size_t last_dot = rlz_parse.find_last_of('.');
    std::string Rf = rlz_parse.substr(0, last_dot) + ".R"; 
    std::ofstream R(Rf, std::ios::binary);
    if (!R) {
        spdlog::error("Error opening {}", Rf);
        std::exit(EXIT_FAILURE);
    }

    // Create C file
    std::string Cf = rlz_parse.substr(0, last_dot) + ".C";
    std::ofstream C(Cf, std::ios::binary);
    if (!C) {
        spdlog::error("Error opening {}", Cf);
        std::exit(EXIT_FAILURE);
    }

    // Process Ref
    auto prepare_ref_start = std::chrono::high_resolution_clock::now();
    prepareRef(rtext);
    auto prepare_ref_end = std::chrono::high_resolution_clock::now();
    prepare_ref_time += prepare_ref_end - prepare_ref_start;

    // Clear char representation of Ref
    rtext.clear();
    rtext.shrink_to_fit();

    // Calculate the size of the original sequence file using the len field of the pairs
    auto calculate_size_start = std::chrono::high_resolution_clock::now();
    psize = calculateParseBytes(pfile);
    auto calculate_size_end = std::chrono::high_resolution_clock::now();
    calculate_size_time += calculate_size_end - calculate_size_start;

    // Call populatePhrases function
    auto parse_phrase_start = std::chrono::high_resolution_clock::now();
    populatePhrases(pfile);
    auto parse_phrase_end = std::chrono::high_resolution_clock::now();
    populate_phrase_time += parse_phrase_end - parse_phrase_start;

    // Call createMaxHeap function
    auto create_heap_start = std::chrono::high_resolution_clock::now();
    createMaxHeap(pfile);
    auto create_heap_end = std::chrono::high_resolution_clock::now();
    create_heap_time += create_heap_end - create_heap_start;

    // RePair
    repair(R, C);

    auto total_time_end = std::chrono::high_resolution_clock::now();
    spdlog::debug("Total Prepare Reference Time (s): {:.6f}", std::chrono::duration<double>(prepare_ref_time).count());
    spdlog::debug("Total Calculate Parse Size Time (s): {:.6f}", std::chrono::duration<double>(calculate_size_time).count());
    spdlog::debug("Total Populate Phrase Time (s): {:.6f}", std::chrono::duration<double>(populate_phrase_time).count());
    spdlog::debug("Total Create Max Heap Time (s): {:.6f}", std::chrono::duration<double>(create_heap_time).count());
    spdlog::debug("Total Phrase Boundary Time (s): {:.6f}", std::chrono::duration<double>(phrase_boundary_time).count());
    spdlog::debug("Total Source Boundary Time (s): {:.6f}", std::chrono::duration<double>(source_boundary_time).count());
    spdlog::debug("Total Build Interval Tree Time (s): {:.6f}", std::chrono::duration<double>(build_interval_time).count());
    spdlog::debug("Total Update Interval Tree Time (s): {:.6f}", std::chrono::duration<double>(update_interval_time).count());
    spdlog::debug("Total Non-explicit Phrase Time (s): {:.6f}", std::chrono::duration<double>(nonexplicit_phrase_time).count());
    spdlog::debug("Total Explicit Phrase Time (s): {:.6f}", std::chrono::duration<double>(explicit_phrase_time).count());
    spdlog::debug("Total Hash Range Update Time (s): {:.6f}", std::chrono::duration<double>(hash_range_time).count());
    total_time = total_time_end - total_time_start;
    spdlog::debug("*********************************************");
    spdlog::debug("Total Time (s): {:.6f}", std::chrono::duration<double>(total_time).count());
    
    R.close();
    C.close();
    rfile.close();
    pfile.close();

    // Statistics
    spdlog::info("File Statistics");
    spdlog::info("Parse Encoded Chars: {}", psize);
    spdlog::info("Number of Rules: {}", n - alpha);
    spdlog::info("Final sequence length: {}", c);

    return 0;
}