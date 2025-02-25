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
 * @brief Calculate number of invalid consecutive same chars in the phrases
 * For example, eeee -> 1 (the 2nd ee)
 * For example, aaa -> 1 (the 2nd aa)
 * For example, iiiii -> 2 (2nd and 4th ii)
 */

int invalidSameCharPair(int letter)
{
    PhraseNode* curr_phrase = plist.getHead();
    int invalidCount = 0;
    int count = 0;
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
                //spdlog::debug("Pair: ({},{})", printSymbol(left_elem), printSymbol(right_elem));
                if (left_elem == letter && right_elem == letter){
                    count++;
                }
                else{
                    count = 0;
                }
                if (count > 0 && count % 2 == 0){
                    invalidCount++;
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
                //spdlog::debug("Pair: ({},{})", printSymbol(left_elem), printSymbol(right_elem));
                if (left_elem == letter && right_elem == letter){
                    count++;
                }
                else{
                    count = 0;
                }
                if (count > 0 && count % 2 == 0){
                    invalidCount++;
                }
                left_elem = right_elem;
                it++;
            }
        } 
        curr_phrase = curr_phrase->next;
    }
    return invalidCount;
}

/**
 * @brief 
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
                    if (exp_pairs[{*pit,*it}].find(ExpPair(curr_phrase, pit, it)) == exp_pairs[{*pit,*it}].end()){
                        printPhrase(curr_phrase);
                        spdlog::debug("Pair missing: ({},{})", printSymbol(*pit), printSymbol(*it));
                        // std::cout << "Address of phrase: " << &(*curr_phrase) << std::endl;
                        // std::cout << "Address of left element pointed to: " << &(*pit) << std::endl;
                        // std::cout << "Address of right element pointed to: " << &(*it) << std::endl;
                        return false;
                    }
                }
                else
                {
                    if (exp_pairs[{*pit,*it}].find(ExpPair(curr_phrase, pit, it)) != exp_pairs[{*pit,*it}].end()){
                        printPhrase(curr_phrase);
                        spdlog::debug("Pair exists that should not exist: ({},{})", printSymbol(*pit), printSymbol(*it));
                        return false;
                    }
                }
                pit = it;
                it = std::next(it);
            }
        }        
        curr_phrase = curr_phrase->next;
    }
    spdlog::debug("Everything is correct");
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
        }
    }

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
        auto prevIt = std::prev(it);
        while (prevIt != p->content.begin())
        {
            if (*prevIt != letter){
                prevIt = std::next(prevIt);
                break;
            }
            prevIt = std::prev(prevIt);
        }
        if (prevIt == p->content.begin() && *prevIt != letter){
            prevIt = std::next(prevIt);
        }
        auto nextIt = std::next(prevIt);
        while(nextIt != p->content.end())
        {
            if (*nextIt != letter){
                break;
            }
            else{
                sameCharCount++;
                if (sameCharCount % 2 == 0){
                    exp_pairs[{letter,letter}].insert(ExpPair(p, prevIt, nextIt));
                }
                else{
                    exp_pairs[{letter,letter}].erase(ExpPair(p, prevIt, nextIt));
                }
            }
            prevIt = nextIt;
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

void phraseBoundaries(int left_elem, int right_elem)
{
    PhraseNode* curr_phrase = plist.getHead();

    // Iterate through the phrases in the phrase list
    while (curr_phrase != nullptr) 
    {
        curr_phrase->ltmp = -1;
        curr_phrase->rtmp = -1;

        PhraseNode* next_phrase = curr_phrase->next;
        // If there is a next phrase, check the phrase boundaries.
        if (next_phrase != nullptr)
        {
            // Both phrases not explicit
            if (!(curr_phrase->exp) && !(next_phrase->exp))
            {
                if (rlist.findNearestRef(curr_phrase->rnode)->val == left_elem && rlist.findNearestRef(next_phrase->lnode)->val == right_elem)
                {
                    // Indicates whether a non-explicit phrase got deleted.
                    bool deleteCurr = false;
                    bool deleteNext = false;

                    std::list<int> content;
                    content.push_back(left_elem);
                    content.push_back(right_elem);

                    // First remove the offending entries in the tree.
                    spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
                    spdlog::trace("Removing ({},{}) from tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);
                    next_phrase->lnode = rlist.findNearestRef(next_phrase->lnode);
                    next_phrase->rnode = rlist.findNearestRef(next_phrase->rnode);

                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);
                    next_phrase->lnode = rlist.findForwardRef(next_phrase->lnode);
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
                        spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                        spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                    }
                    else if (deleteCurr && !deleteNext){
                        spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                    }
                    else if (!deleteCurr && deleteNext){
                        spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
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
                    spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

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
                        spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
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
                    spdlog::trace("Removing ({},{}) from tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                    auto update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.remove({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                    auto update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;

                    // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                    next_phrase->lnode = rlist.findNearestRef(next_phrase->lnode);
                    next_phrase->rnode = rlist.findNearestRef(next_phrase->rnode);

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
                        spdlog::trace("Adding ({},{}) to the tree", next_phrase->lnode->pos, next_phrase->rnode->pos);
                        update_interval_start = std::chrono::high_resolution_clock::now();
                        phrase_tree.insert({next_phrase->lnode->pos, next_phrase->rnode->pos}, next_phrase);
                        update_interval_end = std::chrono::high_resolution_clock::now();
                        update_interval_time += update_interval_end - update_interval_start;
                    }
                    continue;
                }
            }
            // Both explicit phrases, move content of next phrase into first phrase and delete next phrase
            else{
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
            curr_phrase = curr_phrase->next;
        }
        else{ // No more phrase boundaries left
            break;
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
                }

                // Indicates whether a non-explicit phrase got deleted.
                bool deleteCurr = false;
                // First remove the offending entry from the tree
                spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                auto update_interval_start = std::chrono::high_resolution_clock::now();
                phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                auto update_interval_end = std::chrono::high_resolution_clock::now();
                update_interval_time += update_interval_end - update_interval_start;

                // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

                // Update the information of the current phrase (left node pointer & existence).
                curr_phrase->lnode = rlist.findForwardRef(curr_phrase->lnode);                    
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->lnode == nullptr || curr_phrase->lnode->pos > curr_phrase->rnode->pos){ 
                    deleteCurr = true;
                    curr_phrase = plist.remove(curr_phrase); // Removes the current phrase and sets the curr phrase to be the previous phrase
                }

                // Update the tree with the new entry
                if (!deleteCurr){
                    spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
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
                    }
                }

                // Indicates whether a non-explicit phrase got deleted.
                bool deleteCurr = false;
                // First remove the offending entry from the tree
                spdlog::trace("Removing ({},{}) from tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                auto update_interval_start = std::chrono::high_resolution_clock::now();
                phrase_tree.remove({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                auto update_interval_end = std::chrono::high_resolution_clock::now();
                update_interval_time += update_interval_end - update_interval_start;

                // Only when modifying the tree can the lnode and rnode pointers of the phrases be updated.
                curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);

                // Update the information of the current phrase (right node pointer & existence).
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);                       
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->rnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){
                    deleteCurr = true; 
                    curr_phrase = plist.remove(curr_phrase); // Removes the current phrase and sets the curr phrase to be the previous phrase
                }

                // Update the tree with the new entry
                if (!deleteCurr){
                    spdlog::trace("Adding ({},{}) to the tree", curr_phrase->lnode->pos, curr_phrase->rnode->pos);
                    update_interval_start = std::chrono::high_resolution_clock::now();
                    phrase_tree.insert({curr_phrase->lnode->pos, curr_phrase->rnode->pos}, curr_phrase);
                    update_interval_end = std::chrono::high_resolution_clock::now();
                    update_interval_time += update_interval_end - update_interval_start;
                }

                continue; // If it gets to here then the next iteration will have the same curr phrase.
            }
        }
        // If explicit phrase then check surrounding phrases if they are explicit and if so merge.
        else{
            PhraseNode* prev_phrase = curr_phrase->prev;
            // If prev phrase is also explicit merge into the prev phrase
            if (prev_phrase != nullptr && prev_phrase->exp){
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
    if (verbosity == 2){
        spdlog::trace("Decrease frequency: ({},{})", printSymbol(left), printSymbol(right));
    }
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
    if (verbosity == 2){
        spdlog::trace("Increase frequency: ({},{})", printSymbol(left), printSymbol(right));
    }
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

        // Calculate number of invalid consecutive pairs of chars
        int invalidFreq = 0;
        if (orec->pair.left == orec->pair.right){
            invalidFreq = invalidSameCharPair(orec->pair.left);
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
        auto pbound_end = std::chrono::high_resolution_clock::now();
        phrase_boundary_time += pbound_end - pbound_start;

        auto sbound_start = std::chrono::high_resolution_clock::now();
        sourceBoundaries(left_elem, right_elem);
        auto sbound_end = std::chrono::high_resolution_clock::now();
        source_boundary_time += sbound_end - sbound_start;
        
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
                auto curr_elem = it->left;
                auto next_elem = it->right;
                bool rightAddConsecutive = false; // Adding from the right causes consecutive same chars 
                bool leftAddConsecutive = false; // Adding from the left causes consecutive same chars
                bool rightDeleteConsecutive = false;  // Deleting from the right causes shift of consecutive same chars 
                bool leftDeleteConsecutive = false; // Deleting from the left causes shift of consecutive same chars

                if (*curr_elem == left_elem && *next_elem == right_elem)
                {
                    // Have to change the frequencies in max heap.
                    // If the two elements are not at the beginning or end of the phrase then replace happens fully within the phrase
                    if (curr_elem != curr_phrase->content.begin() && next_elem != std::prev(curr_phrase->content.end()))
                    {
                        // Decrease frequency of left pair effected by merge.
                        decreaseFrequency(*(std::prev(curr_elem)), *(curr_elem));
                        exp_pairs[{*(std::prev(curr_elem)), *(curr_elem)}].erase(ExpPair(curr_phrase, std::prev(curr_elem), curr_elem));
                        // Increase frequency of new pair.
                        increaseFrequency(*(std::prev(curr_elem)), n);
                        exp_pairs[{*(std::prev(curr_elem)),n}].insert(ExpPair(curr_phrase, std::prev(curr_elem), curr_elem));
                        if (*(std::prev(curr_elem)) == n){
                            rightAddConsecutive = true;
                        }
                        else if (*(std::prev(curr_elem)) == *(curr_elem))
                        {
                            rightDeleteConsecutive = true;
                        }
                        // Decrease frequency of right pair effected by merge.
                        decreaseFrequency(*(next_elem), *(std::next(next_elem)));
                        exp_pairs[{*(next_elem), *(std::next(next_elem))}].erase(ExpPair(curr_phrase, next_elem, std::next(next_elem)));
                        // Increase frequency of new pair.
                        increaseFrequency(n, *(std::next(next_elem)));
                        exp_pairs[{n, *(std::next(next_elem))}].insert(ExpPair(curr_phrase, curr_elem, std::next(next_elem)));
                        if (*(std::next(next_elem)) == n){
                            leftAddConsecutive = true;
                        }
                        else if (*(std::next(next_elem)) == *(next_elem)){
                            leftDeleteConsecutive = true;
                        }
                    }
                    // If both elments are at the beginning and end of the phrases, then have to look at the end of the prev and start of the next
                    else if (curr_elem == curr_phrase->content.begin() && next_elem == std::prev(curr_phrase->content.end()))
                    {
                        // For the left pair effected have to look at previous phrase
                        if (curr_phrase != plist.getHead())
                        {
                            PhraseNode* prev_phrase = curr_phrase->prev;
                            if (prev_phrase->exp)
                            {
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(prev_phrase->content.back(), *(curr_elem));
                                // Increase frequency of new pair.
                                increaseFrequency(prev_phrase->content.back(), n);
                            }
                            else
                            {
                                int leftleftElem = rlist.findNearestRef(prev_phrase->rnode)->val;
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(leftleftElem, *(curr_elem));
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
                                decreaseFrequency(*(next_elem), next_phrase->content.front());
                                // Increase frequency of new pair
                                increaseFrequency(n, next_phrase->content.front());
                            }
                            else{
                                int rightrightElem = rlist.findNearestRef(next_phrase->lnode)->val;
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*(next_elem), rightrightElem);
                                // Increase frequency of new pair
                                increaseFrequency(n, rightrightElem);
                            }
                        }
                    }
                    // If the left elem is at the beginning of phrase and the right elem is in middle of phrase
                    else if (curr_elem == curr_phrase->content.begin() && next_elem != std::prev(curr_phrase->content.end()))
                    {
                        // For the left pair effected have to look at previous phrase
                        if (curr_phrase != plist.getHead())
                        {
                            PhraseNode* prev_phrase = curr_phrase->prev;
                            if (prev_phrase->exp)
                            {
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(prev_phrase->content.back(), *(curr_elem));
                                // Increase frequency of new pair.
                                increaseFrequency(prev_phrase->content.back(), n);
                            }
                            else
                            {
                                int leftleftElem = rlist.findNearestRef(prev_phrase->rnode)->val;
                                // Decrease frequency of left pair effected by merge.
                                decreaseFrequency(leftleftElem, *(curr_elem));
                                // Increase frequency of new pair.
                                increaseFrequency(leftleftElem, n);
                            }
                        }
                        // Decrease frequency of right pair effected by merge.
                        decreaseFrequency(*(next_elem), *(std::next(next_elem)));
                        exp_pairs[{*(next_elem), *(std::next(next_elem))}].erase(ExpPair(curr_phrase, next_elem, std::next(next_elem)));
                        // Increase frequency of new pair.
                        increaseFrequency(n, *(std::next(next_elem)));
                        exp_pairs[{n, *(std::next(next_elem))}].insert(ExpPair(curr_phrase, curr_elem, std::next(next_elem)));
                        if (*(std::next(next_elem)) == n){
                            leftAddConsecutive = true;
                        }
                        else if (*(std::next(next_elem)) == *(next_elem)){
                            leftDeleteConsecutive = true;
                        }
                    }
                    // If the left elem is in the middle of the phrase and the right elem is at the end of the phrase
                    else if (curr_elem != curr_phrase->content.begin() && next_elem == std::prev(curr_phrase->content.end()))
                    {
                        // Decrease frequency of left pair effected by merge.
                        decreaseFrequency(*(std::prev(curr_elem)), *(curr_elem));
                        exp_pairs[{*(std::prev(curr_elem)), *(curr_elem)}].erase(ExpPair(curr_phrase, std::prev(curr_elem), curr_elem));
                        // Increase frequency of new pair.
                        increaseFrequency(*(std::prev(curr_elem)), n);
                        exp_pairs[{*(std::prev(curr_elem)), n}].insert(ExpPair(curr_phrase, std::prev(curr_elem), curr_elem));
                        if (*(std::prev(curr_elem)) == n){
                            rightAddConsecutive = true;
                        }
                        else if (*(std::prev(curr_elem)) == *(curr_elem)){
                            rightDeleteConsecutive = true;
                        }

                        // For the right pair effected have to look at next phrase
                        if (curr_phrase != plist.getTail())
                        {
                            PhraseNode* next_phrase = curr_phrase->next;
                            if (next_phrase->exp)
                            {
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*(next_elem), next_phrase->content.front());
                                // Increase frequency of new pair
                                increaseFrequency(n, next_phrase->content.front());
                            }
                            else
                            {
                                int rightrightElem = rlist.findNearestRef(next_phrase->lnode)->val;
                                // Decrease frequency of right pair effected by merge.
                                decreaseFrequency(*(next_elem), rightrightElem);
                                // Increase frequency of new pair
                                increaseFrequency(n, rightrightElem);
                            }
                        }
                    } 
                    // Replace current elem val with n and remove next elem 
                    *curr_elem = n;
                    curr_phrase->content.erase(next_elem);

                    // Potentially fix consecutive same char pairs in exp_pairs
                    if (leftAddConsecutive){
                        updateExpPairs(curr_phrase, curr_elem, true);
                    }
                    if (leftDeleteConsecutive){
                        updateExpPairs(curr_phrase, std::next(curr_elem), true);
                    }
                    if (rightAddConsecutive){
                        updateExpPairs(curr_phrase, curr_elem, false);
                    }
                    if (rightDeleteConsecutive){
                        updateExpPairs(curr_phrase, std::prev(curr_elem), false);
                    }
                }
            }
            exp_pairs.erase(max_pair);
        }
        auto exp_end = std::chrono::high_resolution_clock::now();
        explicit_phrase_time += exp_end - exp_start;

        // Check the number of chars replaced is correct
        spdlog::debug("----------------------------");
        printMaxPair(n, orec);
        int phrase_length = checkPhraseSizes();
        int num_pairs_replaced = start_size - phrase_length;
        int max_freq = orec->freq - invalidFreq;
        spdlog::debug("Number of occurences: {}", max_freq);
        spdlog::debug("Number replaced: {}", num_pairs_replaced);
        if (num_pairs_replaced != max_freq){
            spdlog::error("Something is wrong");
        } 
        start_size = phrase_length;
        spdlog::debug("----------------------------");

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
            printAllRecords();
            //checkExpPairs();
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