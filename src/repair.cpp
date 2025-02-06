#include <CLI11.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <unordered_map>
#include <list>
#include <functional>
#include <utility>
#include <stack> 
#include "repair.h"
#include "phrase.h"
#include "IntervalTree.h"
#include "doublelinkedlist.h"

extern "C" {
    #include "heap.h"
    #include "array.h"
    #include "basics.h"
    #include "records.h"
    #include "hash.h"
}

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
std::unordered_map<std::pair<int, int>, std::vector<int>, pair_int_hash> hash_ranges; 

// List of explicit and non explicit phrases
PhraseLinkedList plist;

// Interval Tree for the non-explicit phrases
typedef IntervalTree<int, PhraseNode*> ITree;
ITree phrase_tree;

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
            for (unsigned int i : curr_phrase->content){
                content += printSymbol(i) + " ";
            }
            spdlog::trace("Phrase (explicit): \033[1;31m{}\033[0m", content);
        }
        curr_phrase = curr_phrase->next;
    }
    spdlog::trace("");
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
    for (int i = 0; i < rtext.size(); i++){
        unsigned char x = rtext[i];
        if (chars[x] == -1){
            chars[x] = alpha++;
        }
        curr_elem = rlist.push_back(chars[x]);
        ref_pair.second = curr_elem->val;
        if (i != 0){
            hash_ranges[ref_pair].emplace_back(i-1);
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

    spdlog::trace("Reference at the start");
    spdlog::trace("Size: {}", rlist.getSize());
    printRef();
    spdlog::trace("Hash ranges of reference bi-grams at the start");
    printHashRanges();
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
                    rnode = rlist.findPos(pos);
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
    printAllRecords();

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
            plist.push_back(rlist.findPos(pos), rlist.findPos(pos+len-1));
        }
    }

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);

    // Debug
    spdlog::trace("The non-explicit phrases at the start");
    printPhraseList();
}

/**
 * @brief Builds interval tree from the non-explicit phrases.
 */

void buildIntervalTree()
{
    ITree::interval_vector phrase_intervals;
    PhraseNode* curr_phrase = plist.getHead();
    while(curr_phrase != nullptr)
    {
        if (!curr_phrase->exp){
            phrase_intervals.push_back(ITree::interval(rlist.findNearestRef(curr_phrase->lnode)->pos, 
                                                        rlist.findNearestRef(curr_phrase->rnode)->pos, 
                                                        curr_phrase));
        }
        curr_phrase = curr_phrase->next;
    }
    
    ITree temp_tree(std::move(phrase_intervals));
    phrase_tree = temp_tree;

    // Some debug test
    // ITree::interval_vector phrase_results;
    // phrase_results = phrase_tree.findContained(0,4);
    // spdlog::info("Results of finding (0,4): {}", phrase_results.size());
    
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
        PhraseNode* next_phrase = curr_phrase->next;
        // If there is a next phrase, check the phrase boundaries.
        if (next_phrase != nullptr)
        {
            // Update the rnode and lnode of both the current and next phrase before doing anything.
            if (!curr_phrase->exp){
                curr_phrase->lnode = rlist.findNearestRef(curr_phrase->lnode);
                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode);
            }
            if (!next_phrase->exp){
                next_phrase->lnode = rlist.findNearestRef(next_phrase->lnode);
                next_phrase->rnode = rlist.findNearestRef(next_phrase->rnode);
            }
            
            // Both phrases not explicit
            if (!(curr_phrase->exp) && !(next_phrase->exp))
            {
                if (curr_phrase->rnode->val == left_elem && next_phrase->lnode->val == right_elem)
                {
                    std::list<unsigned int> content;
                    content.push_back(left_elem);
                    content.push_back(right_elem);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);
                    next_phrase->lnode = rlist.findForwardRef(next_phrase->lnode);
                    plist.insert(next_phrase, content); 
                    // If the current or next phrases are empty we delete
                    if (curr_phrase->rnode == nullptr || curr_phrase->lnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){ 
                        curr_phrase = plist.remove(curr_phrase);
                    }
                    if (next_phrase->rnode == nullptr || next_phrase->lnode == nullptr || next_phrase->lnode->pos > next_phrase->rnode->pos){
                        next_phrase = plist.remove(next_phrase);  
                    }
                    continue;
                }
            }
            // Current phrase not explicit and next phrase explicit
            else if (!(curr_phrase->exp) && (next_phrase->exp))
            {
                if ((curr_phrase->rnode->val == left_elem) && (next_phrase->content.front() == right_elem))
                {
                    next_phrase->content.push_front(left_elem);
                    curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);
                    // If the non-explicit phrase is empty we delete it.
                    if (curr_phrase->rnode == nullptr || curr_phrase->lnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){ 
                        curr_phrase = plist.remove(curr_phrase);
                    }
                    continue;
                }
            }
            // Current phrase explicit and next phrase not explicit
            else if ((curr_phrase->exp) && !(next_phrase->exp))
            {
                if ((curr_phrase->content.back() == left_elem) && (next_phrase->lnode->val == right_elem))
                {
                    curr_phrase->content.push_back(right_elem);
                    next_phrase->lnode = rlist.findForwardRef(next_phrase->lnode);
                    if (next_phrase->rnode == nullptr || next_phrase->lnode == nullptr || next_phrase->lnode->pos > next_phrase->rnode->pos){
                        next_phrase = plist.remove(next_phrase);
                    }
                    continue;
                }
            }
            // Both explicit phrases, move content of next phrase into first phrase and delete next phrase
            else{
                curr_phrase->content.splice(curr_phrase->content.end(), next_phrase->content);
                next_phrase = plist.remove(next_phrase);
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
    spdlog::trace("Phrase list after phrase boundary condition.");
    printPhraseList();
}

/**
 * @brief Process the phrase list for the source boundary condition.
 * 
 * If the rightmost elem of a non-explicit phrase and the next elem in the reference form the provided bi-gram
 * or the leftmost elem of a non-explicit phrase and the previous elem in the reference form the provided bi-gram
 * then remove the offending elem from the non-explicit phrase and create an explicit phrase.
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

            // If the bi-gram can be formed then make the leftmost element make an explicit phrase of its own
            if (leftElem != nullptr && rightElem != nullptr && leftElem->val == left_elem && rightElem->val == right_elem){
                PhraseNode* prev_phrase = curr_phrase->prev;
                // If the previous phrase is explicit, we can add directly to it.
                if (prev_phrase != nullptr && prev_phrase->exp){
                    prev_phrase->content.push_back(right_elem);
                }
                // We have to create new explicit phrase anyways
                else{
                    std::list<unsigned int> content;
                    content.push_back(right_elem);
                    plist.insert(curr_phrase, content);
                }
                
                curr_phrase->lnode = rlist.findForwardRef(curr_phrase->lnode);                    
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->lnode == nullptr || curr_phrase->lnode->pos > curr_phrase->rnode->pos){ 
                    curr_phrase = plist.remove(curr_phrase);
                    continue;
                }
            }
            
            // Check if the current phrase with its rightmost elem can form bi-gram with right elem on reference
            leftElem = rlist.findNearestRef(curr_phrase->rnode);
            rightElem = rlist.findForwardRef(leftElem);
            
            // If the bi-gram can be formed then make the rightmost element make an explicit phrase of its own
            if (leftElem != nullptr && rightElem != nullptr && leftElem->val == left_elem && rightElem->val == right_elem){
                PhraseNode* next_phrase = curr_phrase->next;                    
                // If the next phrase is explicit, we can add directly to it.
                if (next_phrase != nullptr && next_phrase->exp){
                    next_phrase->content.push_front(left_elem);
                }
                // We have to create new explicit phrase anyways
                else{
                    std::list<unsigned int> content;
                    content.push_back(left_elem);
                    if (next_phrase == nullptr){
                        plist.push_back(content);
                    }
                    else{
                        plist.insert(next_phrase, content);
                    }
                }

                curr_phrase->rnode = rlist.findNearestRef(curr_phrase->rnode->prev);                       
                // If non-explicit phrase is empty we delete it
                if (curr_phrase->rnode == nullptr || curr_phrase->rnode->pos < curr_phrase->lnode->pos){ 
                    curr_phrase = plist.remove(curr_phrase);
                    continue;
                }
            }
        }
        curr_phrase = curr_phrase->next;
    }
    // Debug
    spdlog::trace("Phrase list after source boundary condition.");
    printPhraseList();
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
    spdlog::trace("Decrease frequency: ({},{})", printSymbol(left), printSymbol(right));
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
    spdlog::trace("Increase frequency: ({},{})", printSymbol(left), printSymbol(right));
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

    oid = extractMax(&Heap);
    while (oid != -1)
    {
        Trecord* orec = &Rec.records[oid];
        // When max frequency is 1, RePair ends.
        if (orec->freq == 1){
            break;
        }

        printMaxPair(n, orec);

        // Write pair to R file 
        R.write(reinterpret_cast<const char*>(&(orec->pair)), sizeof(Tpair));
        if (!R) {
            spdlog::error("Error writing max occurrence pair to R file");
        }

        int left_elem = orec->pair.left;
        int right_elem = orec->pair.right;
        
        phraseBoundaries(left_elem, right_elem);
        sourceBoundaries(left_elem, right_elem);
        
        // Check if the current max pair is in hash ranges, if is then have to go through non-explicit phrases.
        if (hash_ranges.find(std::pair<int,int>{left_elem,right_elem}) != hash_ranges.end())
        {
            std::vector<int> ranges = hash_ranges[std::pair<int,int>(left_elem,right_elem)];
            bool firstRange = true;
            int prev_range;
            buildIntervalTree(); 
            for (int curr_range : ranges)
            {
                if (firstRange){
                    prev_range = curr_range;
                    firstRange = false;
                }
                else if (curr_range - prev_range == 1){
                    continue;
                }  
                // Replace in Ref
                RefNode* lref = rlist.findPos(curr_range);
                RefNode* rref = rlist.findForwardRef(lref);
                ITree::interval_vector phrase_results;
                spdlog::trace("Pair to replace: ({},{})", lref->pos, rref->pos);
                phrase_results = phrase_tree.findContained(lref->pos,rref->pos); // Fully contained within interval
                spdlog::trace("{} non-explicit phrases contain the pair to replace.", phrase_results.size());
                for (int i = 0; i < phrase_results.size(); i++)
                {
                    PhraseNode* nexp_phrase = phrase_results[i].value;
                    int leftElem =  lref->val;
                    int rightElem = rref->val;
                    int leftleftElem;
                    int rightrightElem;
                    // If range fully contained within range then we can do the decrease and increase frequencies
                    if (phrase_results[i].start != lref->pos && phrase_results[i].stop != rref->pos)
                    {
                        leftleftElem = rlist.findNearestRef(lref->prev)->val;
                        rightrightElem = rlist.findForwardRef(rref)->val;
                        decreaseFrequency(leftleftElem, leftElem);
                        increaseFrequency(leftleftElem, n);
                        decreaseFrequency(rightElem, rightrightElem);
                        increaseFrequency(n, rightrightElem);

                    }
                    // If range touches the edges
                    else if (phrase_results[i].start == lref->pos && phrase_results[i].stop == rref->pos)
                    {
                        if (nexp_phrase->prev != nullptr && nexp_phrase->prev->exp){
                            leftleftElem = nexp_phrase->prev->content.back();
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        } 
                        else if (nexp_phrase->prev != nullptr && !nexp_phrase->prev->exp) {
                            leftleftElem = rlist.findNearestRef(nexp_phrase->prev->rnode)->val;
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        }

                        if (nexp_phrase->next != nullptr && nexp_phrase->next->exp){
                            rightrightElem = nexp_phrase->next->content.front();
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        } 
                        else if (nexp_phrase->next != nullptr && !nexp_phrase->next->exp) {
                            rightrightElem = rlist.findNearestRef(nexp_phrase->next->lnode)->val;
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        }
                    }
                    // If range touches left edge, we can only decrease the left pair at the moment.
                    else if (phrase_results[i].start == lref->pos &&  phrase_results[i].stop != rref->pos)
                    {
                        if (nexp_phrase->prev != nullptr && nexp_phrase->prev->exp){
                            leftleftElem = nexp_phrase->prev->content.back();
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        } 
                        else if (nexp_phrase->prev != nullptr && !nexp_phrase->prev->exp) {
                            leftleftElem = rlist.findNearestRef(nexp_phrase->prev->rnode)->val;
                            decreaseFrequency(leftleftElem, leftElem);
                            increaseFrequency(leftleftElem, n);
                        }

                        rightrightElem = rlist.findForwardRef(rref)->val;
                        decreaseFrequency(rightElem, rightrightElem);
                        increaseFrequency(n, rightrightElem);
                    }
                    // If range touches right edge
                    else if (phrase_results[i].start != lref->pos &&  phrase_results[i].stop == rref->pos)
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
                            rightrightElem = rlist.findNearestRef(nexp_phrase->next->lnode)->val;
                            decreaseFrequency(rightElem, rightrightElem);
                            increaseFrequency(n, rightrightElem);
                        }
                    }
                    else{
                        spdlog::error("Should not logically happen!");
                    }
                }
                rlist.replacePair(n, lref, rref);
                prev_range = curr_range;
            }
        }
        
        // Have to process the explicit phrases and little bit of the non explicit phrases
        PhraseNode* curr_phrase = plist.getHead();
        while(curr_phrase != nullptr)
        {
            if (curr_phrase->exp)
            {
                auto curr_elem = curr_phrase->content.begin();
                while (curr_elem != curr_phrase->content.end())
                {
                    auto next_elem = std::next(curr_elem);
                    if (next_elem != curr_phrase->content.end() && (*curr_elem == left_elem && *next_elem == right_elem))
                    {
                        // Have to change the frequencies in max heap.
                        // If the two elements are not at the beginning or end of the phrase then replace happens fully within the phrase
                        if (curr_elem != curr_phrase->content.begin() && next_elem != std::prev(curr_phrase->content.end()))
                        {
                            // Decrease frequency of left pair effected by merge.
                            decreaseFrequency(*(std::prev(curr_elem)), *(curr_elem));
                            // Increase frequency of new pair.
                            increaseFrequency(*(std::prev(curr_elem)), n);
                            // Decrease frequency of right pair effected by merge.
                            decreaseFrequency(*(next_elem), *(std::next(next_elem)));
                            // Increase frequency of new pair.
                            increaseFrequency(n, *(std::next(next_elem)));
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
                            // Increase frequency of new pair.
                            increaseFrequency(n, *(std::next(next_elem)));
                        }
                        // If the left elem is in the middle of the phrase and the right elem is at the end of the phrase
                        else if (curr_elem != curr_phrase->content.begin() && next_elem == std::prev(curr_phrase->content.end()))
                        {
                            // Decrease frequency of left pair effected by merge.
                            decreaseFrequency(*(std::prev(curr_elem)), *(curr_elem));
                            // Increase frequency of new pair.
                            increaseFrequency(*(std::prev(curr_elem)), n);

                        
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
                        // Delete twice for both elements in the bi-gram
                        curr_elem = curr_phrase->content.erase(curr_elem);
                        curr_elem = curr_phrase->content.erase(curr_elem);
                        // Replace with n
                        curr_phrase->content.insert(curr_elem, n);
                    }
                    else{
                        curr_elem = std::next(curr_elem);                        
                    } 
                }
            }
            curr_phrase = curr_phrase->next;
        }

        // Think of better way. But going to clear and then repopulate due to laziness
        // We will go in reverse since prev pointers are updated. We store left pointer of pair so do not want the end.
        hash_ranges.clear();
        RefNode* end_elem = rlist.findNearestRef(rlist.getTail());
        RefNode* begin_elem = rlist.findNearestRef(rlist.getHead());
        std::pair<int,int> hash_pair;
        hash_pair.second = end_elem->val;

        while (end_elem != begin_elem)
        {
            end_elem = rlist.findNearestRef(end_elem->prev);
            hash_pair.first = end_elem->val;
            hash_ranges[hash_pair].emplace_back(end_elem->pos);
            hash_pair.second = hash_pair.first;
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
        spdlog::trace("");
        spdlog::trace("*** Information after bi-gram replacement ***");
        printRef();
        printPhraseList();
        printAllRecords();
        spdlog::trace("*********************************************");
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
    int verbosity = 0;
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
    prepareRef(rtext);

    // Clear char representation of Ref
    rtext.clear();
    rtext.shrink_to_fit();

    // Calculate the size of the original sequence file using the len field of the pairs
    psize = calculateParseBytes(pfile);

    // Call populatePhrases function
    populatePhrases(pfile);

    // Call createMaxHeap function
    createMaxHeap(pfile);

    // RePair
    repair(R, C);
    
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