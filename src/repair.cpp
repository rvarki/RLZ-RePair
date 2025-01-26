#include <CLI11.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <unordered_map>
#include <list>
#include "repair.h"
#include "phrase.h"

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
int alpha; // Number of characters used prior to RePair (alpha should actually start at 0, but currently set at 1 to avoid error in building max heap)
int n; // Technically |R| n - alpha - 1 gives number of rules

// Hash table containing the reference ranges of pairs (bi-grams). 
// The vector contains the left endpoint of the range corresponding to the pair in the bi-gram since the range is left endpoint to left endpoint + 1
std::unordered_map<std::string, std::vector<int>> hash_ranges; 

// List of explicit and non explicit phrases
std::list<Phrase> phrase_lst;

/**
 * @brief Calculates the number of bytes encoded in the RLZ parse.
 * 
 * @param[in] pfile [ifstream] the input file stream of the RLZ parse file.
 * @return the number of bytes encoded by the RLZ parse file.
 */
uint64_t calculate_parse_bytes(std::ifstream& pfile)
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
 * @brief Prints all records in the max heap (thanks BigRePair). Debug purposes only.
 * @return void
 */

void print_all_records()
{
    spdlog::debug("Current records in the heap");
    for (int i = 0; i < Rec.size; i++)
    {
        spdlog::debug("({},{}) {} occs", (unsigned int) Rec.records[i].pair.left, (unsigned int) Rec.records[i].pair.right, Rec.records[i].freq);
    }
}

/**
 * @brief Prints specific record in the max heap with a message. Debug purposes only.
 * @param[in] message [string] the message to be printed with the record.
 * @param[in] orec [Trecord*] the record content to be printed.
 * @return void
 */

void print_record(const std::string message, const Trecord* orec)
{
    spdlog::debug("{}",message);
    spdlog::debug("({},{}) {} occs", (unsigned int) orec->pair.left, (unsigned int) orec->pair.right, orec->freq);
}

/**
 * @brief Prints the reference during RePair. Debug purposes only.
 * @param[in] rvec [std::vector<unsigned char>&] reference vector
 * @return void
 */
void print_ref(const std::vector<unsigned char>& rvec)
{
    std::ostringstream oss;
    
    // Use range-based for loop for clarity
    for (unsigned char c : rvec) {
        oss << static_cast<int>(c) << " ";  // Use cast to int for numeric representation
    }

    spdlog::debug("Reference string (in numeric form): {}", oss.str());
}

/**
 * @brief Prints the hash table of ranges contents
 * @return void
 */
void print_hash_ranges()
{
    std::string values = "";
    for (const auto& [key, value] : hash_ranges){
        for (int i = 0; i < value.size(); i++){
            values += std::to_string(value[i]);
            values += " ";
        }
        spdlog::debug("Key: {}, Values: {}", key, values);
        values = "";
    }
}

/**
 * @brief Prints phrases in list
 * @return void
 */

void print_phrase_lst(const std::vector<unsigned char>& rvec)
{
    std::string content;
    for(Phrase phrase : phrase_lst){
        content = "";
        if (!(phrase.exp)){
            for(int i = phrase.lrange; i <= phrase.rrange; i++){
                content += std::to_string(rvec[i]);
            }
        }
        else{
            for (unsigned char i : phrase.content){
                content += std::to_string(i);
            }
        }
        spdlog::debug("Phrase: {}", content);
    }
    spdlog::debug("#################################");
}

/**
 * @brief Remaps the chars in the Ref vector and updates the chars and map arrays (logic from BigRePair).
 * 
 * Reasoning(?): UTF-8 chars are stored in 1 byte (8 bits, 0-255). However the Ref likely does not 
 * contain all chars which means wasted values between 0-255. We remap the chars in Ref to be 
 * between 0-|Ref| and write those values in the apporpriate cell in the chars array. We write the inverse 
 * of this operation in the map array. Doing this will allow us to know how many rules are created at the
 * end of RePair as well.
 * 
 * Additionally create hash table of the reference ranges bi-grams.
 * 
 * @param[in] rbuffer [std::vector<unsigned char>] The unsigned int representation of the reference.
 * @param[in] chars [int*] pointer to the head of the chars array.
 * @param[in] map [char*] pointer to the head of map array.
 * @param[in] size [int] size of the chars and map arrays.
 * 
 * @return void
 */

void prepareRef(std::vector<unsigned char>& rvec, int* chars, char* map, int size)
{
    alpha = 1;
    // Remaps the chars in the ref vector and upate chars array
    for (int i = 0; i < rvec.size(); i++){
        unsigned int x = rvec[i];
        if (chars[x] == -1){
            chars[x] = alpha++;
        }
        rvec[i] = chars[x];          
    }

    // Updates the map array to undo remapping.
    for (int i = 0; i < size; i++){
        if (chars[i] != -1){
            map[chars[i]] = i;
        }
    }
    print_ref(rvec);

    // Populate of hash table of ranges of reference bi-grams 
    for (int i = 1; i < rvec.size(); i++)
    {
        std::string pair_key = std::to_string(rvec[i-1]) + "-" + std::to_string(rvec[i]);
        int lrange = i - 1;
        auto it = hash_ranges.find(pair_key);
        if (it == hash_ranges.end())
            hash_ranges[pair_key].emplace_back(lrange);
        else{
            auto itv = std::find((it->second).begin(), (it->second).end(), lrange);
            if (itv == (it->second).end()){
                hash_ranges[pair_key].emplace_back(lrange);
            }
        }
    }
    print_hash_ranges();
}

/**
 * @brief Creates max heap of the bi-gram within the sequence file via the RLZ parse
 * 
 * We directly are using Heap/Hash/Record data structures in BigRePair to create the max heap
 * Max heap allows us to quickly find the bi-gram with the highest occurence.
 * We will also store the range that the bi-gram spans within a separate hash table.
 * 
 * @param[in] rvec [std::vector<unsigned char>&] the reference text 
 * @param[in] pfile [std::infstream&] the RLZ parse file stream
 * 
 * @return void?
 */
void createMaxHeap(const std::vector<unsigned char>& rvec, std::ifstream& pfile)
{
    // Create pair-freq heap (Thanks BigRePair)
    Rec = createRecords(factor, minsize);
    Heap = createHeap(psize, &Rec, factor, minsize);
    Hash = createHash(256 * 256, &Rec);
    assocRecords(&Rec, &Hash, &Heap, NULL);

    uint64_t num_pairs, pos, len;
    pair.left = '\0';
    pair.right = '\0';
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
                if (pair.left == '\0'){
                    pair.left = rvec[pos];
                }
                else{
                    pair.right = rvec[pos + j];
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
                }
            }
        }
    }
    // Print the heap in debug mode
    print_all_records();

    // Check the max freq pair from heap
    // id = extractMax(&Heap);
    // Trecord* orec = &Rec.records[id];
    // print_record("Maximum freq pair", orec);

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);
}

/**
 * @brief populates the phrases in RLZ parse. Range for each phrase is [left,right]
 * 
 * @param[in] pfile [std::ifstream&] the RLZ parse filestream
 * @return void 
 */

void populatePhrases(const std::vector<unsigned char>& rvec, std::ifstream& pfile)
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
            Phrase phrase(pos, pos + len - 1);
            phrase_lst.push_back(phrase);
        }
    }

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);

    // Debug
    print_phrase_lst(rvec);
}

/**
 * @brief Process the phrase list for the phrase boundary condition
 * 
 * @param[in] rvec
 * @return void
 */

void phraseBoundaries(std::vector<unsigned char>& rvec, int left_elem, int right_elem)
{
    auto curr_phrase = phrase_lst.begin();

    // Iterate through the phrases in the phrase list
    while (curr_phrase != phrase_lst.end()) 
    {
        auto next_phrase = std::next(curr_phrase);
        // If there is a next phrase, check the phrase boundaries.
        if (next_phrase != phrase_lst.end())
        {
            // Both phrases not explicit
            if (!(curr_phrase->exp) && !(next_phrase->exp))
            {
                if ((rvec[curr_phrase->rrange] == left_elem) && (rvec[next_phrase->lrange] == right_elem))
                {
                    std::list<unsigned char> content;
                    content.push_back(left_elem);
                    content.push_back(right_elem);
                    Phrase new_phrase(content);
                    phrase_lst.insert(next_phrase, new_phrase); // Insert before the next record
                    if (curr_phrase->rrange - 1 < curr_phrase->lrange){ // If the non-explicit phrase is empty we delete it.
                        curr_phrase = phrase_lst.erase(curr_phrase);
                    }
                    else{
                        curr_phrase->rrange = curr_phrase->rrange - 1;
                    }
                    if (next_phrase->lrange + 1 > next_phrase->rrange){
                        next_phrase = phrase_lst.erase(next_phrase);   
                    }
                    else{
                        next_phrase->lrange = next_phrase->lrange + 1;
                    }
                    continue;
                }
            }
            // Current phrase not explicit and next phrase explicit
            else if (!(curr_phrase->exp) && (next_phrase->exp))
            {
                if ((rvec[curr_phrase->rrange] == left_elem) && (next_phrase->content.front() == right_elem))
                {
                    next_phrase->content.push_front(left_elem);
                    if (curr_phrase->rrange - 1 < curr_phrase->lrange){ // If the non-explicit phrase is empty we delete it.
                        curr_phrase = phrase_lst.erase(curr_phrase);
                    }
                    else{
                        curr_phrase->rrange = curr_phrase->rrange - 1;
                    }
                    continue;
                }
            }
            // Current phrase explicit and next phrase not explicit
            else if ((curr_phrase->exp) && !(next_phrase->exp))
            {
                if ((curr_phrase->content.back() == left_elem) && (rvec[next_phrase->lrange] == right_elem))
                {
                    curr_phrase->content.push_back(right_elem);
                    if (next_phrase->lrange + 1 > next_phrase->rrange){
                        next_phrase = phrase_lst.erase(next_phrase);
                    }
                    else{
                        next_phrase->lrange = next_phrase->lrange + 1;
                    }
                    continue;
                }
            }
            // Both explicit phrases, move content of next phrase into first phrase and delete next phrase
            else{
                curr_phrase->content.splice(curr_phrase->content.end(), next_phrase->content);
                next_phrase = phrase_lst.erase(next_phrase);
            }

            // Update the phrase if it gets to this point
            curr_phrase = std::next(curr_phrase);
        }
        else{ // No more phrase boundaries left
            break;
        }
    }

    // Debug
    print_phrase_lst(rvec);
}

/**
 * @brief Process the phrase list for the source boundary condition
 */
void sourceBoundaries(std::vector<unsigned char>& rvec, int left_elem, int right_elem)
{
    auto curr_phrase = phrase_lst.begin();
    // Iterate through the phrases in the phrase list
    while (curr_phrase != phrase_lst.end()) 
    {
        // Check if the curr phrase is explicit or not
        if (!(curr_phrase->exp)){
            // Check if the curr phrase can form the bi-gram with the reference using its leftmost element
            if (curr_phrase->lrange != 0){
                // If the bi-gram can be formed then make the leftmost element make an explicit phrase of its own
                if (rvec[curr_phrase->lrange - 1] == left_elem && rvec[curr_phrase->lrange] == right_elem){
                    // If at the beginning of phrases, then no previous phrase
                    if (curr_phrase == phrase_lst.begin()){
                        std::list<unsigned char> content;
                        content.push_back(rvec[curr_phrase->lrange]);
                        Phrase new_phrase(content);
                        phrase_lst.insert(curr_phrase, new_phrase);
                    }
                    // Previous phrase exists.
                    else{
                        auto prev_phrase = std::prev(curr_phrase);
                        // If the previous phrase is explicit, we can add directly to it.
                        if (prev_phrase->exp){
                            prev_phrase->content.push_back(rvec[curr_phrase->lrange]);
                        }
                        // We have to create new explicit phrase anyways
                        else{
                            std::list<unsigned char> content;
                            content.push_back(rvec[curr_phrase->lrange]);
                            Phrase new_phrase(content);
                            phrase_lst.insert(curr_phrase, new_phrase);
                        }                        
                    }
                    // If non-explicit phrase is empty we delete it
                    if (curr_phrase->lrange + 1 > curr_phrase->rrange){ // If the non-explicit phrase is empty we delete it.
                        curr_phrase = phrase_lst.erase(curr_phrase);
                        continue;
                    }
                    else{
                        curr_phrase->lrange = curr_phrase->lrange + 1;
                    }
                }
            }
            if (curr_phrase->rrange != rvec.size()-1){
                // If the bi-gram can be formed then make the rightmost element make an explicit phrase of its own
                if (rvec[curr_phrase->rrange] == left_elem && rvec[curr_phrase->rrange + 1] == right_elem){
                    // If at the end of phrases, then no next phrase
                    if (curr_phrase == std::prev(phrase_lst.end())){
                        std::list<unsigned char> content;
                        content.push_back(rvec[curr_phrase->rrange]);
                        Phrase new_phrase(content);
                        auto next_phrase = std::next(curr_phrase);
                        phrase_lst.insert(next_phrase, new_phrase);
                    }
                    // Next phrase exists.
                    else{
                        auto next_phrase = std::next(curr_phrase);
                        // If the next phrase is explicit, we can add directly to it.
                        if (next_phrase->exp){
                            next_phrase->content.push_front(rvec[curr_phrase->rrange]);
                        }
                        // We have to create new explicit phrase anyways
                        else{
                            spdlog::debug("Gets to here!!!!!\n");
                            std::list<unsigned char> content;
                            content.push_back(rvec[curr_phrase->rrange]);
                            Phrase new_phrase(content);
                            phrase_lst.insert(next_phrase, new_phrase);
                        }                        
                    }
                    // If non-explicit phrase is empty we delete it
                    if (curr_phrase->rrange - 1 < curr_phrase->lrange){ // If the non-explicit phrase is empty we delete it.
                        curr_phrase = phrase_lst.erase(curr_phrase);
                        continue;
                    }
                    else{
                        curr_phrase->rrange = curr_phrase->rrange - 1;
                    }
                }
            }
        }
        curr_phrase = std::next(curr_phrase);
    }

    // Debug
    print_phrase_lst(rvec);
}

/**
 * @brief RePair
 * @return void
 */
void repair(std::vector<unsigned char>& rvec)
{
    int id = extractMax(&Heap);
    Trecord* orec = &Rec.records[id];
    print_record("Maximum freq pair", orec);
    phraseBoundaries(rvec, orec->pair.left, orec->pair.right);
    sourceBoundaries(rvec, orec->pair.left, orec->pair.right);
}

int main(int argc, char *argv[])
{
    CLI::App app("rlz - Run RePair with the RLZ parse.\n\nImplemented by Rahul Varki");
    std::string ref_file;
    std::string rlz_parse;
    bool verbose = false;
    std::string version = "Version: 1.0.0";

    app.add_option("-r,--ref", ref_file, "The reference file used to create the RLZ parse")->configurable()->required();
    app.add_option("-p,--parse", rlz_parse, "The RLZ parse of the sequence file (.rlz)")->configurable()->required();
    app.add_flag("--verbose", verbose, "Verbose output")->configurable();
    app.set_version_flag("-v,--version", version);
    app.footer("Example usage:\n"
           "  Compress: ./repair --ref reference.fasta --parse sequence.rlz\n");
    app.description("Run RePair on RLZ parse");
    CLI11_PARSE(app, argc, argv);
    if (verbose) { spdlog::set_level(spdlog::level::debug); }

    spdlog::info("Starting to RePair the RLZ parse");
    spdlog::info("The reference file provided: {}", ref_file);
    spdlog::info("The RLZ parse file provided: {}", rlz_parse);

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
    std::vector<unsigned char> rvec(rsize);
    rfile.read(reinterpret_cast<char*>(rvec.data()), rsize);

    // Opening RLZ parse
    std::ifstream pfile(rlz_parse, std::ios::binary | std::ios::in);
    if (!pfile) {
        spdlog::error("Error opening {}", rlz_parse);
        std::exit(EXIT_FAILURE);
    }

    // Initialize the chars array to be -1 to indicate that no characters have been mapped yet.
    for (int i = 0; i < 256; i++){
        chars[i] = -1;
    }

    // Process Ref and update the chars and map arrays.
    prepareRef(rvec, chars, map, 256);

    // Calculate the size of the original sequence file using the len field of the pairs
    psize = calculate_parse_bytes(pfile);
    spdlog::debug("The parse encodes for {} bytes", psize);

    // Call populatePhrases function
    populatePhrases(rvec, pfile);

    // Call createMaxHeap function
    createMaxHeap(rvec, pfile);

    // RePair
    repair(rvec);
    
    rfile.close();
    pfile.close();
    return 0;
}