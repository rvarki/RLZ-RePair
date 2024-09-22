#include <CLI11.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include "repair.h"

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

// Calculates the number of bytes the RLZ parse file encodes for
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
    return seq_orig_size;
}

// Logic from BigRePair repair. Print the records in the heap.
void print_all_records()
{
    spdlog::debug("Current records in the heap");
    for (int i = 0; i < Rec.size; i++)
    {
        spdlog::debug("({},{}) {} occs", (char) Rec.records[i].pair.left, (char) Rec.records[i].pair.right, Rec.records[i].freq);
    }
}

void print_record(std::string message, Trecord* orec)
{
    spdlog::debug("{}",message);
    spdlog::debug("({},{}) {} occs", (char) orec->pair.left, (char) orec->pair.right, orec->freq);
}

// Runs repair on the RLZ parse
void repair(const std::string ref_str, std::ifstream& pfile)
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
                    pair.left = ref_str[pos];
                }
                else{
                    pair.right = ref_str[pos + j];
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
    id = extractMax(&Heap);
    Trecord* orec = &Rec.records[id];
    print_record("Maximum freq pair", orec);
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

    // Store the Ref file content in string variable
    std::ostringstream ref_stream;
    ref_stream << rfile.rdbuf();
    std::string ref_str = ref_stream.str();

    // Opening RLZ parse
    std::ifstream pfile(rlz_parse, std::ios::binary | std::ios::in);
    if (!pfile) {
        spdlog::error("Error opening {}", rlz_parse);
        std::exit(EXIT_FAILURE);
    }

    // Calculate the size of the original sequence file using the len field of the pairs
    psize = calculate_parse_bytes(pfile);
    spdlog::debug("The parse encodes for {} bytes", psize);

    // Reset the file pointer to the beginning of the file
    pfile.clear();
    pfile.seekg(0, std::ios::beg);

    // Call repair function
    repair(ref_str, pfile);
    
    rfile.close();
    pfile.close();
    return 0;
}