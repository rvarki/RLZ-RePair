#include "repair_c.h"

typedef struct {
    Trarray Rec;
    Theap Heap;
    Tlist *L;
    BitArray *b;
    Thash Hash;
    float factor;
    int alph;
    int minsize;
    int *seqFileContInt;
    int initial_len;
    int n;
    int MB;
    int *refFileContInt;
    int refLen;
    BitArray *explicit_symbols;
} Metadata;

void init_metadata(Metadata *meta, float factor) {
    meta->factor = factor;
    meta->minsize = 256;
}

/**
 * @brief Creates a BitArray from the contents of a file.
 *
 * This function reads a file containing the size of the BitArray followed by
 * a series of 64-bit integers. It initializes a BitArray of the specified size
 * and sets bits based on the read integers. The length of the BitArray is
 * provided as a parameter.
 *
 * @param fileName The name of the file to read from.
 * @param len The length of the BitArray to be created.
 * @return BitArray* Pointer to the created BitArray, or NULL on failure.
 */
BitArray *createBitArray(const char *fileName, int len) {
    FILE *file = fopen(fileName, "rb");
    if (file == NULL) {
        fprintf(stderr, "Error: cannot open file %s for reading\n", fileName);
        return NULL;
    }

    size_t n;
    if (fread(&n, sizeof(size_t), 1, file) != 1) {
        fprintf(stderr, "Error: failed to read size from file %s\n", fileName);
        fclose(file);
        return NULL;
    }

    BitArray *b = ba_create(len);
    if (b == NULL) {
        fprintf(stderr, "Error: failed to create BitArray\n");
        fclose(file);
        return NULL;
    }

    u_int64_t x;
    u_int64_t curr = 0;

    for (size_t i = 0; i < n; i++) {
        if (fread(&x, sizeof(u_int64_t), 1, file) != 1) {
            fprintf(stderr, "Error: failed to read data from file %s\n",
                    fileName);
            ba_destroy(b);
            fclose(file);
            return NULL;
        }

        if (fread(&x, sizeof(u_int64_t), 1, file) != 1) {
            fprintf(stderr, "Error: failed to read data from file %s\n",
                    fileName);
            ba_destroy(b);
            fclose(file);
            return NULL;
        }

        curr += x;
        ba_assign(b, curr - 1, true);
        if (i != n - 1) {
            ba_assign(b, curr, true);
        }
    }

    fclose(file);
    // ba_print(b, len, stdout);
    return b;
}

/**
 * @brief Populates an integer array with the contents of a file.
 *
 * This function reads the contents of a file into a character array,
 * then converts each character to an integer and stores it in the provided
 * integer array. It ensures that the entire file is read and converted.
 *
 * @param filePtr Pointer to the file to be read.
 * @param fileContInt Pointer to the integer array to be populated.
 * @param len Length of the file and the integer array.
 * @return int Returns 0 on success, 1 on failure.
 */
int populate_int_array(FILE *filePtr, int *fileContInt, relong len) {
    char *fileCont = (void *)malloc(len * sizeof(char));
    if (fread(fileCont, sizeof(char), len, filePtr) != len) {
        fprintf(stderr, "Error: failed to read file\n");
        return 1;
    }
    for (int i = 0; i < len; i++) {
        fileContInt[i] = (int)fileCont[i];
    }
    free(fileCont);
    return 0;
}

/**
 * @brief Populates various data structures within the Metadata struct.
 *
 * This function initializes and populates several data structures within the
 * provided Metadata struct based on the contents of an integer array and other
 * parameters. It calculates the alphabet size, creates records, heaps, and
 * hashes, and associates records with these structures. It then processes
 * pairs of integers from the array, updating the heap and records accordingly.
 *
 * @param len The length of the integer array to be processed.
 * @param meta Pointer to the Metadata struct to be populated.
 */
void populateStructures(int len, Metadata *meta) {
    int i, id;
    Tpair pair;
    meta->alph = 0;

    for (i = 0; i < meta->initial_len; i++) {
        if (meta->seqFileContInt[i] > meta->alph) {
            meta->alph = meta->seqFileContInt[i];
        }
    }
    meta->n = ++meta->alph;

    meta->Rec = createRecords(meta->factor, meta->minsize);
    meta->Heap =
        createHeap(meta->initial_len, &meta->Rec, meta->factor, meta->minsize);
    meta->Hash = createHash(256 * 256, &meta->Rec);
    assocRecords(&meta->Rec, &meta->Hash, &meta->Heap, meta->L);

    // check what this does
    // if((len/1024/1024) * 3 * sizeof(int) <= meta->MB) return;

    for (i = 0; i < len - 1; i++) {
        pair.left = meta->seqFileContInt[i];
        pair.right = meta->seqFileContInt[i + 1];
        id = searchHash(meta->Hash, pair);
        if (id == -1) {
            id = insertRecord(&meta->Rec, pair, meta->factor);
        } else {
            incFreq(&meta->Heap, id);
        }
    }
    purgeHeap(&meta->Heap);
}

/**
 * @brief Replaces a pair of integers in the reference sequence with a new
 * value.
 *
 * This function iterates through the reference sequence and replaces each
 * occurrence of the specified pair of integers (left and right) with a new
 * integer value (repVal). It updates the reference sequence in place and
 * adjusts its length accordingly.
 *
 * @param meta Pointer to the Metadata struct containing the reference sequence.
 * @param left The left integer of the pair to be replaced.
 * @param right The right integer of the pair to be replaced.
 * @param repVal The new integer value to replace the pair with.
 */
void replace_pair(Metadata *meta, int left, int right, int repVal) {
    // replace the pair with the repVal
    int i = 0, j = 0;
    while (j < meta->refLen - 1) {
        if (meta->refFileContInt[j] == left &&
            meta->refFileContInt[j + 1] == right) {
            meta->refFileContInt[i++] = repVal;
        }
        j++;
    }
    meta->refLen = i;
    meta->refFileContInt = (void *)realloc((void *)meta->refFileContInt,
                                           meta->refLen * sizeof(int));
}

/**
 * @brief Checks if a pair of integers in the sequence is replaceable.
 *
 * This function determines if a pair of integers (left and right) in the
 * sequence can be replaced based on the current position and a reference
 * sequence. It first checks if the current integer in the sequence matches the
 * left or right integer based on the checkRight flag. If there is a mismatch,
 * the pair is considered replaceable. It then iterates through the reference
 * sequence to check if the pair exists. If the pair is found in the reference
 * sequence, it is not replaceable.
 *
 * @param meta Pointer to the Metadata struct containing the sequences.
 * @param left The left integer of the pair to be checked.
 * @param right The right integer of the pair to be checked.
 * @param curr The current position in the sequence.
 * @param checkRight Flag indicating whether to check the right integer (true)
 * or the left integer (false).
 * @return bool Returns true if the pair is replaceable, false otherwise.
 */
bool is_replaceable(Metadata *meta, int left, int right, int curr,
                    bool checkRight) {
    if (checkRight && meta->seqFileContInt[curr] != left) {
        return true;
    } else if (!checkRight && meta->seqFileContInt[curr] != right) {
        return true;
    }
    for (int i = 0; i < meta->refLen - 1; i++) {
        if (checkRight) {
            if (meta->refFileContInt[i] == left &&
                meta->refFileContInt[i + 1] == right) {
                return false;
            }
        } else {
            if (meta->refFileContInt[i] == left &&
                meta->refFileContInt[i + 1] == right) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Checks if a given position in the sequence is a boundary.
 *
 * This function checks if the given position and the next position in the
 * BitArray are both set to true, indicating a boundary in the sequence.
 *
 * @param meta Pointer to the Metadata struct containing the BitArray.
 * @param pos The position in the sequence to check.
 * @return bool Returns true if both the given position and the next position
 *              in the BitArray are set to true, indicating a boundary.
 */
bool is_boundary(Metadata *meta, int pos) {
    return ba_value(meta->b, pos) && ba_value(meta->b, pos + 1);
}

/**
 * @brief Prints the contents of the sequence for testing purposes.
 *
 * This function iterates through the integer array in the Metadata struct
 * and prints each integer as a character followed by a space. The output
 * is terminated with a newline character. This function is intended for
 * testing purposes only.
 *
 * @param meta Pointer to the Metadata struct containing the sequence.
 * @param len The number of elements to print from the sequence.
 */
void check(Metadata *meta, int len) {
    for (int i = 0; i < len; i++) {
        printf("%c ", meta->seqFileContInt[i]);
    }
    printf("\n");
}

/**
 * @brief Prints the contents of a BitArray.
 *
 * This function iterates through the BitArray and prints each bit value
 * as an integer (0 or 1) followed by a space. The output is terminated
 * with a newline character.
 *
 * @param b Pointer to the BitArray to be printed.
 * @param len The number of bits to print from the BitArray.
 */
void print_bit_array(BitArray *b, int len) {
    for (int i = 0; i < len; i++) {
        printf("%d ", ba_value(b, i));
    }
    printf("\n");
}

/**
 * @brief Repairs a sequence by processing pairs and updating metadata.
 *
 * This function processes a sequence of integers, identifies pairs, and updates
 * the sequence and metadata accordingly. It writes the alphabet size and pairs
 * to a rule file, manages the heap and hash structures, and performs
 * substitutions in the sequence. The function iterates until no more pairs can
 * be processed.
 *
 * @param len The length of the sequence to be repaired.
 * @param RuleFileptr Pointer to the file where rules will be written.
 * @param meta Pointer to the Metadata struct containing sequence and auxiliary
 * data.
 * @return int The new length of the repaired sequence.
 */
int repair(int len, FILE *RuleFileptr, Metadata *meta) {
    int oid, id, cpos, pos;
    Trecord *orec;
    Tpair pair;
    int left, right;
    if (fwrite(&meta->alph, sizeof(int), 1, RuleFileptr) != 1) {
        fprintf(stderr, "Error: failed to write to file\n");
        return 1;
    }
    BitArray *new_b = ba_create(len);
    BitArray *new_explicit = ba_create(len);
    printf("printing boundary array: \n");
    print_bit_array(meta->b, len);
    printf("printing explicit array: \n");
    print_bit_array(meta->explicit_symbols, len);
    printf("printing sequence array: \n");
    check(meta, len);
    while (meta->n + 1 > 0) {
        // not sure what this is for
        // if ((len/1024/1024) * 3 * sizeof(int) <= meta->MB) return 0;
        oid = extractMax(&meta->Heap);
        if (oid == -1) {
            break;
        }
        orec = &meta->Rec.records[oid];
        if (fwrite(&orec->pair, sizeof(Tpair), 1, RuleFileptr) != 1) {
            fprintf(stderr, "Error: failed to write to file\n");
            return 1;
        }
        left = orec->pair.left;
        right = orec->pair.right;
        pos = 0;
        for (cpos = 0; cpos < len - 1; cpos++) {
            if ((cpos != 0) &&
                (is_boundary(meta, cpos - 1) &&
                 !is_replaceable(meta, left, right, cpos, false))) {
                // if the pair is in the beginning of the phrase
                ba_assign(new_explicit, pos, true);
                ba_assign(meta->explicit_symbols, cpos, true);
            } else if ((cpos != len - 1) && is_boundary(meta, cpos) &&
                       !is_replaceable(meta, left, right, cpos, true)) {
                // if the pair is in the end of the phrase
                ba_assign(new_explicit, pos, true);
                ba_assign(meta->explicit_symbols, cpos, true);
            }

            // if the pair is not found, copy the pair to the new sequence
            if ((meta->seqFileContInt[cpos] != left) ||
                (meta->seqFileContInt[cpos + 1] != right) ||
                ba_value(meta->explicit_symbols, cpos) ||
                ((cpos != len - 1) &&
                 ba_value(meta->explicit_symbols, cpos + 1))) {
                meta->seqFileContInt[pos] = meta->seqFileContInt[cpos];
                if (ba_value(meta->b, cpos)) {
                    ba_assign(new_b, pos, true);
                }
                if (ba_value(meta->explicit_symbols, cpos)) {
                    ba_assign(new_explicit, pos, true);
                }
            } else if (ba_value(meta->b, cpos) && ba_value(meta->b, cpos + 1)) {
                // if the pair is found on boundary
                meta->seqFileContInt[pos] = meta->seqFileContInt[cpos];
                if (ba_value(meta->b, cpos)) {
                    ba_assign(new_b, pos, true);
                }
                if (ba_value(meta->explicit_symbols, cpos)) {
                    ba_assign(new_explicit, pos, true);
                }
            } else {
                // decrement the frequency of the disappearing pair
                if (pos > 0) {
                    pair.left = meta->seqFileContInt[pos - 1];
                    pair.right = meta->seqFileContInt[cpos];
                    id = searchHash(meta->Hash, pair);
                    if (id != -1) {
                        if (id != oid) {
                            decFreq(&meta->Heap, id);
                        }
                    }
                }
                if (cpos < len - 2) {
                    pair.left = meta->seqFileContInt[cpos + 1];
                    pair.right = meta->seqFileContInt[cpos + 2];
                    id = searchHash(meta->Hash, pair);
                    if (id != -1) {
                        if (id != oid) decFreq(&meta->Heap, id);
                    }
                }

                // increment the frequency of the appearing pair
                if (pos > 0) {
                    pair.left = meta->seqFileContInt[pos - 1];
                    pair.right = meta->n;
                    id = searchHash(meta->Hash, pair);
                    if (id == -1) {
                        id = insertRecord(&meta->Rec, pair, meta->factor);
                    } else {
                        incFreq(&meta->Heap, id);
                    }
                }

                if (cpos < len - 2) {
                    pair.left = meta->n;
                    pair.right = meta->seqFileContInt[cpos + 2];
                    id = searchHash(meta->Hash, pair);
                    if (id == -1) {
                        id = insertRecord(&meta->Rec, pair, meta->factor);
                    } else {
                        incFreq(&meta->Heap, id);
                    }
                }

                // do the substitution
                meta->seqFileContInt[pos] = meta->n;
                if (ba_value(meta->b, cpos) || ba_value(meta->b, cpos + 1)) {
                    ba_assign(new_b, pos, true);
                }
                if (ba_value(meta->explicit_symbols, cpos)) {
                    ba_assign(new_explicit, pos, true);
                }
                cpos++;
            }
            pos++;
        }

        if (cpos == len - 1) {
            meta->seqFileContInt[pos++] = meta->seqFileContInt[cpos];
            if (ba_value(meta->b, cpos)) {
                ba_assign(new_b, pos, true);
            }
            if (ba_value(meta->explicit_symbols, cpos)) {
                ba_assign(new_explicit, pos, true);
            }
        }
        printf("printing boundary array: \n");
        print_bit_array(new_b, pos);
        printf("printing explicit array: \n");
        print_bit_array(new_explicit, pos);
        printf("printing sequence array: \n");
        check(meta, pos);
        len = pos;
        removeRecord(&meta->Rec, oid);
        meta->n++;
        purgeHeap(&meta->Heap);
        meta->seqFileContInt =
            (void *)realloc((void *)meta->seqFileContInt, len * sizeof(int));

        ba_destroy(meta->b);
        meta->b = new_b;
        new_b = ba_create(len);
        ba_destroy(meta->explicit_symbols);
        meta->explicit_symbols = new_explicit;
        new_explicit = ba_create(len);
    }
    purgeHeap(&meta->Heap);
    return len;
}

/**
 * @brief Main function for running repair on a sequence file.
 *
 * This function performs the main operations for repairing a sequence file.
 * It initializes metadata, reads the sequence file into an integer array,
 * creates a BitArray from another file, and processes the sequence to generate
 * rules and a repaired sequence. The results are written to output files.
 *
 * @param seqfile The name of the sequence file to be repaired.
 * @param rlzFile The name of the file containing the BitArray data.
 * @param refFile The name of the reference file used in the repair process.
 * @param MB The memory budget in megabytes.
 * @return int Returns 0 on success, 1 on failure.
 */
int repair_main(const char *seqfile, const char *rlzFile, const char *refFile, float factor,
                int MB) {
    int olen, len;
    Metadata *meta = (void *)malloc(sizeof(Metadata));
    init_metadata(meta, factor);
    meta->MB = MB;

    struct stat stFileInfo;
    if (stat(seqfile, &stFileInfo) != 0) {
        fprintf(stderr, "Error: cannot stat file %s\n", seqfile);
        return 1;
    }

    olen = len = stFileInfo.st_size;
    meta->initial_len = len;
    meta->explicit_symbols = ba_create(meta->initial_len);

    FILE *seqFileptr = fopen(seqfile, "r");
    if (seqFileptr == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", seqfile);
        return 1;
    }

    meta->seqFileContInt = (void *)malloc(len * sizeof(int));
    if (populate_int_array(seqFileptr, meta->seqFileContInt, len) != 0) {
        fprintf(stderr, "Error: failed to read file %s\n", seqfile);
        fclose(seqFileptr);
        return 1;
    }
    fclose(seqFileptr);

    FILE *refFileptr = fopen(refFile, "r");
    if (refFileptr == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", refFile);
        return 1;
    }

    if (stat(refFile, &stFileInfo) != 0) {
        fprintf(stderr, "Error: cannot stat file %s\n", refFile);
        return 1;
    }
    meta->refLen = stFileInfo.st_size;
    meta->refFileContInt = (void *)malloc(meta->refLen * sizeof(int));
    if (populate_int_array(refFileptr, meta->refFileContInt, meta->refLen) !=
        0) {
        fprintf(stderr, "Error: failed to read file %s\n", refFile);
        fclose(refFileptr);
        return 1;
    }

    meta->b = createBitArray(rlzFile, len);
    if (meta->b == NULL) {
        fprintf(stderr, "Error: failed to create BitArray\n");
        return 1;
    }

    char RuleFileName[1024];
    strcpy(RuleFileName, seqfile);
    strcat(RuleFileName, ".R");
    FILE *RuleFileptr = fopen(RuleFileName, "w");
    if (RuleFileptr == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", RuleFileName);
        return 1;
    }

    populateStructures(len, meta);

    len = repair(len, RuleFileptr, meta);
    char startSeqFileName[1024];
    strcpy(startSeqFileName, seqfile);
    strcat(startSeqFileName, ".S");
    FILE *startSeqFileptr = fopen(startSeqFileName, "wb");
    if (startSeqFileptr == NULL) {
        fprintf(stderr, "Error: cannot open file %s\n", startSeqFileName);
        return 1;
    }
    if (fwrite(&len, sizeof(int), 1, startSeqFileptr) != 1) {
        fprintf(stderr, "Error: failed to write to file\n");
        return 1;
    }
    if (fwrite(meta->seqFileContInt, sizeof(int), len, startSeqFileptr) !=
        len) {
        fprintf(stderr, "Error: failed to write to file\n");
        return 1;
    }

    check(meta, len);

    fclose(RuleFileptr);

    free(meta->seqFileContInt);
    ba_destroy(meta->b);
    free(meta);

    printf("Repair done\n");
    return 0;
}
