#ifndef BITARR_H
#define BITARR_H

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef size_t elem_t;
typedef unsigned char bit;

typedef struct {
    elem_t size;
    bit* vector;
} BitArray;

/**
 * Initializes the bit array module.
 * @return The size of a bit in the current environment.
 */
elem_t ba_init();

/**
 * Creates a new bit array with the specified number of elements.
 * @param nelems The number of elements in the bit array.
 * @return A pointer to the created BitArray, or NULL if creation failed.
 */
BitArray* ba_create(elem_t nelems);

/**
 * Destroys the given bit array, freeing its memory.
 * @param b The bit array to destroy.
 */
void ba_destroy(BitArray* b);

/**
 * Copies the contents of one bit array to another.
 * @param dst The destination bit array.
 * @param src The source bit array.
 * @param size The number of elements to copy.
 */
void ba_copy(BitArray* dst, const BitArray* src, elem_t size);

/**
 * Assigns a value to a specific element in the bit array.
 * @param b The bit array.
 * @param elem The element index.
 * @param value The value to assign (true for 1, false for 0).
 */
void ba_assign(BitArray* b, elem_t elem, bool value);

/**
 * Retrieves the value of a specific element in the bit array.
 * @param b The bit array.
 * @param elem The element index.
 * @return The value of the element (true for 1, false for 0).
 */
bool ba_value(const BitArray* b, elem_t elem);

/**
 * Toggles the value of a specific element in the bit array.
 * @param b The bit array.
 * @param elem The element index.
 */
void ba_toggle(BitArray* b, elem_t elem);

/**
 * Assigns a value to all elements in the bit array.
 * @param b The bit array.
 * @param size The number of elements.
 * @param value The value to assign (true for 1, false for 0).
 */
void ba_all_assign(BitArray* b, elem_t size, bool value);

/**
 * Converts an unsigned long number to a bit array.
 * @param num The number to convert.
 * @param b The bit array to store the result.
 * @param size The size of the bit array.
 * @return A pointer to the bit array, or NULL if conversion failed.
 */
BitArray* ba_ul2b(unsigned long num, BitArray* b, elem_t* size);

/**
 * Counts the number of set bits (1s) in the bit array.
 * @param b The bit array.
 * @param size The number of elements.
 * @return The number of set bits.
 */
unsigned long ba_count(const BitArray* b, elem_t size);

/**
 * Computes the intersection of two bit arrays.
 * @param first The first bit array.
 * @param second The second bit array.
 * @param result The resulting bit array.
 * @return True if the operation was successful, false otherwise.
 */
bool ba_intersection(const BitArray* first, const BitArray* second,
                     BitArray** result);

/**
 * Computes the union of two bit arrays.
 * @param first The first bit array.
 * @param second The second bit array.
 * @param result The resulting bit array.
 * @return True if the operation was successful, false otherwise.
 */
bool ba_union(const BitArray* first, const BitArray* second, BitArray** result);

/**
 * Computes the difference of two bit arrays.
 * @param first The first bit array.
 * @param second The second bit array.
 * @param result The resulting bit array.
 * @return True if the operation was successful, false otherwise.
 */
bool ba_diff(const BitArray* first, const BitArray* second, BitArray** result);

/**
 * Computes the complement of a bit array.
 * @param b The bit array.
 * @param size The number of elements.
 */
void ba_complement(BitArray* b, elem_t size);

/**
 * Computes the dot product of two bit arrays.
 * @param first The first bit array.
 * @param second The second bit array.
 * @param size_first The number of elements in the first bit array.
 * @param size_second The number of elements in the second bit array.
 * @return The dot product.
 */
unsigned long ba_dotprod(const BitArray* first, const BitArray* second,
                         elem_t size_first, elem_t size_second);

/**
 * Converts a bit array to a string representation.
 * @param b The bit array.
 * @param size The number of elements.
 * @param dest The destination string buffer.
 * @return The string representation of the bit array.
 */
char* ba_b2str(const BitArray* b, elem_t size, char* dest);

/**
 * Prints the bit array to the specified file stream.
 * @param b The bit array.
 * @param size The number of elements.
 * @param dest The file stream to print to.
 * @return True if the operation was successful, false otherwise.
 */
bool ba_print(const BitArray* b, elem_t size, FILE* dest);

#endif  // BITARR_H