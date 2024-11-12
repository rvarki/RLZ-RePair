#include "bitarr.h"

#define BITS_SZ (CHAR_BIT)  // Number of bits in a char

/**
 * Calculates the number of elements needed to store n bits, given the number of
 * bits per element.
 * @param n The number of bits.
 * @param elem_per The number of bits per element.
 * @return The number of elements needed.
 */
static inline elem_t num_elements(elem_t n, elem_t elem_per) {
    return (n + elem_per - 1) / elem_per;
}

/**
 * Canonizes the bit array by masking out the unused bits in the last element.
 * @param array The bit array.
 * @param num_ints The number of elements in the array.
 * @param num_elem The total number of bits.
 */
static inline void canonize(bit* array, elem_t num_ints, elem_t num_elem) {
    array[num_ints - 1] &=
        (bit)~0 >>
        (BITS_SZ - ((num_elem % BITS_SZ) ? (num_elem % BITS_SZ) : BITS_SZ));
}

elem_t ba_init() { return BITS_SZ; }

BitArray* ba_create(elem_t nelems) {
    BitArray* b = (BitArray*)malloc(
        sizeof(BitArray));  // Allocate memory for the BitArray structure
    if (b) {
        b->size = nelems;
        size_t how_many = num_elements(
            nelems, BITS_SZ);  // Calculate the number of elements needed
        b->vector = (bit*)calloc(
            how_many,
            sizeof(bit));  // Allocate and zero-initialize the bit array
        if (!b->vector) {
            free(b);  // Free the BitArray structure if allocation failed
            return NULL;
        }
    }
    return b;
}

void ba_destroy(BitArray* b) {
    if (b) {
        free(b->vector);  // Free the bit array
        free(b);          // Free the BitArray structure
    }
}

void ba_copy(BitArray* dst, const BitArray* src, elem_t size) {
    elem_t nelem = num_elements(
        size, BITS_SZ);  // Calculate the number of elements to copy
    for (elem_t i = 0; i < nelem; ++i) {
        dst->vector[i] = src->vector[i];  // Copy each element
    }
}

void ba_assign(BitArray* b, elem_t elem, bool value) {
    if (value) {
        b->vector[elem / BITS_SZ] |=
            (1 << (elem % BITS_SZ));  // Set the bit to 1
    } else {
        b->vector[elem / BITS_SZ] &=
            ~(1 << (elem % BITS_SZ));  // Set the bit to 0
    }
}

bool ba_value(const BitArray* b, elem_t elem) {
    return (b->vector[elem / BITS_SZ] & (1 << (elem % BITS_SZ))) !=
           0;  // Check if the bit is set
}

void ba_toggle(BitArray* b, elem_t elem) {
    b->vector[elem / BITS_SZ] ^= (1 << (elem % BITS_SZ));  // Toggle the bit
}

void ba_all_assign(BitArray* b, elem_t size, bool value) {
    elem_t nelem =
        num_elements(size, BITS_SZ);  // Calculate the number of elements
    bit setval =
        (value) ? ~0
                : 0;  // Determine the value to set (all bits 1 or all bits 0)

    for (elem_t i = 0; i < nelem; ++i) {
        b->vector[i] = setval;  // Set each element
    }

    canonize(b->vector, nelem,
             size);  // Mask out the unused bits in the last element
}

unsigned long ba_count(const BitArray* b, elem_t size) {
    static const unsigned char bitcount[256] = {
        0, 1, 1, 2, 1, 2, 2, 3,  // ... (initialize all 256 values)
        // Complete the bitcount array
    };
    unsigned long count = 0;
    elem_t nelem =
        num_elements(size, BITS_SZ);  // Calculate the number of elements

    for (elem_t i = 0; i < nelem; ++i) {
        count += bitcount[b->vector[i]];  // Count the bits in each element
    }

    return count;
}

char* ba_b2str(const BitArray* b, elem_t size, char* dest) {
    if (dest == NULL) {
        dest = (char*)malloc(size + 1);  // Allocate memory for the string
        if (dest == NULL) {
            return NULL;  // Return NULL if allocation failed
        }
    }

    for (elem_t i = 0; i < size; ++i) {
        dest[i] = ba_value(b, i) ? '1' : '0';  // Convert each bit to '1' or '0'
    }

    dest[size] = '\0';  // Null-terminate the string
    return dest;
}

bool ba_print(const BitArray* b, elem_t size, FILE* dest) {
    char* to_print =
        ba_b2str(b, size, NULL);  // Convert the bit array to a string

    if (to_print != NULL) {
        bool status =
            (EOF !=
             fputs(to_print, dest));  // Print the string to the file stream
        free(to_print);               // Free the allocated string
        return status;
    }

    return false;
}