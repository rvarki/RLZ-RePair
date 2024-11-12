#include "despair_c.h"

#define MAX_FILENAME 1024
#define MAX_RECURSION_DEPTH 1000

typedef struct {
    Tpair* rules;
    int n;
    int alph;
    FILE* output_file;
    const char* output_filename;
    int max_depth;
} DecompressionContext;

static int expand(DecompressionContext* ctx, int symbol, int depth) {
    if (depth > ctx->max_depth) {
        ctx->max_depth = depth;
    }

    if (depth > MAX_RECURSION_DEPTH) {
        fprintf(stderr, "Error: Maximum recursion depth exceeded\n");
        return -1;
    }

    if (symbol < ctx->alph) {
        char c = (char)symbol;
        if (fwrite(&c, sizeof(char), 1, ctx->output_file) != 1) {
            fprintf(stderr, "Error: Cannot write to file %s\n", ctx->output_filename);
            return -1;
        }
        return 1;
    }

    int rule_index = symbol - ctx->alph;
    if (rule_index >= ctx->n) {
        fprintf(stderr, "Error: Invalid rule index\n");
        return -1;
    }

    int left_count = expand(ctx, ctx->rules[rule_index].left, depth + 1);
    if (left_count == -1) return -1;

    int right_count = expand(ctx, ctx->rules[rule_index].right, depth + 1);
    if (right_count == -1) return -1;

    return left_count + right_count;
}

int despair(const char* basename) {
    char rules_filename[MAX_FILENAME];
    char compressed_filename[MAX_FILENAME];
    char output_filename[MAX_FILENAME];
    FILE *rules_file = NULL, *compressed_file = NULL, *output_file = NULL;
    DecompressionContext ctx = {0};
    int result = 0;
    long long total_symbols = 0;

    snprintf(rules_filename, sizeof(rules_filename), "%s.R", basename);
    snprintf(compressed_filename, sizeof(compressed_filename), "%s.S", basename);
    snprintf(output_filename, sizeof(output_filename), "%s.decompress", basename);

    // Open files
    rules_file = fopen(rules_filename, "rb");
    compressed_file = fopen(compressed_filename, "rb");
    output_file = fopen(output_filename, "w");  // Changed to "w" for text mode

    if (!rules_file || !compressed_file || !output_file) {
        fprintf(stderr, "Error: Cannot open input or output files\n");
        result = 1;
        goto cleanup;
    }

    // Read alphabet size
    if (fread(&ctx.alph, sizeof(int), 1, rules_file) != 1) {
        fprintf(stderr, "Error: Cannot read alphabet size\n");
        result = 1;
        goto cleanup;
    }

    // Read rules
    struct stat s;
    if (stat(rules_filename, &s) != 0) {
        fprintf(stderr, "Error: Cannot stat rules file\n");
        result = 1;
        goto cleanup;
    }
    ctx.n = (s.st_size - sizeof(int)) / sizeof(Tpair);
    ctx.rules = malloc(ctx.n * sizeof(Tpair));
    if (!ctx.rules) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        result = 1;
        goto cleanup;
    }
    if (fread(ctx.rules, sizeof(Tpair), ctx.n, rules_file) != ctx.n) {
        fprintf(stderr, "Error: Cannot read rules\n");
        result = 1;
        goto cleanup;
    }

    ctx.output_file = output_file;
    ctx.output_filename = output_filename;

    // Process compressed file
    int symbol;
    int first_symbol = 1;
    while (fread(&symbol, sizeof(int), 1, compressed_file) == 1) {
        if (first_symbol) {
            first_symbol = 0;
            continue;
        }
        int expanded_count = expand(&ctx, symbol, 0);
        if (expanded_count == -1) {
            result = 1;
            goto cleanup;
        }
        if (total_symbols > LLONG_MAX - expanded_count) {
            fprintf(stderr, "Error: Integer overflow in total symbol count\n");
            result = 1;
            goto cleanup;
        }
        total_symbols += expanded_count;
    }

    // Print statistics
    fprintf(stderr, "DesPair succeeded\n\n");
    fprintf(stderr, "   Original chars: %lld\n", total_symbols);
    fprintf(stderr, "   Number of rules: %d\n", ctx.n);
    fprintf(stderr, "   Maximum rule depth: %d\n", ctx.max_depth);

cleanup:
    if (rules_file) fclose(rules_file);
    if (compressed_file) fclose(compressed_file);
    if (output_file) fclose(output_file);
    free(ctx.rules);

    return result;
}