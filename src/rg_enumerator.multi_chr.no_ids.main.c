#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib/glib.h>
#include "rg_enumerator_classes.c"
#include "rg_enumerator.multi_chr.no_ids.c"

int N_CHRS, IS_DIPLOID, MAX_DEPTH_DUP, MAX_DEPTH_NONDUP;
GHashTable *seen_somatic_genomes;

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Need input parameters n_chrs, diploid, max_dup_depth, max_overall_depth. Exiting.\n");
        fprintf(stderr, "Usage: /nfs/users/nfs_y/yl3/programs/rg_library_c/rg_enumerator.multi_chr <n_chrs> <diploid> <max_dup_depth> <max_overall_depth>\n");
        exit(1);
    }
    
    sscanf(argv[1], "%d", &N_CHRS);
    sscanf(argv[2], "%d", &IS_DIPLOID);
    sscanf(argv[3], "%d", &MAX_DEPTH_DUP);
    sscanf(argv[4], "%d", &MAX_DEPTH_NONDUP);
    fprintf(stderr, "Using %d chromosomes (%s)...\n", N_CHRS, (IS_DIPLOID == 0 ? "haploid" : "diploid"));
    fprintf(stderr, "Enumerating down to maximum of %d duplicative and %d overall rearrangements...\n", MAX_DEPTH_DUP, MAX_DEPTH_NONDUP);
    seen_somatic_genomes = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    struct genome *g_ptr = create_genome(N_CHRS, IS_DIPLOID);
    bridge(g_ptr);

    return(0);
}
