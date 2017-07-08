/* defines the basic classes and functions needed for enumerating
   rearrangements. */

/*
    A chromosome is coded as a list of segments.

    Segments have following attributes:
    * parent, which was broken down to generate the current segment
    * child_index, index differentiating segment from other children of the parent
    * is_plus, indicating whether the segment is currently sitting in reference or antisense orientation
    * what is the next segment on the chromosome?

    A genome is coded as a list of root segments.
*/
struct seg {
    int *seg_indexes;  /* Array of indexes for unique determination of unique segments */
    int times_divided;  /* How many times _genome_ has been divided to generate this segment. Intact chromosome = minimum = 1. */
    int is_plus;
    int is_maternal;
};
struct chromosome {
    struct seg **root_seg;  /* A dynamic array of segment structs */
    int n_segs;
};
enum rg_type {
    DEL,
    TD,
    INV_DUP,
    INV,
    TEL_BREAK,
    FOLD_BACK,
    BAL_TRANSLOC,
    UNBAL_TRANSLOC,
    WC_DUP,
    WC_DEL,
    WG_DUP
};
char* rg_type_to_txt(enum rg_type rg) {
    switch(rg) {
        case DEL : return("del");
        case TD : return("td");
        case INV_DUP : return("id");
        case INV : return("inv");
        case TEL_BREAK : return("tb");
        case FOLD_BACK : return("fb");
        case BAL_TRANSLOC : return("bt");
        case UNBAL_TRANSLOC : return("ut");
        case WC_DUP : return("wcg");
        case WC_DEL : return("wcl");
        case WG_DUP : return("wgd");
    }
}
struct genome {
    struct chromosome **root_chr;  /* A dynamic array of chromosome structs */
    int n_chrs;
    struct seg **genome_segs;  /* This refers to the list of segments present in this genome */
    int n_genome_segs;         /* This refers to the number of different types of segments in this genome */
    enum rg_type *history;  /* What events led to this genome? */
    int *history_idx;  /* Which of the possible applications of the events ever applied to get this genome? */
    int depth;  /* how many events have happened in this genome so far? */
    int dup_depth;  /* Number of duplicative events happened in this genome so far? */
    int wgd_depth;
};

extern int MAX_DEPTH_DUP, MAX_DEPTH_NONDUP;

#define MAX_TIMES_DIVIDED 256

/*
    Function prototypes
*/
struct seg* create_seg(int name, int is_maternal);
struct seg* copy_seg(struct seg* s_ptr);
void delete_seg(struct seg* s_ptr);
struct chromosome* create_chromosome(int name, int is_maternal);
struct chromosome* copy_chromosome(struct chromosome* c_ptr);
void delete_chromosome(struct chromosome* c_ptr);
struct genome* create_genome(int n_chrs, int paired);
struct genome* copy_genome(struct genome* g_ptr);
void delete_genome(struct genome* g_ptr);
void lose_chromosome_in_genome(struct genome* g_ptr, int c_idx);

void splice_one_seg(struct chromosome *c_ptr, int seg_idx, int split_into);
void splice_all_segs(struct genome *g_ptr, int *seg_indexes, int times_divided, int split_into);
void delete_segs_from_chr(struct chromosome *c_ptr, int from, int to);
void _validate_genome(struct genome *g_ptr, char *source);
void _validate_chromosome(struct chromosome *c_ptr, char *source);
void _validate_seg(struct seg *s_ptr, char *source);

void int_arr_to_string(int* int_arr, int arr_size, char *dest);
void print_genome(struct genome* g_ptr, char *unique_genome_string);

void get_unique_genome_string(struct genome *g_ptr, char *out_string_ptr);
/*
    End function prototypes
*/



/*
    Functions for creating, copying and deleting genomes
*/
struct genome* create_genome(int n_chrs, int paired) {
    struct genome *g_ptr = malloc(sizeof(struct genome));
    if (g_ptr == NULL) {
        fprintf(stderr, "\nCreation of genome node failed. Exiting.\n");
        exit(1);
    }
    g_ptr->n_chrs = (paired ? 2 * n_chrs : n_chrs);
    if (n_chrs == 0) {
        g_ptr->root_chr = NULL;
        return(g_ptr);
    }

    g_ptr->root_chr = malloc((paired ? 2*n_chrs : n_chrs)*sizeof(struct chromosome*));
    if (g_ptr->root_chr == NULL) {
        fprintf(stderr, "\nCreation of root_chr node failed. Exiting.\n");
        exit(1);
    }

    g_ptr->n_genome_segs = n_chrs;  // Because each chromosome gets its own segment
    g_ptr->genome_segs = malloc(n_chrs * sizeof(struct seg*));
    if (g_ptr->genome_segs == NULL) {
        fprintf(stderr, "\nCreation of g_ptr->genome_segs failed. Exiting.\n");
        exit(1);
    }

    g_ptr->history = NULL;
    g_ptr->history_idx = NULL;
    g_ptr->depth = 0;
    g_ptr->dup_depth = 0;
    g_ptr->wgd_depth = 0;

    /* Create the chromosomes */
    int i;
    for (i = 0; i < n_chrs; i++) {
        if (paired) {
            *(g_ptr->root_chr+2*i)   = create_chromosome(i, 0);
            *(g_ptr->root_chr+2*i+1) = create_chromosome(i, 1);
        }
        else {
            *(g_ptr->root_chr+i) = create_chromosome(i, 0);
        }
        *(g_ptr->genome_segs+i) = create_seg(i, 0);
    }

    return(g_ptr);
}

struct genome* copy_genome(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "copy_genome()");
    int i;

    struct genome *new_g_ptr = malloc(sizeof(struct genome));
    if (new_g_ptr == NULL) {
        fprintf(stderr, "\nCreation of genome node failed. Exiting.\n");
        exit(1);
    }
    new_g_ptr->n_chrs = g_ptr->n_chrs;

    // Copy rearrangement history over
    if (g_ptr->depth > 0) {
        int d = g_ptr->depth;
        new_g_ptr->history = malloc(d * sizeof(enum rg_type));
        if (new_g_ptr->history == NULL) {
            fprintf(stderr, "\nCreation of new_g_ptr->history failed. Exiting.\n");
            exit(1);
        }
        memcpy(new_g_ptr->history, g_ptr->history, d * sizeof(enum rg_type));

        new_g_ptr->history_idx = malloc(d * sizeof(int));
        if (new_g_ptr->history_idx == NULL) {
            fprintf(stderr, "\nCreation of g_ptr->history_idx failed. Exiting.\n");
            exit(1);
        }
        memcpy(new_g_ptr->history_idx, g_ptr->history_idx, d * sizeof(int));
    }
    else {
        new_g_ptr->history = NULL;
        new_g_ptr->history_idx = NULL;
    }

    new_g_ptr->depth = g_ptr->depth;
    new_g_ptr->dup_depth = g_ptr->dup_depth;
    new_g_ptr->wgd_depth = g_ptr->wgd_depth;

    /* Copy the genome_segs */
    new_g_ptr->n_genome_segs = g_ptr->n_genome_segs;
    new_g_ptr->genome_segs = malloc(new_g_ptr->n_genome_segs * sizeof(struct seg*));
    if (new_g_ptr->genome_segs == NULL) {
        fprintf(stderr, "\nCreation of new_g_ptr->genome_segs failed. Exiting.\n");
        exit(1);
    }
    for (i=0; i<new_g_ptr->n_genome_segs; i++) {
        *(new_g_ptr->genome_segs+i) = copy_seg(*(g_ptr->genome_segs+i));
    }

    new_g_ptr->root_chr = malloc(g_ptr->n_chrs*sizeof(struct chromosome*));
    if (new_g_ptr->root_chr == NULL) {
        fprintf(stderr, "\nCreation of root_chr node failed. Exiting.\n");
        exit(1);
    }

    /* Copy the chromosomes */
    if (g_ptr->n_chrs == 0) {
        new_g_ptr->root_chr = NULL;
        return(new_g_ptr);
    }
    else {
        for (i = 0; i < g_ptr->n_chrs; i++) {
            *(new_g_ptr->root_chr+i) = copy_chromosome(*(g_ptr->root_chr+i));
        }
    }

    return(new_g_ptr);
}

void delete_genome(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "delete_genome()");
    int i;
    for (i=0; i<g_ptr->n_chrs; i++) {
        delete_chromosome(*(g_ptr->root_chr+i));
    }
    free(g_ptr->root_chr);
    for (i=0; i<g_ptr->n_genome_segs; i++) {
        delete_seg(*(g_ptr->genome_segs+i));
    }
    free(g_ptr->genome_segs);
    free(g_ptr->history);
    free(g_ptr->history_idx);
    free(g_ptr);
    return;
}

void lose_chromosome_in_genome(struct genome* g_ptr, int c_idx) {
    // _validate_genome(g_ptr, "lose_chromosome_in_genome()");
    
    // Delete the chromosome and shift the pointers in g_ptr->root_chr
    delete_chromosome(*(g_ptr->root_chr+c_idx));
    while (c_idx <= g_ptr->n_chrs - 2) {
        *(g_ptr->root_chr+c_idx) = *(g_ptr->root_chr+c_idx+1);
        c_idx++;
    }
    g_ptr->n_chrs -= 1;
    g_ptr->root_chr = realloc(g_ptr->root_chr, (g_ptr->n_chrs) * sizeof(struct chromosome*));
    if (g_ptr->root_chr == NULL) {
        fprintf(stderr, "\nRealloc of g_ptr->root_chr failed. Exiting.\n");
        exit(1);
    }

    return;
}

void make_history(struct genome *g_ptr, enum rg_type rg, int idx) {  /* Make history */
    // _validate_genome(g_ptr, "make_history()");
    g_ptr->depth += 1;
    g_ptr->history = realloc(g_ptr->history, g_ptr->depth * sizeof(enum rg_type));
    if (g_ptr->history == NULL) {
            fprintf(stderr, "\nrealloc of g_ptr->history failed. Exiting.\n");
            exit(1);
    }
    *(g_ptr->history+g_ptr->depth-1) = rg;
    g_ptr->history_idx = realloc(g_ptr->history_idx, g_ptr->depth * sizeof(int));
    if (g_ptr->history_idx == NULL) {
            fprintf(stderr, "\nrealloc of g_ptr->history_idx failed. Exiting.\n");
            exit(1);
    }
    *(g_ptr->history_idx+g_ptr->depth-1) = idx;

    switch(rg) {
        case TD        : g_ptr->dup_depth++; break;
        case INV_DUP   : g_ptr->dup_depth++; break;
        case FOLD_BACK : g_ptr->dup_depth++; break;
        case WC_DUP    : g_ptr->dup_depth++; break;
        case WG_DUP    : g_ptr->wgd_depth++; g_ptr->dup_depth++; break;
    }

    return;
}
/*
    End genome functions
*/


/*
    Functions for creating, copying and deleting chromosomes
*/
struct chromosome* create_chromosome(int name, int is_maternal) {
    struct chromosome *c_ptr = malloc(sizeof(struct chromosome));
    if (c_ptr == NULL) {
        fprintf(stderr, "\nCreation of chromosome node failed. Exiting.\n");
        exit(1);
    }
    
    c_ptr->n_segs = 1;
    c_ptr->root_seg = malloc(sizeof(struct seg*));
    if (c_ptr->root_seg == NULL) {
        fprintf(stderr, "\nCreation of root seg node failed. Exiting.\n");
        exit(1);
    }
    *(c_ptr->root_seg+0) = create_seg(name, is_maternal);

    return(c_ptr);
}

struct chromosome* copy_chromosome(struct chromosome* c_ptr) {
    if (c_ptr == NULL) {
        fprintf(stderr, "\nNULL pointer passed to copy_chromosome(). Exiting.\n");
        exit(1);
    }

    if (c_ptr->root_seg == NULL || c_ptr->n_segs == 0) {
        fprintf(stderr, "\nInput chromosome does not have segments at copy_chromosome(). Exiting.\n");
        exit(1);
    }

    /* Create memory for the new chromosome */
    struct chromosome *new_c_ptr = malloc(sizeof(struct chromosome));
    new_c_ptr->n_segs = c_ptr->n_segs;
    new_c_ptr->root_seg = malloc(c_ptr->n_segs*sizeof(struct seg*));
    if (new_c_ptr->root_seg == NULL) {
        fprintf(stderr, "\nCreation of root seg node failed. Exiting.\n");
        exit(1);
    }

    /* Copy the segments */
    int i;
    for (i = 0; i < c_ptr->n_segs; i++) {
        *(new_c_ptr->root_seg + i) = copy_seg(*(c_ptr->root_seg + i));
    }

    return(new_c_ptr);
}

void delete_chromosome(struct chromosome* c_ptr) {
    char diagnostic_info[256] = "delete_chromosome()";
    // _validate_chromosome(c_ptr, diagnostic_info);
    int i;
    for (i=0; i<c_ptr->n_segs; i++) {
        delete_seg(*(c_ptr->root_seg+i));
    }
    free(c_ptr->root_seg);
    free(c_ptr);
    return;
}
/*
    End chromosome functions
*/


/*
    Functions for creating, copying and deleting segments
*/
struct seg* create_seg(int name, int is_maternal) {
    struct seg *s_ptr = malloc(sizeof(struct seg));
    if (s_ptr == NULL) {
        fprintf(stderr, "\nCreation of seg node failed. Exiting.\n");
        exit(1);
    }
    s_ptr->times_divided = 1;  // A chromosome is divided "once" since it has one seg_index. 
    s_ptr->seg_indexes = malloc(sizeof(int));
    if (s_ptr->seg_indexes == NULL) {
        fprintf(stderr, "\nFailed to malloc s_ptr->seg_indexes. Exiting\n");
        exit(1);
    }
    *(s_ptr->seg_indexes+0) = name;
    s_ptr->is_plus = 1;
    s_ptr->is_maternal = is_maternal;

    return(s_ptr);
}

struct seg* copy_seg(struct seg* s_ptr) {
    if (s_ptr == NULL) {
        fprintf(stderr, "\nNULL pointer passed to s_ptr. Exiting.\n");
        exit(1);
    }

    struct seg *new_s_ptr = malloc(sizeof(struct seg));
    if (new_s_ptr == NULL) {
        fprintf(stderr, "\nCreation of seg node failed. Exiting.\n");
        exit(1);
    }

    new_s_ptr->times_divided = s_ptr->times_divided;
    new_s_ptr->seg_indexes = malloc(s_ptr->times_divided * sizeof(int));
    if (new_s_ptr->seg_indexes == NULL) {
        fprintf(stderr, "\nCreation of seg indexes nodes failed. Exiting.\n");
        exit(1);
    }
    memcpy(new_s_ptr->seg_indexes, s_ptr->seg_indexes, new_s_ptr->times_divided * sizeof(int));

    new_s_ptr->is_plus = s_ptr->is_plus;
    new_s_ptr->is_maternal = s_ptr->is_maternal;

    return(new_s_ptr);
}

void delete_seg (struct seg* s_ptr) {
    char diagnostic_info[256] = "delete_seg()";
    // _validate_seg(s_ptr, diagnostic_info);
    free(s_ptr->seg_indexes);
    free(s_ptr);
    return;
}
/*
    End segment functions
*/

int int_array_cmp(int len1, int *val1, int len2, int *val2) {
    if (len1 != len2) { return(0); }
    int i;
    for (i=0; i<len1; i++) {
        if (*(val1+i) != *(val2+i)) { return(0); }
    }
    return(1);
}

/* Below function splices only one segment */
void splice_one_seg(struct chromosome *c_ptr, int seg_idx, int split_into) {
    char diagnostic_info[256] = "splice_one_seg()";
    // _validate_chromosome(c_ptr, diagnostic_info);

    c_ptr->n_segs += split_into - 1;
    c_ptr->root_seg = realloc(c_ptr->root_seg, c_ptr->n_segs * sizeof(struct seg*));
    if (c_ptr->root_seg == NULL) {
        fprintf(stderr, "NULL segment pointer at c_ptr->root_seg in function splice_one_seg(). Exiting.\n");
        exit(1);
    }

    int i;
    for (i=c_ptr->n_segs-1; i>=seg_idx+split_into; i--) {  /* Shift all segments appearing after the segment to be spliced */
        *(c_ptr->root_seg+i) = *(c_ptr->root_seg+i-split_into+1);
    }

    struct seg *tmp_s_ptr;
    for (i=split_into-1; i>=0; i--) {  /* Splice the segment into split_into new segments */
        if (i > 0) {
            *(c_ptr->root_seg+seg_idx+i) = copy_seg(*(c_ptr->root_seg+seg_idx));
        }
        tmp_s_ptr = *(c_ptr->root_seg+seg_idx+i);
        tmp_s_ptr->times_divided += 1;
        tmp_s_ptr->seg_indexes = realloc(tmp_s_ptr->seg_indexes, tmp_s_ptr->times_divided * sizeof(int));
        if (tmp_s_ptr->seg_indexes == NULL) {
            fprintf(stderr, "Failed to realloc seg_indexes in function splice_one_seg(). Exiting.\n");
            exit(1);
        }
        *(tmp_s_ptr->seg_indexes + tmp_s_ptr->times_divided - 1) = tmp_s_ptr->is_plus ? i : split_into - 1 - i;
    }

    return;
}

/* Below function splices all segments in *g_ptr that have segment indexes same as *seg_indexes */
/* Also splices g_ptr->genome_segs */
void splice_all_segs(struct genome *g_ptr, int *seg_indexes, int times_divided, int split_into) {
    /* In the input genome, look for segments that have an index given by *seg_indexes
       (length: times_divided), and split these segments into split_into child segments. */

    char diagnostic_info[256] = "splice_all_segs()";
    // _validate_genome(g_ptr, diagnostic_info);

    struct chromosome *cur_c_ptr;
    struct seg *cur_s_ptr, *tmp_s_ptr;
    int c_idx, s_idx;
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {  /* Loop through all chromosomes */
        cur_c_ptr = *(g_ptr->root_chr+c_idx);
        sprintf(diagnostic_info, "splice_all_segs(), c_idx %d", c_idx);
        // _validate_chromosome(cur_c_ptr, diagnostic_info);

        for (s_idx=0; s_idx<cur_c_ptr->n_segs; s_idx++) {  /* Loop through all segments in current chromosome */
            sprintf(diagnostic_info, "splice_all_segs(), c_idx %d, s_idx %d", c_idx, s_idx);
            cur_s_ptr = *(cur_c_ptr->root_seg+s_idx);
            // _validate_seg(cur_s_ptr, diagnostic_info);

            /* This segment has to be spliced? */
            if (int_array_cmp(times_divided, seg_indexes, cur_s_ptr->times_divided, cur_s_ptr->seg_indexes)) {
                splice_one_seg(cur_c_ptr, s_idx, split_into);
            }
        }
    }


    // Splice genome_segs
    g_ptr->n_genome_segs += split_into - 1;
    g_ptr->genome_segs = realloc(g_ptr->genome_segs, g_ptr->n_genome_segs * sizeof(struct seg*));
    if (g_ptr->genome_segs == NULL) {
        fprintf(stderr, "\nFailed to realloc g_ptr->genome_segs in function splice_all_segs(). Exiting.\n");
        exit(1);
    }
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        cur_s_ptr = *(g_ptr->genome_segs+s_idx);
        sprintf(diagnostic_info, "splice_all_segs(), genome_segs s_idx %d", s_idx);
        // _validate_seg(cur_s_ptr, diagnostic_info);

        /* This segment has to be spliced? */
        if (int_array_cmp(times_divided, seg_indexes, cur_s_ptr->times_divided, cur_s_ptr->seg_indexes)) {
            int i;
            for (i=g_ptr->n_genome_segs-1; i>=s_idx+split_into; i--) {  /* Shift all segments appearing after the segment to be spliced */
                *(g_ptr->genome_segs+i) = *(g_ptr->genome_segs+i-split_into+1);
            }   

            for (i=split_into-1; i>=0; i--) {  /* Splice the segment into split_into new segments */
                if (i > 0) {
                    *(g_ptr->genome_segs+s_idx+i) = copy_seg(*(g_ptr->genome_segs+s_idx));
                }   
                tmp_s_ptr = *(g_ptr->genome_segs+s_idx+i);
                tmp_s_ptr->times_divided += 1;
                tmp_s_ptr->seg_indexes = realloc(tmp_s_ptr->seg_indexes, tmp_s_ptr->times_divided * sizeof(int));
                if (tmp_s_ptr->seg_indexes == NULL) {
                    fprintf(stderr, "Failed to realloc seg_indexes in function splice_all_segs(). Exiting.\n");
                    exit(1);
                }   
                *(tmp_s_ptr->seg_indexes + tmp_s_ptr->times_divided - 1) = i;
            }   

            break;
        }   
    }

    return;
}


void delete_segs_from_chr(struct chromosome *c_ptr, int from, int to) {
    /* Splices out segments from from to to */
    // _validate_chromosome(c_ptr, "delete_segs_from_chr()");

    if (to < from) {
        /* Nothing to delete */
        return;
    }

    struct seg **segs_ptr = c_ptr->root_seg;
    int del_len = to - from + 1;
    int i;
    for (i=from; i<=to; i++) {
        delete_seg(*(segs_ptr+i));
    }
    for (i=from; i+del_len < c_ptr->n_segs; i++) {
        *(segs_ptr+i) = *(segs_ptr+i+del_len);
    }
    c_ptr->n_segs -= del_len;
    c_ptr->root_seg = realloc(c_ptr->root_seg, c_ptr->n_segs * sizeof(struct seg*));
    if (c_ptr->root_seg == NULL) {
        fprintf(stderr, "Failed to realloc c_ptr->root_seg in delete_segs(). Exiting.\n");
        exit(1);
    }

    return;
}

struct chromosome* yank_segments(struct chromosome *c_ptr, int from, int to) {
    // _validate_chromosome(c_ptr, "yank_segments()");
    struct chromosome *new_c_ptr = copy_chromosome(c_ptr);
    delete_segs_from_chr(new_c_ptr, to+1, new_c_ptr->n_segs - 1);
    delete_segs_from_chr(new_c_ptr, 0, from-1);
    new_c_ptr->n_segs = to - from + 1;
    return new_c_ptr;
}

void insert_segs_into_chr(struct chromosome *c_ptr, struct chromosome *segs_to_insert, int insert_before) {
    // _validate_chromosome(c_ptr, "insert_segs_int_chr(), c_ptr");
    // _validate_chromosome(segs_to_insert, "insert_segs_int_chr(), segs_to_insert");
    c_ptr->n_segs += segs_to_insert->n_segs;
    c_ptr->root_seg = realloc(c_ptr->root_seg, c_ptr->n_segs * sizeof(struct seg*));
    if (c_ptr->root_seg == NULL) {
        fprintf(stderr, "\nFailed to realloc c_ptr->root_seg in function insert_segs_into_chr. Exiting.\n");
        exit(1);
    }

    // Shift all the segments towards the end
    int i;
    for (i=c_ptr->n_segs-1; i-segs_to_insert->n_segs >= insert_before; i--) {
        *(c_ptr->root_seg+i) = *(c_ptr->root_seg+i-(segs_to_insert->n_segs));
    }

    // Add in the new segments
    for (i=0; i<segs_to_insert->n_segs; i++) {
        *(c_ptr->root_seg+insert_before+i) = copy_seg(*(segs_to_insert->root_seg+i));
    }

    return;
}

void invert_segs_in_chr(struct chromosome *c_ptr, int from, int to) {
    // _validate_chromosome(c_ptr, "invert_segs_in_chr()");
    struct seg **segs = malloc((to - from + 1) * sizeof(struct seg*));
    if (segs == NULL) {
        fprintf(stderr, "\nFailed to malloc segs in invert_segs_in_chr(). Exiting.\n");
        exit(1);
    }
    memcpy(segs, c_ptr->root_seg+from, (to - from + 1) * sizeof(struct seg*));

    int i;
    for (i=0; i<to - from + 1; i++) {
        *(c_ptr->root_seg+from+i) = *(segs + (to-from) - i);
        (*(c_ptr->root_seg+from+i))->is_plus = ((*(c_ptr->root_seg+from+i))->is_plus == 1 ? 0 : 1);
    }

    free(segs);
    return;
}


/*
    Validation functions
*/
void _validate_genome(struct genome* g_ptr, char *source) {
    if (g_ptr == NULL) {
        fprintf(stderr, "Validation error: NULL pointer g_ptr passed to function.\n%s\n", source);
        exit(1);
    }
    if (g_ptr->n_chrs == 0 || g_ptr->root_chr == NULL || *(g_ptr->root_chr+0) == NULL) {
        fprintf(stderr, "Validation error: input genome pointer has no chromosome.\n%s\n", source);
        exit(1);
    }
}
void _validate_chromosome(struct chromosome* c_ptr, char *source) {
    if (c_ptr == NULL) {
        fprintf(stderr, "Validation error: NULL chromosome passed to function with following message.\n%s\n", source);
        exit(1);
    }
    if (c_ptr->n_segs == 0 || c_ptr->root_seg == NULL || *(c_ptr->root_seg+0) == NULL) {
        fprintf(stderr, "Validation error: input chromosome has no segments.\n%s\n", source);
        exit(1);
    }
}
void _validate_seg(struct seg *s_ptr, char *source) {
    if (s_ptr == NULL) {
        fprintf(stderr, "\nValidation error: NULL segment passed to function with following message.\n%s\n", source);
        exit(1);
    }
}
/*
    End validation functions
*/

/*
    Function for simplifying genomes by
    removing unused segment breakpoints.
*/
void key_destroyed(gpointer data) {
    free(data);
}

void simplify_genome(struct genome *g_ptr) {
    // Strategy:
    // 1. Find out which segment breakpoints are not used anymore
    // 2. Join segments at such breakpoints
    // 3. Rename segments

    int c_idx, s_idx, *s_idx_ptr;

    // Get the indexes of each genome_segs member
    char *string_to_be_stored;
    char seg_idx_string[MAX_TIMES_DIVIDED];  // Acts as a temporary string holder for the function
    GHashTable* idx_of_seg = g_hash_table_new_full(g_str_hash, g_str_equal, (GDestroyNotify)key_destroyed, g_free);
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        int_arr_to_string(
            (*(g_ptr->genome_segs+s_idx))->seg_indexes,
            (*(g_ptr->genome_segs+s_idx))->times_divided,
            seg_idx_string
        );

        // Initiate hash for segment indexes
        string_to_be_stored = g_strdup(seg_idx_string);
        s_idx_ptr = malloc(sizeof(int));
        if (s_idx_ptr == NULL) {
            fprintf(stderr, "\nFailed to malloc s_idx_ptr in simplify_genome(). Exiting.\n");
            exit(1);
        }
        *s_idx_ptr = s_idx;
        g_hash_table_insert(
            idx_of_seg,
            string_to_be_stored,
            s_idx_ptr
        );
    }

    // Arrays for remembering which joins are not used, and which segments are completely deleted. 
    int *has_only_natural_joins_with_next, *has_only_natural_joins_with_prev;
    has_only_natural_joins_with_next = malloc(g_ptr->n_genome_segs * sizeof(int));
    if (has_only_natural_joins_with_next == NULL) {
        fprintf(stderr, "\nFailed to malloc has_only_natural_joins_with_next in simplify_genome(). Exiting.\n");
        exit(1);
    }
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        *(has_only_natural_joins_with_next+s_idx) = 1;
    }
    has_only_natural_joins_with_prev = malloc(g_ptr->n_genome_segs * sizeof(int));
    if (has_only_natural_joins_with_prev == NULL) {
        fprintf(stderr, "\nFailed to malloc has_only_natural_joins_with_prev in simplify_genome(). Exiting.\n");
        exit(1);
    }
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        *(has_only_natural_joins_with_prev+s_idx) = 1;
    }

    // First and last segment of each chromosome do not have respective natural joins
    struct chromosome *c_ptr;
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        c_ptr = *(g_ptr->root_chr+c_idx);

        s_idx = 0;
        int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
        s_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);

        if ( (*(c_ptr->root_seg+s_idx))->is_plus == 1 ) {
            *(has_only_natural_joins_with_prev+*s_idx_ptr) = 0;
            if (*s_idx_ptr > 0) {
                *(has_only_natural_joins_with_next+*s_idx_ptr-1) = 0;
            }
        }
        else {
            *(has_only_natural_joins_with_next+*s_idx_ptr) = 0;
            if (*s_idx_ptr < g_ptr->n_genome_segs-1) {
                *(has_only_natural_joins_with_prev+*s_idx_ptr+1) = 0;
            }
        }

        s_idx = c_ptr->n_segs - 1;
        int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
        s_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);

        if ( (*(c_ptr->root_seg+s_idx))->is_plus == 1 ) {
            *(has_only_natural_joins_with_next+*s_idx_ptr) = 0;
            if (*s_idx_ptr < g_ptr->n_genome_segs-1) {
                *(has_only_natural_joins_with_prev+*s_idx_ptr+1) = 0;
            }
        }
        else {
            *(has_only_natural_joins_with_prev+*s_idx_ptr) = 0;
            if (*s_idx_ptr > 0) {
                *(has_only_natural_joins_with_next+*s_idx_ptr-1) = 0;
            }
        }
    }
    

    // Find out which unnatural breakpoints are there in the dataset
    int seg1_idx, seg2_idx;
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        c_ptr = *(g_ptr->root_chr+c_idx);

        for (s_idx=0; s_idx<c_ptr->n_segs-1; s_idx++) {
            int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
            s_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);
            seg1_idx = *s_idx_ptr;
            int_arr_to_string((*(c_ptr->root_seg+s_idx+1))->seg_indexes, (*(c_ptr->root_seg+s_idx+1))->times_divided, seg_idx_string);
            s_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);
            seg2_idx = *s_idx_ptr;

            if (!(
                (
                    (seg1_idx + 1 == seg2_idx && (*(c_ptr->root_seg+s_idx))->is_plus == 1 && (*(c_ptr->root_seg+s_idx+1))->is_plus == 1) ||
                    (seg1_idx - 1 == seg2_idx && (*(c_ptr->root_seg+s_idx))->is_plus == 0 && (*(c_ptr->root_seg+s_idx+1))->is_plus == 0)
                ) &&
                (*(c_ptr->root_seg+s_idx))->is_maternal == (*(c_ptr->root_seg+s_idx+1))->is_maternal  // Both segments must be of the same parental origin
            )) {
                if ((*(c_ptr->root_seg+s_idx))->is_plus == 1) {
                    *(has_only_natural_joins_with_next+seg1_idx) = 0;
                    if (seg1_idx < g_ptr->n_genome_segs - 1) {
                        *(has_only_natural_joins_with_prev+seg1_idx+1) = 0;  // Added bugfix 2014-09-12
                    }
                }
                else {
                    *(has_only_natural_joins_with_prev+seg1_idx) = 0;
                    if (seg1_idx > 0) {
                        *(has_only_natural_joins_with_next+seg1_idx-1) = 0;  // Added bugfix 2014-09-12
                    }
                }
                if ((*(c_ptr->root_seg+s_idx+1))->is_plus == 1) {
                    *(has_only_natural_joins_with_prev+seg2_idx) = 0;
                    if (seg2_idx > 0) {
                        *(has_only_natural_joins_with_next+seg2_idx-1) = 0;  // Added bugfix 2014-09-12
                    }
                }
                else {
                    *(has_only_natural_joins_with_next+seg2_idx) = 0;
                    if (seg2_idx < g_ptr->n_genome_segs - 1) {
                        *(has_only_natural_joins_with_prev+seg2_idx+1) = 0;  // Added bugfix 2014-09-12
                    }
                }
            }
        }
    }

    // Merge segments, first in g_ptr->n_genome_segs and then in each chromosome.
    // In each run of segments to be merged, only the first one is kept. 
    int *index_to_be_removed = malloc(g_ptr->n_genome_segs * sizeof(int));
    if (index_to_be_removed == NULL) {
        fprintf(stderr, "\nFailed to malloc index_to_be_removed in simplify_genome(). Exiting.\n");
        exit(1);
    }
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        *(index_to_be_removed+s_idx) = 0;
    }
    for (s_idx=g_ptr->n_genome_segs-2; s_idx>=0; s_idx--) {
        if (
            (*(has_only_natural_joins_with_next+s_idx)==1 && *(has_only_natural_joins_with_prev+s_idx+1)==1) &&
            *((*(g_ptr->genome_segs+s_idx))->seg_indexes) == *((*(g_ptr->genome_segs+s_idx+1))->seg_indexes)
        ) {
            delete_seg(*(g_ptr->genome_segs+s_idx+1));
            int i;
            for (i=s_idx+1; i<g_ptr->n_genome_segs-1; i++) {
                *(g_ptr->genome_segs+i) = *(g_ptr->genome_segs+i+1);
            }
            g_ptr->n_genome_segs--;
            *(index_to_be_removed+s_idx+1) = 1;
        }
    }
    g_ptr->genome_segs = realloc(g_ptr->genome_segs, g_ptr->n_genome_segs * sizeof(struct seg*));
    if (g_ptr->genome_segs == NULL) {
        fprintf(stderr, "Failed to realloc g_ptr->genome_segs in simplify_genome(). Exiting.\n");
        exit(1);
    }

    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        c_ptr = *(g_ptr->root_chr+c_idx);
    
        for (s_idx=c_ptr->n_segs-1; s_idx>=0; s_idx--) {
            int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
            s_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);
            if (!(*(index_to_be_removed+*s_idx_ptr))) {
                continue;
            }

            delete_seg(*(c_ptr->root_seg+s_idx));
            int i;
            for (i=s_idx; i<c_ptr->n_segs-1; i++) {
                *(c_ptr->root_seg+i) = *(c_ptr->root_seg+i+1);
            }

            c_ptr->n_segs--;
        }

        c_ptr->root_seg = realloc(c_ptr->root_seg, c_ptr->n_segs * sizeof(struct seg*));
        if (c_ptr->root_seg == NULL) {
            fprintf(stderr, "Failed to realloc c_ptr->root_seg in simplify_genome(). Exiting.\n");
            exit(1);
        }
    }


    // Free memory
    g_hash_table_destroy(idx_of_seg);
    free(has_only_natural_joins_with_prev);
    free(has_only_natural_joins_with_next);
    free(index_to_be_removed);
}
/*
    End function for simplifying genomes by removing unused segment breakpoints.
*/


/*
    Functions for printing genomes
*/
void int_arr_to_string(int* int_arr, int arr_size, char *dest) {
    if (arr_size > MAX_TIMES_DIVIDED) {
        fprintf(stderr, "\nin int_arr_to_string() int arr_size > MAX_TIMES_DIVIDED. Exiting.\n");
        exit(1);
    }

    // OBS: char *dest assumed to be lenth MAX_TIMES_DIVIDED

    int i;
    for (i=0; i<arr_size; i++) {
        if (*(int_arr+i) > 61) {
            fprintf(stderr, "Segment index >61 at *(int_arr+%d) in function int_arr_to_string(). Exiting.\n", i);
            exit(1);
        }
        *(dest+i) = 65 + *(int_arr+i);
    }
    *(dest+arr_size) = '\0';

    return;
}

void print_genome(struct genome* g_ptr, char *unique_genome_string) {
    // _validate_genome(g_ptr, "print_genome()");

    int c_idx, s_idx, *cn_ptr;
    char seg_idx_string[MAX_TIMES_DIVIDED];  // Acts as a temporary string holder for the function
    char *string_to_be_stored;
    char diagnostic_info[256];

    g_ptr = copy_genome(g_ptr);
    simplify_genome(g_ptr);

    // Initiate copy numbers for all the segments in current genome
    // Also populate segment index table simultaneously.
    GHashTable* cn_of_seg = g_hash_table_new_full(g_str_hash, g_str_equal, (GDestroyNotify)key_destroyed, g_free);
    GHashTable* idx_of_seg = g_hash_table_new_full(g_str_hash, g_str_equal, (GDestroyNotify)key_destroyed, g_free);
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        sprintf(diagnostic_info, "print_genome genome_segs s_idx %d", s_idx);
        // _validate_seg(*(g_ptr->genome_segs+s_idx), diagnostic_info);
        int_arr_to_string(
            (*(g_ptr->genome_segs+s_idx))->seg_indexes,
            (*(g_ptr->genome_segs+s_idx))->times_divided,
            seg_idx_string
        );
        string_to_be_stored = g_strdup(seg_idx_string);

        // Initiate hash for CNs
        cn_ptr = malloc(2 * sizeof(int));
        if (cn_ptr == NULL) {
            fprintf(stderr, "\nFailed to malloc cn_ptr in print_genome(). Exiting.\n");
            exit(1);
        }
        *(cn_ptr+0) = *(cn_ptr+1) = 0;

        g_hash_table_insert(
            cn_of_seg,
            string_to_be_stored,
            cn_ptr
        );

        // Initiate hash for segment indexes
        string_to_be_stored = g_strdup(seg_idx_string);
        cn_ptr = malloc(sizeof(int));  // Use cn_ptr because too lazy to create another pointer
        if (cn_ptr == NULL) {
            fprintf(stderr, "\nFailed to malloc cn_ptr in print_genome(). Exiting.\n");
            exit(1);
        }
        *cn_ptr = s_idx;
        g_hash_table_insert(
            idx_of_seg,
            string_to_be_stored,
            cn_ptr
        );
    }

    struct chromosome *c_ptr;
    struct seg *s_ptr;
    int *seg_2_idx_ptr;
    int prev_seg_idx, prev_seg_is_plus;
    struct rg {
        int seg1_idx;
        int seg1_is_plus;
        int seg2_idx;
        int seg2_is_plus;
    } cur_rg;

    // Count the copy number for each segment.
    // At the same time, collect rearrangements.
    GArray* rgs = g_array_new(0, 0, sizeof(struct rg));
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        sprintf(diagnostic_info, "print_genome(), c_idx = %d", c_idx);
        c_ptr = *(g_ptr->root_chr+c_idx);
        // _validate_chromosome(c_ptr, diagnostic_info); 
        for (s_idx=0; s_idx<c_ptr->n_segs; s_idx++) {
            // Update copy number for segment
            sprintf(diagnostic_info, "print_genome(), c_idx = %d, s_idx = %d", c_idx, s_idx);
            s_ptr = *(c_ptr->root_seg+s_idx);
            // _validate_seg(s_ptr, diagnostic_info);

            int_arr_to_string(s_ptr->seg_indexes, s_ptr->times_divided, seg_idx_string);
            cn_ptr = (int*)g_hash_table_lookup(cn_of_seg, seg_idx_string);
            *(cn_ptr + s_ptr->is_maternal) += 1;

            // Save the current transition (potential rearrangement) between segments
            seg_2_idx_ptr = (int*)g_hash_table_lookup(idx_of_seg, seg_idx_string);
            if (s_idx != 0) {
                cur_rg.seg1_idx = prev_seg_idx;
                cur_rg.seg1_is_plus = prev_seg_is_plus;
                cur_rg.seg2_idx = *seg_2_idx_ptr;
                cur_rg.seg2_is_plus = s_ptr->is_plus;

                // Reverse the rearrangement paired orientation if needed
                if (
                        cur_rg.seg1_idx > cur_rg.seg2_idx ||
                        (cur_rg.seg1_idx == cur_rg.seg2_idx && cur_rg.seg1_is_plus)
                ) {
                    int temp_int = cur_rg.seg1_idx;
                    cur_rg.seg1_idx = cur_rg.seg2_idx;
                    cur_rg.seg2_idx = temp_int;
                    temp_int = cur_rg.seg1_is_plus;
                    cur_rg.seg1_is_plus = cur_rg.seg2_is_plus ^ 1;
                    cur_rg.seg2_is_plus = temp_int ^ 1;
                }
                g_array_append_val(rgs, cur_rg);
            }
            prev_seg_idx = *seg_2_idx_ptr;
            prev_seg_is_plus = s_ptr->is_plus;
        }
    }

    // Print out current detailed history
    int i;
    for (i=0; i<g_ptr->depth; i++) {
        printf(
            "%s%d%s",
            rg_type_to_txt(*(g_ptr->history+i)),
            *(g_ptr->history_idx+i),
            (i == g_ptr->depth - 1 ? " " : "-")
        );
    }

    // Print out current history
    for (i=0; i<g_ptr->depth; i++) {
        printf(
            "%s%s",
            rg_type_to_txt(*(g_ptr->history+i)),
            (i == g_ptr->depth - 1 ? " " : "-")
        );
    }
    
    // Print out copy numbers for each segment
    char separator;
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        int_arr_to_string(
            (*(g_ptr->genome_segs+s_idx))->seg_indexes,
            (*(g_ptr->genome_segs+s_idx))->times_divided,
            seg_idx_string
        );
        cn_ptr = (int*)g_hash_table_lookup(cn_of_seg, seg_idx_string);
        if (s_idx == g_ptr->n_genome_segs - 1) {
            separator = ' ';
        }
        else if (
            *((*(g_ptr->genome_segs+s_idx  ))->seg_indexes+0) ==
            *((*(g_ptr->genome_segs+s_idx+1))->seg_indexes+0)
        ) {
            separator = '/';
        }
        else {
            separator = ';';
        }
        printf(
            "%d,%d%c",
            *(cn_ptr+0),  // Paternal allele CN
            *(cn_ptr+1),  // Maternal allele CN
            separator
        );
    }

    // Print out current rearrangements
    gint rgs_sort_fun(gconstpointer a, gconstpointer b) {
        struct rg *a_rg = (struct rg*)a;
        struct rg *b_rg = (struct rg*)b;

        if (a_rg->seg1_idx == b_rg->seg1_idx) {
            if (a_rg->seg1_is_plus == 0 && b_rg->seg1_is_plus == 1) {
                return -1;
            }
            else if (b_rg->seg1_is_plus == 0 && a_rg->seg1_is_plus == 1) {
                return 1;
            }
        }
        else if (a_rg->seg1_idx < b_rg->seg1_idx) {
            return -1;
        }
        else {
            return 1;
        }

        // If got this far, then a_rg and b_rg tied for low end of the rearrangement.
        // Try to use second end as tie breaker. 
        if (a_rg->seg2_idx == b_rg->seg2_idx) {
            if (a_rg->seg2_is_plus == 0 && b_rg->seg2_is_plus == 1) {
                return -1;
            }
            else if (b_rg->seg2_is_plus == 0 && a_rg->seg2_is_plus == 1) {
                return 1;
            }
        }
        else if (a_rg->seg2_idx < b_rg->seg2_idx) {
            return -1;
        }
        else {
            return 1;
        }

        return 0;
    }
    int is_rg(struct rg *rg) {
        if (
            (rg->seg1_idx == rg->seg2_idx-1 && rg->seg1_is_plus == 1 && rg->seg2_is_plus == 1) ||
            (rg->seg1_idx == rg->seg2_idx+1 && rg->seg1_is_plus == 0 && rg->seg2_is_plus == 0)
        ) {
            return 0;
        }
        return 1;
    }
    int rg_eq(struct rg *rg1, struct rg *rg2) {
        if (
            rg1->seg1_idx     == rg2->seg1_idx  &&
            rg1->seg1_is_plus == rg2->seg1_is_plus  &&
            rg1->seg2_idx     == rg2->seg2_idx  &&
            rg1->seg2_is_plus == rg2->seg2_is_plus
        ) {
            return 1;
        }

        return 0;
    }
    g_array_sort(rgs, (GCompareFunc)rgs_sort_fun);
    struct rg prev_rg;
    prev_rg.seg1_idx = prev_rg.seg2_idx = prev_rg.seg1_is_plus = prev_rg.seg2_is_plus = -1;
    for (i=0; i<rgs->len; i++) {
        cur_rg = g_array_index(rgs, struct rg, i);
        // if (!is_rg(&cur_rg) || rg_eq(&cur_rg, &prev_rg)) {
        if (!is_rg(&cur_rg)) {
            continue;
        }
        if (rg_eq(&prev_rg, &cur_rg)) {
            continue;
        }
        if (prev_rg.seg1_idx != -1) { printf("/"); }
        printf("%d%s,%d%s", cur_rg.seg1_idx, (cur_rg.seg1_is_plus ? "+" : "-"), cur_rg.seg2_idx, (cur_rg.seg2_is_plus ? "-" : "+"));
        prev_rg = cur_rg;
    }

    // Print unique somatic genome string
    // char out_genome_str[1000];
    // get_unique_genome_string(g_ptr, out_genome_str);
    // printf(" %s", out_genome_str);
    printf(" %s", unique_genome_string);

    printf("\n");

    g_hash_table_destroy(cn_of_seg);
    g_hash_table_destroy(idx_of_seg);
    if (rgs->len > 0) {
        g_array_remove_range(rgs, 0, rgs->len);
    }
    g_array_free(rgs, 1);
    delete_genome(g_ptr);  // Because added simplify_genome() into this function

    return;
}

void print_seg_full(struct seg *s_ptr) {
    printf("      times_divided: %d\n", s_ptr->times_divided);
    printf("      seg_indexes: ");
    int i;
    for (i=0; i<s_ptr->times_divided; i++) {
        printf("%d", *(s_ptr->seg_indexes+i));
    }
    printf("\n");
    printf("      is_plus: %d\n", s_ptr->is_plus);
    printf("      is_maternal: %d\n", s_ptr->is_maternal);

    return;
}
void print_chromosome_full(struct chromosome *c_ptr) {
    int s_idx;
    struct seg *s_ptr;
    printf("  n_segs: %d\n", c_ptr->n_segs);
    for (s_idx=0; s_idx<c_ptr->n_segs; s_idx++) {
        s_ptr = *(c_ptr->root_seg+s_idx);
        printf("    Segment %d\n", s_idx);
        print_seg_full(s_ptr);
    }

    return;
}
void print_genome_full(struct genome *g_ptr) {
    int c_idx, s_idx, i;
    struct chromosome *c_ptr;

    printf("Within genome\n");
    printf("History:\n");
    for (i=0; i<g_ptr->depth; i++) {
        printf("%s%d ", rg_type_to_txt(*(g_ptr->history+i)), *(g_ptr->history_idx+i));
    }
    printf("\n");
    printf("n_genome_segs: %d\n", g_ptr->n_genome_segs);
    for (s_idx=0; s_idx<g_ptr->n_genome_segs; s_idx++) {
        print_seg_full(*(g_ptr->genome_segs+s_idx));
    }

    printf("n_chrs: %d\n", g_ptr->n_chrs);
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        c_ptr = *(g_ptr->root_chr + c_idx);
        printf("  Chromosome %d\n", c_idx);
        print_chromosome_full(c_ptr);
    }

    return;
}
/*
    End functions for printing genomes
*/



/*
    Begin function to create unique representations of rearranged genomes
*/

/*
    Example representation of a genomic string:

    {0,0,0;1,0,0;4,0,1}[3,2]
          |            ^Lengths of reference chromosomes are represented at the end by square brackets.
          |             The numbers here correspond to numbers of segments in each WT chr, i.e. chr1 has
          |             segs 0-2 and chr2 has segs 3-4.
          ^Somatic chromosomes are shown inside curly brackets. Segments are separated by semicolons.
           Segments are described by three numbers. 1st one corresponds to segment ID, second one to
           is_maternal and third one to is_minus.

    Every somatic genome is to be described by the single string that has the smallest lexicographic order. 

    Algorithm. 
    1. Try putting each somatic chromosome as the first one in the genomic string. Build up the
       reference chromosome ordering and orientation as the segments of the reference chromosomes
       are encountered in the somatic chromosomes. 
    2. Prune the tree by removing strings that are not lexicographically smallest.
    3. Extend the branches by adding another somatic chromosome, and again building up the reference
       chromosome orderings and orientations as new reference chromosome segments get used.
    4. Continue until all the somatic chromosome have been used. 
    5. At this point only fully deleted chromosomes have not been encountered. Add these + other
       used reference chromosome segments and finish the reference chromosomes section of the
       string. 
*/

struct genome_string_branch {
    char cur_string[10000];
    int n_somatic_chrs;
    int* somatic_chr_is_used;
    GHashTable *map_of_encountered_segments;
    int next_seg_id;
    int wt_chr_lens[100];
};

struct map_seg_elements {
    int seg_id;
    int reversed;
    int maternal_is_paternal;
};

GHashTable* g_hash_table_copy(GHashTable* hash) {
    GHashTable* new_hash = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    struct map_seg_elements *new_val;
    char *new_key;
    void iterator(gpointer key, gpointer val, gpointer h) {
        new_key = g_strdup(key);
        new_val = malloc(sizeof(struct map_seg_elements));
        if (new_val == NULL) {
            fprintf(stderr, "Failed to malloc new_val in g_hash_table_copy(). Exiting. \n");
            exit(1);
        }
        *new_val = *((struct map_seg_elements*)val);
        g_hash_table_insert(h, new_key, new_val);
    }
    g_hash_table_foreach(hash, (GHFunc)iterator, new_hash);

    return(new_hash);
}

struct genome_string_branch* create_genome_string_branch(struct genome *g_ptr) {
    struct genome_string_branch *gsb_ptr = malloc(sizeof(struct genome_string_branch));

    if (gsb_ptr == NULL) {
        fprintf(stderr, "malloc of gsb_ptr failed in create_genome_string_branch(). Exiting. \n");
        exit(1);
    }
    gsb_ptr->n_somatic_chrs = g_ptr->n_chrs;
    gsb_ptr->somatic_chr_is_used = malloc(gsb_ptr->n_somatic_chrs * sizeof(int));
    if (gsb_ptr->somatic_chr_is_used == NULL) {
        fprintf(stderr, "malloc of gsb_ptr->somatic_chr_is_used failed in create_genome_string_branch(). Exiting. \n");
        exit(1);
    }
    int i;
    for (i=0; i<g_ptr->n_chrs; i++) {
        *(gsb_ptr->somatic_chr_is_used+i) = 0;
    }
    gsb_ptr->map_of_encountered_segments = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    gsb_ptr->cur_string[0] = '\0';
    gsb_ptr->next_seg_id = 0;

    // Initiate these values as well
    for (i=0; i<100; i++) {
        *(gsb_ptr->wt_chr_lens+i) = 0;
    }

    return(gsb_ptr);
}

struct genome_string_branch* copy_genome_string_branch(struct genome_string_branch* gsb_ptr) {
    struct genome_string_branch* new_gsb_ptr = malloc(sizeof(struct genome_string_branch));
    if (new_gsb_ptr == NULL) {
        fprintf(stderr, "malloc of new_gsb_ptr failed in copy_genome_string_branch(). Exiting. \n");
        exit(1);
    }

    strcpy(new_gsb_ptr->cur_string, gsb_ptr->cur_string);

    new_gsb_ptr->n_somatic_chrs = gsb_ptr->n_somatic_chrs;
    new_gsb_ptr->somatic_chr_is_used = malloc(new_gsb_ptr->n_somatic_chrs * sizeof(int));
    if (new_gsb_ptr->somatic_chr_is_used == NULL) {
        fprintf(stderr, "malloc of new_gsb_ptr->somatic_chr_is_used failed in copy_genome_string_branch(). Exiting. \n");
        exit(1);
    }

    // Below doesn't work. WHY???
    // memcpy(new_gsb_ptr->somatic_chr_is_used, gsb_ptr->somatic_chr_is_used, sizeof(int) * new_gsb_ptr->n_somatic_chrs);
    int i;
    for (i=0; i<new_gsb_ptr->n_somatic_chrs; i++) {
        *(new_gsb_ptr->somatic_chr_is_used+i) = *(gsb_ptr->somatic_chr_is_used+i);
    }

    new_gsb_ptr->next_seg_id = gsb_ptr->next_seg_id;
    memcpy(new_gsb_ptr->wt_chr_lens, gsb_ptr->wt_chr_lens, 100 * sizeof(int));

    new_gsb_ptr->map_of_encountered_segments = g_hash_table_copy(gsb_ptr->map_of_encountered_segments);

    return(new_gsb_ptr);
}

void delete_genome_string_branch(struct genome_string_branch *gsb_ptr) {
    free(gsb_ptr->somatic_chr_is_used);
    g_hash_table_destroy(gsb_ptr->map_of_encountered_segments);
    free(gsb_ptr);

    return;
}

void iterate_segs_and_update_map(struct genome_string_branch *gsb_ptr, int s_idx, struct chromosome *c_ptr, GArray *nodes_array, struct genome *g_ptr) {
    char seg_idx_string[MAX_TIMES_DIVIDED];  // Acts as a temporary string holder for the function
    char new_string_to_be_cat[MAX_TIMES_DIVIDED];
    int chr_name, maternal_is_paternal, cur_wt_chr_len, i;
    struct genome_string_branch *new_gsb_ptr;
    char *orig_key = NULL;
    char **orig_key_ptr = &orig_key;
    struct map_seg_elements *orig_val = NULL;
    struct map_seg_elements **orig_val_ptr = &orig_val;
    gboolean key_exists;
    struct map_seg_elements *map_seg;

    int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
    key_exists = g_hash_table_lookup_extended(gsb_ptr->map_of_encountered_segments, seg_idx_string, (gpointer*)orig_key_ptr, (gpointer*)orig_val_ptr);
    map_seg = *orig_val_ptr;
    if (key_exists) {
        // Just calmly proceed to the next s_idx
        sprintf(
            new_string_to_be_cat,
            "%d,%d,%d",
            map_seg->seg_id,
            (map_seg->maternal_is_paternal + (*(c_ptr->root_seg+s_idx))->is_maternal) % 2 == 1 ? 1 : 0,
            (map_seg->reversed + (*(c_ptr->root_seg+s_idx))->is_plus) % 2 == 1 ? 0 : 1
        );
        strncat(gsb_ptr->cur_string, new_string_to_be_cat, strlen(new_string_to_be_cat));

        // Was this last segment?
        if (s_idx == c_ptr->n_segs - 1) {
            strncat(gsb_ptr->cur_string, "}", 1);  // Close up somatic chromosome
            g_array_append_val(nodes_array, gsb_ptr);  // This is the only place where new elements are added to nodes_array
        }
        else {
            strncat(gsb_ptr->cur_string, ";", 1);  // Close up current segment
            iterate_segs_and_update_map(gsb_ptr, s_idx+1, c_ptr, nodes_array, g_ptr);  // Pass the structure to handling the next segment
        }
    }
    else {
        // If the current segment has not been included as one of the included WT chromosomes
        // yet, this has to be done first and creates two new WT chr branches
        chr_name = *((*(c_ptr->root_seg+s_idx))->seg_indexes+0);
        maternal_is_paternal = (*(c_ptr->root_seg+s_idx))->is_maternal;
        
        //
        // In the first case, the WT chromosome is presented in its default orientation
        new_gsb_ptr = copy_genome_string_branch(gsb_ptr);
        cur_wt_chr_len = 0;

        for (i=0; i<g_ptr->n_genome_segs; i++) {
            if ( *((*(g_ptr->genome_segs+i))->seg_indexes+0) == chr_name ) {
                int_arr_to_string((*(g_ptr->genome_segs+i))->seg_indexes, (*(g_ptr->genome_segs+i))->times_divided, seg_idx_string);
                map_seg = malloc(sizeof(struct map_seg_elements));
                if (map_seg == NULL) {
                    fprintf(stderr, "Failed to malloc map_seg in iterate_segs_and_update_map(). Exiting.\n");
                    exit(1);
                }
                map_seg->seg_id = new_gsb_ptr->next_seg_id++;
                map_seg->maternal_is_paternal = maternal_is_paternal;
                map_seg->reversed = 0;
                g_hash_table_insert(new_gsb_ptr->map_of_encountered_segments, g_strdup(seg_idx_string), map_seg);
                cur_wt_chr_len++;
            }
        }

        i = 0;
        while (*(new_gsb_ptr->wt_chr_lens+i) != 0) { i++; }
        *(new_gsb_ptr->wt_chr_lens+i) = cur_wt_chr_len;

        // Now finally we can add the current segment to the string
        int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
        key_exists = g_hash_table_lookup_extended(new_gsb_ptr->map_of_encountered_segments, seg_idx_string, (gpointer*)orig_key_ptr, (gpointer*)orig_val_ptr);
        map_seg = (struct map_seg_elements*)orig_val;
        if (!key_exists) {
            fprintf(stderr, "key_exists is FALSE unexpectedly. Exiting.\n");
            exit(1);
        }
        sprintf(
            new_string_to_be_cat,
            "%d,%d,%d",
            map_seg->seg_id,
            (map_seg->maternal_is_paternal + (*(c_ptr->root_seg+s_idx))->is_maternal) % 2 == 1 ? 1 : 0,
            (map_seg->reversed + (*(c_ptr->root_seg+s_idx))->is_plus) % 2 == 1 ? 0 : 1
        );
        strncat(new_gsb_ptr->cur_string, new_string_to_be_cat, strlen(new_string_to_be_cat));

        // Was this last segment?
        if (s_idx == c_ptr->n_segs - 1) {
            strncat(new_gsb_ptr->cur_string, "}", 1);  // Close up somatic chromosome
            g_array_append_val(nodes_array, new_gsb_ptr);  // This is the only place where new elements are added to nodes_array
        }
        else {
            strncat(new_gsb_ptr->cur_string, ";", 1);  // Close up current segment
            iterate_segs_and_update_map(new_gsb_ptr, s_idx+1, c_ptr, nodes_array, g_ptr);  // Pass the structure to handling the next segment
        }
        // End scenario 1
        //


        //
        // In the second scenario, the WT chromosome is presented in reverse order
        new_gsb_ptr = copy_genome_string_branch(gsb_ptr);
        cur_wt_chr_len = 0;

        for (i=g_ptr->n_genome_segs-1; i>=0; i--) {
            if ( *((*(g_ptr->genome_segs+i))->seg_indexes+0) == chr_name ) {
                int_arr_to_string((*(g_ptr->genome_segs+i))->seg_indexes, (*(g_ptr->genome_segs+i))->times_divided, seg_idx_string);
                map_seg = malloc(sizeof(struct map_seg_elements));
                if (map_seg == NULL) {
                    fprintf(stderr, "Failed to malloc map_seg in iterate_segs_and_update_map(). Exiting.\n");
                    exit(1);
                }
                map_seg->seg_id = new_gsb_ptr->next_seg_id++;
                map_seg->maternal_is_paternal = maternal_is_paternal;
                map_seg->reversed = 1;
                g_hash_table_insert(new_gsb_ptr->map_of_encountered_segments, g_strdup(seg_idx_string), map_seg);
                cur_wt_chr_len++;
            }
        }

        i = 0;
        while (*(new_gsb_ptr->wt_chr_lens+i) != 0) { i++; }
        *(new_gsb_ptr->wt_chr_lens+i) = cur_wt_chr_len;

        // Now finally we can add the current segment to the string
        int_arr_to_string((*(c_ptr->root_seg+s_idx))->seg_indexes, (*(c_ptr->root_seg+s_idx))->times_divided, seg_idx_string);
        key_exists = g_hash_table_lookup_extended(new_gsb_ptr->map_of_encountered_segments, seg_idx_string, (gpointer*)orig_key_ptr, (gpointer*)orig_val_ptr);
        if (!key_exists) {
            fprintf(stderr, "key_exists is FALSE unexpectedly. Exiting.\n");
            exit(1);
        }
        map_seg = (struct map_seg_elements*)orig_val;
        sprintf(
            new_string_to_be_cat,
            "%d,%d,%d",
            map_seg->seg_id,
            (map_seg->maternal_is_paternal + (*(c_ptr->root_seg+s_idx))->is_maternal) % 2 == 1 ? 1 : 0,
            (map_seg->reversed + (*(c_ptr->root_seg+s_idx))->is_plus) % 2 == 1 ? 0 : 1
        );
        strncat(new_gsb_ptr->cur_string, new_string_to_be_cat, strlen(new_string_to_be_cat));

        // Was this last segment?
        if (s_idx == c_ptr->n_segs - 1) {
            strncat(new_gsb_ptr->cur_string, "}", 1);  // Close up somatic chromosome
            g_array_append_val(nodes_array, new_gsb_ptr);  // This is the only place where new elements are added to nodes_array
        }
        else {
            strncat(new_gsb_ptr->cur_string, ";", 1);  // Close up current segment
            iterate_segs_and_update_map(new_gsb_ptr, s_idx+1, c_ptr, nodes_array, g_ptr);  // Pass the structure to handling the next segment
        }
        // End scenario 2
        //


        // In this part we branched so delete this internal node
        delete_genome_string_branch(gsb_ptr);
    }
}

void remove_non_smallest_members_from_genome_branch_array(GArray *nodes) {
    char min_genome_string[10000];
    struct genome_string_branch *debug_gsb_ptr;
    int i;
    
    strcpy(min_genome_string, (g_array_index(nodes, struct genome_string_branch*, 0))->cur_string);
    for (i=1; i<nodes->len; i++) {
        if (strcmp(min_genome_string, (g_array_index(nodes, struct genome_string_branch*, i))->cur_string) > 0) {
            strcpy(min_genome_string, (g_array_index(nodes, struct genome_string_branch*, i))->cur_string);
        }
    }
    i=0;
    while (i < nodes->len) {
        if (strcmp((g_array_index(nodes, struct genome_string_branch*, i))->cur_string, min_genome_string) > 0) {
            debug_gsb_ptr = g_array_index(nodes, struct genome_string_branch*, i);
            g_array_remove_index(nodes, i);
        }
        else {
            i++;
        }
    }
}

// void nodes_array_element_free(struct genome_string_branch* gsb_ptr) {
void nodes_array_element_free(gpointer gsb_ptr_ptr) {
    struct genome_string_branch *gsb_ptr = *((struct genome_string_branch **)gsb_ptr_ptr);
    free(gsb_ptr->somatic_chr_is_used);
    gsb_ptr->somatic_chr_is_used = NULL;
    g_hash_table_destroy(gsb_ptr->map_of_encountered_segments);
    gsb_ptr->map_of_encountered_segments = NULL;
    free(gsb_ptr);

    return;
}

void get_unique_genome_string(struct genome* g_ptr, char* out_string_ptr) {
    int c_idx, s_idx, i, j, n_chrs_used, cur_node_count;
    char seg_idx_string[MAX_TIMES_DIVIDED];  // Acts as a temporary string holder for the function
    struct chromosome *c_ptr;
    struct genome_string_branch *gsb_ptr;

    // Use a GArray to store all the terminal nodes
    GArray *nodes = g_array_new(0, 0, sizeof(struct genome_string_branch*));
    g_array_set_clear_func(nodes, (GDestroyNotify)nodes_array_element_free);
    
    // Initiate by taking all somatic chromosomes and both orientations as the root somatic chromosome
    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        // Present this chr in forward orientation
        c_ptr = *(g_ptr->root_chr+c_idx);
        gsb_ptr = create_genome_string_branch(g_ptr);
        strncat(gsb_ptr->cur_string, "{", 1);  // Start of first somatic chromosome
        *(gsb_ptr->somatic_chr_is_used+c_idx) = 1; // Mark current somatic chromosome as used
        iterate_segs_and_update_map(gsb_ptr, 0, c_ptr, nodes, g_ptr);  // After this step, genome_string_branch have been extended by an additional level of depth

        // Present this chr in reverse orientation
        c_ptr = copy_chromosome(c_ptr);
        invert_segs_in_chr(c_ptr, 0, c_ptr->n_segs - 1);
        gsb_ptr = create_genome_string_branch(g_ptr);
        strncat(gsb_ptr->cur_string, "{", 1);
        *(gsb_ptr->somatic_chr_is_used+c_idx) = 1; // Mark current somatic chromosome as used
        iterate_segs_and_update_map(gsb_ptr, 0, c_ptr, nodes, g_ptr);
        delete_chromosome(c_ptr);
    }

    // Compare and remove bad starts
    remove_non_smallest_members_from_genome_branch_array(nodes);
    n_chrs_used = 1;

    // Repeat until all chromosomes have been used up
    while (n_chrs_used < g_ptr->n_chrs) {
        cur_node_count = nodes->len;

        for (i=0; i<cur_node_count; i++) {
            for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
                gsb_ptr = copy_genome_string_branch(g_array_index(nodes, struct genome_string_branch*, i));

                if (*(gsb_ptr->somatic_chr_is_used+c_idx)) {
                    delete_genome_string_branch(gsb_ptr);
                    continue;
                }
                c_ptr = *(g_ptr->root_chr+c_idx);

                // Forward orientation of current somatic chromosome
                *(gsb_ptr->somatic_chr_is_used+c_idx) = 1;
                strncat(gsb_ptr->cur_string, "{", 1);  // Start of the next somatic chromosome
                iterate_segs_and_update_map(gsb_ptr, 0, c_ptr, nodes, g_ptr);


                // Reverse orientation of current somatic chromosome
                gsb_ptr = copy_genome_string_branch(g_array_index(nodes, struct genome_string_branch*, i));
                *(gsb_ptr->somatic_chr_is_used+c_idx) = 1;
                strncat(gsb_ptr->cur_string, "{", 1);
                c_ptr = copy_chromosome(c_ptr);
                invert_segs_in_chr(c_ptr, 0, c_ptr->n_segs - 1);
                iterate_segs_and_update_map(gsb_ptr, 0, c_ptr, nodes, g_ptr);
                delete_chromosome(c_ptr);
            }
        }

        // Remove the nodes created at the previous round
        g_array_remove_range(nodes, 0, cur_node_count);
        cur_node_count = nodes->len;
        n_chrs_used++;

        // Remove nodes with genome strings larger than minimum genome string
        remove_non_smallest_members_from_genome_branch_array(nodes);
    }

    // Now all chromosomes have been added to nodes->genome_string_branch variables, and all the
    // remaining variables have the lexicographically smallest representation. Now simply add
    // the WT chr lengths to the strings.
    for (i=0; i<nodes->len; i++) {
        gsb_ptr = g_array_index(nodes, struct genome_string_branch*, i);
        strncat(gsb_ptr->cur_string, "[", 1);

        j=0;
        while (*(gsb_ptr->wt_chr_lens+j) != 0) {
            if (j != 0) {
                strncat(gsb_ptr->cur_string, ",", 1);
            }
            sprintf(seg_idx_string, "%d", *(gsb_ptr->wt_chr_lens+j));
            strncat(gsb_ptr->cur_string, seg_idx_string, strlen(seg_idx_string));

            j++;
        }
        strncat(gsb_ptr->cur_string, "]", 1);
    }

    remove_non_smallest_members_from_genome_branch_array(nodes);
    gsb_ptr = g_array_index(nodes, struct genome_string_branch*, 0);
    strcpy(out_string_ptr, gsb_ptr->cur_string);
    g_array_free(nodes, 1);

    return;
}

/*
    End function to create unique representation of rearranged genomes
*/


