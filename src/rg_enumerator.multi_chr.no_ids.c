/*
    These functions act on genomes consisting of a single chromosome. 
    Functions take genomes as input and output nothing. 
    
    Key is the bridge function that applies a given genome to all the possible downstream
    possible rearrangement functions.

    Decision whether to bridge or terminate is dependent on a predefined depth.
*/

extern int MAX_DEPTH_DUP, MAX_DEPTH_NONDUP;
extern GHashTable *seen_somatic_genomes;

void bridge(struct genome *g_ptr);
void enum_dels(struct genome *g_ptr);
void enum_tds(struct genome *g_ptr);
void enum_invs(struct genome *g_ptr);
void enum_inv_dups(struct genome *g_ptr);
void enum_tel_break(struct genome *g_ptr);
void enum_bal_transloc(struct genome *g_ptr);
void enum_unbal_transloc(struct genome *g_ptr);
void enum_wc_dup(struct genome *g_ptr);
void enum_wc_del(struct genome *g_ptr);
void enum_wg_dup(struct genome *g_ptr);
void enum_fbs(struct genome *g_ptr);

void bridge(struct genome *g_ptr) {
    if (g_ptr->depth < MAX_DEPTH_NONDUP) {
        enum_dels(g_ptr);
        enum_invs(g_ptr);
        enum_tel_break(g_ptr);
        enum_bal_transloc(g_ptr);
        enum_unbal_transloc(g_ptr);
        if (g_ptr->n_chrs > 1) {
            enum_wc_del(g_ptr);
        }

        if (g_ptr->dup_depth < MAX_DEPTH_DUP) {
            enum_tds(g_ptr);
            // enum_inv_dups(g_ptr);
            enum_fbs(g_ptr);
            enum_wc_dup(g_ptr);

            if (g_ptr->wgd_depth == 0) {
                enum_wg_dup(g_ptr);
            }
        }
    }

    delete_genome(g_ptr);

    return;
}


/*
    Little helper function
*/
char *get_detailed_history(struct genome *g_ptr) {
    char *bfr = malloc(100 * sizeof(char));  // Just allow 100 characters should be plenty
    char bfr2[100];
    if (bfr == NULL) {
        fprintf(stderr, "Failed to malloc bfb in get_detailed_history(). Exiting. \n");
        exit(1);
    }
    *bfr = '\0';

    int i;
    for (i=0; i<g_ptr->depth; i++) {
        sprintf(
            bfr2,
            "%s%d%s",
            rg_type_to_txt(*(g_ptr->history+i)),
            *(g_ptr->history_idx+i),
            (i == g_ptr->depth - 1 ? " " : "-")
        );
        strncat(bfr, bfr2, strlen(bfr2));
    }

    return(bfr);
}

int get_dup_depth_from_genome_history_string(char *genome_hist) {
    int depth = 0;
    char bfr[4];
    char c;
    int i = 0, j = 0;
    while (*(genome_hist+i) != '\0') {
        c = *(genome_hist+i);
        if (isdigit(c)) {
            i++;
            continue;
        }
        else if (c == '-') {
            bfr[j] = '\0';
            if      (strcmp(bfr, "td")  == 0) { depth++; }
            else if (strcmp(bfr, "id")  == 0) { depth++; }
            else if (strcmp(bfr, "fb")  == 0) { depth++; }
            else if (strcmp(bfr, "wcg") == 0) { depth++; }
            else if (strcmp(bfr, "wgd") == 0) { depth++; }
            j = 0;
        }
        else {
            bfr[j++] = c;
        }

        i++;
    }

    return(depth);
}

int get_overall_depth_from_genome_history_string(char *genome_hist) {
    int depth = 1;
    char c;
    int i = 0;
    while (*(genome_hist+i) != '\0') {
        c = *(genome_hist+i);
        if (c == '-') {
            depth++;
        }
        i++;
    }

    return(depth);
}

void handle_next_step(struct genome *g_ptr) {
    char *previous_somatic_genome;
    char unique_genome_string[2000];
    int prev_depth, prev_dup_depth;

    simplify_genome(g_ptr);
    get_unique_genome_string(g_ptr, unique_genome_string);  // Get the unique genome string of *g_ptr and assign to genome_string_ptr

    if (g_hash_table_contains(seen_somatic_genomes, unique_genome_string)) {
        previous_somatic_genome = (char*)g_hash_table_lookup(seen_somatic_genomes, unique_genome_string);
        prev_depth = get_overall_depth_from_genome_history_string(previous_somatic_genome);
        prev_dup_depth = get_dup_depth_from_genome_history_string(previous_somatic_genome);

        // Previous genome with the same configuration as the current one was reached with fewer events?
        if (prev_depth <= g_ptr->depth && prev_dup_depth <= g_ptr->dup_depth) {
            print_genome(g_ptr, previous_somatic_genome);
            delete_genome(g_ptr);  // Need to delete here since not passing to bridge().
        }
        else {
            print_genome(g_ptr, unique_genome_string);
            g_hash_table_replace(seen_somatic_genomes, g_strdup(unique_genome_string), get_detailed_history(g_ptr));  // No need to free this since need to keep in memory
            bridge(g_ptr);
        }
    }
    else {
        print_genome(g_ptr, unique_genome_string);
        g_hash_table_insert(seen_somatic_genomes, g_strdup(unique_genome_string), get_detailed_history(g_ptr));  // No need to free this since need to keep in memory
        bridge(g_ptr);
    }
}

void handle_next_step_after_fold_back(struct genome *g_ptr) {
    char *previous_somatic_genome;
    char unique_genome_string[2000];
    int prev_depth, prev_dup_depth;

    simplify_genome(g_ptr);
    get_unique_genome_string(g_ptr, unique_genome_string);  // Get the unique genome string of *g_ptr and assign to genome_string_ptr

    if (g_hash_table_contains(seen_somatic_genomes, unique_genome_string)) {
        previous_somatic_genome = (char*)g_hash_table_lookup(seen_somatic_genomes, unique_genome_string);
        prev_depth = get_overall_depth_from_genome_history_string(previous_somatic_genome);
        prev_dup_depth = get_dup_depth_from_genome_history_string(previous_somatic_genome);

        if (prev_depth <= g_ptr->depth && prev_dup_depth <= g_ptr->dup_depth) {
            print_genome(g_ptr, previous_somatic_genome);
        }
        else {
            print_genome(g_ptr, unique_genome_string);
            g_hash_table_replace(seen_somatic_genomes, g_strdup(unique_genome_string), get_detailed_history(g_ptr));  // No need to free this since need to keep in memory
            if (g_ptr->depth < MAX_DEPTH_NONDUP) {
                enum_tel_break(g_ptr);
                if (g_ptr->dup_depth < MAX_DEPTH_DUP) {
                    enum_fbs(g_ptr);
                }
            }
        }
    }
    else {
        print_genome(g_ptr, unique_genome_string);
        g_hash_table_insert(seen_somatic_genomes, g_strdup(unique_genome_string), get_detailed_history(g_ptr));  // No need to free this since need to keep in memory
        if (g_ptr->depth < MAX_DEPTH_NONDUP) {
            enum_tel_break(g_ptr);
            if (g_ptr->dup_depth < MAX_DEPTH_DUP) {
                enum_fbs(g_ptr);
            }
        }
    }

    delete_genome(g_ptr);
}
/*
    End helper functions
*/


/*
    Rg enumeration functions
*/
void enum_dels(struct genome *g_ptr) {
    // _validate_genome(g_ptr, "enum_dels()");
    int hist_idx = 0;


    /* Declare reusable variables */
    int b1, b2;
    struct genome *new_g_ptr;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c_idx;

    /*
        Two rules for selecting deletion breakpoints:
        1. b1 < b2 in chromosome coordinates.
        2. b1 and b2 can happen in same or different segments.
        3. b1 and b2 denote segments at which the breaks happen.
    */

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            /* First case: both breakpoints at exactly same segment
               In this case the affected segment is broken into three pieces. */
            b2 = b1;
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, DEL, hist_idx++);

            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));


            /*
                Splice affected segment into three.
                Then delete the centre segment.
                Then splice the remaining similar segments.
            */
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), b1+1, b1+1);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );
            free(b1_seg_indexes);

            handle_next_step(new_g_ptr);


            /* Next, when the two breakpoints occur at two physically different DNA segments */
            int two_segments_look_identical;
            int delete_from, delete_to;
            for (b2 = b1+1; b2 < (*(g_ptr->root_chr+c_idx))->n_segs; b2++) {
                /* Do the two affected segments have the same seg_indexes? */
                two_segments_look_identical = int_array_cmp(
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes
                );
                if (two_segments_look_identical) {
                    /* If the two segments have the same indexes,
                       there are two ways the two breakpoints can occur
                       in the segment. */

                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    /* Option 1: looking at plus strand, s1 < s2 at the segment.
                       Segment in question is broken into following:
                          b1  b2
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, DEL, hist_idx++);

                    /*
                        Splice the two affected segments.
                        Then deleted everything in between.
                        Then splice the remaining similar segments.
                    */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    delete_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 1 : b1 + 2;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    delete_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 1 : 0);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), delete_from, delete_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    
                    /* Option 2: looking at plus strand, s1 > s2 at the segment.
                       Segment in question is broken into following:
                          b2  b1
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, DEL, hist_idx++);

                    /* Same story as above. */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    delete_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 2 : b1 + 1;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    delete_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 0 : 1);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), delete_from, delete_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                }
                else {
                    /* If the two broken segments have different indexes,
                       simply break them and splice intervening segments out. */

                    /* But first, let's store the identities of the affected segments first */
                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided;
                    b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
                    if (b2_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));

                    
                    /* Finally, let's break segments */
                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, DEL, hist_idx++);

                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+1, 2);
                    delete_from = b1 + 1;
                    delete_to = 1 + b2;
                    delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), delete_from, delete_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 2);
                    splice_all_segs(new_g_ptr, b2_seg_indexes, b2_seg_indexes_ln, 2);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                    free(b2_seg_indexes);
                }
            }
        }  // for b1
    }  // for each c_idx

    return;
}


void enum_tds(struct genome *g_ptr) {
    // _validate_genome(g_ptr, "enum_tds()");

    /* Declare reusable variables */
    int hist_idx = 0;
    int b1, b2;
    struct genome *new_g_ptr;
    struct chromosome *segs_to_be_dup;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c_idx;

    /*
        Two rules for selecting duplication breakpoints:
        1. b1 < b2 in chromosome coordinates.
        2. b1 and b2 can happen in same or different segments.
        3. b1 and b2 denote segments at which the breaks happen.
    */

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            /* First case: both breakpoints at exactly same segment
               In this case the affected segment is broken into three pieces. */
            b2 = b1;
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, TD, hist_idx++);

            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            /*
                Splice affected segment into three.
                Then yank the centre segment and insert it after the yanked position.
                Then splice the remaining similar segments.
            */
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
            segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), b1+1, b1+1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, b1+2);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );
            delete_chromosome(segs_to_be_dup);
            free(b1_seg_indexes);

            handle_next_step(new_g_ptr);


            // Next, when the two breakpoints occur at two physically different DNA segments
            int two_segments_look_identical;
            int yank_from, yank_to;
            for (b2 = b1+1; b2 < (*(g_ptr->root_chr+c_idx))->n_segs; b2++) {
                /* Do the two affected segments have the same seg_indexes? */
                two_segments_look_identical = int_array_cmp(
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes
                );
                if (two_segments_look_identical) {
                    /* If the two segments have the same indexes,
                       there are two ways the two breakpoints can occur
                       in the segment. */

                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    /* Option 1: looking at plus strand, s1 < s2 at the segment.
                       Segment in question is broken into following:
                          b1  b2
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, TD, hist_idx++);

                    /*
                        Splice the two affected segments.
                        Then deleted everything in between.
                        Then splice the remaining similar segments.
                    */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    yank_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 1 : b1 + 2;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    yank_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 1 : 0);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);
                    delete_chromosome(segs_to_be_dup);

                    handle_next_step(new_g_ptr);

                    
                    /* Option 2: looking at plus strand, s1 > s2 at the segment.
                       Segment in question is broken into following:
                          b2  b1
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, TD, hist_idx++);

                    /* Same story as above. */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    yank_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 2 : b1 + 1;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    yank_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 0 : 1);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);
                    delete_chromosome(segs_to_be_dup);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                }
                else {
                    /* If the two broken segments have different indexes,
                       simply break them and splice intervening segments out. */

                    /* But first, let's store the identities of the affected segments first */
                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided;
                    b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
                    if (b2_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));

                    
                    /* Finally, let's break segments */
                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, TD, hist_idx++);

                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+1, 2);
                    yank_from = b1 + 1;
                    yank_to = 1 + b2;
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 2);
                    splice_all_segs(new_g_ptr, b2_seg_indexes, b2_seg_indexes_ln, 2);
                    delete_chromosome(segs_to_be_dup);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                    free(b2_seg_indexes);
                }
            }
        }  // for each b1
    }  // for each c_idx

    return;
}

void enum_invs(struct genome *g_ptr) {
    // _validate_genome(g_ptr, "enum_invs()");

    /* Declare reusable variables */
    int hist_idx = 0;
    int b1, b2;
    struct genome *new_g_ptr;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c_idx;

    /*
        Two rules for selecting inversion breakpoints:
        1. b1 < b2 in chromosome coordinates.
        2. b1 and b2 can happen in same or different segments.
        3. b1 and b2 denote segments at which the breaks happen.
    */

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            /* First case: both breakpoints at exactly same segment
               In this case the affected segment is broken into three pieces. */
            b2 = b1;
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, INV, hist_idx++);

            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            /*
                Splice affected segment into three.
                Then yank the centre segment and insert it after the yanked position.
                Then splice the remaining similar segments.
            */
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
            invert_segs_in_chr(*(new_g_ptr->root_chr+c_idx), b1+1, b1+1);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );
            free(b1_seg_indexes);

            handle_next_step(new_g_ptr);


            // Next, when the two breakpoints occur at two physically different DNA segments
            int two_segments_look_identical;
            int inv_from, inv_to;
            for (b2 = b1+1; b2 < (*(g_ptr->root_chr+c_idx))->n_segs; b2++) {
                /* Do the two affected segments have the same seg_indexes? */
                two_segments_look_identical = int_array_cmp(
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes
                );
                if (two_segments_look_identical) {
                    /* If the two segments have the same indexes,
                       there are two ways the two breakpoints can occur
                       in the segment. */

                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    /* Option 1: looking at plus strand, s1 < s2 at the segment.
                       Segment in question is broken into following:
                          b1  b2
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV, hist_idx++);

                    /*
                        Splice the two affected segments.
                        Then deleted everything in between.
                        Then splice the remaining similar segments.
                    */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    inv_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 1 : b1 + 2;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    inv_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 1 : 0);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    invert_segs_in_chr(*(new_g_ptr->root_chr+c_idx), inv_from, inv_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    
                    /* Option 2: looking at plus strand, s1 > s2 at the segment.
                       Segment in question is broken into following:
                          b2  b1
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV, hist_idx++);

                    /* Same story as above. */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    inv_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 2 : b1 + 1;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    inv_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 0 : 1);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    invert_segs_in_chr(*(new_g_ptr->root_chr+c_idx), inv_from, inv_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                }
                else {
                    /* If the two broken segments have different indexes,
                       simply break them and splice intervening segments out. */

                    /* But first, let's store the identities of the affected segments first */
                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided;
                    b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
                    if (b2_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));

                    
                    /* Finally, let's break segments */
                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV, hist_idx++);

                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+1, 2);
                    inv_from = b1 + 1;
                    inv_to = 1 + b2;
                    invert_segs_in_chr(*(new_g_ptr->root_chr+c_idx), inv_from, inv_to);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 2);
                    splice_all_segs(new_g_ptr, b2_seg_indexes, b2_seg_indexes_ln, 2);

                    handle_next_step(new_g_ptr);

                    free(b1_seg_indexes);
                    free(b2_seg_indexes);
                }
            }
        }  // for each b1
    }  // for each c_idx

    return;
}


void enum_inv_dups(struct genome *g_ptr) {
    // _validate_genome(g_ptr, "enum_inv_dups()");

    /* Declare reusable variables */
    int hist_idx = 0;
    int b1, b2;
    struct genome *new_g_ptr, *new_g_ptr2;
    struct chromosome *segs_to_be_dup;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c_idx;

    /*
        Two rules for selecting duplication breakpoints:
        1. b1 < b2 in chromosome coordinates.
        2. b1 and b2 can happen in same or different segments.
        3. b1 and b2 denote segments at which the breaks happen.
    */

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            /* First case: both breakpoints at exactly same segment
               In this case the affected segment is broken into three pieces. */
            b2 = b1;
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, INV_DUP, hist_idx++);
            new_g_ptr2 = copy_genome(g_ptr);
            make_history(new_g_ptr2, INV_DUP, hist_idx++);

            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            /*
                Splice affected segment into three.
                Then yank the centre segment and insert it after the yanked position.
                Then splice the remaining similar segments.
            */
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
            segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), b1+1, b1+1);
            invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs-1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, b1+2);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b1, 3);
            insert_segs_into_chr(*(new_g_ptr2->root_chr+c_idx), segs_to_be_dup, b1+1);
            splice_all_segs(
                new_g_ptr2,         /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr2);

            delete_chromosome(segs_to_be_dup);
            free(b1_seg_indexes);


            // Next, when the two breakpoints occur at two physically different DNA segments
            int two_segments_look_identical;
            int yank_from, yank_to;
            for (b2 = b1+1; b2 < (*(g_ptr->root_chr+c_idx))->n_segs; b2++) {
                /* Do the two affected segments have the same seg_indexes? */
                two_segments_look_identical = int_array_cmp(
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided,
                    (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes
                );
                if (two_segments_look_identical) {
                    /* If the two segments have the same indexes,
                       there are two ways the two breakpoints can occur
                       in the segment. */

                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    /* Option 1: looking at plus strand, s1 < s2 at the segment.
                       Segment in question is broken into following:
                          b1  b2
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV_DUP, hist_idx++);
                    new_g_ptr2 = copy_genome(g_ptr);
                    make_history(new_g_ptr2, INV_DUP, hist_idx++);

                    /*
                        Splice the two affected segments.
                        Then deleted everything in between.
                        Then splice the remaining similar segments.
                    */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    yank_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 1 : b1 + 2;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    yank_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 1 : 0);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs-1);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b2+2, 3);
                    insert_segs_into_chr(*(new_g_ptr2->root_chr+c_idx), segs_to_be_dup, yank_from);
                    splice_all_segs(new_g_ptr2, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr2);

                    delete_chromosome(segs_to_be_dup);
                    

                    /* Option 2: looking at plus strand, s1 > s2 at the segment.
                       Segment in question is broken into following:
                          b2  b1
                        A ^ B ^ C
                       === === ===>
                    */

                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV_DUP, hist_idx++);
                    new_g_ptr2 = copy_genome(g_ptr);
                    make_history(new_g_ptr2, INV_DUP, hist_idx++);

                    /* Same story as above. */
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+2, 3);
                    yank_from = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->is_plus ? b1 + 2 : b1 + 1;  /* After splicing, segment at b1 becomes
                                                                                                               three segments. */
                    yank_to = 2 + b2 + ( (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b2+2))->is_plus ? 0 : 1);  /* 2 + b2 comes from the fact that segment at
                                                                                                               b1 has now been split from one to 3 segments. */
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs-1);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr);

                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b1, 3);
                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b2+2, 3);
                    insert_segs_into_chr(*(new_g_ptr2->root_chr+c_idx), segs_to_be_dup, yank_from);
                    splice_all_segs(new_g_ptr2, b1_seg_indexes, b1_seg_indexes_ln, 3);

                    handle_next_step(new_g_ptr2);

                    delete_chromosome(segs_to_be_dup);
                    free(b1_seg_indexes);
                }
                else {
                    /* If the two broken segments have different indexes,
                       simply break them and splice intervening segments out. */

                    /* But first, let's store the identities of the affected segments first */
                    b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
                    b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
                    if (b1_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

                    b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->times_divided;
                    b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
                    if (b2_seg_indexes == NULL) {
                        fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                        exit(1);
                    }
                    memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));

                    
                    /* Finally, let's break segments */
                    new_g_ptr = copy_genome(g_ptr);
                    make_history(new_g_ptr, INV_DUP, hist_idx++);
                    new_g_ptr2 = copy_genome(g_ptr);
                    make_history(new_g_ptr2, INV_DUP, hist_idx++);

                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
                    splice_one_seg(*(new_g_ptr->root_chr+c_idx), b2+1, 2);
                    yank_from = b1 + 1;
                    yank_to = 1 + b2;
                    segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), yank_from, yank_to);
                    invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs-1);
                    insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, yank_to + 1);
                    splice_all_segs(new_g_ptr, b1_seg_indexes, b1_seg_indexes_ln, 2);
                    splice_all_segs(new_g_ptr, b2_seg_indexes, b2_seg_indexes_ln, 2);

                    handle_next_step(new_g_ptr);

                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b1, 2);
                    splice_one_seg(*(new_g_ptr2->root_chr+c_idx), b2+1, 2);
                    insert_segs_into_chr(*(new_g_ptr2->root_chr+c_idx), segs_to_be_dup, yank_from);
                    splice_all_segs(new_g_ptr2, b1_seg_indexes, b1_seg_indexes_ln, 2);
                    splice_all_segs(new_g_ptr2, b2_seg_indexes, b2_seg_indexes_ln, 2);

                    handle_next_step(new_g_ptr2);

                    delete_chromosome(segs_to_be_dup);
                    free(b1_seg_indexes);
                    free(b2_seg_indexes);
                }
            }
        }  // for each b1
    }  // for each c_idx

    return;
}


//
// BFB related enumerators
//
void enum_tel_break(struct genome *g_ptr) {  // Telomeric break without fold-back rearrangement
    // _validate_genome(g_ptr, "enum_dels()");
    int hist_idx = 0;

    /* Declare reusable variables */
    int b1;
    struct genome *new_g_ptr;
    int *b1_seg_indexes, b1_seg_indexes_ln;
    int c_idx;

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            //
            // Left telomere, no fusion, but neotelomerization
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, TEL_BREAK, hist_idx++);

            // Below has to be done only once
            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), 0, b1);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            //
            // Right telomere, no fusion, but neotelomerization
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, TEL_BREAK, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), b1+1, (*(new_g_ptr->root_chr+c_idx))->n_segs - 1);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            free(b1_seg_indexes);
        }  // for each b1
    }  // for each c_idx
}

void enum_fbs(struct genome *g_ptr) {  // Fold-back rearrangement
    // _validate_genome(g_ptr, "enum_dels()");
    int hist_idx = 0;

    /* Declare reusable variables */
    int b1;
    struct genome *new_g_ptr;
    struct chromosome *segs_to_be_dup;
    int *b1_seg_indexes, b1_seg_indexes_ln;
    int c_idx;

    for (c_idx = 0; c_idx < g_ptr->n_chrs; c_idx++) {
        for (b1 = 0; b1 < (*(g_ptr->root_chr+c_idx))->n_segs; b1++) {
            //
            // Left telomere, telomeric fusion
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, FOLD_BACK, hist_idx++);

            // Below has to be done only once
            b1_seg_indexes_ln = (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "\nFailed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(new_g_ptr->root_chr+c_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), 0, b1);
            segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), 0, (*(new_g_ptr->root_chr+c_idx))->n_segs - 1);
            invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs - 1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, 0);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step_after_fold_back(new_g_ptr);

            delete_chromosome(segs_to_be_dup);

            //
            // Right telomere, telomeric fusion
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, FOLD_BACK, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c_idx), b1, 2);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c_idx), b1+1, (*(new_g_ptr->root_chr+c_idx))->n_segs - 1);
            segs_to_be_dup = yank_segments(*(new_g_ptr->root_chr+c_idx), 0, (*(new_g_ptr->root_chr+c_idx))->n_segs - 1);
            invert_segs_in_chr(segs_to_be_dup, 0, segs_to_be_dup->n_segs - 1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c_idx), segs_to_be_dup, (*(new_g_ptr->root_chr+c_idx))->n_segs);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step_after_fold_back(new_g_ptr);

            delete_chromosome(segs_to_be_dup);

            free(b1_seg_indexes);
        }  // for each b1
    }  // for each c_idx
}
//
// End BFB related enumerators
//



void enum_bal_transloc(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "enum_transloc()");
    int hist_idx = 0;

    /* Declare reusable variables */
    int b1, b2;
    struct genome *new_g_ptr;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c1_idx, c2_idx, two_segments_look_identical;
    struct chromosome *seg_holder1, *seg_holder2;
    
    // Go through all chromosomes and all segments
    for (c1_idx=0; c1_idx<g_ptr->n_chrs; c1_idx++) {
    for (c2_idx=c1_idx+1; c2_idx<g_ptr->n_chrs; c2_idx++) {
    for (b1=0; b1 < (*(g_ptr->root_chr+c1_idx))->n_segs; b1++) {
    for (b2=0; b2 < (*(g_ptr->root_chr+c2_idx))->n_segs; b2++) {
        two_segments_look_identical = int_array_cmp(
            (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided,
            (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes,
            (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->times_divided,
            (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->seg_indexes
        );

        // Are the two affected segments the same segment?
        if (two_segments_look_identical) {
            b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));


            /* Option 1: looking at plus strand, b1 < b2 at the segment.
               Segment in question is broken into following:
                  b1  b2
                A ^ B ^ C
               === === ===>
            */

            // Case 1A, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(  // Get the q-telomeric pieces
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1
            );
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(  // Remove the translocated piece from between
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);
            

            // Case 1B: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2)
            );
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c1_idx),
                seg_holder2,
                (*(new_g_ptr->root_chr+c1_idx))->n_segs
            );
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c2_idx),
                seg_holder1,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1)
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2)
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            /* Option 2: looking at plus strand, b1 > b2 at the segment.
               Segment in question is broken into following:
                  b2  b1
                A ^ B ^ C
               === === ===>
            */

            // Case 2A, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(  // Get the q-telomeric pieces
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1
            );
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(  // Remove the translocated piece from between
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);
            

            // Case 2B: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2 : b2+1)
            );
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c1_idx),
                seg_holder2,
                (*(new_g_ptr->root_chr+c1_idx))->n_segs
            );
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c2_idx),
                seg_holder1,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2)
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2 : b2+1)
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            free(b1_seg_indexes);
        }
        else {
            // We are here because the two affected segments are not the same

            b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->times_divided;
            b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
            if (b2_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));


            // First case, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 2);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 2);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1);  // Get the q-telomeric pieces
            seg_holder2 = yank_segments(*(new_g_ptr->root_chr+c2_idx), b2+1, (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs);  // Remove the translocated piece from between
            delete_segs_from_chr(*(new_g_ptr->root_chr+c2_idx), b2+1, (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b2_seg_indexes,     /* The index of the segment to be split */
                b2_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            // Second case: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            make_history(new_g_ptr, BAL_TRANSLOC, hist_idx++);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 2);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 2);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1);  // Get the q-telomeric pieces
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(*(new_g_ptr->root_chr+c2_idx), 0, b2);
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, b2+1);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c2_idx), 0, b2);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b2_seg_indexes,     /* The index of the segment to be split */
                b2_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            handle_next_step(new_g_ptr);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            free(b1_seg_indexes);
            free(b2_seg_indexes);
        }
    }
    }
    }
    }
}


void enum_unbal_transloc(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "enum_unbal_transloc()");
    int hist_idx = 0;

    /* Declare reusable variables */
    int b1, b2;
    struct genome *new_g_ptr, *new_g_ptr2;
    int *b1_seg_indexes, *b2_seg_indexes, b1_seg_indexes_ln, b2_seg_indexes_ln;
    int c1_idx, c2_idx, two_segments_look_identical;
    struct chromosome *seg_holder1, *seg_holder2;
    
    // Go through all chromosomes and all segments
    for (c1_idx=0; c1_idx<g_ptr->n_chrs; c1_idx++) {
    for (c2_idx=c1_idx+1; c2_idx<g_ptr->n_chrs; c2_idx++) {
    for (b1=0; b1 < (*(g_ptr->root_chr+c1_idx))->n_segs; b1++) {
    for (b2=0; b2 < (*(g_ptr->root_chr+c2_idx))->n_segs; b2++) {
        two_segments_look_identical = int_array_cmp(
            (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided,
            (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes,
            (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->times_divided,
            (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->seg_indexes
        );

        // Are the two affected segments the same segment?
        if (two_segments_look_identical) {
            b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));


            /* Option 1: looking at plus strand, b1 < b2 at the segment.
               Segment in question is broken into following:
                  b1  b2
                A ^ B ^ C
               === === ===>
            */

            // Case 1A, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(  // Get the q-telomeric pieces
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1
            );
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(  // Remove the translocated piece from between
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            // Either c1_idx or c2_idx gets lost
            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);
            

            // Case 1B: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2)
            );
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c1_idx),
                seg_holder2,
                (*(new_g_ptr->root_chr+c1_idx))->n_segs
            );
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c2_idx),
                seg_holder1,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+2 : b2+1)
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+1 : b1+2),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2)
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            /* Option 2: looking at plus strand, b1 > b2 at the segment.
               Segment in question is broken into following:
                  b2  b1
                A ^ B ^ C
               === === ===>
            */

            // Case 2A, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(  // Get the q-telomeric pieces
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1
            );
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(  // Remove the translocated piece from between
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2),
                (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);
            

            // Case 2B: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 3);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 3);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1
            );
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2 : b2+1)
            );
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c1_idx),
                seg_holder2,
                (*(new_g_ptr->root_chr+c1_idx))->n_segs
            );
            insert_segs_into_chr(
                *(new_g_ptr->root_chr+c2_idx),
                seg_holder1,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2+1 : b2+2)
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c1_idx),
                ( (*((*(new_g_ptr->root_chr+c1_idx))->root_seg+b1))->is_plus ? b1+2 : b1+1),
                (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs
            );
            delete_segs_from_chr(
                *(new_g_ptr->root_chr+c2_idx),
                0,
                ( (*((*(new_g_ptr->root_chr+c2_idx))->root_seg+b2))->is_plus ? b2 : b2+1)
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                3                   /* Split the segments into three */
            );

            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            free(b1_seg_indexes);
        }
        else {
            // We are here because the two affected segments are not the same

            b1_seg_indexes_ln = (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->times_divided;
            b1_seg_indexes = malloc(b1_seg_indexes_ln * sizeof(int));
            if (b1_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b1_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b1_seg_indexes, (*((*(g_ptr->root_chr+c1_idx))->root_seg+b1))->seg_indexes, b1_seg_indexes_ln * sizeof(int));

            b2_seg_indexes_ln = (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->times_divided;
            b2_seg_indexes = malloc(b2_seg_indexes_ln * sizeof(int));
            if (b2_seg_indexes == NULL) {
                fprintf(stderr, "Failed malloc for b2_seg_indexes. Exiting\n");
                exit(1);
            }
            memcpy(b2_seg_indexes, (*((*(g_ptr->root_chr+c2_idx))->root_seg+b2))->seg_indexes, b2_seg_indexes_ln * sizeof(int));


            // First case, two +- rearrangements
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 2);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 2);

            // Swap the q-telomeric pieces
            seg_holder1 = yank_segments(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1);  // Get the q-telomeric pieces
            seg_holder2 = yank_segments(*(new_g_ptr->root_chr+c2_idx), b2+1, (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);  // Insert the pieces to the ends of the chromosomes
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, (*(new_g_ptr->root_chr+c2_idx))->n_segs);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs);  // Remove the translocated piece from between
            delete_segs_from_chr(*(new_g_ptr->root_chr+c2_idx), b2+1, (*(new_g_ptr->root_chr+c2_idx))->n_segs - 1 - seg_holder1->n_segs);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b2_seg_indexes,     /* The index of the segment to be split */
                b2_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            // Second case: a ++ and a -- rearrangement
            new_g_ptr = copy_genome(g_ptr);
            splice_one_seg(*(new_g_ptr->root_chr+c1_idx), b1, 2);
            splice_one_seg(*(new_g_ptr->root_chr+c2_idx), b2, 2);

            // c1_idx gets the p-telomeric pieces, c2_idx gets the q-telomeric pieces
            seg_holder1 = yank_segments(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1);  // Get the q-telomeric pieces
            invert_segs_in_chr(seg_holder1, 0, seg_holder1->n_segs-1);
            seg_holder2 = yank_segments(*(new_g_ptr->root_chr+c2_idx), 0, b2);
            invert_segs_in_chr(seg_holder2, 0, seg_holder2->n_segs-1);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c1_idx), seg_holder2, (*(new_g_ptr->root_chr+c1_idx))->n_segs);
            insert_segs_into_chr(*(new_g_ptr->root_chr+c2_idx), seg_holder1, b2+1);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c1_idx), b1+1, (*(new_g_ptr->root_chr+c1_idx))->n_segs - 1 - seg_holder2->n_segs);
            delete_segs_from_chr(*(new_g_ptr->root_chr+c2_idx), 0, b2);
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b1_seg_indexes,     /* The index of the segment to be split */
                b1_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );
            splice_all_segs(
                new_g_ptr,          /* Genome where the segments are to be split */
                b2_seg_indexes,     /* The index of the segment to be split */
                b2_seg_indexes_ln,  /* How many times this segment has been split */
                2                   /* Split the segments into three */
            );

            new_g_ptr2 = copy_genome(new_g_ptr);
            make_history(new_g_ptr, UNBAL_TRANSLOC, hist_idx++);
            make_history(new_g_ptr2, UNBAL_TRANSLOC, hist_idx++);
            lose_chromosome_in_genome(new_g_ptr, c1_idx);
            lose_chromosome_in_genome(new_g_ptr2, c2_idx);
            handle_next_step(new_g_ptr);
            handle_next_step(new_g_ptr2);

            delete_chromosome(seg_holder1);
            delete_chromosome(seg_holder2);


            free(b1_seg_indexes);
            free(b2_seg_indexes);
        }
    }
    }
    }
    }
}


void enum_wc_dup(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "enum_wc_dup()");
    int hist_idx = 0;

    /* Declare reusable variables */
    struct genome *new_g_ptr;
    int c_idx;

    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        new_g_ptr = copy_genome(g_ptr);
        make_history(new_g_ptr, WC_DUP, hist_idx++);
        new_g_ptr->n_chrs += 1;
        new_g_ptr->root_chr = realloc(new_g_ptr->root_chr, new_g_ptr->n_chrs * (sizeof(struct chromosome*)));
        if (new_g_ptr->root_chr == NULL) {
            fprintf(stderr, "\nrealloc of new_g_ptr->root_chr failed. Exiting.\n");
            exit(1);
        }
        *(new_g_ptr->root_chr + new_g_ptr->n_chrs - 1) = copy_chromosome(*(new_g_ptr->root_chr + c_idx));
        handle_next_step(new_g_ptr);
    }

    return;
}


void enum_wg_dup(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "enum_wg_dup()");
    struct genome *new_g_ptr;
    new_g_ptr = copy_genome(g_ptr);
    make_history(new_g_ptr, WG_DUP, 0);
    int c_idx, n_chrs;
    n_chrs = new_g_ptr->n_chrs;
    new_g_ptr->root_chr = realloc(new_g_ptr->root_chr, 2 * n_chrs * (sizeof(struct chromosome*)));
    if (new_g_ptr->root_chr == NULL) {
        fprintf(stderr, "\nrealloc of new_g_ptr->root_chr failed. Exiting.\n");
        exit(1);
    }
    for (c_idx=0; c_idx<n_chrs; c_idx++) {
        *(new_g_ptr->root_chr + n_chrs + c_idx) = copy_chromosome(*(new_g_ptr->root_chr + c_idx));
    }
    new_g_ptr->n_chrs = 2 * n_chrs;
    handle_next_step(new_g_ptr);
    return;
}


void enum_wc_del(struct genome* g_ptr) {
    // _validate_genome(g_ptr, "enum_wc_del()");
    int hist_idx = 0;

    /* Declare reusable variables */
    struct genome *new_g_ptr;
    int c_idx, i;

    for (c_idx=0; c_idx<g_ptr->n_chrs; c_idx++) {
        new_g_ptr = copy_genome(g_ptr);
        make_history(new_g_ptr, WC_DEL, hist_idx++);
        delete_chromosome(*(new_g_ptr->root_chr + c_idx));

        for (i=c_idx; i<=new_g_ptr->n_chrs-2; i++) {
            *(new_g_ptr->root_chr+i) = *(new_g_ptr->root_chr+i+1);
        }

        new_g_ptr->n_chrs -= 1;
        new_g_ptr->root_chr = realloc(new_g_ptr->root_chr, new_g_ptr->n_chrs * (sizeof(struct chromosome*)));
        if (new_g_ptr->root_chr == NULL) {
            fprintf(stderr, "\nrealloc of new_g_ptr->root_chr failed. Exiting.\n");
            exit(1);
        }
        handle_next_step(new_g_ptr);
    }

    return;
}

