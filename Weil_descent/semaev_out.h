//
//  semaev_out.h
//  weil
//
//  Created by Monika Trimoska on 27/09/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#ifndef semaev_out_h
#define semaev_out_h

#include <stdio.h>

void init_vars_grobner_out(FILE *out);
void init_to_unsym_out(FILE *out);
void e_X_correspondence_grobner_out(FILE *out);
void e_X_correspondence_to_unsym_out(FILE *out);
int e_X_correspondence_dimacs_out(FILE *out);
int e_X_correspondence_xorDand_out(FILE *out);
int create_cnf_dimacs_out(FILE *out, int *offsets);

#endif /* semaev_out_h */
