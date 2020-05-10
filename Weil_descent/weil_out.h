//
//  weil_out.h
//  weil
//
//  Created by Monika Trimoska on 24/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#ifndef weil_out_h
#define weil_out_h

#include <stdio.h>
#include "vect_bin.h"
#include "params.h"

void print_coefs(_vect_bin_t *v);
void system_grobner_out(FILE *out, _vect_bin_t *f, int *constant);
void system_to_unsym_out(FILE *out, _vect_bin_t *f, int *constant);
int system_dimacs_out(FILE *out, _vect_bin_t *f, int *constant);
int system_xorDand_out(FILE *out, _vect_bin_t *f, int *constant);
void header_dimacs_out(FILE *out, int nb_lines);
void header_xorDand_out(FILE *out);
#endif /* weil_out_h */
