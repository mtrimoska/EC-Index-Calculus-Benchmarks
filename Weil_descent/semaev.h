//
//  semaev.h
//  weil
//
//  Created by Monika Trimoska on 24/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#ifndef semaev_h
#define semaev_h

#include <stdio.h>
#include "semaev_masks.h"
#include "constants.h"

extern int X_offset[__MAX_VARS_X__];
extern int e_offset[__MAX_VARS_E__];

void create_semaev(void);
void compute_offsets_X(void);
void compute_offsets_e(int nb_terms);
#endif /* semaev_h */
