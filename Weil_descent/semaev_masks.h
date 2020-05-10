//
//  semaev_masks.h
//  weil
//
//  Created by Monika Trimoska on 25/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#ifndef semaev_masks_h
#define semaev_masks_h

#include <stdio.h>

#include "vect_bin.h"
#include "params.h"

extern _vect_bin_t X[__MAX_M__][__ARRAY_SIZE__]; //starts at 1
extern _vect_bin_t e[__MAX_M__][__ARRAY_SIZE__]; //starts at 1
extern _vect_bin_t e_X_correspondence[__MAX_M__][__ARRAY_SIZE__]; // starts at 1

int degree(_vect_bin_t *v);
void create_mask_coef(_vect_bin_t *v, int d);
void multiply(_vect_bin_t *product, _vect_bin_t *factor1, _vect_bin_t *factor2);
void power(_vect_bin_t *result, _vect_bin_t *v, int power);
void create_e_masks(void);
void create_X_masks(void);
void create_e_X_correspondence(void);
#endif /* semaev_masks_h */
