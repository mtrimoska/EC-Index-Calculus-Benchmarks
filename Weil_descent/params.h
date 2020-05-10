#ifndef params_h
#define params_h

#include <stdio.h>
#include "vect_bin.h"

extern int n;
extern int l;
extern int m;
extern int _n;
extern int m_vars;
extern int m_deg;
extern int T;
extern int L_e[__MAX_M__]; //starts at 1
extern int P_offset[__MAX_M__][__MAX_C__]; //starts at 1
extern int indices_buffer[__MAX_M__][2]; //starts at 0

int comb_norep(int n, int k);
void set_params(int new_n, int new_m, int new_l, int new__n, int new_deg);
int get_num(int i, ...);
int num_to_indices(int num);
int indices_to_num(int i);
#endif /* params_h */
