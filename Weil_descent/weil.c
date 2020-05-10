#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdint.h>
#include "vect_bin.h"
#include "params.h"
#include "weil.h"
#include "semaev.h"
#include "terms.h"

int check_solution(_vect_bin_t *v, char *solution)
{
	int i, j, it, k, sum, t, k_T, one_term_val;
	int values[__MAX_M__][__MAX_E_LEN__] = {0};
	it = 0;
	for(i = 1; i <= m_vars; i++)
	{
		for(j = 0; j < L_e[i]; j++)
		{
			values[i][j] = solution[it++] - 48;
			printf("values[%d][%d] = %d\n", i, j , values[i][j]);
		}
	}
	for(k = 0; k < n; k++)
	{
		sum = 0;
		k_T = k * T;
		for(t = 0; t < T; t++)
		{
			if(vect_bin_get_bit(v, k_T + t))
			{
				i = num_to_indices(t);
				one_term_val = 1;
				for(it = 0; it < i; it++)
				{
					if(values[indices_buffer[it][0]][indices_buffer[it][1]] == 0)
					{
						one_term_val = 0;
						break;
					}
					//printf("J'ajoute à %d a_%d_%d = %d\n", k, indices_buffer[it][0], indices_buffer[it][1], values[indices_buffer[it][0]][indices_buffer[it][1]]);
				}
				sum += one_term_val;
			}
		}
		printf("deg %d : %d\n\n", k, sum % 2);
	}
	return 1;
}

void multiply_with_const(_vect_bin_t *v, _vect_bin_t *term, char *Xr, int power)
{
	int cpt, shift, p;
	init(cpy_term);
	init(cpy_term_impair);
	
	if(power % 2 == 0)
	{
		p = power;
		vect_bin_or(cpy_term, term);
	}
	else
	{
		p = power - 1;
		vect_bin_or(cpy_term_impair, term);
		shift = -1;
		for(cpt = 0; cpt < n; cpt++)
		{
			shift++;
			if(Xr[cpt] == 49)
			{
				vect_bin_left_shift(cpy_term_impair, shift * T);
				vect_bin_xor(cpy_term, cpy_term_impair);
				shift = 0;
			}
		}
	}
	if(power == 1 || power == 0)
	{
		vect_bin_xor(v, cpy_term);
	}
	else
	{
		shift = -1;
		for(cpt = 0; cpt < n; cpt++)
		{
			shift++;
			if(Xr[cpt] == 49)
			{
				vect_bin_left_shift(cpy_term, p * shift * T);
				vect_bin_xor(v, cpy_term);
				shift = 0;
			}
		}
	}
}


/// do v modulo irreductible (Euclidean)
void modulo(_vect_bin_t *v, char *irr)
{
	int dividendD, shift, i;
	init(reduce);
	init(coef);
	dividendD = degree(v);
	while(dividendD > n - 1)
	{
		create_mask_coef(coef, dividendD);
		vect_bin_and(coef, v);
		shift = 0;
		for(i = n; i >= 0; i--)
		{
			if(irr[i] == 49)
			{
				vect_bin_right_shift(coef, shift * T);
				vect_bin_or(reduce, coef);
				shift = 0;
			}
			shift++;
		}
		vect_bin_xor(v, reduce);
		vect_bin_t_reset(reduce);
		vect_bin_t_reset(coef);
		dividendD = degree(v);
	}
}

void compute_constant(int *constant, char *Xr, char *irr)
{
	///n'est pas générique, juste pour 4ème pol.
	int d;
	init(const_vect)
	for(d = 0; d < n; d++)
	{
		if(Xr[d] == 49)
		{
			vect_bin_set_1(const_vect, d * 4 * T);
		}
	}
	modulo(const_vect, irr);
	for(d = 0; d < n; d++)
	{
		if(vect_bin_get_bit(const_vect, d * T))
		{
			constant[d] = 1;
		}
	}
}

///Compute Weil descent
void create_weil(_vect_bin_t *f, int *constant, char *irr, char *Xr, char *solution, int nb_terms)
{
	int i;
	for(i = 0; i < nb_terms; i++)
	{
		multiply_with_const(f, term_mask[i], Xr, term_const_power[i]);
	}
	modulo(f, irr);
	//print_coefs(f);
	//check_solution(f, solution);
	compute_constant(constant, Xr, irr);
}

