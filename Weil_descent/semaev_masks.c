//
//  semaev_masks.c
//  weil
//
//  Created by Monika Trimoska on 25/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#include <string.h>
#include <assert.h>
#include "semaev_masks.h"
#include "vect_bin.h"
#include "params.h"
#include "constants.h"

_vect_bin_t X[__MAX_M__][__ARRAY_SIZE__];
_vect_bin_t e[__MAX_M__][__ARRAY_SIZE__];
_vect_bin_t e_X_correspondence[__MAX_M__][__ARRAY_SIZE__];

/// find highest degree in vector
int degree(_vect_bin_t *v)
{
	int i;
	for(i = _vect_bin_size; i > 0; i--)
	{
		if(vect_bin_get_bit(v, i))
		{
			return (i/T);
		}
	}
	return i;
}

/// create mask to extract coefficient on degree d
void create_mask_coef(_vect_bin_t *v, int d)
{
	int i;
	for(i = d * T; i < (d+1) * T; i++)
	{
		vect_bin_set_1(v, i);
	}
}

_bool_t is_zero_degree(_vect_bin_t *v, int d)
{
	int i;
	for(i = d * T; i < (d+1) * T; i++)
	{
		if(vect_bin_get_bit(v, i))
		{
			return _false;
		}
	}
	return _true;
}

void create_X_masks()
{
	int i, d;
	for(i = 1; i <= m_vars; i++)
	{
		vect_bin_t_reset(X[i]);
		for(d = 0; d < l; d++)
		{
			vect_bin_set_1(X[i], d * T + get_num(1, i, d));
		}
	}
}

void create_e_masks()
{
	int i, d;
	for(i = 1; i <= m_vars; i++)
	{
		vect_bin_t_reset(e[i]);
		for(d = 0; d < L_e[i]; d++)
		{
			vect_bin_set_1(e[i], d * T + get_num(1, i, d));
		}
	}
}

void create_e_X_correspondence()
{
	int i, comb_nb, comb_num, c, ii, ii_val, subtract;
	init(combination);
	init(temp_multiply);
	for(i = 1; i <= m_vars; i++)
	{
		vect_bin_t_reset(e_X_correspondence[i]);
		comb_nb = comb_norep(m_vars, i);
		for(comb_num = 0; comb_num < comb_nb; comb_num++)
		{
			ii = 1;
			ii_val = 1;
			c = comb_num;
			vect_bin_t_reset(combination);
			while(i - ii > 0)
			{
				subtract = comb_norep(m_vars - ii_val, i - ii);
				while(c - subtract >= 0 && ii_val < m_deg)
				{
					c = c - subtract;
					ii_val++;
					subtract = comb_norep(m_vars - ii_val, i - ii);
				}
				if(ii == 1)
				{
					vect_bin_or(combination, X[ii_val]);
				}
				else
				{
					multiply(temp_multiply, combination, X[ii_val]);
					vect_bin_t_reset(combination);
					vect_bin_or(combination, temp_multiply);
					vect_bin_t_reset(temp_multiply);
				}
				ii++;
				ii_val++;
			}
			ii_val += c;
			if(ii == 1)
			{
				vect_bin_or(combination, X[ii_val]);
			}
			else
			{
				multiply(temp_multiply, combination, X[ii_val]);
				vect_bin_t_reset(combination);
				vect_bin_or(combination, temp_multiply);
				vect_bin_t_reset(temp_multiply);
			}
			vect_bin_xor(e_X_correspondence[i], combination);
		}
	}
}

/// multiply two polynomial vectors
void multiply(_vect_bin_t *product, _vect_bin_t *factor1, _vect_bin_t *factor2)
{
	int d1, d2, i1, i2, c1, c2, it, i, num_c1, num_c2, it1, it2;
	int temp_indices1[__MAX_M__][2] = {0};
	int temp_indices2[__MAX_M__][2] = {0};
	d1 = -1;
	d2 = -1;
	for(c1 = 0; c1 < _n * T; c1++)
	{
		if(c1 % T == 0) d1++;
		if(vect_bin_get_bit(factor1, c1))
		{
			num_c1 = c1 - (d1 * T);
			i1 = num_to_indices(num_c1);
			memcpy(temp_indices1, indices_buffer, sizeof(int) * i1 * 2);
			for(c2 = 0; c2 < _n * T; c2++)
			{
				if(c2 % T == 0) d2++;
				if(vect_bin_get_bit(factor2, c2))
				{
					num_c2 = c2 - (d2 * T);
					i2 = num_to_indices(num_c2);
					memcpy(temp_indices2, indices_buffer, sizeof(int) * i2 * 2);
					i = i1 + i2;
					it1 = 0;
					it2 = 0;
					for(it = 0; it < i; it++)
					{
						if((it2 >= i2) || ((it1 < i1) && (temp_indices1[it1][0] < temp_indices2[it2][0])))
						{
							indices_buffer[it][0] = temp_indices1[it1][0];
							indices_buffer[it][1] = temp_indices1[it1][1];
							it1++;
						}
						else
						{
							if((it2 < i2) && (it1 < i1) && (temp_indices1[it1][0] == temp_indices2[it2][0]))
							{
								if(temp_indices1[it1][1] <= temp_indices2[it2][1])
								{
									if(temp_indices1[it1][1] == temp_indices2[it2][1]) //car a^2 == a
									{
										it2++;
										i--;
									}
									indices_buffer[it][0] = temp_indices1[it1][0];
									indices_buffer[it][1] = temp_indices1[it1][1];
									it1++;
								}
								else
								{
									
									indices_buffer[it][0] = temp_indices2[it2][0];
									indices_buffer[it][1] = temp_indices2[it2][1];
									it2++;
								}
							}
							else
							{
								indices_buffer[it][0] = temp_indices2[it2][0];
								indices_buffer[it][1] = temp_indices2[it2][1];
								it2++;
							}
						}
						/*if((indices_buffer[it][0] == indices_buffer[it - 1][0]) && (indices_buffer[it][1] == indices_buffer[it - 1][1])) //car a^2 == a
						{
							it--;
							i--;
						}*/
					}
					vect_bin_set_1(product, (d1 + d2) * T + indices_to_num(i));
				}
			}
			d2 = -1;
		}
	}
	
}

/// raise a polynomial vector to a power
void power(_vect_bin_t *result, _vect_bin_t *v, int power)
{
	if(power == 1)
	{
		vect_bin_or(result, v);
	}
	else
	{
		int d, i, mask_d, p;
		init(cpy);
		init(result_inter);
		init(mask);
		init(add);
		d = degree(v);
		if(power % 2 == 0)
		{
			p = power;
		}
		else
		{
			p = power - 1;
		}
		mask_d = 0;
		vect_bin_or(cpy, v);
		create_mask_coef(mask, mask_d);
		vect_bin_or(add, cpy);
		vect_bin_and(add, mask);
		vect_bin_or(result_inter, add);
		for(i = 1; i <= d; i++)
		{
			vect_bin_t_reset(mask);
			vect_bin_t_reset(add);
			vect_bin_left_shift(cpy, (p - 1) * T);
			mask_d = mask_d + p;
			create_mask_coef(mask, mask_d);
			vect_bin_or(add, cpy);
			vect_bin_and(add, mask);
			vect_bin_or(result_inter, add);
		}
		if(power % 2 == 0)
		{
			vect_bin_or(result, result_inter);
		}
		else
		{
			multiply(result, result_inter, v);
		}
	}
}
