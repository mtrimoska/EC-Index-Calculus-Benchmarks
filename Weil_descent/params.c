#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include "vect_bin.h"
#include "constants.h"
#include "params.h"

int n;
int l;
int m;
int _n;
int m_vars;
int m_deg;
int T;
int L_e[__MAX_M__];
int P_offset[__MAX_M__][__MAX_C__];
int indices_buffer[__MAX_M__][2];

/* Manipulating coefs */

/// Compute number of combinations without repetition
int comb_norep(int n, int k)
{
	int res;
	int i;
	int max;
	if(n - k > k)
	{
		max = n - k;
	}
	else
	{
		max = k;
		k = n - k;
	}
	res = 1;
	for(i = n; i > max; i--)
	{
		res = res * i;
	}
	for(i = k; i > 1; i--)
	{
		res = res / i;
	}
	return res;
}

/// Compute number of combinations with repetition
int comb_rep(int n, int k)
{
	int res;
	int i;
	int max;
	int top = n + k - 1;
	if(n - 1 > k)
	{
		max = n - 1;
	}
	else
	{
		max = k;
		k = n - 1;
	}
	res = 1;
	for(i = top; i > max; i--)
	{
		res = res * i;
	}
	for(i = k; i > 1; i--)
	{
		res = res / i;
	}
	return res;
}

int indices_to_num(int i)
{
	int res = 0;
	int comb_num = 0;
	int it, p, it2;
	comb_num = comb_num + (comb_rep(m_vars, i)
				 - comb_rep(m_vars - indices_buffer[0][0] + 1, i));
	
	for(it = 1; it < i; it++)
	{
		comb_num = comb_num + (comb_rep(m_vars - indices_buffer[it-1][0] + 1, i - it)
					 - comb_rep(m_vars - indices_buffer[it][0] + 1, i - it));
	}
	res = P_offset[i][comb_num];
	
	for(it = 0; it < i; it++)
	{
		p = 1;
		for(it2 = it + 1; it2 < i; it2++)
		{
			p *= L_e[indices_buffer[it2][0]];
		}
		res = res + indices_buffer[it][1] * p;
	}
	
	return res;
}

int get_num(int i, ...)
{
	int ii;
	va_list valist;
	va_start(valist, i);
	for(ii = 0; ii < i; ii++)
	{
		indices_buffer[ii][0] = va_arg(valist, int);
		indices_buffer[ii][1] = va_arg(valist, int);
	}
	va_end(valist);
	for(ii = 0; ii < i - 1; ii++)
	{
		assert(indices_buffer[ii][0] < indices_buffer[ii + 1][0]);
	}
	return indices_to_num(i);
}

//returns the number of the group of vars
int num_to_indices(int num)
{
	int i, comb_num, perm_num, ii, ii_val, subtract, comb_nb, it2, p, c;
	comb_num = 0;
	//find i and combination number
	for(i = 1; i <= m_deg; i++)
	{
		comb_nb = comb_rep(m_vars, i);
		for(comb_num = 0; comb_num < comb_nb - 1; comb_num++)
		{
			if(num < P_offset[i][comb_num + 1])
			{
				break;
			}
		}
		if(i == m_deg || num < P_offset[i + 1][0])
		{
			break;
		}
	}
	
	//find combination indices i
	ii = 1;
	ii_val = 1;
	c = comb_num;
	while(i - ii > 0)
	{
		subtract = comb_rep(m_vars - ii_val + 1, i - ii);
		while(c - subtract >= 0 && ii_val <= m_deg)
		{
			c = c - subtract;
			ii_val++;
			subtract = comb_rep(m_vars - ii_val + 1, i - ii);
		}
		indices_buffer[ii - 1][0] = ii_val;
		ii++;
	}
	indices_buffer[ii - 1][0] = c + ii_val;
	
	//find permutation indices j
	perm_num = num - P_offset[i][comb_num];
	for(ii = 0; ii < i; ii++)
	{
		p = 1;
		for(it2 = ii + 1; it2 < i; it2++)
		{
			p *= L_e[indices_buffer[it2][0]];
		}
		indices_buffer[ii][1] = perm_num / p;
		perm_num = perm_num % p;
	}
	
	return i;
}

/// Setting all parameters
void set_params(int new_n, int new_m, int new_l, int new__n, int new_deg)
{
	int i, comb_nb, c, comb_num, ii, ii_val, subtract, P;

	n = new_n;
	m = new_m;
	m_vars = m - 1;
	m_deg = new_deg;
	if(new_l == -1) l = n/m_vars; else l = new_l;
	_n = new__n;
	T = 0;
	
	for(i = 1; i <= m_vars; i++)
	{
		L_e[i] = i * l - (i - 1);
	}
	
	for(i = 1; i <= m_deg; i++)
	{
		comb_nb = comb_rep(m_vars, i);
		for(c = 0; c < comb_nb; c++)
		{
			P_offset[i][c] = T;
			comb_num = c;
			P = 1;
			
			//find combination indices i
			ii = 1;
			ii_val = 1;
			while(i - ii > 0)
			{
				subtract = comb_rep(m_vars - ii_val + 1, i - ii);
				while(comb_num - subtract >= 0 && ii_val < m_vars)
				{
					comb_num = comb_num - subtract;
					ii_val++;
					subtract = comb_rep(m_vars - ii_val + 1, i - ii);
				}
				P *= L_e[ii_val];
				ii++;
			}
			P *= L_e[comb_num + ii_val];
			T += P;
		}
	}
	set_vect_size(_n * T);
	//printf("T:%d\n", T);
}



