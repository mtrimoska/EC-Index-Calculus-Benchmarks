//
//  weil_out.c
//  weil
//
//  Created by Monika Trimoska on 24/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//
#include "params.h"
#include "semaev.h"
#include "weil_out.h"
#include "e_offset_weil.h"

void print_coefs(_vect_bin_t *v)
{
	int d, k, i, j;
	for(d = _n - 1; d >= 0; d--)
	{
		printf("e^%d: ", d);
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(v, d * T + k))
			{
				i = num_to_indices(k);
				for(j = 0; j < i; j++)
					printf("a[%d,%d]", indices_buffer[j][0], indices_buffer[j][1]);
				printf(" + ");
			}
		}
		printf("\n\n");
	}
}

void system_grobner_out(FILE *out, _vect_bin_t *f, int *constant)
{
	int i, j, d, k;
	for(d = 0; d < n; d++)
	{
		fprintf(out, "\n, %d", constant[d]);
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(f, d * T + k))
			{
				i = num_to_indices(k);
				fprintf(out, " + e_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
				for(j = 1; j < i; j++)
					fprintf(out, "*e_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
			}
		}
	}
	fprintf(out, "];\nI := Ideal(B);\nVariety(I);\n");
}

void system_to_unsym_out(FILE *out, _vect_bin_t *f, int *constant)
{
	int i, j, d, k;
	fprintf(out, "list_eq=[");
	d = 0;
	fprintf(out, "%d", constant[d]);
	for(k = 0; k < T; k++)
	{
		if(vect_bin_get_bit(f, d * T + k))
		{
			i = num_to_indices(k);
			fprintf(out, " + e_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
			for(j = 1; j < i; j++)
				fprintf(out, "*e_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
		}
	}
	for(d = 1; d < n; d++)
	{
		fprintf(out, "\n, %d", constant[d]);
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(f, d * T + k))
			{
				i = num_to_indices(k);
				fprintf(out, " + e_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
				for(j = 1; j < i; j++)
					fprintf(out, "*e_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
			}
		}
	}
	fprintf(out, "]\n");
}

int system_dimacs_out(FILE *out, _vect_bin_t *f, int *constant)
{
	int d, k;
	for(d = 0; d < n; d++)
	{
		fprintf(out, "\nx ");
		if(constant[d] == 0)
			fprintf(out, "-");
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(f, d * T + k))
			{
				fprintf(out, "%d ", k + e_offset_weil[k]);
			}
		}
		fprintf(out, "0");
	}
	return n;
}

int system_xorDand_out(FILE *out, _vect_bin_t *f, int *constant)
{
	int d, k, i, j;
	int factors[__MAX_M__] = {0};
	for(d = 0; d < n; d++)
	{
		fprintf(out, "x ");
		if(constant[d] == 0)
			fprintf(out, "T ");
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(f, d * T + k))
			{
				if(k < P_offset[2][0])
				{
					fprintf(out, "%d ", k + 3*l + 1);
				}
				else
				{
					i = num_to_indices(k);
					if(i == 2) fprintf(out, ".2 ");
					if(i == 3) fprintf(out, ".3 ");
					for(j = 0; j < i; j++)
					{
						indices_buffer[0][0] = indices_buffer[j][0];
						indices_buffer[0][1] = indices_buffer[j][1];
						factors[j] = indices_to_num(1);
						fprintf(out, "%d ", factors[j] + 3*l + 1);
					}
				}
			}
		}
		fprintf(out, "0\n");
	}
	return n;
}

void header_dimacs_out(FILE *out, int nb_lines)
{
	int last = T - 1;
	while(e_offset_weil[last] == (__MAX_VARS_X__ + __MAX_VARS_E__))
	{
		last--;
	}
	fprintf(out, "p cnf %d %d", last + e_offset_weil[last], nb_lines);
}

void header_xorDand_out(FILE *out)
{
	int i, nb_vars, nb_lines;
	nb_vars = m_vars * l;
	nb_lines = 0;
	for(i = 1; i <= m_vars; i++)
	{
		nb_vars += L_e[i];
		nb_lines += L_e[i];
	}
	nb_lines += n;
	fprintf(out, "p cnf %d %d\n", nb_vars, nb_lines);
}

