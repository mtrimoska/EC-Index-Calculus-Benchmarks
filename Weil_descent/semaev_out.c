//
//  semaev_out.c
//  weil
//
//  Created by Monika Trimoska on 27/09/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//

#include "semaev_out.h"
#include "semaev_masks.h"
#include "semaev.h"
#include "params.h"
#include "constants.h"


void init_vars_grobner_out(FILE *out)
{
	int i, j;
	fprintf(out, "R<");
	
	//X coefs
	fprintf(out, "x_1_0");
	for(j = 1; j < l; j++)
	{
		fprintf(out,",x_1_%d", j);
	}
	for(i = 2; i <= m_vars; i++)
	{
		for(j = 0; j < l; j++)
		{
			fprintf(out, ",x_%d_%d", i, j);
		}
	}
	
	//e coefs
	for(i = 1; i <= m_vars; i++)
	{
		for(j = 0; j < L_e[i]; j++)
		{
			fprintf(out, ",e_%d_%d", i, j);
		}
	}
	fprintf(out, "> := BooleanPolynomialRing(%d, \"grevlex\");\n", P_offset[2][0] + m_vars * l);
}

void init_to_unsym_out(FILE *out)
{
	int i, j;
	fprintf(out, "Fpk.<t> = GF(2^%d)\n", n);
	fprintf(out, "var(e_1_0");
	for(j = 1; j < L_e[1]; j++)
	{
			fprintf(out, " e_%d_%d", 1, j);
	}
	for(i = 2; i <= m_vars; i++)
	{
		for(j = 0; j < L_e[i]; j++)
		{
			fprintf(out, " e_%d_%d", i, j);
		}
	}
	fprintf(out, ")\n");
	fprintf(out, "f3.<x_1_0");
	for(j = 1; j < l; j++)
	{
		fprintf(out,",x_1_%d", j);
	}
	for(i = 2; i <= m_vars; i++)
	{
		for(j = 0; j < l; j++)
		{
			fprintf(out, ",x_%d_%d", i, j);
		}
	}
	fprintf(out, ">=PolynomialRing(Fpk,order='invlex')\n");
	
}

void e_X_correspondence_grobner_out(FILE *out)
{
	int i, j, d, k;
	
	// e1
	i = 1;
	fprintf(out, "B := [ e_1_0");
	d = 0;
	for(k = 0; k < T; k++)
	{
		if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
		{
			i = num_to_indices(k);
			fprintf(out, " + x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
			for(j = 1; j < i; j++)
				fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
		}
	}
	for(d = 1; d < L_e[i]; d++)
	{
		fprintf(out, "\n, e_%d_%d", i, d);
		for(k = 0; k < T; k++)
		{
			if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
			{
				i = num_to_indices(k);
				fprintf(out, " + x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
				for(j = 1; j < i; j++)
					fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
			}
		}
	}
	
	//ei
	for(i = 2; i <= m_vars; i++)
	{
		for(d = 0; d < L_e[i]; d++)
		{
			fprintf(out, "\n, e_%d_%d", i, d);
			for(k = 0; k < T; k++)
			{
				if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
				{
					i = num_to_indices(k);
					fprintf(out, " + x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
					for(j = 1; j < i; j++)
						fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
				}
			}
		}
	}
}

void e_X_correspondence_to_unsym_out(FILE *out)
{
	int i, j, d, k, first;
	
	// e1
	i = 1;
	fprintf(out, "for eq in list_eq:\n\tunsym=eq.subs(e_1_0 = ");
	d = 0;
	for(k = 0; k < T; k++)
	{
		first = 1;
		if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
		{
			i = num_to_indices(k);
			if(!first)
				fprintf(out, " + ");
			fprintf(out, "x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
			for(j = 1; j < i; j++)
				fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
			first = 0;
		}
	}
	for(d = 1; d < L_e[i]; d++)
	{
		fprintf(out, "\n, e_%d_%d", i, d);
		for(k = 0; k < T; k++)
		{
			first = 1;
			if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
			{
				i = num_to_indices(k);
				if(!first)
					fprintf(out, " + ");
				fprintf(out, "x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
				for(j = 1; j < i; j++)
					fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
				first = 0;
			}
		}
	}
	
	//ei
	for(i = 2; i <= m_vars; i++)
	{
		for(d = 0; d < L_e[i]; d++)
		{
			fprintf(out, "\n, e_%d_%d", i, d);
			for(k = 0; k < T; k++)
			{
				first = 1;
				if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
				{
					i = num_to_indices(k);
					if(!first)
						fprintf(out, " + ");
					fprintf(out, "x_%d_%d", indices_buffer[0][0], indices_buffer[0][1]);
					for(j = 1; j < i; j++)
						fprintf(out, "*x_%d_%d", indices_buffer[j][0], indices_buffer[j][1]);
					first = 0;
				}
			}
		}
	}
	fprintf(out, ").expand()\n\tstr_sys=str(unsym)\n\tprint \", \", str_sys");
}

int e_X_correspondence_dimacs_out(FILE *out)
{
	int i, d, k, cpt, nb_lines;
	cpt = 0;
	nb_lines = 0;
	for(i = 1; i <= m_vars; i++)
	{
		for(d = 0; d < L_e[i]; d++)
		{
			fprintf(out, "\nx %d ", -(e_offset[0] + cpt));
			for(k = 0; k < T; k++)
			{
				if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
				{
					fprintf(out, "%d ", k + X_offset[k]);
				}
			}
			fprintf(out, "0");
			cpt++;
			nb_lines++;
		}
	}
	return nb_lines;
}

int e_X_correspondence_xorDand_out(FILE *out)
{
	int i, j, d, k, cpt, nb_lines;
	int factors[__MAX_M__] = {0};
	cpt = 1;
	nb_lines = 0;
	//printf("%d,%d,%d,%d\n", L_e[0], L_e[1], L_e[2], L_e[3]);
	for(i = 1; i <= m_vars; i++)
	{
		for(d = 0; d < L_e[i]; d++)
		{
			fprintf(out, "x T %d ", m_vars*l + cpt);
			for(k = 0; k < T; k++)
			{
				if(vect_bin_get_bit(e_X_correspondence[i], d * T + k))
				{
					if(k < m_vars * l)
					{
						fprintf(out, "%d ", k + X_offset[k]);
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
							fprintf(out, "%d ", factors[j] + X_offset[factors[j]]);
							nb_lines++;
						}
					}
					
				}
			}
			fprintf(out, "0\n");
			cpt++;
			nb_lines++;
		}
	}
	return nb_lines;
}

int create_cnf_dimacs_out(FILE *out, int *offsets)
{
	int k, i, j, nb_lines;
	int factors[__MAX_M__] = {0};
	nb_lines = 0;
	
	// X vars
	for(k = P_offset[2][0]; k < T; k++)
	{
		if(offsets[k] != (__MAX_VARS_X__ + __MAX_VARS_E__))
		{
			i = num_to_indices(k);
			for(j = 0; j < i; j++)
			{
				indices_buffer[0][0] = indices_buffer[j][0];
				indices_buffer[0][1] = indices_buffer[j][1];
				factors[j] = indices_to_num(1);
				fprintf(out, "\n%d %d 0", -(k + offsets[k]), factors[j] + offsets[factors[j]]);
				nb_lines++;
			}
			
			fprintf(out, "\n%d ", k + offsets[k]);
			for(j = 0; j < i; j++)
			{
				fprintf(out, "%d ", -(factors[j] + offsets[factors[j]]));
			}
			fprintf(out, "0");
			nb_lines++;
		}
	}
	return nb_lines;
}
