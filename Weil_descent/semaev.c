//
//  semaev.c
//  weil
//
//  Created by Monika Trimoska on 24/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//
#include "semaev.h"
#include "weil_out.h"
#include "offset.h"

_vect_bin_t term_mask_temp[__MAX_TERMS__][__ARRAY_SIZE__];
int term_const_power_temp[__MAX_TERMS__];
int X_offset[__MAX_VARS_X__];
int e_offset[__MAX_VARS_E__];

void create_semaev()
{
	int i, j;
	FILE *terms_file;
	init(temp1);
	init(temp2);

	vect_bin_t_reset(term_mask_temp[0]);
	power(term_mask_temp[0], e[1], 4);
	term_const_power_temp[0] = 0;
	
	vect_bin_t_reset(term_mask_temp[1]);
	//exception
	term_const_power_temp[1] = 4;
	
	vect_bin_t_reset(term_mask_temp[2]);
	power(term_mask_temp[2], e[3], 4);
	term_const_power_temp[2] = 0;
	
	vect_bin_t_reset(term_mask_temp[3]);
	power(term_mask_temp[3], e[2], 4);
	term_const_power_temp[3] = 4;
	
	vect_bin_t_reset(term_mask_temp[4]);
	power(term_mask_temp[4], e[3], 3);
	term_const_power_temp[4] = 1;
	
	vect_bin_t_reset(term_mask_temp[5]);
	power(temp2, e[2], 2);
	multiply(term_mask_temp[5], e[3], temp2);
	vect_bin_t_reset(temp2);
	term_const_power_temp[5] = 3;
	
	vect_bin_t_reset(term_mask_temp[6]);
	power(temp2, e[1], 2);
	multiply(term_mask_temp[6], e[3], temp2);
	vect_bin_t_reset(temp2);
	term_const_power_temp[6] = 1;
	
	vect_bin_t_reset(term_mask_temp[7]);
	power(term_mask_temp[7], e[3], 1);
	term_const_power_temp[7] = 3;
	
	vect_bin_t_reset(term_mask_temp[8]);
	power(temp1, e[1], 2);
	power(temp2, e[3], 2);
	multiply(term_mask_temp[8], temp1, temp2);
	vect_bin_t_reset(temp1);
	vect_bin_t_reset(temp2);
	term_const_power_temp[8] = 2;
	
	vect_bin_t_reset(term_mask_temp[9]);
	power(term_mask_temp[9], e[3], 2);
	term_const_power_temp[9] = 4;
	
	vect_bin_t_reset(term_mask_temp[10]);
	power(term_mask_temp[10], e[3], 2);
	term_const_power_temp[10] = 0;
	
	vect_bin_t_reset(term_mask_temp[11]);
	power(term_mask_temp[11], e[2], 2);
	term_const_power_temp[11] = 2;
	
	/** pol 3*/
	/*vect_bin_t_reset(term_mask[0]);
	power(term_mask[0], e[2], 2);
	term_const_power[0] = 0;
	
	vect_bin_t_reset(term_mask[1]);
	power(term_mask[1], e[1], 2);
	term_const_power[1] = 2;
	
	vect_bin_t_reset(term_mask[2]);
	power(term_mask[2], e[2], 1);
	term_const_power[2] = 1;*/
	
	/*for(int i_temp = 0; i_temp < __MAX_TERMS__; i_temp++)
	{
		printf("\n\n\nterm %d\n",i_temp);
		print_coefs(term_mask[i_temp]);
	}*/
	
	if((terms_file = fopen("terms.h", "w")) == NULL)
	{
		printf("Error file terms.h\n");
	}
	
	fprintf(terms_file, "_vect_bin_t term_mask[__MAX_TERMS__][__ARRAY_SIZE__] = {{%llu", term_mask_temp[0][0]);
	for(j = 1; j < __ARRAY_SIZE__; j++)
	{
		fprintf(terms_file, ",%llu", term_mask_temp[0][j]);
	}
	fprintf(terms_file, "}");
	
	for(int i = 1; i < __MAX_TERMS__; i++)
	{
		fprintf(terms_file, ",{%llu", term_mask_temp[i][0]);
		for(j = 1; j < __ARRAY_SIZE__; j++)
		{
			fprintf(terms_file, ",%llu", term_mask_temp[i][j]);
		}
		fprintf(terms_file, "}");
	}
	fprintf(terms_file, "};\n");
	
	fprintf(terms_file, "int term_const_power[__MAX_TERMS__] = {%d", term_const_power_temp[0]);
	for(int i = 1; i < __MAX_TERMS__; i++)
	{
		fprintf(terms_file, ",%d", term_const_power_temp[i]);
	}
	fprintf(terms_file, "};\n");
	
	fclose(terms_file);
}

void compute_offsets_X()
{
	int i, d, k, missing;
	FILE *offset_file;
	init(occurence)
	init(cpy)
	
	for(i = 1; i <= m_vars; i++)
	{
		vect_bin_or(cpy, e_X_correspondence[i]);
		vect_bin_or(occurence, cpy);
		for(d = 1; d < L_e[i]; d++)
		{
			vect_bin_right_shift(cpy, T);
			vect_bin_or(occurence, cpy);
		}
		vect_bin_t_reset(cpy);
	}
	
	missing = -1; // dimacs vars start at 1
 	for(k = 0; k < T; k++)
	{
		if(vect_bin_get_bit(occurence, k))
		{
			X_offset[k] = (- missing);
		}
		else
		{
			X_offset[k] = (__MAX_VARS_X__ + __MAX_VARS_E__);
			missing++;
		}
	}
	
	e_offset[0] = T - missing;
	if((offset_file = fopen("offset.h", "w")) == NULL)
	{
		printf("Error file terms.h\n");
	}
	fprintf(offset_file, "int new_offset = %d;\n", e_offset[0]);
	fclose(offset_file);

	/*int ii;
	for(int cpt = 0; cpt < T; cpt++)
	{
		ii = num_to_indices(cpt);
		printf("a[%d,%d] ", indices_buffer[0][0], indices_buffer[0][1]);
		if(ii > 1)
		{
			printf("a[%d,%d] ", indices_buffer[1][0], indices_buffer[1][1]);
		}
		if(ii > 2)
		{
			printf("a[%d,%d] ", indices_buffer[2][0], indices_buffer[2][1]);
		}
		printf("--Exists:%d---%d -> %d -> %d\n", vect_bin_get_bit(occurence, cpt), cpt, indices_to_num(ii), X_offset[cpt]);
	}*/
}

void compute_offsets_e(int nb_terms)
{
	int i, d, k, missing;
	FILE *offset_file;
	init(occurence)
	init(cpy)
	
	for(i = 0; i < nb_terms; i++)
	{
		vect_bin_or(cpy, term_mask_temp[i]);
		vect_bin_or(occurence, cpy);
		for(d = 1; d < _n; d++)
		{
			vect_bin_right_shift(cpy, T);
			vect_bin_or(occurence, cpy);
		}
		vect_bin_t_reset(cpy);
	}
	
	missing = - new_offset;
	for(k = 0; k < T; k++)
	{
		if(vect_bin_get_bit(occurence, k))
		{
			e_offset[k] = - missing;
		}
		else
		{
			e_offset[k] = (__MAX_VARS_X__ + __MAX_VARS_E__);
			missing++;
		}
	}
	
	/*int ii;
	for(int cpt = 0; cpt < T; cpt++)
	{
		ii = num_to_indices(cpt);
		printf("a[%d,%d] ", indices_buffer[0][0], indices_buffer[0][1]);
		if(ii > 1)
		{
			printf("a[%d,%d] ", indices_buffer[1][0], indices_buffer[1][1]);
		}
		if(ii > 2)
		{
			printf("a[%d,%d] ", indices_buffer[2][0], indices_buffer[2][1]);
		}
		printf("--Exists:%d---%d -> %d -> %d\n", vect_bin_get_bit(occurence, cpt), cpt, indices_to_num(ii), e_offset[cpt]);
	}*/
	if((offset_file = fopen("e_offset_weil.h", "w")) == NULL)
	{
		printf("Error file terms.h\n");
	}
	fprintf(offset_file, "int e_offset_weil[__MAX_VARS_E__] = {%d", e_offset[0]);
	for(i = 1; i < __MAX_VARS_E__; i++)
	{
		fprintf(offset_file, ",%d", e_offset[i]);
	}
	fprintf(offset_file, "};\n");
	fclose(offset_file);
}
