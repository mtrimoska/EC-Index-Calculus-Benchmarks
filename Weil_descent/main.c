//
//  main.c
//  weil descent
//
//  Created by Monika Trimoska on 23/07/2018.
//  Copyright © 2018 Monika Trimoska. All rights reserved.
//
#define __FILE_PATH__ "."
//#define __FILE_PATH__ "/Users/monika/casspair/Index_calculus/Weil_descent"

#include<stdio.h>
#include<stdlib.h>
#include <getopt.h>
#include <string.h>
#include "vect_bin.h"
#include "weil_out.h"
#include "semaev_out.h"
#include "weil.h"
#include "params.h"
#include "semaev.h"
#include "nb_lines.h"

/*char irr[162] = "1000100001";
char new_Xr[162] = "010001011"; //110011011
char solution[162] = "110001110000000";
int new_n = 9;
int new_l = 3;
int new_m = 4;*/
/*char irr[162] = "1101101";
char new_Xr[162] = "010001";
char solution[162] = "00010001";
int new_n = 6;
int new_l = 3;
int new_m = 3;*/

char irr[162] = "1100000000101000001";
char new_Xr[162] = "100111011010110100";
char solution[162] = "000011101000010100010111111010000";
int new_n = 18;
int new_l = 6;
int new_m = 4;
int action = 1; // 1 2 22 3
int out = 1; // grobner 1, dimacs 2

/**
 * @fn int scan_opt(int argc, char **argv, const char *opt)
 * @brief It scans and checks command line options
 */
int scan_opt(int argc, char **argv, const char *opt) {
	char c;
	while ((c = getopt (argc, argv, opt)) != -1)
		switch (c) {
			case 'n': new_n = atoi(optarg); break;
			case 'l': new_l = atoi(optarg); break;
			case 'm': new_m = atoi(optarg); break;
			case 'a': action = atoi(optarg); break;
			case 'o': strncmp("dimacs", optarg, 6) ? (out = 1) : (out = 2); break;
			case 'r': strcpy(irr, optarg); break;
			case 'x': strcpy(new_Xr, optarg); break;
			case 's': strcpy(solution, optarg); break;
			default: return(-1);
		}
	return(0);
}

int main(int argc, char * argv[])
{
	FILE *out_file;
	int constant[162] = {0};
	scan_opt(argc, argv, "n:l:r:x:s:m:a:o:");
	/** Action 1*/
	if(action == 1){}
	
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
		printf("--OK:%d---%d -> %d\n", indices_to_num(ii) == cpt, cpt, indices_to_num(ii));
	}*/
	
	
	/** Action 2.1*/
	if(action == 2)
	{
		set_params(new_n, new_m, new_l, new_n, new_m - 1);
		create_X_masks();
		create_e_X_correspondence();
		
		/** Action 2.1 - dimacs */
		if(out == 2)
		{
			compute_offsets_X();
			if((out_file = fopen(__FILE_PATH__"/out/2_cnf_x.fordimacs", "w")) == NULL)
			{
				printf("Error file 2_cnf_x.fordimacs\n");
			}
			else
			{
				nb_lines += create_cnf_dimacs_out(out_file, X_offset);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/3_corr.fordimacs", "w")) == NULL)
			{
				printf("Error file 3_corr.fordimacs\n");
			}
			else
			{
				nb_lines += e_X_correspondence_dimacs_out(out_file);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/2_corr.forxorDand", "w")) == NULL)
			{
				printf("Error file 2_corr.forxorDand\n");
			}
			else
			{
				e_X_correspondence_xorDand_out(out_file);
			}
			fclose(out_file);
		}
		
		/** Action 2.1 - grobner*/
		if(out == 1)
		{
			if((out_file = fopen(__FILE_PATH__"/out/1_vars.grobner", "w")) == NULL)
			{
			printf("Error file 1_vars.grobner\n");
			}
			else
			{
				init_vars_grobner_out(out_file);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/2_corr.grobner", "w")) == NULL)
			{
				printf("Error file 2_corr.grobner\n");
			}
			else
			{
				e_X_correspondence_grobner_out(out_file);
			}
			fclose(out_file);
			
			//to unsym
			if((out_file = fopen(__FILE_PATH__"/out/1_vars.sage", "w")) == NULL)
			{
				printf("Error file 1_vars.sage\n");
			}
			else
			{
				init_to_unsym_out(out_file);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/2_corr.sage", "w")) == NULL)
			{
				printf("Error file 2_corr.sage\n");
			}
			else
			{
				e_X_correspondence_to_unsym_out(out_file);
			}
			fclose(out_file);
		}
	}
	/** Action 2.2  */
	if(action == 22)
	{
		set_params(new_n, new_m, new_l, 4 * new_n + 8 * new_l - 11, new_m - 2);
		create_e_masks();
		create_semaev();
		
		/** Action 2.2 - dimacs */
		if(out == 2)
		{
			compute_offsets_e(12);
			if((out_file = fopen(__FILE_PATH__"/out/4_cnf_e.fordimacs", "w")) == NULL)
			{
				printf("Error file 4_cnf_e.fordimacs\n");
			}
			else
			{
				nb_lines += create_cnf_dimacs_out(out_file, e_offset);
			}
			fclose(out_file);
		}
	}
	
	
	
	/** Action 3 */
	if(action == 3)
	{
		set_params(new_n, new_m, new_l, 4 * new_n + 8 * new_l - 11, new_m - 2);
		init(f)
		create_weil(f, constant, irr, new_Xr, solution, 12);
		
		
		/** Action 3 - grobner*/
		if(out == 1)
		{
			if((out_file = fopen(__FILE_PATH__"/out/3_system.grobner", "w")) == NULL)
			{
				printf("Error file 3_system.grobner\n");
			}
			else
			{
				system_grobner_out(out_file, f, constant);
			}
			fclose(out_file);
			
			//to unsym
			if((out_file = fopen(__FILE_PATH__"/out/3_system.sage", "w")) == NULL)
			{
				printf("Error file 3_system.sage\n");
			}
			else
			{
				system_to_unsym_out(out_file, f, constant);
			}
			fclose(out_file);
		}
		
		/** Action 3 - dimacs*/
		if(out == 2)
		{
			if((out_file = fopen(__FILE_PATH__"/out/5_system.fordimacs", "w")) == NULL)
			{
				printf("Error file 5_system.fordimacs\n");
			}
			else
			{
				nb_lines += system_dimacs_out(out_file, f, constant);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/3_system.forxorDand", "w")) == NULL)
			{
				printf("Error file 3_system.forxorDand\n");
			}
			else
			{
				system_xorDand_out(out_file, f, constant);
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/1_header.fordimacs", "w")) == NULL)
			{
				printf("Error file 1_header.fordimacs\n");
			}
			else
			{
				header_dimacs_out(out_file, nb_lines);
				nb_lines = 0;
			}
			fclose(out_file);
			
			if((out_file = fopen(__FILE_PATH__"/out/1_header.forxorDand", "w")) == NULL)
			{
				printf("Error file 1_header.forxorDand\n");
			}
			else
			{
				header_xorDand_out(out_file);
			}
			fclose(out_file);
			printf("*\n");
		}
	}
	
	//update nb of lines
	if((out_file = fopen("nb_lines.h", "w")) == NULL)
	{
		printf("Error file nb_lines.h\n");
	}
	else
	{
		fprintf(out_file, "int nb_lines = %d;", nb_lines);
	}
	fclose(out_file);
	
	return(0);
}
