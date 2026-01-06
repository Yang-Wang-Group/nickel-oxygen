#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <libana.h>
//#include <libana_global_variables.h>

#define CUTOFFBO 0.300
#define DUMPTYPE 1 //0:cobpos-normal 1:cobpos-group

int main() {

	/**************************************************/
	/*         Local Variables Definations            */
	/**************************************************/
	int   i, j, t, id1, id2;
	char  str[512], workdir[64], fname[64];	workdir[0] = '\0';
	FILE* fp;
	int step;

	/********************************************/
	/* Atom Data corresponding to para.rd file  */
	/********************************************/
	set_element(0, 1, "C");
	set_element(1, 2, "H");
	set_element(2, 3, "O");
	set_element(3, 4, "Fe");
	set_element(4, 5, "Al");
	//set_element(5, 6, "Ni");

	/* Read input files for LAMMPS */
	read_lammps_input("in.input");

	//int   asp[MAXATOMSPECIES];
	//float ifb[MAXATOMSPECIES][MAXATOMSPECIES];

	/**************************************************/
	/*                    Main Loop                   */
	/**************************************************/
	step = 960000;
	read_lammps_bond("bondorder", step, CUTOFFBO);
	int sum = 0;

	/* Read dump.pos.X file and get posv, velo, force, and charge */

	read_dump(step, DUMPTYPE);

		
		for (i = 0; i < atomicity; i++)
		{
			nbrNu[i][0];
			//printf("%d\n",nbrNu[i][0]);
			sum += nbrNu[i][0];
			
		}
printf("Average: %f\n", (float)sum / atomicity);

}


