#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>  
#include <libana.h>
//#include <libana_global_variables.h>

#define CUTOFFBO 0.300
#define DUMPTYPE 1 //0:cobpos-normal 1:cobpos-group

int main() {
    /**************************************************/
    /*         Local Variables Definitions            */
    /**************************************************/
    int   i, t, id1, id2;
    char  str[512], fname[256], outfile[256];
    FILE* fp;
    int step;
    int max_step = 1000000; 

    
    mkdir("222", 0777);

    /********************************************/
    /* Atom Data corresponding to para.rd file  */
    /********************************************/
    set_element(0, 1, "C");
    set_element(1, 2, "H");
    set_element(2, 3, "O");
    set_element(3, 4, "Fe");
    set_element(4, 5, "Al");
    set_element(5, 6, "Ni");
    /* Read input files for LAMMPS */
    read_lammps_input("in.input");

    /**************************************************/
    /*                    Main Loop                   */
    /**************************************************/
    for(step = 0; step <= max_step; step += 10000) {
        
        sprintf(fname, "dump.pos.%d", step);
        fp = fopen(fname, "r");
        if (fp == NULL) {
           
            continue;
        }
        fclose(fp);

       
        printf("Processing step: %d\n", step);
        
       
        read_lammps_bond("bondorder", step, CUTOFFBO);
        
       
        read_dump(step, DUMPTYPE);
        
       
        dfs();

        
        int* target = (int*)malloc(sizeof(int) * atomicity);
        memset(target, 0, sizeof(int) * atomicity);

       
        for (i = 0; i < atomicity; i++) {
            if (gnum[gid[i]] < 4) {
                target[i] = 0;
            } else {
                target[i] = 1;
            }
        }
            
        sprintf(outfile, "222/newdump.%dsteps", step);
        write_dump(outfile, target);
        
       
        free(target);
    }

    return 0;
}
