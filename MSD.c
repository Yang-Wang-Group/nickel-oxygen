#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ana.h>

#define CUTOFFBO 0.300
#define DUMPTYPE 1
#define START_STEP 60000
#define END_STEP 960000
#define STEP_INTERVAL 10000


extern int atomicity;
extern int *aSpec;
extern int *gid;
extern int *gnum;
extern float box_x, box_y, box_z;
extern float **pos;


void set_element(int index, int atomic_number, char *symbol);
int read_lammps_input(char *filename);
int read_dump(int step, int dumptype);
int dfs();


typedef struct {
    int index;       
    double x0, y0, z0; 
} SelectedOxygen;


void pbc_correct_displacement(double *dx, double *box) {
    for (int i = 0; i < 3; i++) {
        if (dx[i] > box[i]/2) dx[i] -= box[i];
        else if (dx[i] < -box[i]/2) dx[i] += box[i];
    }
}


double unfold_coord(double current, double initial, double box_size) {
    double diff = current - initial;
    if (diff > box_size/2) return current - box_size;
    else if (diff < -box_size/2) return current + box_size;
    return current;
}

int main() {
    /**************************************************/
    /*         Local Variables Definations            */
    /**************************************************/
    int i, j, step;
    char fname[128];
    FILE* fp_msd = fopen("msd_results.txt", "w"); 
    
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
    /*             MSD Calculation Setup               */
    /**************************************************/
    SelectedOxygen* selectedOxygens = NULL;
    int numSelected = 0;
    double box[3]; 
    
    /**************************************************/
    /*                    Main Loop                   */
    /**************************************************/
    for (step = START_STEP; step <= END_STEP; step += STEP_INTERVAL) {
       
        sprintf(fname, "dump.pos.%d", step);
        
        
        read_dump(step, DUMPTYPE);
        
        
        box[0] = box_x;
        box[1] = box_y;
        box[2] = box_z;
        
       
        if (step == START_STEP) {
            
            read_lammps_bond("bondorder", step, CUTOFFBO);
            dfs();
            
            
            for (i = 0; i < atomicity; i++) {
                if (aSpec[i] == 3 && gnum[gid[i]] > 15) { 
                    numSelected++;
                    selectedOxygens = (SelectedOxygen*)realloc(selectedOxygens, 
                                            numSelected * sizeof(SelectedOxygen));
                    
                    
                    selectedOxygens[numSelected-1].index = i;
                    selectedOxygens[numSelected-1].x0 = pos[i][0];
                    selectedOxygens[numSelected-1].y0 = pos[i][1];
                    selectedOxygens[numSelected-1].z0 = pos[i][2];
                    
                    
                    printf("Selected atom: Index=%d, Pos=(%.6f, %.6f, %.6f)\n", 
                           i, pos[i][0], pos[i][1], pos[i][2]);
                }
            }
            printf("Selected %d oxygen atoms at step %d\n", numSelected, step);
            
            
            fprintf(fp_msd, "%d 0.0\n", step);
            continue;
        }
        
        /**************************************************/
        /*         Calculate MSD for selected atoms       */
        /**************************************************/
        double total_msd = 0.0;
        int found_count = 0;
        
        
        for (j = 0; j < numSelected; j++) {
            int idx = selectedOxygens[j].index;
            
           
            if (idx < 0 || idx >= atomicity) {
                fprintf(stderr, "Warning: Invalid index %d for step %d\n", idx, step);
                continue;
            }
            
           
            printf("Processing atom index %d at step %d: ", idx, step);
            printf("Initial pos=(%.6f, %.6f, %.6f) ", 
                   selectedOxygens[j].x0, selectedOxygens[j].y0, selectedOxygens[j].z0);
            printf("Current pos=(%.6f, %.6f, %.6f)\n", 
                   pos[idx][0], pos[idx][1], pos[idx][2]);
            
           
            double unfolded_x = unfold_coord(pos[idx][0], selectedOxygens[j].x0, box[0]);
            double unfolded_y = unfold_coord(pos[idx][1], selectedOxygens[j].y0, box[1]);
            double unfolded_z = unfold_coord(pos[idx][2], selectedOxygens[j].z0, box[2]);
            
            
            double dx = unfolded_x - selectedOxygens[j].x0;
            double dy = unfolded_y - selectedOxygens[j].y0;
            double dz = unfolded_z - selectedOxygens[j].z0;
            
            
            double dr[3] = {dx, dy, dz};
            pbc_correct_displacement(dr, box);
            
           
            double msd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            total_msd += msd;
            found_count++;
           
            printf("Displacement: dx=%.6f, dy=%.6f, dz=%.6f, MSD=%.6f\n", 
                   dr[0], dr[1], dr[2], msd);
        }
        
       
        double avg_msd = (found_count > 0) ? total_msd / found_count : 0.0;
        fprintf(fp_msd, "%d %.8f\n", step, avg_msd);
        printf("Step %d: Average MSD = %.6f (%d/%d atoms found)\n", 
               step, avg_msd, found_count, numSelected);
    }
    
    /**************************************************/
    /*                   Cleanup                      */
    /**************************************************/
    fclose(fp_msd);
    free(selectedOxygens);
    
    return 0;
}
