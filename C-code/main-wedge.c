/* Event-driven dynamics for 2D disks with thermal walls. */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "event.h"

void collision_rule_disk_disk(long, long);
void collision_rule_disk_wall(long, long);
void check_overlap_other();
void write(double);
void graph(double);
void write_packing();
void get_input();

/* Global variables */
double g, density, Gamma, wb, mu, kn, kt, dampR, meanR;
long nDisks, colN, eventN;
double cTime, runTime, thermalTime;
double timeForGraph, timeForWrite;
double box_w, box_h;

particle *disk, *vdisk, *wall;
double **eTime;
event *cList;
double * temp;

FILE *phase_space_fptr;

long seed;
gsl_rng *rgen;

void init_packing() {
    /* Three radii. */
    double sigma = meanR*0.25;
    double R1 = meanR - sigma, R2 = meanR, R3 = meanR + sigma;

    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, seed);

    /* Allocate memory for 4 walls and nDisks disks.*/
    disk = (particle*)calloc(nDisks, sizeof(particle));
    wall = (particle*)calloc(4, sizeof(particle));
    vdisk = (particle*)calloc(nDisks, sizeof(particle));

    /* Set up box. Other geometries are possible.*/
    wall[0].x = box_w*0.5; wall[0].y = 0; wall[0].nx = 0; wall[0].ny = 1;
    wall[1].x = box_w; wall[1].y = box_h*0.5; wall[1].nx = -1; wall[1].ny = 0;
    wall[2].x = box_w*0.5; wall[2].y = box_h; wall[2].nx = 0;  wall[2].ny = -1;
    wall[3].x = 0; wall[3].y = box_h*0.5;     wall[3].nx = 1;  wall[3].ny = 0;
    for (long i = 0; i < 4; i++) wall[i].type = 1;

    /* Set up disk properties. Consider 3 different radii.*/
    long j;
    for (long i = 0; i < nDisks; i++) {
        if       (i < (long)(nDisks/3))   disk[i].R = R1;
        else if  (i < (long)(2*nDisks/3)) disk[i].R = R2;
        else                              disk[i].R = R3;

        disk[i].mass = M_PI*(disk[i].R*disk[i].R)*density;
        disk[i].iMoment = 0.5*disk[i].mass*(disk[i].R*disk[i].R);
        disk[i].type = 0;
    }
    /*Now make a random permutation (Knuth permutation).*/
    particle tempDisk;
    for (long i = 0; i < nDisks; i++) {
        j = gsl_rng_uniform_int(rgen, i + 1);
        tempDisk = disk[i];
        disk[i] = disk[j];
        disk[j] = tempDisk;
    }

    /* Setup a dense packing using an advancing front algorithm. */
    pack_disks();
}

int main() {
    get_input();

    phase_space_fptr = fopen("phase_space.out", "w");

    init_packing();
    eventN = init_time_array();
    write_packing(); /* Write packing file. */

    cTime = -thermalTime;
    double nextWrite = 0;       /* Do not write thermalization. */
#ifdef GRAPHICS
    double nextGraph = 0;
#endif
    double checkOverlapStep = 1;
    double checkOverlapTime = cTime + checkOverlapStep;

    printf("Running: 0%%"); fflush(stdout);
    while (cTime < runTime) {
        colN++;

        /* Get next event participants. */
        long nextEvent = get_closest_event();
        long i = cList[nextEvent].idx1;
        long j =cList[nextEvent].idx2;
        int type = cList[nextEvent].type;

        /* Recalculate time to avoid accumulation of numerical
           errors. */
        double dt =
            !type ? get_dt_event_disk_disk(i, j) : get_dt_event_disk_wall(i, j);
        double nTime = cTime + dt;
        if (nTime <= cTime) {
            printf("Error: Simulation not advancing. dt=%g\n", dt);
            exit(1);
        }

        while (nTime > nextWrite) {
            write(nextWrite);
            nextWrite += timeForWrite;
            printf("\rRunning: %0.lf%%", cTime/runTime*100); fflush(stdout);
        }

#ifdef GRAPHICS
        while (nTime > nextGraph) {
            graph(nextGraph);
            nextGraph += timeForGraph;
        }
#endif

#ifdef GRAPHICS
        // graphics(i, j, type);
#endif
        /* Advance all particles to t + dt. */
        update_pos_vel(dt);
        cTime = nTime;

        // graphics(i, j, type);

        if (cTime > checkOverlapTime) {
            check_overlap_other(i, j);
        }
        while (cTime > checkOverlapTime) checkOverlapTime += checkOverlapStep;

        /* Apply collision rule to colliding pair. */
        if (type == 0) collision_rule_disk_disk(i, j);
        else collision_rule_disk_wall(i, j);

        // graphics(i, j, type);

        /* Update the events that involve i and j. */
        for ( long l = 0; l < i; l++)
            eTime[l][i] = cTime + get_dt_event_disk_disk(i, l);
        for ( long l = i + 1; l < nDisks; l++)
            eTime[i][l] = cTime + get_dt_event_disk_disk(i, l);
        for ( long l = nDisks; l < nDisks + 4; l++)
            eTime[i][l] = cTime + get_dt_event_disk_wall(i, l - nDisks);

        if (type == 0) {
            for ( long l = 0; l < j; l++)
                eTime[l][j] = cTime + get_dt_event_disk_disk(j, l);
            for ( long l = j + 1; l < nDisks; l++)
                eTime[j][l] = cTime + get_dt_event_disk_disk(j, l);
            for ( long l = nDisks; l < nDisks + 4; l++)
                eTime[j][l] = cTime + get_dt_event_disk_wall(j, l - nDisks);
        } else {
            for (long l = 0; l < nDisks; l++)
                eTime[l][j + nDisks] = cTime + get_dt_event_disk_wall(l, j);
        }
    }
    printf("\rDone: 100%%\n");
    printf("Total collisions: %ld\n", colN);

    write_packing(); /* Write packing file. */

    gsl_rng_free(rgen);
    free(disk); free(wall); free(vdisk); free(eTime); free(temp); free(cList);
    fclose(phase_space_fptr);

    return 0;
}

void write_packing() {
    FILE *packing_fptr;
    packing_fptr = fopen("input_phase_space.out", "w");
    fprintf(packing_fptr, "%ld %ld %e %e\n", nDisks, 4L, box_w, box_h);
    for (long i = 0; i < 4; i++) {
        fprintf(packing_fptr, "%e %e %e %e\n",
                wall[i].x, wall[i].y, wall[i].nx, wall[i].ny);
    }
    for (long i = 0; i < nDisks; i++) {
        fprintf(packing_fptr, "%16.12e %16.12e %16.12e %16.12e %16.12e\n",
                disk[i].x, disk[i].y, disk[i].R, disk[i].mass, disk[i].iMoment);
    }
    fclose(packing_fptr);
#ifdef GRAPHICS
    graphics(-1, -1, 0);
#endif
}

void get_input() {
    char input[200];
    char type[200];
    char value[200];
    // int inputs = 25, countInputs = 0;
    FILE *fp;
    double pr, freq;

    if ( (fp = fopen("input_file", "r")) == NULL ) {
        printf("Error opening input file.\n");
        exit(1);
    }

    while (fgets(input, sizeof(input), fp) != NULL) {
        sscanf(input,"%s %s",type, value);
        // countInputs++;
        if (strcmp(type,"#nParticles") == 0) {
            sscanf(value, "%ld", &nDisks);
            printf("Particles = %ld.\n", nDisks);
        } else if (strcmp(type,"#seed") == 0) {
            sscanf(value, "%ld", &seed);
            printf("Seed for rng = %ld.\n", seed);
        } else if (strcmp(type,"#box_w") == 0) {
            sscanf(value, "%lf", &box_w);
            printf("box_w = %g disks.\n", box_w);
        } else if (strcmp(type,"#box_h") == 0) {
            sscanf(value, "%lf", &box_h);
            printf("box_h = %g disks.\n", box_h);
        } else if (strcmp(type,"#freq") == 0) {
            sscanf(value, "%lf", &freq);
            printf("Bottom frequency = %g Hz.\n", freq);
        } else if (strcmp(type,"#dimensionlessAc") == 0) {
            sscanf(value, "%lf", &Gamma);
            printf("Bottom acceleration = %g g's.\n", Gamma);
        } else if (strcmp(type,"#gravity") == 0) {
            sscanf(value, "%lf", &g);
            printf("Gravity = %g m/s^2.\n", g);
        } else if (strcmp(type,"#runTime") == 0) {
            sscanf(value, "%lf", &runTime);
            printf("Simulation time = %g s. \n", runTime);
        } else if (strcmp(type,"#timeForGraph") == 0) {
            sscanf(value, "%lf", &timeForGraph);
            printf("Time between graphics = %g s. \n", timeForGraph);
        } else if (strcmp(type,"#timeForWriteRun") == 0) {
            sscanf(value, "%lf", &timeForWrite);
            printf("Time between writes during main simulation= %g s. \n",
                   timeForWrite);
        } else if (strcmp(type,"#meanR") == 0) {
            sscanf(value, "%lf", &meanR);
            printf("Mean disk radius = %g m. \n", meanR);
        } else if (strcmp(type,"#density") == 0) {
            sscanf(value, "%lf", &density);
            printf("Disk's density = %g kg/m^3. \n", density);
        } else if (strcmp(type,"#kn") == 0) {
            sscanf(value, "%lf", &kn);
            printf("Normal elastic constant = %g. \n", kn);
        } else if (strcmp(type,"#pr") == 0) {
            sscanf(value, "%lf", &pr);
            printf("Disk's Poisson's ratio = %g. \n", pr);
        } else if (strcmp(type,"#mu") == 0) {
            sscanf(value, "%lf", &mu);
            printf("Disk's friction coefficient = %g. \n", mu);
        } else if (strcmp(type,"#vGamma") == 0) {
            sscanf(value, "%lf", &dampR);
            printf("Viscoelastic damping ratio = %g. \n", dampR);
        } else if (strcmp(type,"#thermalTime") == 0) {
            sscanf(value, "%lf", &thermalTime);
            printf("Thermalization time = %g. \n", thermalTime);
        }
    }
    fclose(fp);

    wb = 2*M_PI*freq;
    kt = (1 - pr)/(1 - 0.5*pr)*kn;
    printf("Calculated tangential elastic constant = %g. \n", kt);
    box_w *= 2*meanR;
    box_h *= 2*meanR;

    return;
}

void write(double nextWrite) {
    /* Virtually updates the system and writes some data. */
    double dt = nextWrite - cTime;
    v_update(dt);
    fprintf(phase_space_fptr, "%e ", cTime);
    for (long i = 0; i < nDisks; i++) {
        fprintf(phase_space_fptr, "%e %e %e %e %e %e %e %e %e ",
                disk[i].R, disk[i].mass, disk[i].iMoment,
                vdisk[i].x, vdisk[i].y, vdisk[i].theta,
                disk[i].vx, vdisk[i].vy, disk[i].w);
    }
    fprintf(phase_space_fptr, "\n");
}
