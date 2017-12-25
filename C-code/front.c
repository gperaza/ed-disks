#include "structures.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

int front_add();
void draw_frame(long i);
void place_disk_dd(long i);
void place_disk_dw(long i, particle*, particle*);
int check_overlap(particle*, long);
void check_overlap_all();

Front front;
int layer;
extern particle *disk;
extern particle *wall;
extern long nDisks;

void pack_disks() {
    /* Builds a disk packing using an advancing front approach. */
    /* Let's define the front using a linked list. */

    /* Event-drive needs disks spaced out, so let's increase each
       radius a bit to ensure no disk touches. */
    for (long i = 0; i < nDisks; i++) {
        disk[i].R *= 1.01;
    }

    /* Build the initial front including walls. */
    /* Front is Wall_L(a) -> Wall_B(b) -> Wall_R. */
    front.a = front.first = &wall[3];  // left wall
    front.a->next = front.b = &wall[0];  // bottom wall
    front.a->prev = NULL;
    front.last = &wall[1];  // right wall
    front.b->next = front.last;
    front.b->prev = front.a;
    front.last->next = NULL;
    front.last->prev = front.b;
    draw_frame(-1);

    /* Add first disk manually and update front. */
    /* Front is Wall_L -> Disk[0](a) -> Wall_B(b) -> Wall_R. */
    disk[0].x = disk[0].R;
    disk[0].y = disk[0].R;
    disk[0].prev = front.a;
    disk[0].next = front.b;
    front.a->next = front.b->prev = &disk[0];
    front.a = &disk[0];
    draw_frame(0);

    /* Build layers. */
    int success;
    for (long i = 1; i < nDisks; i++) {
        success = 0;
        while (!success) {
            success = front_add(i);
            draw_frame(i);
        }
        if (front.b == front.last) {
            /* Restart front. */
            front.a = front.first;
            front.b = front.a->next;
            layer++;
        }
    }

    /* Recover original radii and move disks a bit to ensure
       collisions with walls occur at different times.*/
    gsl_rng *rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, 138491);
    for (long i = 0; i < nDisks; i++) {
        double random = gsl_rng_uniform_pos(rgen);
        disk[i].R /= 1.01;
        /* disk[i].R /= 2; */
        disk[i].x += disk[i].R*0.004*random;
        disk[i].y += disk[i].R*0.004*random;
    }
    gsl_rng_free(rgen);
    /* Check for overlaps. */
    check_overlap_all();
}

int front_add(long i) {
    particle *curr;

    assert(front.a->type + front.b->type <= 1);

    /* Try to place disk. */
    if (front.a->type == 0 && front.b->type == 0) place_disk_dd(i);
    else if (front.a->type == 1) place_disk_dw(i, front.b, front.a);
    else place_disk_dw(i, front.a, front.b);


    /* Check overlap backward. */
    int overlapb = 0;
    particle *backColl = NULL;
    curr = front.a->prev;
    while (curr != NULL && !overlapb) {
        if (check_overlap(curr, i)) {
            overlapb++;
            backColl = curr;
        }
        curr = curr->prev;
    }

    /* Check overlap forward. */
    int overlapf = 0;
    particle *forwardColl = NULL;
    curr = front.b->next;
    while (curr != NULL && !overlapf) {
        if (check_overlap(curr, i)) {
            overlapf++;
            forwardColl = curr;
        }
        curr = curr->next;
    }

    /* If no overlap advance front. */
    if (overlapf + overlapb == 0) {
        front.a->next = &disk[i];
        disk[i].next = front.b;
        disk[i].prev = front.a;
        front.b->prev = &disk[i];
        if (layer == 0) front.a = &disk[i];
        else {
            front.a = front.b;
            front.b = (front.a->next != NULL) ? front.a->next : front.last;
        }
        return 1;
    }

    /* If just forward overlap adjust front. */
    if (overlapf) {
        forwardColl->prev = front.a;
        front.a->next = forwardColl;
        front.b = forwardColl;
    }

    /* If just backward overlap. */
    if (overlapb) {
        backColl->next = front.b;
        front.b->prev = backColl;
        front.a = backColl;
    }

    return 0;
}

int check_overlap(particle* curr, long i) {
    int overlap = 0;
    double dx, dy, d2, d;

    if (curr->type == 0 /*disk*/) {
        dx = disk[i].x - curr->x;
        dy = disk[i].y - curr->y;
        d2 = dx*dx + dy*dy;
        if (d2 < (disk[i].R + curr->R)*(disk[i].R + curr->R)) overlap = 1;
    } else /*wall*/ {
        dx = curr->nx*(disk[i].x - curr->x);
        dy = curr->ny*(disk[i].y - curr->y);
        d = dx+ dy;
        if (fabs(d) < disk[i].R) overlap = 1;
    }

    return overlap;
}

void check_overlap_all() {
    for (long i = 0; i < nDisks; i++) {
        for (long j = i +1; j < nDisks; j++) {
            if (check_overlap(&disk[i], j)) {
                printf("Overlap between disk %ld and disk %ld\n", i+1, j+1);
                // exit(1);
            };
        }
    }
    for (long i = 0; i < nDisks; i++) {
        for (long j = 0; j < 4; j++) {
            if (check_overlap(&wall[j], i)) {
                printf("Overlap between disk %ld and wall %ld\n", i+1, j+1);
                // exit(1);
            }
        }
    }
}

void check_overlap_other(long e1, long e2) {
    for (long i = 0; i < nDisks; i++) {
        for (long j = i +1; j < nDisks; j++) {
            if (i == e1 && j == e2) continue;
            if (check_overlap(&disk[i], j)) {
                printf("Overlap between disk %ld and disk %ld\n", i+1, j+1);
                exit(1);
            };
        }
    }
    for (long i = 0; i < nDisks; i++) {
        for (long j = 0; j < 4; j++) {
            if (i == e1 && j == e2) continue;
            if (check_overlap(&wall[j], i)) {
                printf("Overlap between disk %ld and wall %ld\n", i+1, j+1);
                exit(1);
            }
        }
    }
}

void place_disk_dd(long i) {
    /* Place new disk touching other two.*/
    double Dx = front.b->x - front.a->x;
    double Dy = front.b->y - front.a->y;
    double L2 = Dx*Dx + Dy*Dy;
    double Ria = disk[i].R + front.a->R;
    double Rib = disk[i].R + front.b->R;
    double A = L2 + Ria*Ria - Rib*Rib;
    double B = sqrt((L2 - (Ria - Rib)*(Ria - Rib))*
                    ((Ria + Rib)*(Ria + Rib) - L2));

    disk[i].y = front.a->y + (A*Dy + B*Dx)*0.5/L2;
    disk[i].x = front.a->x + (A*Dx - B*Dy)*0.5/L2;
}

void place_disk_dw(long i, particle *d, particle *w) {
    /* Place new disk touching a wall and a disk.*/
    /* Signs in this function take into account the fact that walls
       are orthogonal and we always want the grater solution. It is
       not a general solution.*/

    double Dx = w->x - d->x;
    double Dy = w->y - d->y;
    double Ria = disk[i].R + d->R;
    double A = w->nx*Dx + w->ny*Dy + disk[i].R;
    double B = sqrt(Ria*Ria - A*A);

    disk[i].y = d->y + A*w->ny + B*fabs(w->nx);
    disk[i].x = d->x + A*w->nx + B*w->ny;
}
