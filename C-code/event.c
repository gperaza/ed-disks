#include "event.h"

extern double g, wb, Ax, Ay;
extern double t;

double get_dt_event_disk_disk(particle *i, particle *j) {
    /* Get the time (dt) of the next event between disks.

       i: Pointer to a disk structure.
       j: pointer to a disk structure.

       return: time until the event involving i and j.

       Code assumes j>i. Need to consider the case if disks are
       overlapping. If they are overlapping and separating conditional 1
       takes care of this. If the case when they are overlapping and getting
       closer should occur we abort. */

    double rijx = i->x - j->x;
    double rijy = i->y - j->y;
    double vijx = i->vx - j->vx;
    double vijy = i->vy - j->vy;

    double vDotR = (rijx*vijx + rijy*vijy);

    if (vDotR > 0) {
        /* Disks are separating. */
        return INFINITY;
    }

    double Rij = i->R + j->R;
    double vijSqrd = vijx*vijx + vijy*vijy;
    double q = rijx*rijx + rijy*rijy - Rij*Rij;
    if (q < 0) {
        printf("Error: Overlapping disks are getting closer.\n");
        exit(1);
    }

    double det = vDotR*vDotR - vijSqrd*q;
    if (det < 0) {
        /* Imaginary time. */
        return INFINITY;
    }

    /* Keep the shortest time, calculate in this manner to avoid
       precision loss for small initial separation. */
    double dt = q/(-vDotR + sqrt(det));

    assert(dt >= 0);
    return dt;
}

double get_dt_event_disk_box(particle *i /*Disk*/, particle *j /*Wall*/) {
    /* Get the time (dt) of the next event between a disk and a wall
       from a square box containind the disk. This function is faster
       than calculating the event time for a wall at an arbitrary
       angle but only works for a disk in a box. The box is fixed.

       i: Pointer to a disk structure.
       j: pointer to a wall structure.
       g: gravity

       return: time until the event involving i and j.

       We need to consider the case when disk and wall are
       overlapping. For no overlap the function is correct. If walls
       are overlapping and going apart we need to chose the other
       solution to the quadratic equation (last conditional); if they
       overlap and are getting closer we abort. */

    double nx = j->nx;
    double ny = j->ny;
    double xj = j->x;
    double yj = j->y;
    /* Walls are static. */
    double vxj = 0, vyj = 0;
    double xi = i->x;
    double yi = i->y;
    double vxi = i->vx;
    double vyi = i->vy;

    double rijDotn = (xi - xj)*nx + (yi - yj)*ny;
    if (rijDotn < 0) {
        nx = -nx;
        ny = -ny;
        rijDotn = -rijDotn;
    }
    double d = rijDotn - i->R;

    double dd = 2*d;

    double vn = (vxi - vxj)*nx + (vyi - vyj)*ny;

    if (ny == 0) {
        if (vn >= 0) return INFINITY;
        else return -0.5*dd/vn;
    }

    double gn = g*ny;
    double det = vn*vn + dd*gn;

    if (det < 0) return INFINITY;

    double detrt = sqrt(det);
    if (detrt == vn) d = dd = 0;

    if (vn <= 0) {
        if (d >= 0) {
            /* The only positive solution is ny > 0, the shortest one
               if ny < 0. */
            return dd/(-vn + detrt);
        }
        else {
            /* No positive solutions. */
            printf("Something went wrong\n");
            exit(1);
        }
    }
    else {
        if (ny > 0) {
            if (d >= 0) return (vn + detrt)/gn; // only positive solutions
            else {
                /* two positive solutions, choose greater one, the
                   other is nonphysical. */
                return (vn + detrt)/gn;
            }
        }
        else return INFINITY; /* should never collide */
    }
}

double get_dt_event_disk_wall(particle *i /*Disk*/, particle *j /*Wall*/) {
    /* Get the time (dt) of the next event between a disk and a wall
       at an arbitrary angle with arbitrary boundary conditions. If
       walls are those of a static box containing the simulation use
       the faster function for boxes.

       i: Pointer to a disk structure.
       j: pointer to a wall structure.

       g: gravity
       t: current time
       wb: frequancy bottom vibration
       Ax: x component of vibration amplitude
       Ay: y component of vibration amplitude

       return: time until the event involving i and j. */

    double dt = 0;

    double dist_disk_wall(double dt) {

        double xi = i->x + i->vx*dt;
        double yi = i->y + i->vy*dt - 0.5*g*dt*dt;

        double xj = j->x + Ax*sin(wb*(t + dt));
        double yj = j->y + Ay*sin(wb*(t + dt));

        double nx = j->nx;
        double ny = j->ny;

        return fabs((xi-xj)*nx + (yi-yj)*ny) - i->R;
    }

    /* Make sure wall and disk are not touching, if they are we
       missed a collision and need to abort. */
    assert(dist_disk_wall(0) > 0);

    /* Numerically solve for dt. We take the bracket interval from 0
       to the first time the function changes sign. If we can't find
       and interval up to 10 seconds we abort. 10 seconds is a long
       time for systems considered in our simulations. */
    double dt1 = 0, dt2 = 1e-6;
    while (dist_disk_wall(dt1)*dist_disk_wall(dt2) >= 0) {
        dt2 += 1e-6;
    }


    return dt;
}

void update_pos_vel(double dt, particle *disk, long nDisks) {
    /* Updates positions and velocities of all disks from current time
       up to current time + dt.

       dt: Time interval. The simulation will be updated from t to t + dt.
       disk: pointer to particle array.
       nDisks: number of disks in the simulation.
       g: gravity*/

    long i;
    for (i=0; i < nDisks; i++) {
        disk[i].x = disk[i].x + disk[i].vx*dt;
        disk[i].y = disk[i].y + disk[i].vy*dt - 0.5*g*dt*dt;
        disk[i].theta += disk[i].w*dt;
        disk[i].vy = disk[i].vy - g*dt;
    }
}

long get_closest_event(event *cList, long eventN) {
    /* Gets the closest event from the event list cList with eventN
       events.

       cList: pointer to the list containing the events.
       eventN: total number of events, the length of cList.

       return: the index of the closest event */

    long min = 0;
    for (long i = 1; i < eventN; i++) {
        if (*cList[i].time < *cList[min].time) min = i;
    }

    return min;
}

long init_time_array(long nDisks, long nWalls,
                     double **eTime, event *cList,
                     particle *disk, particle*wall) {
    /* Allocates the arrays to store times for the next event between
       elements i and j. Calculates and returns the total number of events.

       nDisks: total number of disks
       nWalls: total number of walls
       eTime: two dimensional array to store event times for i j
       cList: on dimensional array to store all events
       disk: array of disk elements
       wall: array of wall elemnts
       g: gravity

       return: total number of events

       Event times are stored in a two dimensional array for fast
       access and fast update. Uses the function to calculate event
       times between a disk an arbitrary wall for
       inizialization. Since it is only used once the overhead when
       using a box is negligible. */

    /* Total number of possible events */
    long eventN = nDisks*(nDisks - 1)*0.5 + nWalls*nDisks;

    /* Allocate array on continuous memory. */
    eTime = (double**)calloc(nDisks, sizeof(double*));
    double *temp = (double*)calloc((nDisks)*(nDisks + nWalls), sizeof(double));
    for (int i = 0; i < nDisks; i++) {
        eTime[i] = temp + (i * (nDisks + nWalls));
    }
    /* Allocate collision list. */
    cList = (event*)calloc(eventN, sizeof(event));

    /* Initialize values for disks. */
    long k = 0;
    for(long i = 0; i < nDisks; i++) {
        for (long j = i + 1; j < nDisks; j++) {
            eTime[i][j] = get_dt_event_disk_disk(&disk[i], &disk[j]);
            cList[k].idx1 = i;
            cList[k].idx2 = j;
            cList[k].type = 0;
            cList[k].time = &eTime[i][j];
            k++;
        }
    }

    /* Initialize values for walls. */
    for (long i = 0; i < nDisks; i++) {
        for (long j = 0; j < nWalls; j++) {
            eTime[i][nDisks + j] = get_dt_event_disk_wall(&disk[i], &wall[j]);
            cList[k].idx1 = i;
            cList[k].idx2 = j;
            cList[k].type = 1;
            cList[k].time = &eTime[i][nDisks + j];
            k++;
        }
    }

    assert(k == eventN);

    return eventN;
}

void v_update(double dt, long nDisks, particle *vdisk, particle *disk) {
    /* Updates an virtual copy of the simulation up to t + dt. The
       virtual copy is not used for any calculations. It is useful to
       perform data storage, graphics, etc. at arbitrary times
       between events without disrupting the simulation flow.

       dt: time interval. Updates virtual copy up to t + dt.
       nDisks: total number of disks
       vdisk: the virtual copy, array of particle elements
       disk: the actual system, array of particle elements */

    for (long i = 0; i < nDisks; i++) {
        vdisk[i].x = disk[i].x + disk[i].vx*dt;
        vdisk[i].y = disk[i].y + disk[i].vy*dt - 0.5*g*dt*dt;
        vdisk[i].theta = disk[i].theta + disk[i].w*dt;
        vdisk[i].vy = disk[i].vy - g*dt;
    }
}
