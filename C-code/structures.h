typedef struct particle
{
    int type; // 0: disk, 1:wall
    double x; double y; double theta;
    double vx; double vy; double w;
    double R;
    double mass;
    double iMoment;
    /* Extra parameters for lines. */
    /* Lines are defined by a point and a normal. */
    double nx, ny;
    /* Linked list next element. Used for packing.*/
    struct particle *next;
    struct particle *prev;
} particle;

typedef struct event {
    long idx1; long idx2;
    double *time;
    int type;
} event;

/* Front structure*/
typedef struct Front {
    particle *first;
    particle *last;
    /* Active front, from a to b. */
    particle *a;
    particle *b;
} Front;

void Front_push();
void Front_pop();
