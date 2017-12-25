#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "structures.h"

double get_dt_event_disk_disk(particle*, particle*);
double get_dt_event_disk_box(particle*, particle*);
double get_dt_event_disk_wall(particle*, particle*);
void update_pos_vel(double, particle*, long);
long get_closest_event(event*, long);
long init_time_array(long, long, double**, event*, particle*, particle*);
void v_update(double, long, particle*, particle*);
