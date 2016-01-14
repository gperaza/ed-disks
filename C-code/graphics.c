#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <math.h>
#include <stdio.h>
#include "structures.h"

#define DSP_W 1000
#define DSP_H ((int)(DSP_W*(box_h/box_w)))
#define SHRINK (0.99)
extern double box_w, box_h;
extern long nDisks;
extern particle *disk;
extern particle *wall;

void draw_arrow(cairo_t *cr, double start_x, double start_y,
                double end_x, double end_y, double arrow_lenght) {
    double angle = atan2(end_y - start_y, end_x - start_x) + M_PI;
    double arrow_degrees = M_PI/6;
    double x1 = end_x + arrow_lenght * cos(angle - arrow_degrees);
    double y1 = end_y + arrow_lenght * sin(angle - arrow_degrees);
    double x2 = end_x + arrow_lenght * cos(angle + arrow_degrees);
    double y2 = end_y + arrow_lenght * sin(angle + arrow_degrees);
    cairo_move_to(cr, start_x, start_y);
    cairo_line_to(cr, end_x, end_y);
    /* cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); */
    cairo_stroke(cr);
    cairo_move_to(cr, end_x, end_y);
    cairo_line_to(cr, x1, y1);
    cairo_line_to(cr, x2, y2);
    cairo_line_to(cr, end_x, end_y);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
}

void draw_disk_label(cairo_t *cr, long id, double x, double y, double r,
                     double angle, int color) {

    double maxrot = angle;
    double intensity = angle/maxrot;
    char label[50];
    intensity = 0.5;
    /*draw and color the disk*/
    cairo_arc(cr, x, y, r, 0, 2*M_PI);
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    if (color) cairo_set_source_rgb(cr, 0.0, 1.0, 0.0);
    else cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_fill_preserve(cr);

    if (angle < 0) {
        cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, intensity);
        cairo_fill_preserve(cr);
    } else if (angle > 0) {
        cairo_set_source_rgba(cr, 0.0, 0.0, 1.0, intensity);
        cairo_fill_preserve(cr);
    }

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_stroke(cr);

    /*draw a line to indicate rotation*/
    cairo_arc(cr, x, y, r, angle, angle);
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_line_to(cr, x, y);
    cairo_stroke(cr);

    /* draw an arrow to indicate direction */
    /* if (disk[id].vx != 0 || disk[id].vy != 0) { */
    /*     double l = sqrt(disk[id].vx*disk[id].vx + disk[id].vy*disk[id].vy); */
    /*     draw_arrow(cr, x, y, x + disk[id].vx/l*r, y + disk[id].vy/l*r, 0.4*r); */
    /* } */

    //draw label
    sprintf(label,"%-3ld",id+1);
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, label);
    cairo_new_path(cr);

    return;
}

void graphics(long e1, long e2, int type) {
    cairo_surface_t *surface;
    cairo_t *cr;
    double scale = SHRINK*fmin(DSP_W/box_w, DSP_H/box_h);
    char out_file[50];
    static long framenumber = 0;
    cairo_matrix_t font_matrix;

    framenumber++;
    sprintf(out_file, "frame%06ld.png", framenumber);
    //surface = cairo_svg_surface_create(out_file, DSP_W, DSP_H);
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, DSP_W, DSP_H);
    cr = cairo_create(surface);
    cairo_translate(cr, 0, DSP_H);
    cairo_scale( cr, 1.0, -1.0 ); //flip y axis
    cairo_translate(cr, (1-SHRINK)/2*DSP_W, (1-SHRINK)/2*DSP_H);
    cairo_scale(cr, scale, scale);
    cairo_set_font_size(cr, fmin(box_h, box_w)*0.02);
    cairo_get_font_matrix(cr, &font_matrix);
    font_matrix.yy = -font_matrix.yy; // negative size to vertically flip text
    cairo_set_font_matrix(cr, &font_matrix);

    /* Draw disks */
    long i;
    cairo_set_line_width(cr,fmin(box_w, box_h)/800);
    for (i = 0;  i < nDisks; i++) {
        if (i == e1 || (i == e2 && type == 0))
            draw_disk_label(cr, i, disk[i].x,
                            disk[i].y, disk[i].R, disk[i].theta, 1);
        else
            draw_disk_label(cr, i, disk[i].x,
                            disk[i].y, disk[i].R, disk[i].theta, 0);
    }

    /* Draw box */
    cairo_set_line_width(cr,fmin(box_w, box_h)/400);
    cairo_rectangle(cr, 0, 0, box_w, box_h);
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_stroke(cr);
    for (long i = 0; i < 4; i++) {
        cairo_set_source_rgb(cr, 0, 0, 0);
        if (type == 1 && i == e2) cairo_set_source_rgb(cr, 0, 1, 0);
        draw_arrow(cr, wall[i].x, wall[i].y,
                   wall[i].nx*box_w*0.02 + wall[i].x,
                   wall[i].ny*box_h*0.02 + wall[i].y, 0.4*0.02);
    }

    cairo_surface_write_to_png(surface, out_file);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    return;
}
