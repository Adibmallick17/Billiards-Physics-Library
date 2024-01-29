#ifndef PHYLIB_H
#define PHYLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
#define PHYLIB_TABLE_WIDTH  1350.0  //mm
#define PHYLIB_TABLE_LENGTH 2700.0  //mm
#define PHYLIB_BALL_RADIUS  28.5    //mm
#define PHYLIB_BALL_DIAMETER  57.0  //mm
#define PHYLIB_HOLE_RADIUS  114.0   //mm
#define PHYLIB_DRAG         150.0   //mm/s^2
#define PHYLIB_VEL_EPSILON  0.01    //mm/s
#define PHYLIB_SIM_RATE     0.0001  //S
#define PHYLIB_MAX_TIME     600.0   //s
#define PHYLIB_MAX_OBJECTS  26

typedef struct {
    double x;
    double y;
} phylib_coord;

typedef enum {
    PHYLIB_STILL_BALL,
    PHYLIB_ROLLING_BALL,
    PHYLIB_HOLE,
    PHYLIB_HCUSHION,
    PHYLIB_VCUSHION
} phylib_object_type;

typedef struct {
    unsigned char number;
    phylib_coord pos;
} phylib_still_ball_data;

typedef struct {
    unsigned char number;
    phylib_coord pos;
    phylib_coord vel;
    phylib_coord acc;
} phylib_rolling_ball_data;

typedef struct {
    phylib_coord pos;
} phylib_hole_data;

typedef struct {
    double y;
} phylib_hcushion_data;

typedef struct {
    double x;
} phylib_vcushion_data;

typedef struct {
    phylib_object_type type;
    union {
        phylib_still_ball_data still_ball;
        phylib_rolling_ball_data rolling_ball;
        phylib_hole_data hole;
        phylib_hcushion_data hcushion;
        phylib_vcushion_data vcushion;
    } obj;
} phylib_object;

typedef struct {
    double time;
    phylib_object *object[PHYLIB_MAX_OBJECTS];
} phylib_table;

// Function prototypes
phylib_object *phylib_new_still_ball(unsigned char number, phylib_coord *pos);
phylib_object *phylib_new_rolling_ball(unsigned char number, phylib_coord *pos, phylib_coord *vel, phylib_coord *acc);
phylib_object *phylib_new_hole(phylib_coord *pos);
phylib_object *phylib_new_hcushion(double y);
phylib_object *phylib_new_vcushion(double x);
phylib_table *phylib_new_table(void);
void phylib_copy_object(phylib_object **dest, phylib_object **src);
phylib_table *phylib_copy_table(phylib_table *table); 
void phylib_add_object(phylib_table *table, phylib_object *object);
void phylib_free_table(phylib_table *table);
phylib_coord phylib_sub(phylib_coord c1, phylib_coord c2);
double phylib_length(phylib_coord c);
double phylib_dot_product(phylib_coord a, phylib_coord b);
double phylib_distance(phylib_object *obj1, phylib_object *obj2);
void phylib_roll(phylib_object *new, phylib_object *old, double time);
unsigned char phylib_stopped(phylib_object *object);
void phylib_bounce(phylib_object **a, phylib_object **b);
unsigned char phylib_rolling(phylib_table *t);
phylib_table *phylib_segment(phylib_table *table);
void phylib_print_table(phylib_table *table);
void phylib_free_object(phylib_object *object);

#endif // PHYLIB_H
