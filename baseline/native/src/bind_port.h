#pragma once
#include "geometry/unit_pattern.h"

// Verbatim port of bind.cpp's max_opening_angle(): the kinematic bound from
// kirigami::max_opening_angle(), optionally reduced by their 1%-step forward
// collision scan.  This is theta_max "as THEY compute it".
double max_opening_angle_with_collisions(const kirigami::UnitPattern &pattern,
                                         bool detect_collisions);

// Diagnostic: same predicate and same deployed configuration as their scan,
// but reports which face pair overlaps.  Faces are indices into the
// merge_close_verts()/remove_unused_verts() mesh their scan builds.
bool first_collision_at(const kirigami::UnitPattern &pattern, double angle,
                        int *f0out, int *f1out);
