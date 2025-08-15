#---------- 0.0 ----------#
import math

def intersect_ray_sphere(O, D):
    # O, D: each a 3‐tuple (x,y,z)
    # Sphere center C=(0,0,0), radius r=1.
    # Solve ‖O + t D‖^2 = 1  ⇒  a t^2 + b t + c = 0

    # Quadratic coefficients:
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2.0*(O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = (O[0]*O[0] + O[1]*O[1] + O[2]*O[2]) - 1.0

    # Discriminant
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0   # no real roots ⇒ no intersection

    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2.0*a)
    t2 = (-b + sqrt_disc) / (2.0*a)

    # We want the smallest t ≥ 0 (ray starts at t=0 and goes forward)
    t_min = min(t1, t2)
    t_max = max(t1, t2)

    if t_min >= 0.0:
        return t_min
    elif t_max >= 0.0:
        return t_max
    else:
        return -1.0   # both intersections are "behind" the ray origin


#---------- 0.1 ----------#
import math

def intersect_sphere(O, D):
    """
    O: ray origin, tuple of 3 floats (ox, oy, oz)
    D: ray direction, tuple of 3 floats (dx, dy, dz)
    returns: smallest positive t such that |O + t D| = 1,
             or -1 if no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    # Quadratic coefficients for |O + t D|^2 = 1
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0

    disc = b*b - 4*a*c
    if disc < 0:
        # no real roots: ray misses sphere
        return -1.0

    sqrt_disc = math.sqrt(disc)
    # two possible solutions
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)

    # we want the smallest positive t
    t_min = min(t0, t1)
    t_max = max(t0, t1)

    if t_min > 0:
        return t_min
    elif t_max > 0:
        return t_max
    else:
        # both intersections are behind the ray origin
        return -1.0


#---------- 0.2 ----------#
import math

def intersect_sphere(O, D, radius=1.0):
    # O, D are 3-tuples (ox, oy, oz), (dx, dy, dz)
    # sphere is centered at (0,0,0) with given radius
    # returns smallest t>=0 or -1 if no hit

    # unpack
    ox, oy, oz = O
    dx, dy, dz = D

    # coefficients of quadratic at^2 + bt + c = 0
    a = dx*dx + dy*dy + dz*dz
    b = 2.0 * (ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - radius*radius

    # discriminant
    disc = b*b - 4.0*a*c
    if disc < 0:
        return -1.0  # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2.0 * a)
    t1 = (-b + sqrt_disc) / (2.0 * a)

    # pick the smallest non‐negative t
    if t0 >= 0:
        return t0
    if t1 >= 0:
        return t1
    return -1.0


#---------- 0.3 ----------#
import math

def intersect_ray_sphere(O, D):
    """
    O: tuple of 3 floats = ray origin (ox, oy, oz)
    D: tuple of 3 floats = ray direction (dx, dy, dz)
    Returns:
      t >= 0 : distance along the ray to the first intersection
      -1     : if there is no intersection
    """
    ox, oy, oz = O
    dx, dy, dz = D

    # Quadratic coefficients
    a = dx*dx + dy*dy + dz*dz
    b = 2.0 * (ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0

    # Discriminant
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0   # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)

    # Two roots of the quadratic
    t1 = (-b - sqrt_disc) / (2.0*a)
    t2 = (-b + sqrt_disc) / (2.0*a)

    # We want the smallest non‐negative t
    t_min = min(t1, t2)
    t_max = max(t1, t2)

    if t_min >= 0.0:
        return t_min
    elif t_max >= 0.0:
        return t_max
    else:
        return -1.0


#---------- 0.4 ----------#
import math

def intersect_sphere(O, D):
    """
    O: ray origin, tuple of 3 floats (ox, oy, oz)
    D: ray direction, tuple of 3 floats (dx, dy, dz)
    Returns the smallest positive t at which O + t D hits the unit
    sphere at the origin, or -1 if there is no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D
    # Quadratic coefficients
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0

    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0

    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)

    # We want the smallest t > 0
    tmin = min(t for t in (t1, t2) if t > 0, default=None)
    return tmin if tmin is not None else -1.0


#---------- 0.5 ----------#
import math

def intersect_ray_sphere(O, D):
    """
    O: ray origin, tuple of 3 floats (Ox, Oy, Oz)
    D: ray direction, tuple of 3 floats (Dx, Dy, Dz)
    returns: the smallest t >= 0 such that O + t D lies on the sphere of radius 1 at (0,0,0),
             or -1 if there is no intersection.
    """
    # coefficients of a t^2 + b t + c = 0
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2*(O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = O[0]*O[0] + O[1]*O[1] + O[2]*O[2] - 1.0

    # discriminant
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0   # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)
    # two solutions of the quadratic
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)

    # we want the smallest non-negative t
    if t0 >= 0:
        return t0
    if t1 >= 0:
        return t1
    return -1.0       # both intersections are “behind” the ray origin


#---------- 0.6 ----------#
import math

def intersect_ray_sphere(O, D):
    """
    Compute the intersection of a ray R(t) = O + t*D with the sphere
    centered at (0,0,0) of radius 1.
    O, D: tuples of 3 floats
    returns: the smallest t >= 0 at which the ray hits the sphere,
             or -1.0 if there is no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Quadratic coefficients a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2.0 * (Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    # Discriminant
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0   # no real roots, no intersection

    sqrt_disc = math.sqrt(disc)

    # Two solutions to the quadratic
    t0 = (-b - sqrt_disc) / (2.0*a)
    t1 = (-b + sqrt_disc) / (2.0*a)

    # We want the smallest non-negative t
    if t0 > t1:
        t0, t1 = t1, t0

    if t0 >= 0.0:
        return t0
    if t1 >= 0.0:
        return t1
    return -1.0


#---------- 0.7 ----------#
import math

def intersect_sphere(O, D):
    """
    Ray-sphere intersection.
    Sphere: center = (0,0,0), radius = 1.
    Ray: R(t) = O + t*D, t >= 0
    O, D: 3-tuples or lists of floats.
    Returns the smallest t >= 0 at which the ray hits the sphere,
    or -1.0 if there is no intersection.
    """
    # Quadratic coefficients for |O + tD|^2 = 1
    # a t^2 + b t + c = 0
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2.0 * (O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = (O[0]*O[0] + O[1]*O[1] + O[2]*O[2]) - 1.0

    disc = b*b - 4.0*a*c
    if disc < 0.0:
        # no real roots → no intersection
        return -1.0

    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2.0 * a)
    t1 = (-b + sqrt_disc) / (2.0 * a)

    # We want the smallest non-negative t
    if t0 >= 0.0:
        return t0
    if t1 >= 0.0:
        return t1
    return -1.0


#---------- 0.8 ----------#
import math

def intersect_sphere(O, D):
    """
    O: tuple of 3 floats, the ray origin (Ox, Oy, Oz)
    D: tuple of 3 floats, the ray direction (Dx, Dy, Dz)
       (does not have to be normalized)
    Returns:
      t >= 0 (float) = distance along the ray to the first intersection with the
                       unit sphere centered at (0,0,0),
      or -1 if there is no intersection.
    """
    # Unpack for convenience
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Compute quadratic coefficients
    # a = D·D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    # b = 2 O·D
    b = 2.0 * (Ox*Dx + Oy*Dy + Oz*Dz)
    # c = O·O - R^2, here R=1
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    # Discriminant
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0   # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)

    # Two solutions of the quadratic
    t1 = (-b - sqrt_disc) / (2.0*a)
    t2 = (-b + sqrt_disc) / (2.0*a)

    # We want the smallest non‐negative t
    t_min = float('inf')
    if t1 >= 0.0:
        t_min = t1
    if t2 >= 0.0 and t2 < t_min:
        t_min = t2

    if t_min == float('inf'):
        return -1.0  # both intersections are behind the ray origin

    return t_min


#---------- 0.9 ----------#
import math

def intersect_sphere(O, D):
    """
    O: tuple of 3 floats, ray origin
    D: tuple of 3 floats, ray direction
    returns: the smallest positive t such that |O + t D| = 1,
             or -1 if there is no intersection.
    """
    # Coefficients of the quadratic a t^2 + b t + c = 0
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2 * (O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = O[0]*O[0] + O[1]*O[1] + O[2]*O[2] - 1.0   # radius^2 = 1

    # Discriminant
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0     # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)

    # the two roots of the quadratic
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)

    # We assume the ray starts outside the sphere, so c>0.
    # We want the smallest positive t:
    if t0 > 0 and t1 > 0:
        return min(t0, t1)
    elif t0 > 0:
        return t0
    elif t1 > 0:
        return t1
    else:
        return -1.0     # both intersections are “behind” the origin


#---------- 1.0 ----------#
def intersect_ray_unit_cube(O, D):
    """
    Compute the intersection of a ray R(t) = O + t D with the unit cube [0,1]^3.
    O, D are length‐3 tuples or lists.  Returns the smallest t >= 0 at which
    the ray hits the cube, or -1 if there is no intersection.
    """
    # Initialize entry and exit t’s
    t_min = -float('inf')
    t_max =  float('inf')
    
    # For each axis, compute intersection t’s with the two planes
    for i in range(3):
        origin = O[i]
        direction = D[i]
        slab_min = 0.0
        slab_max = 1.0

        if abs(direction) < 1e-12:
            # Ray is parallel to this pair of planes.
            # If origin is not between the planes, no intersection.
            if origin < slab_min or origin > slab_max:
                return -1.0
            # otherwise it imposes no new constraint on t
        else:
            # Compute t where the ray intersects the two planes
            t1 = (slab_min - origin) / direction
            t2 = (slab_max - origin) / direction
            # Order them
            t_near = min(t1, t2)
            t_far  = max(t1, t2)
            # Intersect this interval with the running [t_min, t_max]
            t_min = max(t_min, t_near)
            t_max = min(t_max, t_far)
            # If your entry is after your exit, no intersection
            if t_min > t_max:
                return -1.0

    # At this point [t_min, t_max] is the intersection interval with the cube.
    # We assumed the ray starts outside, so we want the first non‐negative t.
    if t_max < 0:
        # Entire intersection is “behind” the ray origin
        return -1.0

    # If t_min is negative, the origin is inside the cube,
    # so the ray exits at t_max.  Otherwise it enters at t_min.
    t_hit = t_min if t_min >= 0 else t_max
    return t_hit


#---------- 1.1 ----------#
def intersect_unit_cube(O, D, eps=1e-8):
    """
    Intersect a ray R(t) = O + t*D against the axis-aligned unit cube
    defined by 0 <= x,y,z <= 1.
    O, D : sequences of length 3
    Returns: the smallest t >= 0 at which the ray hits the cube,
             or -1.0 if there is no intersection.
    """
    t_near = -float('inf')
    t_far  =  float('inf')

    for i in range(3):
        o = O[i]
        d = D[i]

        if abs(d) < eps:
            # Ray is parallel to this pair of planes.
            # If origin is not between the slabs, no hit.
            if o < 0.0 or o > 1.0:
                return -1.0
            # otherwise, it “always” overlaps this slab, so we do nothing
        else:
            # Compute intersection with the two planes perpendicular to axis i
            t1 = (0.0 - o) / d
            t2 = (1.0 - o) / d
            # Order them so t1 is entry, t2 is exit
            t_entry = min(t1, t2)
            t_exit  = max(t1, t2)
            # Intersect this slab’s t-interval [t_entry, t_exit]
            t_near = max(t_near, t_entry)
            t_far  = min(t_far,  t_exit)
            # Early out: if we’ve already missed
            if t_near > t_far:
                return -1.0

    # If the whole intersection interval is behind the ray origin → no hit
    if t_far < 0.0:
        return -1.0

    # Otherwise the first valid intersection is max(t_near, 0).
    # (But the problem statement guarantees the ray starts outside,
    #  so t_near will be >= 0 if there is a hit.)
    return t_near if t_near >= 0.0 else t_far


# Example usage:
O = (2.0, 0.5, 0.5)
D = (-1.0, 0.0, 0.0)
t = intersect_unit_cube(O, D)
if t > 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit cube at t =", t, "point =", P)
else:
    print("No intersection")


#---------- 1.2 ----------#
def ray_unit_cube_intersect(O, D, eps=1e-8):
    """
    Intersect ray R(t) = O + t D with the unit cube [0,1]^3.
    O, D: 3‐tuples or lists of floats.
    Returns the smallest t >= 0 at which the ray hits or grazes the cube,
    or -1 if there is no intersection.
    """
    t_min = float('-inf')
    t_max = float('inf')

    for i in range(3):
        origin = O[i]
        direction = D[i]

        if abs(direction) < eps:
            # Ray is parallel to this axis-slab.
            # If origin not within the slab [0,1], no hit:
            if origin < 0.0 or origin > 1.0:
                return -1.0
            # else, it stays inside this slab for all t, so skip.
        else:
            # Compute intersection t’s with the two planes perpendicular to axis i
            t1 = (0.0 - origin) / direction
            t2 = (1.0 - origin) / direction
            t_near = min(t1, t2)
            t_far  = max(t1, t2)

            # Shrink our global intersection interval [t_min, t_max]
            t_min = max(t_min, t_near)
            t_max = min(t_max, t_far)

            # If the interval is empty, no intersection
            if t_min > t_max:
                return -1.0

    # At this point [t_min, t_max] is the intersection interval.
    # We assume the ray starts outside the cube, so the entry is t_min.
    # If t_max < 0, the whole interval is behind the ray origin.
    if t_max < 0.0:
        return -1.0

    # If t_min is negative but t_max positive, the origin was inside—
    # but the problem says the ray always starts outside, so normally
    # t_min >= 0. Return t_min in all valid cases.
    return t_min if t_min >= 0.0 else t_max


# Example tests

print(ray_unit_cube_intersect((2,0.5,0.5), (-1,0,0)))
# should print 1.0, because the ray hits x=1 at t=(1-2)/(-1)=1

print(ray_unit_cube_intersect((0.5,0.5,2), (0,0,-1)))
# should also print 1.0, because it hits z=1 at t=(1-2)/(-1)=1

print(ray_unit_cube_intersect((2,2,2), (1,1,1)))
# should print -1.0 (ray goes away from the cube)


#---------- 1.3 ----------#
def intersect_ray_unit_cube(O, D):
    """
    O: origin of the ray, tuple/list of 3 floats
    D: direction of the ray, tuple/list of 3 floats
    returns: smallest t >= 0 such that O + t D lies in [0,1]^3,
             or –1 if no intersection.
    """
    t_min = float('-inf')
    t_max = float('inf')

    # for each axis, intersect the ray with the pair of planes {x=0, x=1}, {y=0, y=1}, {z=0, z=1}
    for i in range(3):
        o = O[i]
        d = D[i]
        if abs(d) < 1e-9:
            # Ray is parallel to the planes perpendicular to this axis.
            # If the origin is outside the slab [0,1] on this axis, no intersection.
            if o < 0.0 or o > 1.0:
                return -1.0
            # otherwise, it imposes no new constraints on t
        else:
            # Compute the t's where the ray crosses the two planes
            t1 = (0.0 - o) / d
            t2 = (1.0 - o) / d
            t_near = min(t1, t2)
            t_far  = max(t1, t2)

            # Intersect this slab’s interval [t_near, t_far] with our running [t_min, t_max]
            t_min = max(t_min, t_near)
            t_max = min(t_max, t_far)

            # If the intervals ever become disjoint, there is no intersection
            if t_min > t_max:
                return -1.0

    # At this point [t_min, t_max] is the interval on the ray that lies inside the cube.
    # We want the closest intersection in front of the origin, so t >= 0.
    if t_max < 0.0:
        # the whole intersection interval is behind the ray origin
        return -1.0

    # If t_min >= 0, that’s our first hit.  Otherwise the ray origin was inside the cube
    # (but the problem guarantees the ray starts outside), so we’d take t_max.
    return t_min if t_min >= 0.0 else t_max


#---------- 1.4 ----------#
def intersect_unit_cube(O, D):
    """
    O: tuple of floats (Ox, Oy, Oz) – ray origin
    D: tuple of floats (Dx, Dy, Dz) – ray direction
    Returns:
      t >= 0 : the smallest distance along the ray to the cube
      -1     : if no intersection
    """
    t_min = float('-inf')
    t_max = float('inf')
    bounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    
    for i in range(3):
        Oi = O[i]
        Di = D[i]
        lo, hi = bounds[i]
        
        if abs(Di) < 1e-12:
            # Ray is parallel to this slab. If origin not within slab, no intersection.
            if Oi < lo or Oi > hi:
                return -1.0
            # else it imposes no new constraints on t
        else:
            # Compute intersection t with the two planes of this slab:
            t1 = (lo - Oi) / Di
            t2 = (hi - Oi) / Di
            # Order them so t1 is entry, t2 is exit
            t_entry_i = min(t1, t2)
            t_exit_i  = max(t1, t2)
            
            # Merge with global interval [t_min, t_max]
            if t_entry_i > t_min:
                t_min = t_entry_i
            if t_exit_i < t_max:
                t_max = t_exit_i
            
            # If at any point the intervals become disjoint, no hit:
            if t_min > t_max:
                return -1.0
    
    # At this point [t_min, t_max] is the intersection interval
    # We want the nearest non-negative t:
    if t_max < 0:
        # Intersection is “behind” the ray origin
        return -1.0
    
    # If t_min < 0 < t_max, the ray origin is inside the box;
    # but the problem guarantees the ray starts outside, so t_min>=0.
    return t_min if t_min >= 0 else t_max


# Example usage:
O = (2.0, 0.5, 0.5)
D = (-1.0, 0.0, 0.0)
t = intersect_unit_cube(O, D)
if t >= 0:
    hit_point = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", hit_point)
else:
    print("No intersection")


#---------- 1.5 ----------#
def intersect_ray_unit_cube(O, D):
    """
    O: tuple of 3 floats, ray origin (Ox,Oy,Oz)
    D: tuple of 3 floats, ray direction (Dx,Dy,Dz)
    returns: smallest t >= 0 such that O + t*D hits the cube [0,1]^3, or -1 if no hit
    """
    t_near = -float('inf')
    t_far  =  float('inf')

    # for each slab X, Y, Z
    for i in range(3):
        origin_i = O[i]
        dir_i    = D[i]
        if abs(dir_i) < 1e-12:
            # Ray is parallel to this axis
            # If origin is outside the slab, no hit
            if origin_i < 0.0 or origin_i > 1.0:
                return -1.0
            # else we just skip updating t_near/t_far for this axis
        else:
            # Compute intersection t's with the two planes of the slab
            t1 = (0.0 - origin_i) / dir_i
            t2 = (1.0 - origin_i) / dir_i
            t_min_i = min(t1, t2)
            t_max_i = max(t1, t2)
            # Shrink our global [t_near,t_far] interval
            t_near = max(t_near, t_min_i)
            t_far  = min(t_far,  t_max_i)
            # If the intervals don’t overlap, there is no hit
            if t_near > t_far:
                return -1.0

    # If the whole box is “behind” the ray origin, no hit
    if t_far < 0.0:
        return -1.0

    # If t_near is negative, the origin is inside the box, so we exit at t_far
    return t_near if t_near >= 0.0 else t_far


#---------- 1.6 ----------#
def intersect_ray_unit_cube(O, D):
    """
    O, D: 3‐tuples (ox, oy, oz), (dx, dy, dz)
    Returns the smallest t >= 0 such that O + t*D lies on or in the cube [0,1]^3,
    or -1 if there is no intersection.
    """
    t_min = float('-inf')
    t_max = float('inf')
    
    # for each axis i = 0,1,2 (x,y,z)
    for i in range(3):
        o = O[i]
        d = D[i]
        
        if abs(d) > 1e-12:
            # compute intersection parameters with the two planes for this axis
            t1 = (0.0 - o) / d
            t2 = (1.0 - o) / d
            t_near = min(t1, t2)
            t_far  = max(t1, t2)
        else:
            # ray is parallel to this axis
            # if origin not between slab planes, no hit
            if o < 0.0 or o > 1.0:
                return -1.0
            # otherwise it imposes no new constraint
            t_near = float('-inf')
            t_far  = float('inf')
        
        # tighten our global t‐interval
        t_min = max(t_min, t_near)
        t_max = min(t_max, t_far)
        
        # if the intervals do not overlap, there is no intersection
        if t_max < t_min:
            return -1.0
    
    # at this point [t_min, t_max] is the interval of intersection
    # we want the first point t >= 0
    if t_max < 0.0:
        # whole intersection is behind the ray origin
        return -1.0
    
    # if the “entry” t is negative but the exit is positive,
    # it means the ray origin is inside the cube—but the
    # problem guarantees the origin is outside, so t_min should
    # already be >= 0.  We’ll clamp anyway:
    return t_min if t_min >= 0.0 else t_max


#---------- 1.7 ----------#
def intersect_unit_cube(O, D):
    """
    O: tuple of 3 floats = ray origin
    D: tuple of 3 floats = ray direction
    returns: smallest t>=0 so that O + t*D lies on/in the cube [0,1]^3,
             or -1 if no intersection.
    """
    t_min = -float('inf')
    t_max =  float('inf')
    # for each axis x=0, y=1, z=2
    for i in range(3):
        Oi = O[i]
        Di = D[i]
        if abs(Di) < 1e-12:
            # Ray parallel to this slab; must lie between the planes
            if Oi < 0 or Oi > 1:
                return -1.0
        else:
            # compute intersection t's with the two planes of the slab
            t1 = (0.0 - Oi) / Di
            t2 = (1.0 - Oi) / Di
            t_near = min(t1, t2)
            t_far  = max(t1, t2)
            # tighten overall intersection interval [t_min, t_max]
            if t_near > t_min:
                t_min = t_near
            if t_far  < t_max:
                t_max = t_far
            # if ever empty, no hit
            if t_min > t_max:
                return -1.0

    # at this point [t_min, t_max] is the interval along the ray
    # where it is inside the cube
    if t_max < 0:
        # whole interval is behind the ray origin
        return -1.0

    # t_min is the first entry point; if that’s behind the origin
    # then the ray must start inside the cube (not allowed here by assumption),
    # but we can return t_max as the exit if you want that.
    return t_min if t_min >= 0 else t_max


#---------- 1.8 ----------#
def intersect_unit_cube(O, D):
    """
    O = (Ox, Oy, Oz)  ray origin
    D = (Dx, Dy, Dz)  ray direction
    returns the smallest t >= 0 such that O + t*D ∈ [0,1]^3,
    or -1 if no intersection.
    """
    t_enter = float('-inf')
    t_exit  = float('inf')

    for i in range(3):
        Oi = O[i]
        Di = D[i]
        if abs(Di) < 1e-12:
            # Ray is parallel to the slab.  If origin is not inside, no hit.
            if Oi < 0.0 or Oi > 1.0:
                return -1.0
            # else the ray is inside these two planes for all t, so no update
        else:
            t1 = (0.0 - Oi) / Di
            t2 = (1.0 - Oi) / Di
            t_near = min(t1, t2)
            t_far  = max(t1, t2)
            # tighten the global interval
            if t_near > t_enter:
                t_enter = t_near
            if t_far < t_exit:
                t_exit = t_far
            # if the interval is empty, no hit
            if t_enter > t_exit:
                return -1.0

    # at this point [t_enter, t_exit] is the intersection interval
    if t_enter >= 0.0:
        return t_enter
    else:
        # the cube is “behind” the ray origin
        return -1.0


#---------- 1.9 ----------#
def intersect_ray_unit_cube(O, D):
    """
    O: tuple (Ox, Oy, Oz) ray origin
    D: tuple (Dx, Dy, Dz) ray direction
    returns: smallest t >= 0 such that O + t D lies on the cube [0,1]^3,
             or -1 if there is no intersection.
    """
    import math

    # Initialize running entry and exit ts
    t_enter = -math.inf
    t_exit  =  math.inf

    # For each axis do the slab test
    for i in range(3):
        o = O[i]
        d = D[i]
        slab_min = 0.0
        slab_max = 1.0

        if abs(d) < 1e-12:
            # Ray is parallel to this axis. If origin not within the slab, no hit.
            if o < slab_min or o > slab_max:
                return -1.0
            # Otherwise the ray is “inside” this pair of planes for all t, so skip
        else:
            # Compute intersection t’s with the two planes
            t0 = (slab_min - o) / d
            t1 = (slab_max - o) / d
            # t0 might be the “far” intersection if d<0
            t_near = min(t0, t1)
            t_far  = max(t0, t1)
            # Intersect the intervals
            t_enter = max(t_enter, t_near)
            t_exit  = min(t_exit,  t_far)
            # If they’ve crossed, there’s no intersection
            if t_enter > t_exit:
                return -1.0

    # At this point [t_enter, t_exit] is the intersection interval with the cube.
    # We only want t >= 0, and the first hit:
    if t_exit < 0:
        # The entire hit interval is “behind” the ray origin
        return -1.0
    if t_enter >= 0:
        return t_enter
    else:
        # Ray origin is inside the cube; the first positive hit is t_exit
        return t_exit


# Example usage:
O = (-1.0, 0.5, 0.5)   # some ray origin
D = ( 1.0,    0.0,  0.0)   # some ray direction
t = intersect_ray_unit_cube(O, D)
if t >= 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", P)
else:
    print("No intersection")


#---------- 2.0 ----------#
def intersect_square(O, D):
    """
    O, D : 3‐tuples or lists (ox,oy,oz), (dx,dy,dz)
    returns the smallest t >= 0 so that O + t D hits the unit
    square in the plane z=0 with corners (0,0,0)..(1,1,0),
    or -1 if there is no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    # 1) Ray parallel to the plane?
    if abs(dz) < 1e-8:
        return -1.0

    # 2) Solve for intersection with z=0: oz + t*dz = 0 → t = -oz/dz
    t = -oz / dz

    # 3) We only want intersections in front of the ray origin
    if t < 0:
        return -1.0

    # 4) Compute intersection point and test if it's inside the square
    x = ox + t*dx
    y = oy + t*dy
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.1 ----------#
def intersect_unit_square(O, D):
    """
    O, D: 3‐tuples or lists of floats giving the ray origin and direction.
    Returns the smallest t >= 0 such that O + t D hits the square
    with corners (0,0,0),(1,0,0),(1,1,0),(0,1,0), or -1 if no hit.
    """
    # Plane z=0 has normal (0,0,1).  Ray: O + t D.
    # solve O.z + t D.z = 0  =>  t = -O.z / D.z
    Dz = D[2]
    # parallel to the plane?
    if abs(Dz) < 1e-9:
        return -1.0

    t = -O[2] / Dz
    # we only care about intersections in front of the ray origin
    if t < 0:
        return -1.0

    # compute the hit point
    x = O[0] + t * D[0]
    y = O[1] + t * D[1]
    # check if it’s inside the square [0,1]×[0,1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.2 ----------#
def intersect_ray_unit_square(O, D):
    """
    Intersect a ray R(t) = O + t D with the unit square in the z=0 plane
    whose corners (in CCW order) are
         p0=(0,0,0), p1=(1,0,0), p2=(1,1,0), p3=(0,1,0).
    Returns the smallest t>=0 at which the ray hits the square,
    or -1 if there is no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect the plane z=0
    #    solve Oz + t Dz = 0  =>  t = -Oz/Dz
    if abs(Dz) < 1e-9:
        # Ray is parallel to the plane (or nearly so) -> no intersection
        return -1.0

    t = -Oz / Dz
    if t < 0:
        # Intersection is "behind" the ray origin
        return -1.0

    # 2) compute the hit point P = O + t D
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    # Pz would be zero by construction

    # 3) test whether P lies inside the unit square [0,1]×[0,1]
    #    (we include the boundary so grazing counts as an intersection)
    if 0.0 <= Px <= 1.0 and 0.0 <= Py <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.3 ----------#
def intersect_ray_square(O, D, eps=1e-8):
    """
    O: ray origin, tuple of 3 floats (Ox, Oy, Oz)
    D: ray direction, tuple of 3 floats (Dx, Dy, Dz)
    returns: smallest t >= 0 such that O + t*D hits the square,
             or -1 if no intersection
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) Intersect with plane z = 0
    #    plane equation: z = 0  =>  Oz + t*Dz = 0  =>  t = -Oz / Dz
    if abs(Dz) < eps:
        # Ray is parallel to the plane, no hit (we assume origin is "outside")
        return -1.0

    t = -Oz / Dz
    if t < 0:
        # Intersection is "behind" the ray origin
        return -1.0

    # 2) Compute the hit point in the plane
    x = Ox + t * Dx
    y = Oy + t * Dy

    # 3) Test whether (x,y) is inside the unit square [0,1]×[0,1]
    #    we include boundary hits as intersections (grazing counts)
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.4 ----------#
def intersect_unit_square(O, D, eps=1e-8):
    """
    O, D: 3‐tuples (Ox,Oy,Oz), (Dx,Dy,Dz)
    returns smallest t>=0 so that O + t D hits the square,
    or -1 if no hit.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with the plane z=0
    #    plane eqn: z=0 -> 0 = Oz + t*Dz  ==>  t = -Oz/Dz
    if abs(Dz) < eps:
        # Ray is (nearly) parallel to the plane; no intersection
        return -1.0

    t = -Oz / Dz
    if t < 0:
        # intersection is “behind” the ray origin
        return -1.0

    # 2) compute the intersection point P = O + t D
    x = Ox + t * Dx
    y = Oy + t * Dy

    # 3) test whether (x,y) lies inside [0,1]×[0,1]
    #    (grazing on an edge is allowed so we use <= and >=)
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.5 ----------#
def intersect_square(O, D):
    """
    O, D: each a tuple (x, y, z) of floats
    returns: t >= 0 of closest intersection, or -1 if none
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) if ray is (nearly) parallel to the plane z=0, no intersection
    if abs(Dz) < 1e-8:
        return -1.0

    # 2) solve Oz + t*Dz = 0  =>  t = -Oz / Dz
    t = -Oz / Dz

    # 3) we only care about intersections in front of the ray origin
    if t < 0:
        return -1.0

    # 4) compute the hit-point
    x = Ox + t * Dx
    y = Oy + t * Dy
    # z = 0 by construction

    # 5) check if (x,y) lies within [0,1]×[0,1] (inclusive)
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


#---------- 2.6 ----------#
def intersect_square(O, D, eps=1e-8):
    """
    O: tuple of 3 floats = ray origin (Ox,Oy,Oz)
    D: tuple of 3 floats = ray direction (Dx,Dy,Dz)
    returns: float t >= 0 or -1.0 if no intersection
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with the plane z=0
    #    plane equation: Z = 0
    #    solve Oz + t*Dz = 0  ==>  t = -Oz/Dz
    if abs(Dz) < eps:
        # ray is parallel to the plane
        return -1.0

    t = -Oz / Dz
    if t < 0:
        # intersection is behind the ray origin
        return -1.0

    # 2) compute intersection point
    x = Ox + t*Dx
    y = Oy + t*Dy

    # 3) check if (x,y) lies inside [0,1]×[0,1]
    #    (we include the boundary so grazing counts)
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t

    return -1.0


#---------- 2.7 ----------#
def ray_square_intersection(O, D):
    """
    Ray–square intersection.
    O : tuple of 3 floats = ray origin  (Ox,Oy,Oz)
    D : tuple of 3 floats = ray direction (Dx,Dy,Dz)
    
    Square is assumed given by the corners (in CCW order):
      p0 = (0,0,0)
      p1 = (1,0,0)
      p2 = (1,1,0)
      p3 = (0,1,0)
    
    Returns:
      t >= 0  = distance along R(t)=O+tD to the first intersection,
      -1      = no intersection
    """
    # square corners
    p0 = (0.0, 0.0, 0.0)
    p1 = (1.0, 0.0, 0.0)
    p3 = (0.0, 1.0, 0.0)
    
    # helper vector ops
    def sub(a,b):    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a,b):    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a,b):  return ( a[1]*b[2] - a[2]*b[1],
                              a[2]*b[0] - a[0]*b[2],
                              a[0]*b[1] - a[1]*b[0] )
    
    # 1) form plane from p0,p1,p3
    u = sub(p1, p0)
    v = sub(p3, p0)
    n = cross(u, v)        # plane normal
    
    # 2) intersect ray with plane: t = dot(p0 - O, n) / dot(D, n)
    denom = dot(D, n)
    if abs(denom) < 1e-8:
        # parallel => no intersection
        return -1.0
    
    numer = dot(sub(p0, O), n)
    t = numer / denom
    if t < 0:
        # intersection is “behind” the ray origin
        return -1.0
    
    # 3) compute the intersection point P
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])
    
    # 4) express P in the (u,v) basis of the square:
    #    P = p0 + a*u + b*v  -->  solve for (a,b)
    w = sub(P, p0)
    uu = dot(u,u)
    uv = dot(u,v)
    vv = dot(v,v)
    wu = dot(w,u)
    wv = dot(w,v)
    
    # matrix M = [[uu,uv],[uv,vv]] ; rhs = [wu, wv]
    det = uu*vv - uv*uv
    if abs(det) < 1e-8:
        # degenerate square?
        return -1.0
    
    a = ( vv*wu - uv*wv ) / det
    b = ( uu*wv - uv*wu ) / det
    
    # inside if 0 <= a <= 1 and 0 <= b <= 1
    if 0.0 <= a <= 1.0 and 0.0 <= b <= 1.0:
        return t
    else:
        return -1.0


# Example usage:

# Ray straight down onto the square:
O = (0.3, 0.4, 2.0)
D = (0.0, 0.0, -1.0)
print(ray_square_intersection(O, D))   # should be 2.0

# Ray missing the square:
O2 = (1.5, 0.5, 2.0)
D2 = (0.0, 0.0, -1.0)
print(ray_square_intersection(O2, D2))  # should be -1.0


#---------- 2.8 ----------#
def intersect_square(O, D, eps=1e-8):
    """
    Ray-square intersection for the square with corners
      (0,0,0), (1,0,0), (1,1,0), (0,1,0).
    O, D: 3‐tuples or lists (ox,oy,oz), (dx,dy,dz)
    Returns the smallest t>=0 such that O + t D hits the square,
    or -1 if no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    # 1) intersect with the plane z = 0
    #    plane normal = (0,0,1), plane eqn: z = 0
    #    denom = D·N = dz
    if abs(dz) < eps:
        # Ray is parallel to the plane
        return -1.0

    # t such that oz + t*dz = 0  =>  t = -oz / dz
    t = -oz / dz
    if t < 0:
        # Intersection is "behind" the ray origin
        return -1.0

    # 2) compute the hit point
    x = ox + t*dx
    y = oy + t*dy

    # 3) check if inside the unit‐square [0,1]×[0,1]
    #    we include the boundary (grazing counts as hit)
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t

    return -1.0


#---------- 2.9 ----------#
def ray_square_intersection(O, D):
    """
    O : tuple of 3 floats = ray origin (Ox,Oy,Oz)
    D : tuple of 3 floats = ray direction (Dx,Dy,Dz)
    returns the t >= 0 of the intersection R(t)=O+tD with the unit square
    (0,0,0),(1,0,0),(1,1,0),(0,1,0), or -1.0 if no intersection.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with plane z = 0
    #    plane equation: z = 0
    #    Oz + t*Dz = 0   ==>   t = -Oz / Dz
    eps = 1e-8
    if abs(Dz) < eps:
        # ray is parallel to the plane
        return -1.0

    t = -Oz / Dz
    if t < 0.0:
        # intersection is "behind" the ray origin
        return -1.0

    # 2) compute the hit point
    x = Ox + t * Dx
    y = Oy + t * Dy
    # z would be zero (by construction)

    # 3) test if (x,y) is inside the unit square [0,1] x [0,1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    else:
        return -1.0


# --- simple test cases ---

# shoot straight down onto the square
O = (0.5, 0.5, 5.0)
D = (0.0, 0.0, -1.0)
print(ray_square_intersection(O, D))  # should print 5.0

# miss to the side
O = (2.0, 0.5, 5.0)
D = (0.0, 0.0, -1.0)
print(ray_square_intersection(O, D))  # should print -1.0

# grazing the edge at x=1
O = (2.0, 1.0, 5.0)
D = (0.0, 0.0, -1.0)
print(ray_square_intersection(O, D))  # should print 5.0  (grazing counts)


#---------- 3.0 ----------#
import math

def intersect_ray_cylinder(O, D):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: smallest t >= 0 at which the ray O + t D hits the cylinder,
             or -1 if there is no intersection
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # We'll collect all valid intersections here:
    ts = []

    # 1) Side of the infinite cylinder x^2 + z^2 = 1
    #    Solve  (Ox + t Dx)^2 + (Oz + t Dz)^2 = 1
    A = Dx*Dx + Dz*Dz
    B = 2*(Ox*Dx + Oz*Dz)
    C = Ox*Ox + Oz*Oz - 1

    if abs(A) > 1e-12:
        disc = B*B - 4*A*C
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-B - sqrt_disc) / (2*A)
            t2 = (-B + sqrt_disc) / (2*A)
            for t in (t1, t2):
                if t >= 0:
                    y_hit = Oy + t*Dy
                    # check if hit is between y=0 and y=1
                    if 0 <= y_hit <= 1:
                        ts.append(t)

    # 2) Bottom cap: plane y = 0
    if abs(Dy) > 1e-12:
        t_bot = (0 - Oy) / Dy
        if t_bot >= 0:
            x_hit = Ox + t_bot*Dx
            z_hit = Oz + t_bot*Dz
            if x_hit*x_hit + z_hit*z_hit <= 1:
                ts.append(t_bot)

        # 3) Top cap: plane y = 1
        t_top = (1 - Oy) / Dy
        if t_top >= 0:
            x_hit = Ox + t_top*Dx
            z_hit = Oz + t_top*Dz
            if x_hit*x_hit + z_hit*z_hit <= 1:
                ts.append(t_top)

    if not ts:
        return -1.0

    return min(ts)


# Example usage:
if __name__ == "__main__":
    # shoot a ray from (2, 0.5, 0) in direction (-1, 0, 0) 
    # i.e. straight towards the side of the cylinder
    O = (2.0, 0.5, 0.0)
    D = (-1.0, 0.0, 0.0)
    t = intersect_ray_cylinder(O, D)
    print("t =", t)   # Should print t = 1.0, because at t=1 we are at x=1, y=0.5, z=0


#---------- 3.1 ----------#
import math

def intersect_ray_cylinder(O, D):
    """
    O: (Ox, Oy, Oz) ray origin
    D: (Dx, Dy, Dz) ray direction
    returns t >= 0 of first intersection, or -1 if none
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []

    # 1) intersect infinite cylinder x^2 + z^2 = 1
    #    solve a*t^2 + b*t + c = 0
    a = Dx*Dx + Dz*Dz
    b = 2*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1

    eps = 1e-8
    if abs(a) > eps:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b - sqrt_disc) / (2*a)
            t2 = (-b + sqrt_disc) / (2*a)
            for t in (t1, t2):
                if t >= 0:
                    y = Oy + t*Dy
                    # clip to the finite cylinder 0 <= y <= 1
                    if 0 <= y <= 1:
                        ts.append(t)

    # 2) intersect bottom cap at y=0 and top cap at y=1
    #    plane y = ycap  =>  t = (ycap - Oy)/Dy
    #    check if (x,z) lies inside unit disk
    if abs(Dy) > eps:
        for ycap in (0.0, 1.0):
            t = (ycap - Oy) / Dy
            if t >= 0:
                x = Ox + t*Dx
                z = Oz + t*Dz
                if x*x + z*z <= 1.0:
                    ts.append(t)

    if not ts:
        return -1.0

    return min(ts)


#---------- 3.2 ----------#
import math

def intersect_cylinder(O, D):
    """
    O: (Ox, Oy, Oz) ray origin
    D: (Dx, Dy, Dz) ray direction, assumed normalized or not (doesn't matter)
    returns: smallest positive t at which R(t)=O + t D hits the cylinder
             x^2+z^2=1, 0<=y<=1 (including caps), or -1 if no hit.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []

    # 1) SIDE of infinite cylinder x^2 + z^2 = 1
    a = Dx*Dx + Dz*Dz
    if a != 0.0:  # otherwise ray is parallel to cylinder axis, no side‐hit
        b = 2*(Ox*Dx + Oz*Dz)
        c = Ox*Ox + Oz*Oz - 1
        disc = b*b - 4*a*c
        if disc >= 0.0:
            sqrt_disc = math.sqrt(disc)
            t0 = (-b - sqrt_disc) / (2*a)
            t1 = (-b + sqrt_disc) / (2*a)
            for t in (t0, t1):
                if t > 0:
                    y = Oy + t*Dy
                    if 0.0 <= y <= 1.0:
                        ts.append(t)

    # 2) CAP at y=0 plane
    if Dy != 0.0:
        t_cap0 = (0.0 - Oy)/Dy
        if t_cap0 > 0.0:
            x = Ox + t_cap0*Dx
            z = Oz + t_cap0*Dz
            if x*x + z*z <= 1.0:
                ts.append(t_cap0)

        # 3) Cap at y=1 plane
        t_cap1 = (1.0 - Oy)/Dy
        if t_cap1 > 0.0:
            x = Ox + t_cap1*Dx
            z = Oz + t_cap1*Dz
            if x*x + z*z <= 1.0:
                ts.append(t_cap1)

    if not ts:
        return -1.0

    return min(ts)

# Example usage:
if __name__ == '__main__':
    # a ray that starts at (2,0.5,0) pointing towards the cylinder
    O = (2.0, 0.5, 0.0)
    D = (-1.0, 0.0, 0.0)    # straight along −x
    t = intersect_cylinder(O, D)
    if t > 0:
        P = (O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2])
        print("Hit at t =", t, "point =", P)
    else:
        print("Missed!")


#---------- 3.3 ----------#
import math

def intersect_cylinder(O, D):
    """
    O : tuple of 3 floats, the ray origin (Ox, Oy, Oz)
    D : tuple of 3 floats, the ray direction (Dx, Dy, Dz)
    Returns the smallest t>0 such that O + t D hits the cylinder x^2+z^2=1, 0<=y<=1,
    including the caps y=0 and y=1. Returns -1 if no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    t_candidates = []

    # 1) Intersect infinite side (x^2 + z^2 = 1)
    #    Solve (Ox + t Dx)^2 + (Oz + t Dz)^2 = 1
    a = Dx*Dx + Dz*Dz
    b = 2*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1

    if abs(a) > 1e-8:  # ray not parallel to cylinder axis in xz-plane
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t0 = (-b - sqrt_disc) / (2*a)
            t1 = (-b + sqrt_disc) / (2*a)
            for t in (t0, t1):
                if t > 0:
                    y = Oy + t*Dy
                    if 0 <= y <= 1:
                        t_candidates.append(t)

    # 2) Intersect bottom cap y = 0
    if abs(Dy) > 1e-8:
        t_bot = (0 - Oy) / Dy
        if t_bot > 0:
            x_bot = Ox + t_bot*Dx
            z_bot = Oz + t_bot*Dz
            if x_bot*x_bot + z_bot*z_bot <= 1:
                t_candidates.append(t_bot)

    # 3) Intersect top cap y = 1
    if abs(Dy) > 1e-8:
        t_top = (1 - Oy) / Dy
        if t_top > 0:
            x_top = Ox + t_top*Dx
            z_top = Oz + t_top*Dz
            if x_top*x_top + z_top*z_top <= 1:
                t_candidates.append(t_top)

    if not t_candidates:
        return -1.0

    # return the nearest positive intersection
    return min(t_candidates)


#---------- 3.4 ----------#
import math

def intersect_cylinder(O, D, eps=1e-8):
    """
    O: ray origin, tuple of 3 floats (ox, oy, oz)
    D: ray direction, tuple of 3 floats (dx, dy, dz)
    Returns: smallest t >= 0 for which O + t D hits the cylinder,
             or -1 if no hit.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    ts = []

    # 1) Side intersection: solve (ox + t dx)^2 + (oz + t dz)^2 = 1
    a = dx*dx + dz*dz
    b = 2*(ox*dx + oz*dz)
    c = ox*ox + oz*oz - 1

    if abs(a) > eps:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b - sqrt_disc) / (2*a)
            t2 = (-b + sqrt_disc) / (2*a)
            for t in (t1, t2):
                if t >= 0:
                    y = oy + t*dy
                    if 0 <= y <= 1:
                        ts.append(t)

    # 2) Cap y=0
    if abs(dy) > eps:
        t_cap0 = (0 - oy) / dy
        if t_cap0 >= 0:
            x_at = ox + t_cap0*dx
            z_at = oz + t_cap0*dz
            if x_at*x_at + z_at*z_at <= 1 + eps:
                ts.append(t_cap0)

        # 3) Cap y=1
        t_cap1 = (1 - oy) / dy
        if t_cap1 >= 0:
            x_at = ox + t_cap1*dx
            z_at = oz + t_cap1*dz
            if x_at*x_at + z_at*z_at <= 1 + eps:
                ts.append(t_cap1)

    if not ts:
        return -1.0

    return min(ts)

# --- Example usage / quick test ---
if __name__ == "__main__":
    # Ray from (2,0.5,0) pointing towards negative x-axis
    O = (2.0, 0.5, 0.0)
    D = (-1.0, 0.0, 0.0)
    t = intersect_cylinder(O, D)
    if t > 0:
        P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
        print("Hit at t =", t, "point =", P)
    else:
        print("No intersection")


#---------- 3.5 ----------#
import math

def intersect_ray_cylinder(O, D):
    """
    Intersect a ray R(t) = O + t*D with the finite cylinder
      x^2 + z^2 = 1,   0 <= y <= 1
    Returns the smallest t>=0 at which the ray hits the cylinder (side or caps),
    or -1.0 if there is no intersection.
    
    O: (Ox, Oy, Oz) ray origin
    D: (Dx, Dy, Dz) ray direction
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []

    # 1) Intersect with the infinite cylinder x^2 + z^2 = 1
    #    (then we will check the y‐coordinate stay between 0 and 1)
    a = Dx*Dx + Dz*Dz
    b = 2.0*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1.0

    if abs(a) > 1e-9:
        disc = b*b - 4.0*a*c
        if disc >= 0.0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b - sqrt_disc) / (2.0 * a)
            t2 = (-b + sqrt_disc) / (2.0 * a)
            for t in (t1, t2):
                if t >= 0.0:
                    yhit = Oy + t*Dy
                    if 0.0 <= yhit <= 1.0:
                        ts.append(t)

    # 2) Intersect with the caps y=0 and y=1 (if the ray is not parallel to them)
    if abs(Dy) > 1e-9:
        # bottom cap at y=0
        t_bot = (0.0 - Oy) / Dy
        if t_bot >= 0.0:
            x_bot = Ox + t_bot*Dx
            z_bot = Oz + t_bot*Dz
            if x_bot*x_bot + z_bot*z_bot <= 1.0:
                ts.append(t_bot)
        # top cap at y=1
        t_top = (1.0 - Oy) / Dy
        if t_top >= 0.0:
            x_top = Ox + t_top*Dx
            z_top = Oz + t_top*Dz
            if x_top*x_top + z_top*z_top <= 1.0:
                ts.append(t_top)

    if not ts:
        return -1.0

    return min(ts)


#---------- 3.6 ----------#
import math

def intersect_cylinder(O, D):
    """
    O: tuple (ox, oy, oz)  -- ray origin
    D: tuple (dx, dy, dz)  -- ray direction
    Returns: t>0 distance to first intersection, or -1 if none.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    t_min = float('inf')

    # 1) SIDE of cylinder: x^2 + z^2 = 1, 0 <= y <= 1
    #    solve a t^2 + b t + c = 0
    a = dx*dx + dz*dz
    if abs(a) > 1e-8:
        b = 2 * (ox*dx + oz*dz)
        c = ox*ox + oz*oz - 1
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            # two possible solutions
            for t in [(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)]:
                if t > 1e-8:
                    y_at_t = oy + t*dy
                    if 0.0 <= y_at_t <= 1.0:
                        t_min = min(t_min, t)

    # 2) BOTTOM CAP y = 0  (x^2+z^2 <= 1)
    if abs(dy) > 1e-8:
        t_bot = (0.0 - oy) / dy
        if t_bot > 1e-8:
            x_bot = ox + t_bot*dx
            z_bot = oz + t_bot*dz
            if x_bot*x_bot + z_bot*z_bot <= 1.0:
                t_min = min(t_min, t_bot)

        # 3) TOP CAP y = 1
        t_top = (1.0 - oy) / dy
        if t_top > 1e-8:
            x_top = ox + t_top*dx
            z_top = oz + t_top*dz
            if x_top*x_top + z_top*z_top <= 1.0:
                t_min = min(t_min, t_top)

    if t_min < float('inf'):
        return t_min
    else:
        return -1.0


#---------- 3.7 ----------#
def intersect_ray_cylinder(O, D, eps=1e-8):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: smallest positive t at which R(t)=O+t*D hits the cylinder or its caps,
             or -1 if no such positive intersection exists.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []  # will collect all valid intersections

    # 1) Side intersection: solve (Ox + t*Dx)^2 + (Oz + t*Dz)^2 = 1
    a = Dx*Dx + Dz*Dz
    b = 2*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1

    if abs(a) > eps:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = disc**0.5
            t0 = (-b - sqrt_disc) / (2*a)
            t1 = (-b + sqrt_disc) / (2*a)
            # check the smaller t first
            for t in sorted((t0, t1)):
                if t > eps:
                    y = Oy + t*Dy
                    if 0.0 <= y <= 1.0:
                        ts.append(t)
                        break

    # 2) Cap intersections: planes y=0 and y=1
    if abs(Dy) > eps:
        # bottom cap at y=0
        t_bot = (0.0 - Oy) / Dy
        if t_bot > eps:
            x = Ox + t_bot*Dx
            z = Oz + t_bot*Dz
            if x*x + z*z <= 1.0:
                ts.append(t_bot)
        # top cap at y=1
        t_top = (1.0 - Oy) / Dy
        if t_top > eps:
            x = Ox + t_top*Dx
            z = Oz + t_top*Dz
            if x*x + z*z <= 1.0:
                ts.append(t_top)

    if not ts:
        return -1.0

    return min(ts)


#---------- 3.8 ----------#
import math

def intersect_ray_cylinder(O, D):
    """
    O = (ox, oy, oz)  ray origin
    D = (dx, dy, dz)  ray direction
    Cylinder: radius=1, axis = y-axis, y in [0,1]
    Returns smallest t>=0 or -1 if no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    ts = []

    # 1) Side surface: solve (ox + t dx)^2 + (oz + t dz)^2 = 1
    a = dx*dx + dz*dz
    b = 2*(ox*dx + oz*dz)
    c = ox*ox + oz*oz - 1

    if abs(a) > 1e-12:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b - sqrt_disc) / (2*a)
            t2 = (-b + sqrt_disc) / (2*a)
            for t in (t1, t2):
                if t >= 0:
                    y = oy + t*dy
                    if 0 <= y <= 1:
                        ts.append(t)
    # else a==0: ray is exactly parallel to cylinder axis in xz-plane,
    # it never hits the curved side (unless it starts exactly on it, but
    # we assume the ray starts outside)

    # 2) Bottom cap: plane y=0, disk x^2+z^2 <= 1
    if abs(dy) > 1e-12:
        t_bot = (0 - oy)/dy
        if t_bot >= 0:
            x_bot = ox + t_bot*dx
            z_bot = oz + t_bot*dz
            if x_bot*x_bot + z_bot*z_bot <= 1:
                ts.append(t_bot)

        # 3) Top cap: plane y=1
        t_top = (1 - oy)/dy
        if t_top >= 0:
            x_top = ox + t_top*dx
            z_top = oz + t_top*dz
            if x_top*x_top + z_top*z_top <= 1:
                ts.append(t_top)

    # pick the smallest non-negative t, if any
    if not ts:
        return -1.0
    return min(ts)

# Example sanity‐check:
if __name__ == "__main__":
    # A ray from (2,0.5,0) pointing toward the cylinder
    O = (2.0, 0.5, 0.0)
    D = (-1.0, 0.0, 0.0)
    print(intersect_ray_cylinder(O, D))   # should be 1.0, since at t=1: x=1,y=0.5,z=0

    # A ray from above pointing down onto the top cap
    O = (0.2, 2.0, 0.3)
    D = (0.0, -1.0, 0.0)
    print(intersect_ray_cylinder(O, D))   # should be 1.0, since it hits y=1 cap first

    # A ray missing entirely
    O = (2,2,2)
    D = (1,0,0)
    print(intersect_ray_cylinder(O, D))   # should be -1.0


#---------- 3.9 ----------#
import math

def ray_cylinder_intersection(O, D):
    """
    O: Ray origin, tuple (Ox, Oy, Oz)
    D: Ray direction, tuple (Dx, Dy, Dz)
    Returns: smallest t>0 at which R(t)=O+t*D hits the unit cylinder
             aligned with the Y-axis between y=0 and y=1, or -1 if no hit.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []

    # 1) Side‐wall intersection: solve (Ox + Dx t)^2 + (Oz + Dz t)^2 = 1
    a = Dx*Dx + Dz*Dz
    b = 2*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1

    if abs(a) > 1e-12:  # ray not parallel to cylinder's axis in xz‐plane
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b - sqrt_disc) / (2*a)
            t2 = (-b + sqrt_disc) / (2*a)
            for t in (t1, t2):
                if t > 0:
                    y = Oy + Dy*t
                    if 0.0 <= y <= 1.0:
                        ts.append(t)

    # 2) Bottom cap (y = 0) if ray not parallel to Y
    if abs(Dy) > 1e-12:
        t_bot = (0.0 - Oy) / Dy
        if t_bot > 0:
            x_bot = Ox + Dx*t_bot
            z_bot = Oz + Dz*t_bot
            if x_bot*x_bot + z_bot*z_bot <= 1.0:
                ts.append(t_bot)

        # 3) Top cap (y = 1)
        t_top = (1.0 - Oy) / Dy
        if t_top > 0:
            x_top = Ox + Dx*t_top
            z_top = Oz + Dz*t_top
            if x_top*x_top + z_top*z_top <= 1.0:
                ts.append(t_top)

    if not ts:
        return -1.0
    return min(ts)

# Example usage:

# Ray from (2,0.5,0) pointing towards the cylinder center:
O = (2.0, 0.5, 0.0)
D = (-1.0, 0.0, 0.0)  # already normalized in length=1 for simplicity
t = ray_cylinder_intersection(O, D)
if t > 0:
    hit_point = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", hit_point)
else:
    print("No intersection")


#---------- 4.0 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O: ray origin as a tuple (Ox, Oy, Oz)
    D: ray direction as a tuple (Dx, Dy, Dz)
    
    Returns the smallest positive t such that O + t*D
    lies on or inside either sphere, or -1.0 if no hit.
    """
    best_t = float('inf')
    # define the two spheres: (center, radius)
    spheres = [((0.0, 0.0, 0.0), 1.0),
               ((1.0, 0.0, 0.0), 1.0)]
    
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    
    for (Cx, Cy, Cz), R in spheres:
        # move into sphere's local coords
        ocx = Ox - Cx
        ocy = Oy - Cy
        ocz = Oz - Cz
        
        # quadratic coefficients a*t^2 + b*t + c = 0
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2.0 * (Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - R*R
        
        disc = b*b - 4.0*a*c
        if disc < 0.0:
            # no real roots → no intersection with this sphere
            continue
        
        sqrt_disc = math.sqrt(disc)
        # two possible solutions
        t0 = (-b - sqrt_disc) / (2.0*a)
        t1 = (-b + sqrt_disc) / (2.0*a)
        
        # we only care about the first time the ray hits (smallest non‐negative t)
        if t0 >= 0.0 and t0 < best_t:
            best_t = t0
        if t1 >= 0.0 and t1 < best_t:
            best_t = t1

    return best_t if best_t < float('inf') else -1.0


#---------- 4.1 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O: tuple of 3 floats = ray origin
    D: tuple of 3 floats = ray direction (should be normalized)
    returns: smallest t >= 0 so that O + t*D hits either sphere,
             or -1 if no hit.
    """
    centers = ((0.0,0.0,0.0), (1.0,0.0,0.0))
    radius = 1.0

    t_min = float('inf')

    # helper for dot product
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

    for C in centers:
        # compute quadratic coefficients for |O + t D - C|^2 = r^2
        OC = (O[0]-C[0], O[1]-C[1], O[2]-C[2])
        a = dot(D, D)               # = 1.0 if D is normalized
        b = 2.0 * dot(D, OC)
        c = dot(OC, OC) - radius*radius

        disc = b*b - 4*a*c
        if disc < 0.0:
            continue    # no real roots → no intersection with this sphere

        sqrt_disc = math.sqrt(disc)
        # two solutions of the quadratic
        t1 = (-b - sqrt_disc) / (2.0*a)
        t2 = (-b + sqrt_disc) / (2.0*a)

        # we only care about t >= 0
        if t1 >= 0.0 and t1 < t_min:
            t_min = t1
        if t2 >= 0.0 and t2 < t_min:
            t_min = t2

    if t_min == float('inf'):
        return -1.0
    return t_min


#---------- 4.2 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O: tuple of 3 floats, ray origin
    D: tuple of 3 floats, ray direction (need not be normalized)
    Returns the smallest t>=0 such that O + t*D hits either sphere,
    or -1 if no intersection.
    """
    best_t = float('inf')
    
    # centers of the two spheres
    centers = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    r = 1.0
    
    # precompute a = D·D
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    
    for C in centers:
        # vector from sphere center C to ray origin O
        OCx = O[0] - C[0]
        OCy = O[1] - C[1]
        OCz = O[2] - C[2]
        
        # coefficients of quadratic a t^2 + b t + c = 0
        b = 2.0 * (D[0]*OCx + D[1]*OCy + D[2]*OCz)
        c = OCx*OCx + OCy*OCy + OCz*OCz - r*r
        
        disc = b*b - 4*a*c
        if disc < 0:
            # no real intersection with this sphere
            continue
        
        sqrt_disc = math.sqrt(disc)
        # two solutions
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        
        # we only care about intersections in front of the ray origin (t >= 0)
        # and we want the closest one
        if t0 >= 0 and t0 < best_t:
            best_t = t0
        if t1 >= 0 and t1 < best_t:
            best_t = t1
    
    if best_t < float('inf'):
        return best_t
    else:
        return -1.0


#---------- 4.3 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O: tuple of 3 floats, the ray origin (Ox, Oy, Oz)
    D: tuple of 3 floats, the ray direction (Dx, Dy, Dz)
       (does not have to be unit length)
    Returns: the distance t ≥ 0 to the first intersection with
             the union of the two unit spheres centered at
             (0,0,0) and (1,0,0), or -1 if there's no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    t_hits = []

    # List of sphere centers
    centers = [(0.0, 0.0, 0.0),
               (1.0, 0.0, 0.0)]
    radius = 1.0

    # Precompute a = D·D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    if a == 0.0:
        return -1.0    # no valid direction

    for (Cx, Cy, Cz) in centers:
        # Solve   ||O + t D − C||^2 = r^2
        # which gives a t^2 + b t + c = 0
        Lx = Ox - Cx
        Ly = Oy - Cy
        Lz = Oz - Cz

        b = 2.0 * (Dx*Lx + Dy*Ly + Dz*Lz)
        c = (Lx*Lx + Ly*Ly + Lz*Lz) - radius*radius

        disc = b*b - 4.0*a*c
        if disc < 0.0:
            continue   # no real roots → miss this sphere

        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2.0*a)
        t1 = (-b + sqrt_disc) / (2.0*a)

        # t0 ≤ t1 by construction.  We want the first positive hit.
        if t0 >= 0.0:
            t_hits.append(t0)
        elif t1 >= 0.0:
            # t0 < 0 < t1  → ray origin inside sphere but they guaranteed
            #                we start outside the object, so this case
            #                shouldn't normally happen.  We still handle it:
            t_hits.append(t1)

    if not t_hits:
        return -1.0

    # Return the closest non-negative intersection
    return min(t_hits)


#---------- 4.4 ----------#
import math

def dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def sub(u, v):
    return (u[0]-v[0], u[1]-v[1], u[2]-v[2])

def ray_sphere_interval(O, D, C, r):
    """
    Solve (O + t D - C)^2 = r^2 for t.
    Return None if no real intersections,
    otherwise return (t1, t2) with t1 <= t2.
    """
    # Compute quadratic coefficients
    OC = sub(O, C)
    a = dot(D, D)
    b = 2 * dot(D, OC)
    c = dot(OC, OC) - r*r

    disc = b*b - 4*a*c
    if disc < 0:
        return None
    elif abs(disc) < 1e-8:
        t = -b / (2*a)
        return (t, t)
    else:
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        if t1 > t2:
            t1, t2 = t2, t1
        return (t1, t2)

def intersect_dual_sphere(O, D):
    """
    O, D are 3‐tuples.
    Returns the smallest t >= 0 where the ray O + t D
    hits the union of the two unit spheres centered at
    (0,0,0) and (1,0,0).  Returns -1 if no hit.
    """
    # Define our two spheres
    spheres = [
        ((0.0, 0.0, 0.0), 1.0),
        ((1.0, 0.0, 0.0), 1.0),
    ]

    intervals = []
    for C, r in spheres:
        iv = ray_sphere_interval(O, D, C, r)
        if iv is not None:
            t0, t1 = iv
            # discard intervals that lie completely behind the ray
            if t1 < 0:
                continue
            # clamp entry to zero if it straddles the origin
            t0 = max(t0, 0.0)
            intervals.append((t0, t1))

    if not intervals:
        return -1.0

    # sort by entry time
    intervals.sort(key=lambda x: x[0])

    # merge overlapping intervals
    merged = []
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end + 1e-8:    # overlap or touch
            cur_end = max(cur_end, e)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = s, e
    merged.append((cur_start, cur_end))

    # the ray starts outside the union, so the first interval
    # in merged[] is the entry into the object
    entry, exit = merged[0]
    return entry

# --- Example of use ------------------------------------

# Ray origin and direction
O = ( -5.0,  0.0,  0.0 )
D = (  1.0,  0.0,  0.0 )

t = intersect_dual_sphere(O, D)
if t >= 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point", P)
else:
    print("No intersection")


#---------- 4.5 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O: ray origin, tuple of 3 floats
    D: ray direction, tuple of 3 floats
    returns: t >= 0 of first intersection, or -1.0 if none
    """
    def intersect_sphere(C, r=1.0):
        """
        Solve ||O + t D − C||^2 = r^2.
        Returns a list of 0, 1 or 2 t‐values (may be negative).
        """
        # oc = O − C
        ocx = O[0] - C[0]
        ocy = O[1] - C[1]
        ocz = O[2] - C[2]

        # quadratic coefficients a t^2 + b t + c = 0
        a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
        b = 2.0 * (D[0]*ocx + D[1]*ocy + D[2]*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - r*r

        disc = b*b - 4.0*a*c
        if disc < 0.0:
            return []            # no intersection
        if abs(disc) < 1e-8:
            # grazing hit (one repeated root)
            t = -b / (2.0*a)
            return [t]
        # two distinct intersections
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2.0*a)
        t2 = (-b + sqrt_disc) / (2.0*a)
        return [t1, t2]

    # our two sphere centers
    centers = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    ts = []

    # collect all non‐negative intersections
    for C in centers:
        for t in intersect_sphere(C, r=1.0):
            if t >= 0.0:
                ts.append(t)

    if not ts:
        return -1.0   # no hit
    return min(ts)    # first intersection

# Example usage:
O = (-5.0, 0.0, 0.0)
D = (1.0, 0.0, 0.0)
print(intersect_dual_sphere(O, D))   # should print 4.0, since the left‐hand sphere at (0,0,0) is hit at t=5−1=4


#---------- 4.6 ----------#
import math

def _intersect_sphere(O, D, C, r):
    """
    Solve |O + t D - C|^2 = r^2 for t.
    Returns a sorted list of all t >= 0 (possibly empty).
    """
    ox, oy, oz = O
    dx, dy, dz = D
    cx, cy, cz = C

    # Shift ray origin into sphere's local space
    ocx = ox - cx
    ocy = oy - cy
    ocz = oz - cz

    a = dx*dx + dy*dy + dz*dz
    b = 2*(dx*ocx + dy*ocy + dz*ocz)
    c = ocx*ocx + ocy*ocy + ocz*ocz - r*r

    disc = b*b - 4*a*c
    if disc < 0:
        return []        # no real roots → no intersection

    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)

    ts = []
    if t0 >= 0:
        ts.append(t0)
    if t1 >= 0 and disc > 0:
        # If disc==0, t0==t1; we only append once.
        ts.append(t1)

    ts.sort()
    return ts


def intersect_dual_sphere(O, D):
    """
    O: ray origin, tuple of 3 floats
    D: ray direction, tuple of 3 floats (need not be unit length)
    Returns the smallest t>=0 at which R(t)=O+tD hits the union
    of
      Sphere 1: center=(0,0,0), r=1
      Sphere 2: center=(1,0,0), r=1
    or -1 if there is no intersection.
    """
    # sphere #1
    ts1 = _intersect_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    # sphere #2
    ts2 = _intersect_sphere(O, D, (1.0, 0.0, 0.0), 1.0)

    # collect all positive hits
    all_ts = ts1 + ts2
    if not all_ts:
        return -1.0
    return min(all_ts)


# --- Example usage ---
if __name__ == '__main__':
    # shoot a ray from (–3,0,0) toward +x
    O = (-3.0, 0.0, 0.0)
    D = (1.0, 0.0, 0.0)
    t = intersect_dual_sphere(O, D)
    print("first hit at t =", t)   # should print t = 2.0 (first sphere at x=-1)


#---------- 4.7 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O, D: tuples of 3 floats each, ray origin and direction.
    Returns the smallest t >= 0 such that O + t*D hits either sphere,
    or -1 if there is no intersection.
    """
    centers = [(0.0, 0.0, 0.0),
               (1.0, 0.0, 0.0)]
    r = 1.0

    t_hit = float('inf')

    # helper: dot product of 3‐vectors
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

    for C in centers:
        # compute coefficients of a*t^2 + b*t + c = 0
        # where (O + t D - C)·(O + t D - C) = r^2
        OC = (O[0]-C[0], O[1]-C[1], O[2]-C[2])
        a = dot(D, D)
        b = 2.0 * dot(D, OC)
        c = dot(OC, OC) - r*r

        disc = b*b - 4.0*a*c
        if disc < 0.0:
            # no real roots, ray misses this sphere
            continue

        sqrt_disc = math.sqrt(disc)
        # two roots
        t1 = (-b - sqrt_disc) / (2.0*a)
        t2 = (-b + sqrt_disc) / (2.0*a)

        # we want the first time the ray actually hits (the smaller non‐negative t)
        if t1 >= 0.0:
            t_candidate = t1
        elif t2 >= 0.0:
            t_candidate = t2
        else:
            # both intersections are "behind" the ray origin
            continue

        # keep the closest positive hit
        if t_candidate < t_hit:
            t_hit = t_candidate

    if t_hit == float('inf'):
        return -1.0
    return t_hit


#---------- 4.8 ----------#
import math

def dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def subtract(u, v):
    return (u[0]-v[0], u[1]-v[1], u[2]-v[2])

def intersect_dual_sphere(O, D):
    """
    O: ray origin tuple (ox, oy, oz)
    D: ray direction tuple (dx, dy, dz) – need not be normalized
    returns: smallest t>0 such that O + t*D hits either sphere, or -1 if none
    """
    centers = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    radius = 1.0
    best_t = float('inf')

    a = dot(D, D)
    # if D were normalized, a == 1 and you could skip it

    for C in centers:
        # compute coefficients for |O + t D - C|^2 = radius^2
        OC = subtract(O, C)
        b = 2.0 * dot(OC, D)
        c = dot(OC, OC) - radius*radius

        disc = b*b - 4*a*c
        if disc < 0.0:
            # no real roots → no intersection
            continue

        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)

        # we only care about intersections in front of the ray (t > 0)
        if t1 > 0 and t1 < best_t:
            best_t = t1
        if t2 > 0 and t2 < best_t:
            best_t = t2

    return best_t if best_t < float('inf') else -1.0

# Example usage:
O = (0.0, 0.0, -5.0)
D = (0.0, 0.0, 1.0)
t = intersect_dual_sphere(O, D)
if t > 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", P)
else:
    print("No intersection")


#---------- 4.9 ----------#
import math

def intersect_dual_sphere(O, D):
    """
    O, D are 3‐tuples or lists of floats (the ray origin and direction).
    Returns the smallest t >= 0 such that O + t D hits either sphere,
    or -1 if there is no intersection.
    """

    # sphere data: (center_x, center_y, center_z, radius)
    spheres = [
        (0.0, 0.0, 0.0, 1.0),
        (1.0, 0.0, 0.0, 1.0)
    ]

    best_t = float('inf')

    # precompute D·D
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]

    for cx, cy, cz, r in spheres:
        # vector from center to ray origin
        Lx = O[0] - cx
        Ly = O[1] - cy
        Lz = O[2] - cz

        b = 2.0 * (D[0]*Lx + D[1]*Ly + D[2]*Lz)
        c = (Lx*Lx + Ly*Ly + Lz*Lz) - r*r

        disc = b*b - 4.0*a*c
        if disc < 0.0:
            # no real roots, miss
            continue

        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2.0*a)
        t1 = (-b + sqrt_disc) / (2.0*a)

        # we only want t >= 0, and the closer one if both are >= 0
        if t0 >= 0.0 and t0 < best_t:
            best_t = t0
        if t1 >= 0.0 and t1 < best_t:
            best_t = t1

    if best_t == float('inf'):
        return -1.0
    return best_t


#---------- 5.0 ----------#
import math

def intersect_disc(O, D):
    """
    Intersect a ray R(t) = O + t D with the disc of radius 1
    centered at (0,0,0) whose (unnormalized) normal is n = (1,1,1).
    Return the smallest t>0 at which the ray hits the disc, or -1
    if there is no intersection.
    
    O, D : 3‐tuples of floats
    """
    # Disc parameters
    n = (1.0, 1.0, 1.0)   # plane normal (doesn't need to be unit length here)
    radius = 1.0
    
    # 1) Compute denominator of plane‐intersection formula
    denom = n[0]*D[0] + n[1]*D[1] + n[2]*D[2]
    if abs(denom) < 1e-8:
        # Ray is parallel to the disc’s plane
        return -1.0
    
    # 2) Compute t at which ray hits the plane
    #    Plane passes through C=(0,0,0) so plane eq is n·X = 0
    #    We solve n·(O + t D) = 0  =>  t = -n·O  /  (n·D)
    t = - (n[0]*O[0] + n[1]*O[1] + n[2]*O[2]) / denom
    
    if t < 0.0:
        # Intersection is behind the ray origin
        return -1.0
    
    # 3) Compute the hit point
    Px = O[0] + t * D[0]
    Py = O[1] + t * D[1]
    Pz = O[2] + t * D[2]
    
    # 4) Check if that point lies within the disc’s radius
    if (Px*Px + Py*Py + Pz*Pz) <= radius*radius:
        return t
    else:
        return -1.0

# Example usage:
O = (2.0, 2.0, -1.0)
D = (-1.0, -1.0, 1.0)
# (you’d normally normalize D before tracing)
length = math.sqrt(D[0]*D[0] + D[1]*D[1] + D[2]*D[2])
D = (D[0]/length, D[1]/length, D[2]/length)

t_hit = intersect_disc(O, D)
if t_hit > 0:
    P = (O[0] + t_hit*D[0],
         O[1] + t_hit*D[1],
         O[2] + t_hit*D[2])
    print("Hit at t =", t_hit, "point =", P)
else:
    print("No intersection")


#---------- 5.1 ----------#
def intersect_ray_disc(O, D, eps=1e-6):
    """
    Compute intersection of ray R(t) = O + t*D with the disc
    centered at (0,0,0), radius = 1, normal = (1,1,1).
    O, D are 3‐tuples or lists of floats.  Returns the
    nearest t >= 0 if there is an intersection, otherwise -1.
    """
    # Unnormalized normal of the disc
    N = (1.0, 1.0, 1.0)

    # Dot products helper
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    # 1) Ray‐plane intersection: denom = D·N
    denom = dot(D, N)
    if abs(denom) < eps:
        # Ray is parallel to plane
        return -1.0

    # 2) Compute t for intersection with the plane
    # Since plane passes through C=(0,0,0), (C−O)·N = -O·N
    t = -dot(O, N) / denom
    if t < 0:
        # Intersection is behind the ray origin
        return -1.0

    # 3) Compute the intersection point P = O + t*D
    P = (O[0] + t * D[0],
         O[1] + t * D[1],
         O[2] + t * D[2])

    # 4) Check if P lies inside the disc of radius 1
    # (since center is at the origin, we just test P·P <= 1)
    if dot(P, P) <= 1.0 + eps:
        return t
    else:
        return -1.0


#---------- 5.2 ----------#
import math

def intersect_disc(O, D):
    """
    Compute the intersection of the ray R(t)=O + t D with the unit‐radius disc
    centered at (0,0,0) whose (unnormalized) normal is N=(1,1,1).

    Inputs:
      O, D : 3‐tuples or lists of floats  (ray origin and ray direction)
    Returns:
      t    : the smallest non‐negative t for which R(t) hits the disc,
             or –1.0 if there is no intersection.
    """
    # disc plane normal (doesn't need to be unit length for the formula below)
    N = (1.0, 1.0, 1.0)

    # 1) Ray–plane intersection: find t so that N·(O + t D) = 0  ⇔  N·O + t (N·D) = 0
    denom = N[0]*D[0] + N[1]*D[1] + N[2]*D[2]
    if abs(denom) < 1e-8:
        # Ray is (nearly) parallel to the plane → no intersection
        return -1.0

    t = - (N[0]*O[0] + N[1]*O[1] + N[2]*O[2]) / denom
    if t < 0.0:
        # Intersection is “behind” the ray origin
        return -1.0

    # 2) Compute the intersection point
    Px = O[0] + t * D[0]
    Py = O[1] + t * D[1]
    Pz = O[2] + t * D[2]

    # 3) Check if it lies within the unit‐radius disc (centered at 0,0,0)
    if (Px*Px + Py*Py + Pz*Pz) <= 1.0 + 1e-8:
        return t
    else:
        return -1.0


#---------- 5.3 ----------#
def intersect_disc(O, D, 
                   C=(0.0, 0.0, 0.0),      # disc center
                   N=(1.0, 1.0, 1.0),      # disc normal
                   radius=1.0,
                   eps=1e-8):
    """
    O: tuple of 3 floats = ray origin
    D: tuple of 3 floats = ray direction
    C: disc center (default (0,0,0))
    N: disc normal    (default (1,1,1))
    radius: disc radius (default 1)
    returns: t >= 0 at intersection or -1.0 if no hit
    """
    # 1) Compute denom = N·D
    denom = N[0]*D[0] + N[1]*D[1] + N[2]*D[2]
    # if denom == 0, ray is parallel to the plane
    if abs(denom) < eps:
        return -1.0

    # 2) Compute t for intersection with the plane N·(X - C) = 0
    #    plane eq: N·X + d = 0, here d = -N·C.  Since C=(0,0,0), d=0.
    #    t = - (N·O + d) / (N·D) = - (N·O) / denom
    num = N[0]*O[0] + N[1]*O[1] + N[2]*O[2]
    t = - num / denom

    # if intersection is "behind" the ray origin, discard
    if t < 0.0:
        return -1.0

    # 3) Compute the point P = O + t*D
    Px = O[0] + t * D[0]
    Py = O[1] + t * D[1]
    Pz = O[2] + t * D[2]

    # 4) Check if P lies within the disc radius
    #    distance^2 from P to disc center C
    dx = Px - C[0]
    dy = Py - C[1]
    dz = Pz - C[2]
    if dx*dx + dy*dy + dz*dz <= radius*radius + eps:
        return t
    else:
        return -1.0


# Example usage:

# Ray coming straight down the normal, from (0,0,5) toward the disc
O = (0.0, 0.0, 5.0)
D = (1.0, 1.0, -5.0)   # not unit length, that's fine
t = intersect_disc(O, D)
if t >= 0:
    print("hit at t =", t, "point =", (O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2]))
else:
    print("no intersection")


#---------- 5.4 ----------#
import math

def ray_disc_intersect(O, D):
    """
    O: ray origin, 3‐tuple (x0,y0,z0)
    D: ray direction, 3‐tuple (dx,dy,dz)
    disc: centered at (0,0,0), radius = 1, normal = (1,1,1)
    returns: smallest t ≥ 0 so that O + t D hits the disc, or −1 if no hit
    """
    # 1) set up the plane x + y + z = 0  (normal n = (1,1,1), disc center at origin)
    n = (1.0, 1.0, 1.0)

    # 2) compute denominator = n·D
    denom = n[0]*D[0] + n[1]*D[1] + n[2]*D[2]
    if abs(denom) < 1e-9:
        # ray is parallel to the plane
        return -1.0

    # 3) solve for t so that n · (O + t D) = 0  ->  t = - (n·O) / (n·D)
    numer = n[0]*O[0] + n[1]*O[1] + n[2]*O[2]
    t = - numer / denom

    # 4) we only want intersections in front of the ray origin
    if t < 0:
        return -1.0

    # 5) compute the hit point P = O + t D
    Px = O[0] + t*D[0]
    Py = O[1] + t*D[1]
    Pz = O[2] + t*D[2]

    # 6) check whether P lies within radius = 1 of the disc center
    if math.sqrt(Px*Px + Py*Py + Pz*Pz) <= 1.0 + 1e-9:
        return t
    else:
        return -1.0


#---------- 5.5 ----------#
def intersect_ray_disc(O, D):
    """
    O, D: 3‐tuples of floats
    Returns the smallest t >= 0 such that O + t*D hits the disc
    centered at (0,0,0), radius=1, normal=(1,1,1), or -1 if no hit.
    """
    C = (0.0, 0.0, 0.0)
    N = (1.0, 1.0, 1.0)

    # 1) Ray‐plane intersection
    denom = D[0]*N[0] + D[1]*N[1] + D[2]*N[2]
    if abs(denom) < 1e-8:
        return -1.0   # Ray is parallel to the plane

    # N·(C - O)
    num = N[0]*(C[0] - O[0]) + N[1]*(C[1] - O[1]) + N[2]*(C[2] - O[2])
    t = num / denom
    if t < 0.0:
        return -1.0   # Intersection is behind the ray origin

    # 2) Check if inside the disc radius
    Px = O[0] + t*D[0]
    Py = O[1] + t*D[1]
    Pz = O[2] + t*D[2]

    # distance squared from center
    if Px*Px + Py*Py + Pz*Pz <= 1.0:
        return t
    else:
        return -1.0


#---------- 5.6 ----------#
import math

def intersect_disc(O, D, eps=1e-8):
    """
    Ray: R(t) = O + t * D,  t >= 0
    Disc: center C = (0,0,0), radius = 1, normal n = (1,1,1).
    Returns the smallest t >= 0 at which the ray hits the disc, or -1 if none.
    """
    # unpack
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # plane: n·X = 0,  n = (1,1,1)
    # solve for t: n·(O + t D) = 0  =>  t (n·D) = -n·O
    n_dot_D = Dx + Dy + Dz
    n_dot_O = Ox + Oy + Oz

    # if ray is parallel to plane, no (proper) intersection
    if abs(n_dot_D) < eps:
        return -1.0

    t = - n_dot_O / n_dot_D

    # we only want intersections in front of O
    if t < 0:
        return -1.0

    # compute the hit point
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz

    # test if within radius 1 of the disc center
    if (Px*Px + Py*Py + Pz*Pz) <= 1.0 + eps:
        return t
    else:
        return -1.0


# -- quick sanity checks --

# 1) Ray straight down the normal from above:
O1 = (0.0, 0.0, 5.0)
D1 = (1.0, 1.0, 1.0)   # not normalized, but that's OK
# This ray goes *away* from the plane (n·D > 0), so no hit
print(intersect_disc(O1, D1))   # should be -1

# 2) Ray toward the plane:
O2 = (0.0, 0.0, 5.0)
D2 = (-1.0, -1.0, -1.0)
# Should hit at t = 5 / 3  (since n·O=5, n·D=-3 → t= -5/(-3) = 5/3)
print(intersect_disc(O2, D2))   # ≈ 1.666...

# 3) Grazing ray:
# start at (1,0,5), direction toward the plane so that it lands at (1,0,0)
O3 = (1.0, 0.0, 5.0)
D3 = (0.0, 0.0, -5.0)
# hits exactly at the edge of the unit circle
print(intersect_disc(O3, D3))   # = 1.0

# 4) Ray that misses laterally:
O4 = (2.0, 0.0, 5.0)
D4 = (0.0, 0.0, -1.0)
print(intersect_disc(O4, D4))   # = -1


#---------- 5.7 ----------#
def intersect_disc(O, D):
    """
    O: tuple of 3 floats, ray origin
    D: tuple of 3 floats, ray direction (should be normalized)
    returns: the t > 0 of the first intersection, or -1 if none
    """
    # Disc data
    C = (0.0, 0.0, 0.0)     # center
    N = (1.0, 1.0, 1.0)     # normal (not yet unit‐length, but that's OK)
    R = 1.0                # radius

    # 1) Compute denom = D ⋅ N
    denom = D[0]*N[0] + D[1]*N[1] + D[2]*N[2]
    # if denom == 0 the ray is parallel to the plane → no intersection
    if abs(denom) < 1e-8:
        return -1.0

    # 2) Compute t such that O + tD lies in the plane dot(N, X–C) = 0
    #    plane constant d = dot(N,C)=0 since C=(0,0,0)
    #    so t = – dot(N, O–C) / dot(N, D)
    t = - (N[0]*O[0] + N[1]*O[1] + N[2]*O[2]) / denom

    # 3) We only care about intersections in front of the ray origin
    if t < 0.0:
        return -1.0

    # 4) Compute the intersection point P = O + t D
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 5) Check if P lies within the disc radius in‐plane.
    #    Since P is guaranteed on the plane, its full 3D distance from the center
    #    is the same as its in‐plane distance.
    if (P[0]*P[0] + P[1]*P[1] + P[2]*P[2]) <= R*R + 1e-8:
        return t

    # otherwise it hit the infinite plane but outside the disc
    return -1.0


#---------- 5.8 ----------#
import math

def normalize(v):
    l = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return (v[0]/l, v[1]/l, v[2]/l)

D = normalize(D)


#---------- 5.9 ----------#
import math

def intersect_ray_disc(O, D):
    """
    O : tuple of 3 floats = ray origin
    D : tuple of 3 floats = ray direction (doesn't have to be normalized)
    returns t >= 0 of the nearest intersection, or -1 if none.
    """
    # Disc definition
    P0 = (0.0, 0.0, 0.0)
    radius = 1.0
    N = (1.0, 1.0, 1.0)   # plane normal (not yet unit length)

    # helper: dot product
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

    # 1) intersect ray with the plane of the disc
    denom = dot(D, N)
    if abs(denom) < 1e-8:
        # Ray is parallel to the plane (no single‐point intersection)
        return -1.0

    numer = dot((O[0]-P0[0], O[1]-P0[1], O[2]-P0[2]), N)
    t = - numer / denom

    # 2) we only care about intersections in front of the origin
    if t < 0:
        return -1.0

    # 3) compute the intersection point Q = O + t D
    Q = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 4) check if Q lies inside the radius of the disc
    #    i.e. |Q - P0| <= radius
    dist2 = Q[0]**2 + Q[1]**2 + Q[2]**2  # since P0 is (0,0,0)
    if dist2 <= radius*radius:
        return t
    else:
        return -1.0


#---------- 6.0 ----------#
import math

def dot(a,b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def sub(a,b):   return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def cross(a,b):
    return ( a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] )

def intersect_quad_with_hole(O, D):
    # Plane normal for z=x is n = (−1,0,1)
    n = (-1.0, 0.0, 1.0)
    dn = dot(D, n)
    if abs(dn) < 1e-9:
        return -1.0      # Ray parallel to the quad‐plane

    # plane passes through (0,0,0), so plane offset is zero
    t = -dot(O, n) / dn
    if t < 0:
        return -1.0      # Intersection is behind the origin

    # Candidate intersection point
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 1) Test point‐in‐quad via edge‐cross tests
    verts = [(1,1,1), (-1,1,-1), (-2,-1,-2), (2,-1,2)]
    for i in range(4):
        V0 = verts[i]
        V1 = verts[(i+1)%4]
        E  = sub(V1, V0)
        W  = sub(P, V0)
        C  = cross(E, W)
        if dot(C, n) < 0:
            return -1.0

    # 2) Test the circular hole of radius 1 in that same plane
    #    Basis of the plane: u=(1,0,1)/√2, v=(0,1,0)
    #    Coordinates of P in that plane:
    a = (P[0] + P[2]) / math.sqrt(2.0)
    b = P[1]
    if a*a + b*b <= 1.0:
        return -1.0   # Ray went through the hole

    # If we get here, t is the distance to the first visible intersection
    return t


#---------- 6.1 ----------#
def intersect_quad_with_hole(O, D, eps=1e-9):
    """
    O, D : tuples of 3 floats, the ray origin and (normalized or not) direction
    Returns: the smallest t>0 such that O + t*D hits the quad outside the hole,
             or -1 if no such intersection exists.
    """
    # 1) Set up the four corners of the trapezoid (in CCW order as seen from the front):
    v0 = ( 1.0,  1.0,  1.0)
    v1 = (-1.0,  1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = ( 2.0, -1.0,  2.0)
    
    # little helper‐functions for vector ops
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def cross(u, v):
        return (u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0])
    def add(u, v):
        return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
    def mul_scalar(u, s):
        return (u[0]*s, u[1]*s, u[2]*s)
    
    # 2) Compute plane normal (we don't even have to normalize it)
    N = cross(sub(v1, v0), sub(v3, v0))
    
    # 3) Ray‐plane intersection test:
    denom = dot(D, N)
    if abs(denom) < eps:
        return -1.0   # ray is parallel to the plane
    t = dot(sub(v0, O), N) / denom
    if t < eps:
        return -1.0   # intersection is behind the ray origin (or too close)
    
    # 4) Compute the hit‐point
    P = add(O, mul_scalar(D, t))
    
    # 5) Point‐in‐convex‐quad test via edge‐cross‐product checks
    #    For a convex polygon in CCW order, (edge × (P – edge_start))·N must be >= 0
    def on_same_side(a, b):
        # returns True if P is on the same side of edge a→b as the normal N
        return dot(cross(sub(b, a), sub(P, a)), N) >= -eps
    
    if not (on_same_side(v0, v1) and
            on_same_side(v1, v2) and
            on_same_side(v2, v3) and
            on_same_side(v3, v0)):
        return -1.0   # outside the quad
    
    # 6) Circular hole at (0,0,0), radius=1. Reject if strictly inside r<1
    #    Points exactly at r=1 are on the hole boundary, which we treat as an intersection.
    r2 = P[0]*P[0] + P[1]*P[1] + P[2]*P[2]
    if r2 < 1.0 - eps:
        return -1.0
    
    # 7) This is the first valid intersection
    return t


# --- example usage -----------------------------------------------
# Ray origin
O = (0.0, 0.0, -5.0)
# Ray direction (should be normalized if you care about t in world units)
D = (0.0, 0.0, 1.0)

t = intersect_quad_with_hole(O, D)
if t > 0:
    print("Hit at t =", t, "point =", (O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2]))
else:
    print("No hit")


#---------- 6.2 ----------#
import math

# small epsilon for floating‐point comparisons
EPS = 1e-8

# Vector operations on 3‐tuples
def vadd(a,b):     return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def vsub(a,b):     return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def vdot(a,b):     return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def vcross(a,b):   return (
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
)
def vscale(a,s):   return (a[0]*s, a[1]*s, a[2]*s)
def vlen2(a):      return vdot(a,a)


def intersect_quad_with_circular_hole(O, D):
    """
    O: ray origin,  tuple of 3 floats
    D: ray direction (should be normalized for correct t), tuple of 3 floats
    Returns: the smallest t>0 at which the ray hits the QUAD (but not the hole),
             or -1 if no such intersection occurs.
    """
    # --- 1) define your quad vertices in CCW order (as seen from the "front" side) ---
    #    (these are the problem's four corners)
    v0 = ( 1.0,  1.0,  1.0)
    v1 = (-1.0, -1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = ( 2.0, -1.0,  2.0)
    verts = [v0, v1, v2, v3]

    # --- 2) compute plane normal with any two edges of the quad ---
    #    note: CCW winding means the normal points "outward" by right-hand rule
    edge1 = vsub(v1, v0)
    edge2 = vsub(v3, v0)
    N = vcross(edge1, edge2)

    # --- 3) intersect ray with plane: solve (O + t D - v0)·N = 0  =>  t = (v0 - O)·N / (D·N)
    denom = vdot(D, N)
    if abs(denom) < EPS:
        # Ray is parallel to the plane → no intersection
        return -1.0

    t = vdot(vsub(v0, O), N) / denom
    if t < EPS:
        # Intersection is behind the ray origin (or too close) → ignore
        return -1.0

    # Compute the hit point
    P = vadd(O, vscale(D, t))

    # --- 4) point‐in‐quad test by checking sign on all four edges ---
    #    for each edge vi→v(i+1), the vector from vi to P
    #    must lie to the "left" of edge (since winding is CCW),
    #    i.e. cross(edge, P-vi)·N >= 0
    for i in range(4):
        vi   = verts[i]
        vj   = verts[(i+1) % 4]
        edge = vsub(vj, vi)
        vp   = vsub(P, vi)
        if vdot(vcross(edge, vp), N) < -EPS:
            # point is outside this half‐space → outside quad
            return -1.0

    # --- 5) circular hole test: center at C=(0,0,0), radius=1.0 ---
    #    if P lies *inside* that disk, we consider it a “hole” and let the ray pass
    C = (0.0, 0.0, 0.0)
    if vlen2(vsub(P, C)) < (1.0 - EPS)**2:
        # inside the hole → no intersection with the quad
        return -1.0

    # passed all tests → this is our intersection
    return t


# --- example usage ---
if __name__ == "__main__":
    # a test ray
    O = (0.0, 0.0, -5.0)     # origin
    D = (0.0, 0.0,  1.0)     # toward +z
    t_hit = intersect_quad_with_circular_hole(O, D)
    if t_hit > 0:
        print("Hit at t =", t_hit, "point =", vadd(O, vscale(D, t_hit)))
    else:
        print("No intersection")


#---------- 6.3 ----------#
import math

def intersect_quad_with_hole(O, D):
    """
    O, D: each a 3‐tuple of floats, the ray origin and (normalized) direction.
    Returns the smallest t>0 such that O + t*D lies on the quad but outside
    the circular hole of radius 1 at the plane origin.  Returns -1 if none.
    """

    # 1) Define the quad vertices in CCW order (as seen from the "front" side):
    v0 = ( 1.0,  1.0,  1.0)
    v1 = (-1.0,  1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = ( 2.0, -1.0,  2.0)

    # helper vector operations
    def sub(a,b):   return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a,b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a,b): return ( a[1]*b[2]-a[2]*b[1],
                              a[2]*b[0]-a[0]*b[2],
                              a[0]*b[1]-a[1]*b[0] )
    def add(a,b):   return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def mul(a,s):   return (a[0]*s, a[1]*s, a[2]*s)

    # 2) Compute the plane normal (not necessarily unit length)
    e1 = sub(v1, v0)
    e3 = sub(v3, v0)
    N  = cross(e1, e3)

    # 3) Intersect ray with plane: N·(O + t D − v0) = 0  ⇒  t = N·(v0−O) / N·D
    denom = dot(N, D)
    eps   = 1e-8
    if abs(denom) < eps:
        return -1.0    # ray is parallel to the quad's plane

    t = dot(N, sub(v0, O)) / denom
    if t <= 0:
        return -1.0    # intersection is behind the ray origin

    # 4) Compute the intersection point
    P = add(O, mul(D, t))

    # 5) Check that P lies inside the convex quad
    #    We do a standard “same-side” test for each edge, using the sign of
    #      N · ( edge × (P−vertex) )
    verts = (v0, v1, v2, v3)
    inside = True
    for i in range(4):
        vi = verts[i]
        vj = verts[(i+1)%4]
        edge = sub(vj, vi)
        vp   = sub(P, vi)
        if dot(N, cross(edge, vp)) < -eps:
            inside = False
            break
    if not inside:
        return -1.0

    # 6) Finally, test if P is inside the circular hole of radius 1 around (0,0,0)
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] < 1.0*1.0 - eps:
        # P is in the hole → the quad is “transparent” there
        return -1.0

    # Otherwise we have a valid hit on the solid part of the quad
    return t


#---------- 6.4 ----------#
import math

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return ( a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] )

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def mul(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def intersect_quad_with_hole(O, D):
    # 1) quad vertices in CCW order as seen from the “front”
    v0 = ( 1.0,  1.0,  1.0)
    v1 = (-1.0,  1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = ( 2.0, -1.0,  2.0)
    verts = [v0, v1, v2, v3]

    # 2) plane normal
    #    use edge v0->v1 and v0->v3
    N = cross(sub(v1, v0), sub(v3, v0))

    # 3) ray-plane intersection
    denom = dot(N, D)
    EPS = 1e-7
    if abs(denom) < EPS:
        return -1.0    # ray is parallel to the plane

    t = dot(N, sub(v0, O)) / denom
    if t <= 0.0:
        return -1.0    # intersection is behind the origin or grazing at O

    # 4) compute intersection point
    P = add(O, mul(D, t))

    # 5) check the circular hole (radius = 1 at origin)
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] < 1.0:
        return -1.0    # passes through the hole

    # 6) check inside convex quad by edge‐side tests
    #    for each edge vi->v(i+1), we want cross(edge, P-vi) to point in same
    #    general direction as the normal N:
    for i in range(4):
        vi = verts[i]
        vj = verts[(i+1) % 4]
        edge = sub(vj, vi)
        vp   = sub(P,   vi)
        if dot(N, cross(edge, vp)) < 0.0:
            return -1.0  # outside this edge → no intersection

    # if we get here, P is on the plane, outside the hole and inside the quad
    return t


#---------- 6.5 ----------#
import math

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return ( a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] )

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def scale(v, s):
    return (v[0]*s, v[1]*s, v[2]*s)

def length(v):
    return math.sqrt(dot(v,v))

def intersect_quad_with_hole(O, D):
    """
    O: ray origin, tuple of 3 floats
    D: ray direction (should be normalized or not – we just solve for t), tuple of 3 floats
    returns: smallest t>0 at which ray hits the quad (excluding the circular hole),
             or -1 if no such intersection.
    """

    # 1) the four corners of the trapezoid, in CCW order as seen from the "front" side
    verts = [
        ( 1.0,  1.0,  1.0),
        (-1.0,  1.0, -1.0),
        (-2.0, -1.0, -2.0),
        ( 2.0, -1.0,  2.0),
    ]

    # 2) build the plane normal
    v0, v1, v2 = verts[0], verts[1], verts[2]
    n = cross(sub(v1, v0), sub(v2, v0))
    n_len = length(n)
    if n_len == 0:
        return -1.0
    # normalize n
    n = scale(n, 1.0/n_len)

    # 3) intersect ray with plane
    denom = dot(D, n)
    # if denom == 0, ray is parallel to plane
    if abs(denom) < 1e-8:
        return -1.0

    t = dot(sub(v0, O), n) / denom
    # we only care about intersections in front of the ray
    if t <= 0:
        return -1.0

    # 4) compute the hit point
    P = add(O, scale(D, t))

    # 5) reject if inside the circular hole of radius 1 at the origin
    if dot(P, P) <= 1.0:
        return -1.0

    # 6) test whether P lies inside the convex quad
    #    using the “edge‐cross” test for a CCW polygon in 3D
    for i in range(len(verts)):
        vi = verts[i]
        vj = verts[(i+1) % len(verts)]
        edge = sub(vj, vi)
        toP  = sub(P, vi)
        # if cross(edge, toP) ⋅ n < 0 then P is to the right of this edge → outside
        if dot(cross(edge, toP), n) < -1e-8:
            return -1.0

    # if we get here, P is on the plane, inside the quad, and outside the hole
    return t


#---------- 6.6 ----------#
def intersect_quad_with_hole(O, D):
    """
    O, D: 3‐tuples or lists giving ray origin and (normalized) direction.
    Returns the smallest positive t such that O + t D hits the trapezoid ABCD
    except for the circular hole of radius 1 about (0,0,0).  Returns −1 if no hit.
    """
    # small epsilon for numerical robustness
    eps = 1e-8

    # helper vector functions
    def sub(u, v):       return (u[0]-v[0],  u[1]-v[1],  u[2]-v[2])
    def dot(u, v):       return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u, v):     return (u[1]*v[2] - u[2]*v[1],
                                 u[2]*v[0] - u[0]*v[2],
                                 u[0]*v[1] - u[1]*v[0])

    # the four corners of the trapezoid, in CCW order
    A = ( 1.0,  1.0,  1.0)
    B = (-1.0,  1.0, -1.0)
    C = (-2.0, -1.0, -2.0)
    Dpt = ( 2.0, -1.0,  2.0)
    verts = [A, B, C, Dpt]

    # compute plane normal n = (B−A) × (D−A)
    e1 = sub(B, A)
    e2 = sub(Dpt, A)
    n  = cross(e1, e2)

    # find ray‐plane intersection t
    denom = dot(D, n)
    if abs(denom) < eps:
        return -1.0        # ray is parallel to the quad's plane

    t = dot(sub(A, O), n) / denom
    if t <= eps:
        return -1.0        # intersection behind or too close to origin

    # compute the 3D intersection point
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])

    # 1) check if P lies inside the convex quad (using edge‐normal test)
    for i in range(4):
        Vi = verts[i]
        Vj = verts[(i+1) % 4]
        edge = sub(Vj, Vi)
        vp   = sub(P, Vi)
        # if cross(edge, vp)·n < 0 then P is outside this edge
        if dot(cross(edge, vp), n) < -eps:
            return -1.0

    # 2) check the circular hole of radius 1 around the plane origin (0,0,0).
    #    Since the plane passes through (0,0,0), the distance in 3D from the origin
    #    is exactly the in‐plane distance.
    dist2 = P[0]*P[0] + P[1]*P[1] + P[2]*P[2]
    if dist2 < (1.0 - eps)**2:
        # strictly inside the hole → no hit
        return -1.0
    # if exactly on the circle (dist2 ≈ 1), we count it as a grazing intersection

    return t


#---------- 6.7 ----------#
import math

def ray_quad_with_hole(O, D):
    """
    O, D: tuples of length 3, ray origin and (normalized) direction
    Returns the nearest positive t such that O + t D hits the trapezoidal quad
    with a circular hole of radius 1 in its plane; or -1 if no hit.
    """
    # Quad vertices in CCW order
    A = ( 1.0,  1.0,  1.0)
    B = (-1.0,  1.0, -1.0)
    C = (-2.0, -1.0, -2.0)
    Dv= ( 2.0, -1.0,  2.0)  # named Dv to avoid shadowing ray‐direction D

    eps = 1e-8

    # vector utilities
    def sub(u,v): return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def dot(u,v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u,v):
        return (u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0])
    def add(u,v): return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
    def scale(u,s): return (u[0]*s, u[1]*s, u[2]*s)
    def length(u): return math.sqrt(dot(u,u))

    # 1) Compute plane normal N = (B-A) × (D-A)
    BA = sub(B, A)
    DA = sub(Dv, A)
    N = cross(BA, DA)
    Nlen = length(N)
    if Nlen < eps:
        return -1.0
    # normalize N
    N = scale(N, 1.0/Nlen)

    # 2) denom = D·N
    denom = dot(D, N)
    if abs(denom) < eps:
        # ray is parallel to plane
        return -1.0

    # 3) t = (A - O)·N / (D·N)
    t = dot(sub(A, O), N) / denom
    if t <= eps:
        return -1.0

    # 4) Compute intersection point P
    P = add(O, scale(D, t))

    # 5) Test if P is inside quad ABCD by splitting into triangles ABC and ACD
    def point_in_triangle(P, V0, V1, V2, N):
        # for a CCW triangle (V0,V1,V2) in plane normal N,
        # P is inside iff all edge‐tests have the same sign:
        # (V1−V0)×(P−V0) · N ≥ 0, etc.
        e0 = dot(cross(sub(V1, V0), sub(P, V0)), N)
        e1 = dot(cross(sub(V2, V1), sub(P, V1)), N)
        e2 = dot(cross(sub(V0, V2), sub(P, V2)), N)
        # allow grazing (≥ 0)
        return e0 >= -eps and e1 >= -eps and e2 >= -eps

    inside_tri1 = point_in_triangle(P, A, B, C, N)
    inside_tri2 = point_in_triangle(P, A, C, Dv, N)
    if not (inside_tri1 or inside_tri2):
        return -1.0

    # 6) Check circular hole of radius 1 centered at origin (0,0,0)
    # We must measure the 3D distance from P to (0,0,0):
    r = length(P)
    # if strictly inside the hole, no hit; if on or outside, it's a hit
    if r < 1.0 - eps:
        return -1.0

    # Otherwise this is the first valid intersection
    return t


# Example usage:
if __name__ == "__main__":
    # Ray from (0,0,-5) toward (0,0,1):
    t = ray_quad_with_hole((0,0,-5), (0,0,1))
    if t > 0:
        P = (0 + t*0, 0 + t*0, -5 + t*1)
        print("Hit at t =", t, "at P =", P)
    else:
        print("No hit")


#---------- 6.8 ----------#
def intersect_quad_with_hole(O, D):
    """
    Compute intersection of ray R(t) = O + t D with the trapezoidal quad
    having corners (1,1,1), (-1,1,-1), (-2,-1,-2), (2,-1,2) in CCW order,
    but with the circular hole of radius 1 at the origin removed.
    Return the smallest t>0 if the ray hits the quad outside the hole,
    or -1 if there is no such intersection.
    """

    # --- helper vector functions ---
    def sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b):
        return (
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        )
    def mul(a, s):
        return (a[0]*s, a[1]*s, a[2]*s)

    # --- define the quad corners in CCW order ---
    P1 = ( 1.0,  1.0,  1.0)
    P2 = (-1.0,  1.0, -1.0)
    P3 = (-2.0, -1.0, -2.0)
    P4 = ( 2.0, -1.0,  2.0)

    # Compute the plane normal N = (P2-P1) x (P4-P1)
    v21 = sub(P2, P1)
    v41 = sub(P4, P1)
    N   = cross(v21, v41)

    # Check if ray is parallel to the plane
    D_dot_N = dot(D, N)
    eps = 1e-9
    if abs(D_dot_N) < eps:
        return -1

    # Compute t of intersection with the infinite plane
    t = dot(sub(P1, O), N) / D_dot_N
    # we only care about intersections in front of the ray
    if t <= 0:
        return -1

    # Compute the intersection point P = O + t*D
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])

    # 1) Check whether P lies inside the (convex) quad.
    #    We do that by checking that for each edge (Ai -> Bi) the
    #    point P is on the same side as the interior, i.e.
    #    (Bi - Ai) x (P - Ai)  has the same direction as N
    corners = [P1, P2, P3, P4]
    for i in range(4):
        A = corners[i]
        B = corners[(i+1) % 4]
        edge       = sub(B, A)
        toPoint    = sub(P, A)
        c = cross(edge, toPoint)
        if dot(c, N) < -eps:
            # If it's strictly on the "outside" side of any edge,
            # then P is outside the quad.
            return -1

    # 2) Check whether P lies inside the circular hole of radius 1
    #    centered at the origin.  If so, the ray just goes through and
    #    does *not* hit the quad.
    r2 = P[0]*P[0] + P[1]*P[1] + P[2]*P[2]
    if r2 < 1.0*1.0 + eps:
        return -1

    # Passed both tests ⇒ this is a valid intersection
    return t


# Example usage:
if __name__ == "__main__":
    # Ray from (0,0,-5) towards (0,0,1):
    O = (0.0, 0.0, -5.0)
    D = (0.0, 0.0,  1.0)
    t = intersect_quad_with_hole(O, D)
    if t > 0:
        P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
        print("Hit at t =", t, "point =", P)
    else:
        print("No intersection")


#---------- 6.9 ----------#
def intersect_quad_with_hole(O, D):
    """
    O, D: 3‐tuples or lists of floats (ray origin and direction)
    returns: the t > 0 of the first intersection with the quad outside
             the unit‐radius hole at the origin, or -1 if none.
    """
    eps = 1e-8

    # 1) plane of the quad is x - z = 0 (normal = (1,0,-1)), passing through A=(1,1,1)
    n = (1.0, 0.0, -1.0)
    A = (1.0, 1.0, 1.0)

    # compute denom = n·D
    denom = n[0]*D[0] + n[1]*D[1] + n[2]*D[2]
    if abs(denom) < eps:
        # ray is parallel to the plane
        return -1.0

    # compute t = n·(A - O) / (n·D)
    num = n[0]*(A[0] - O[0]) + n[1]*(A[1] - O[1]) + n[2]*(A[2] - O[2])
    t = num / denom
    if t <= eps:
        # intersection is behind the origin or too close to count
        return -1.0

    # 2) compute the intersection point P
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 3) test if P lies inside the convex quad ABCD
    verts = [
        (1.0,  1.0,  1.0),   # A
        (-1.0, 1.0, -1.0),   # B
        (-2.0,-1.0, -2.0),   # C
        (2.0, -1.0,  2.0)    # D
    ]
    # for each edge vi->vnext, require (edge × (P - vi))·n >= 0
    for i in range(4):
        v0 = verts[i]
        v1 = verts[(i+1) % 4]
        # edge vector e = v1 - v0
        e = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
        # w = P - v0
        w = (P[0] - v0[0], P[1] - v0[1], P[2] - v0[2])
        # cross e × w
        cx = e[1]*w[2] - e[2]*w[1]
        cy = e[2]*w[0] - e[0]*w[2]
        cz = e[0]*w[1] - e[1]*w[0]
        # if the point is ever to the “right” of an edge (w.r.t. n),
        # it is outside the quad
        if cx*n[0] + cy*n[1] + cz*n[2] < -eps:
            return -1.0

    # 4) test if P falls inside the circular hole of radius 1 at the origin
    #    (distance in the plane == full 3D distance since P is on the plane x=z)
    if (P[0]*P[0] + P[1]*P[1] + P[2]*P[2]) < 1.0:
        # it passed through the hole, so no intersection with the quad
        return -1.0

    # otherwise this is the first hit of the quad
    return t


#---------- 7.0 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    O: ray origin, tuple of 3 floats (ox, oy, oz)
    D: ray direction, tuple of 3 floats (dx, dy, dz)
    returns: the smallest t > 0 such that O + t*D hits the plane z=y
             but not inside the triangle A=(1,1,1), B=(-1,0,0), C=(0,-1,-1).
             If no such intersection exists, returns -1.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    # Plane is z = y  <=>  (0, -1, 1) · (x,y,z) = 0
    # so denom = D·n = dz - dy
    denom = dz - dy
    if abs(denom) < eps:
        # Ray is (nearly) parallel to the plane
        return -1.0

    # Solve (O + t D)·n = 0  =>  t = (O·n) / (–D·n), but since n·(x,y,z)= -y+z:
    #    O·n = -oy + oz
    t = (oz - oy) / denom

    if t <= 0:
        # Intersection is behind the ray origin
        return -1.0

    # Compute intersection point P
    px = ox + t*dx
    py = oy + t*dy
    pz = oz + t*dz

    # Now do a point‐in‐triangle test in 3D for triangle A,B,C
    A = (1.0,  1.0,  1.0)
    B = (-1.0, 0.0,  0.0)
    C = (0.0, -1.0, -1.0)

    # Build vectors
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])

    v0 = sub(B, A)
    v1 = sub(C, A)
    v2 = (px - A[0], py - A[1], pz - A[2])

    # Compute dot products
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot11 = dot(v1, v1)
    dot20 = dot(v2, v0)
    dot21 = dot(v2, v1)

    # Compute barycentric coordinates (u,v)
    denom2 = dot00*dot11 - dot01*dot01
    if abs(denom2) < eps:
        # Degenerate triangle?  Should not happen here
        return -1.0

    u = ( dot11*dot20 - dot01*dot21 ) / denom2
    v = ( dot00*dot21 - dot01*dot20 ) / denom2

    # If inside triangle, then the ray "passes through" the hole => no intersection
    if u >= 0 and v >= 0 and (u + v) <= 1:
        return -1.0

    # Otherwise the plane outside the hole blocks the ray here
    return t


#---------- 7.1 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    Intersect ray R(t) = O + t*D with the plane y = z, which has a triangular hole
    with vertices (1,1,1), (-1,0,0), (0,-1,-1) (in CCW order).  Return the smallest
    positive t if the ray hits the plane outside the hole; otherwise return -1.
    
    O, D are 3‐tuples or lists of floats (Ox,Oy,Oz), (Dx,Dy,Dz).
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) Find intersection t with plane y = z:
    denom = Dy - Dz
    if abs(denom) < eps:
        # Ray is parallel to the plane (or lies in it) → no valid single intersection
        return -1.0

    t = (Oz - Oy) / denom
    # 2) We only want intersections in front of the ray origin:
    if t <= 0.0:
        return -1.0

    # 3) Compute the intersection point P:
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz
    P = (Px, Py, Pz)

    # 4) Test whether P lies *strictly* inside the triangle
    #    Triangle vertices (in CCW order):
    v0 = (1.0,  1.0,  1.0)
    v1 = (-1.0, 0.0,  0.0)
    v2 = (0.0, -1.0, -1.0)

    # set up for barycentric‐coordinate test
    def sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    u = sub(v1, v0)
    v = sub(v2, v0)
    w = sub(P,  v0)

    uu = dot(u,u)
    uv = dot(u,v)
    vv = dot(v,v)
    wu = dot(w,u)
    wv = dot(w,v)

    # denominator for barycentric coords
    Dden = uu*vv - uv*uv
    if abs(Dden) < eps:
        # Degenerate triangle (should not happen here)
        return t

    # barycentric coordinates
    s = ( vv*wu - uv*wv ) / Dden
    r = ( uu*wv - uv*wu ) / Dden

    # If 0 < s < 1, 0 < r < 1, and s+r < 1, P is strictly inside the open triangle.
    # In that case the ray “passes through the hole,” so we return -1.
    if s > eps and r > eps and (s + r) < 1.0 - eps:
        return -1.0

    # Otherwise P is either outside the triangle or on its boundary → a valid plane hit
    return t


#---------- 7.2 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    O: origin of ray, tuple of 3 floats
    D: direction of ray, tuple of 3 floats
    returns: t >= 0 at intersection, or -1 if no intersection
    """

    # Triangle (hole) vertices in the plane y=z
    A = ( 1.0,  1.0,  1.0)
    B = (-1.0,  0.0,  0.0)
    C = ( 0.0, -1.0, -1.0)

    # 1) Solve for intersection with plane y - z = 0:
    #    O_y + t D_y = O_z + t D_z
    # => t * (D_y - D_z) = -(O_y - O_z)
    denom = D[1] - D[2]
    if abs(denom) < eps:
        # Ray is (nearly) parallel to the plane
        return -1.0

    t = -(O[1] - O[2]) / denom
    if t < 0:
        # Intersection is behind the ray origin
        return -1.0

    # Compute the 3D intersection point
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 2) Test whether P lies inside the triangle ABC (the hole).
    #    We'll use the "same‐side" / barycentric‐sign test in 3D.
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u, v):
        return (
            u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0]
        )

    AB = sub(B, A)
    BC = sub(C, B)
    CA = sub(A, C)

    PA = sub(P, A)
    PB = sub(P, B)
    PC = sub(P, C)

    # triangle normal
    n = cross(AB, sub(C, A))

    # compute the three cross‐products and project onto the normal
    c1 = dot(n, cross(AB, PA))
    c2 = dot(n, cross(BC, PB))
    c3 = dot(n, cross(CA, PC))

    # if P is strictly inside or on edge of ABC, all three c_i
    # have the same sign (or zero).  In that case it's in the hole,
    # so the ray passes through and we report NO intersection -> -1.
    inside_or_on = (
        (c1 >= -eps and c2 >= -eps and c3 >= -eps) or
        (c1 <=  eps and c2 <=  eps and c3 <=  eps)
    )
    if inside_or_on:
        return -1.0

    # otherwise P is on the plane but outside the triangular hole
    return t


#---------- 7.3 ----------#
def intersect_plane_with_tri_hole(O, D):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: the smallest positive t such that O + t*D hits the plane y=z
             but not inside the triangle hole; or -1 if no such intersection.
    """
    # 1) Plane is y = z  <=>  (0,1,-1) · (X - any_plane_point) = 0
    #    we can pick the plane point A = (1,1,1) (one corner of the hole)
    #    normal n = (0,1,-1)
    n = (0.0, 1.0, -1.0)

    # denominator = D·n
    denom = D[0]*n[0] + D[1]*n[1] + D[2]*n[2]
    if abs(denom) < 1e-8:
        # ray is parallel to the plane
        return -1.0

    # numerator = (A - O)·n
    A = (1.0, 1.0, 1.0)
    AO = (A[0] - O[0], A[1] - O[1], A[2] - O[2])
    numer = AO[0]*n[0] + AO[1]*n[1] + AO[2]*n[2]

    t = numer / denom
    if t < 0:
        # intersection is “behind” the ray origin
        return -1.0

    # compute intersection point P
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # 2) Test whether P lies inside the triangle A,B,C
    B = (-1.0, 0.0, 0.0)
    C = ( 0.0,-1.0,-1.0)

    if point_in_triangle(P, A, B, C):
        # it hit the hole → no hit
        return -1.0
    else:
        # it hit the plane outside the hole
        return t


def point_in_triangle(P, A, B, C):
    """
    Barycentric‐coordinate test for P inside triangle ABC.
    returns True if P is inside or on the edge of the triangle.
    """
    # vectors
    v0 = (C[0]-A[0], C[1]-A[1], C[2]-A[2])
    v1 = (B[0]-A[0], B[1]-A[1], B[2]-A[2])
    v2 = (P[0]-A[0], P[1]-A[1], P[2]-A[2])

    # dot products
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    # barycentric coordinates
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) < 1e-8:
        return False

    invDen = 1.0 / denom
    u = ( dot11*dot02 - dot01*dot12) * invDen
    v = (-dot01*dot02 + dot00*dot12) * invDen

    # check if point is in triangle
    return (u >= 0) and (v >= 0) and (u + v <= 1)


# --- Example usage ---
O = (0.0, 2.0, 0.0)
D = (0.0, -1.0, 0.0)
t = intersect_plane_with_tri_hole(O, D)
print("t =", t)


#---------- 7.4 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    O, D: 3‐tuples (Ox,Oy,Oz), (Dx,Dy,Dz)
    returns: t >= 0 if the ray R(t)=O+tD hits the plane outside the hole,
             or -1 if it misses or passes through the hole.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with plane y = z
    denom = Dy - Dz
    if abs(denom) < eps:
        # Ray is (nearly) parallel to the plane
        return -1.0

    t = (Oz - Oy) / denom
    if t < 0:
        # Intersection is behind the origin
        return -1.0

    # 2) compute the intersection point P
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz

    # 3) test against triangle A,B,C
    A = (1.0,  1.0,  1.0)
    B = (-1.0, 0.0,  0.0)
    C = (0.0, -1.0, -1.0)

    # build the two edge vectors of the triangle in 3D
    v0 = (C[0]-A[0], C[1]-A[1], C[2]-A[2])
    v1 = (B[0]-A[0], B[1]-A[1], B[2]-A[2])
    v2 = (Px - A[0], Py - A[1], Pz - A[2])

    # dot‐products
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    # barycentric coords
    invDenom = 1.0 / (dot00*dot11 - dot01*dot01)
    u = (dot11*dot02 - dot01*dot12) * invDenom
    v = (dot00*dot12 - dot01*dot02) * invDenom

    # strictly inside the triangle means: u>0, v>0, u+v<1
    if (u > eps) and (v > eps) and (u + v < 1.0 - eps):
        # it passed through the hole
        return -1.0

    # otherwise it hit the opaque part of the plane
    return t


#---------- 7.5 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    O, D: 3‐tuples or lists of floats
         O = ray origin, D = ray direction
    returns: the smallest t >= 0 so that O + t*D
             hits the plane y−z=0 outside the triangle hole,
             or -1 if no such intersection exists.
    """
    # Plane normal n for y − z = 0 is (0,1,−1)
    n = (0.0, 1.0, -1.0)

    # Dot(D,n)
    denom = D[0]*n[0] + D[1]*n[1] + D[2]*n[2]
    # If denom ≈ 0 the ray is parallel to the plane
    if abs(denom) < eps:
        return -1.0

    # Solve O·n + t (D·n) = 0  ⇒  t = −(O·n)/(D·n)
    Ond = O[0]*n[0] + O[1]*n[1] + O[2]*n[2]
    t = -Ond / denom

    # we only care about intersections in front of the origin
    if t < 0:
        return -1.0

    # Compute the hit point P = O + t D
    P = (O[0] + t*D[0],
         O[1] + t*D[1],
         O[2] + t*D[2])

    # Now test whether P lies inside the triangular hole
    v0 = ( 1.0,  1.0,  1.0)
    v1 = (-1.0,  0.0,  0.0)
    v2 = ( 0.0, -1.0, -1.0)

    # helper: vector subtraction
    def sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

    # helper: dot product
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    # edge vectors
    v0v1 = sub(v1, v0)
    v0v2 = sub(v2, v0)
    v0p  = sub(P,  v0)

    # Compute barycentric coordinates (s, t2) of P in triangle v0,v1,v2
    d00 = dot(v0v1, v0v1)
    d01 = dot(v0v1, v0v2)
    d11 = dot(v0v2, v0v2)
    d20 = dot(v0p,  v0v1)
    d21 = dot(v0p,  v0v2)

    denom_tri = d00*d11 - d01*d01
    # if denom_tri is zero the triangle is degenerate, but it's not in our case
    s = ( d11*d20 - d01*d21 ) / denom_tri
    t2 = ( d00*d21 - d01*d20 ) / denom_tri

    # P is inside the triangle iff s >= 0, t2 >= 0 and s + t2 <= 1
    inside_hole = (s >= -eps) and (t2 >= -eps) and (s + t2 <= 1 + eps)

    # If the ray hit falls _inside_ the hole, we ignore it → return -1
    if inside_hole:
        return -1.0

    # Otherwise this is our first intersection with the plane exterior
    return t


#---------- 7.6 ----------#
def intersect_plane_with_tri_hole(O, D, eps=1e-9):
    """
    O, D: 3‐tuples (Ox,Oy,Oz), (Dx,Dy,Dz)
    returns:
      t >= 0   if the ray hits the plane OUTSIDE (or ON THE EDGE) of the hole,
      -1      if it misses (parallel, behind origin, or goes through the hole interior).
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with plane y - z = 0
    denom = Dy - Dz
    if abs(denom) < eps:
        return -1.0    # parallel → no hit

    t = (Oz - Oy) / denom
    if t < 0:
        return -1.0    # behind the ray origin → no hit

    # 2) compute intersection point P
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz
    P = (Px, Py, Pz)

    # 3) triangle hole vertices
    A = ( 1.0,  1.0,  1.0)
    B = (-1.0,  0.0,  0.0)
    C = ( 0.0, -1.0, -1.0)

    # plane normal
    n = (0.0, 1.0, -1.0)

    # helper functions
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u, v):
        return (u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0])

    # edges and point‐vectors
    AB = sub(B, A);   AP = sub(P, A)
    BC = sub(C, B);   BP = sub(P, B)
    CA = sub(A, C);   CP = sub(P, C)

    # cross each edge with the corresponding P‐vector
    c0 = cross(AB, AP)
    c1 = cross(BC, BP)
    c2 = cross(CA, CP)

    # dot with normal
    s0 = dot(c0, n)
    s1 = dot(c1, n)
    s2 = dot(c2, n)

    # if all three are strictly > 0, we're in the *interior* of the hole → MISS
    if s0 > eps and s1 > eps and s2 > eps:
        return -1.0

    # otherwise we hit the plane (outside or on edge of hole)
    return t


#---------- 7.7 ----------#
def intersect_ray_with_plane_hole(O, D, eps=1e-8):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: the smallest t >= 0 such that O + t D lies on the plane y = z
             but not inside the triangle hole with vertices
             A=(1,1,1), B=(-1,0,0), C=(0,-1,-1).
             If there is no such intersection, returns -1.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with plane y = z
    denom = Dy - Dz
    if abs(denom) < eps:
        # ray is parallel to plane (or lies in it); no valid intersection
        return -1.0

    # solve Oy + t Dy == Oz + t Dz  =>  t = (Oz - Oy) / (Dy - Dz)
    t = (Oz - Oy) / denom

    if t < 0:
        # intersection is "behind" the ray origin
        return -1.0

    # 2) compute the intersection point P
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    # (we know Pz == Py because it's on the plane y = z)

    # 3) test whether P lies inside the triangle A',B',C' in the (x,y)-plane
    #    where A'=(1,1), B'=(-1,0), C'=(0,-1).  If it is inside, that's
    #    the hole, so we return -1.  Otherwise we hit the opaque part of the plane.

    # 2D cross‐product sign function
    def sign(p1, p2, p3):
        # computes (p1.x - p3.x)*(p2.y - p3.y) - (p2.x - p3.x)*(p1.y - p3.y)
        return (p1[0] - p3[0])*(p2[1] - p3[1]) - (p2[0] - p3[0])*(p1[1] - p3[1])

    P2 = (Px, Py)
    A2 = (1.0,  1.0)
    B2 = (-1.0, 0.0)
    C2 = (0.0, -1.0)

    b1 = sign(P2, A2, B2) < -eps
    b2 = sign(P2, B2, C2) < -eps
    b3 = sign(P2, C2, A2) < -eps

    inside_triangle = (b1 == b2 == b3)

    if inside_triangle:
        # P is in the hole ⇒ no intersection
        return -1.0

    # otherwise, this is the correct hit
    return t


#---------- 7.8 ----------#
def intersect_plane_with_triangular_hole(O, D, eps=1e-8):
    """
    O, D : tuples of 3 floats (ox, oy, oz), (dx, dy, dz)
    Returns t >= 0 if the ray R(t)=O + t*D hits the plane y=z *outside* the
    triangular hole, or -1 if there is no such intersection (ray misses or
    passes through the hole).
    """

    # 1) intersect with the plane y = z
    oy, oz = O[1], O[2]
    dy, dz = D[1], D[2]
    denom = dy - dz
    if abs(denom) < eps:
        # Ray is parallel to the plane y=z (or almost so)
        return -1.0

    t = (oz - oy) / denom
    if t < 0.0:
        # Intersection is "behind" the ray origin
        return -1.0

    # Compute the 3D point of intersection
    Px = O[0] + t * D[0]
    Py = O[1] + t * D[1]
    Pz = O[2] + t * D[2]
    P = (Px, Py, Pz)

    # 2) set up the triangle in the plane y=z, given CCW in 3D:
    A = ( 1.0,  1.0,  1.0)
    B = (-1.0,  0.0,  0.0)
    C = ( 0.0, -1.0, -1.0)

    # helper vector ops
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u, v):
        return (u[1]*v[2] - u[2]*v[1],
                u[2]*v[0] - u[0]*v[2],
                u[0]*v[1] - u[1]*v[0])

    # Precompute the triangle's face normal (orientation)
    AB = sub(B, A)
    AC = sub(C, A)
    N  = cross(AB, AC)

    # 3) edge‐tests: for each directed edge, see if P lies to the "inside" side.
    #    We will *only* declare "no intersection" (hole) if P is strictly 
    #    inside *all* three edges.  If P is outside or exactly on any edge 
    #    boundary, we count that as a hit on the plane!
    def edge_test(V0, V1):
        E = sub(V1, V0)
        VP = sub(P,   V0)
        return dot(cross(E, VP), N)

    # compute the three edge‐dot‐products
    dab = edge_test(A, B)
    dbc = edge_test(B, C)
    dca = edge_test(C, A)

    # if all three are > 0, P is strictly inside the triangle => hole => no hit
    if dab > eps and dbc > eps and dca > eps:
        return -1.0

    # otherwise P is outside the hole or on its boundary => we hit the plane
    return t


#---------- 7.9 ----------#
def intersect_plane_with_triangular_hole(O, D):
    """
    O: tuple of 3 floats = ray origin (Ox, Oy, Oz)
    D: tuple of 3 floats = ray direction (Dx, Dy, Dz)
    returns: the smallest t>0 such that O + t D hits the plane y=z 
             but NOT inside the triangle with corners
             v0=(1,1,1), v1=(-1,0,0), v2=(0,-1,-1).
             returns -1 if no valid intersection.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Plane is y − z = 0  =>  n = (0,1,−1);  n·X = 0
    denom = (Dy - Dz)
    if abs(denom) < 1e-8:
        # ray is (nearly) parallel to plane
        return -1.0

    # solve O + t D lies in y=z  =>  t = (Oz - Oy)/(Dy - Dz)
    t = (Oz - Oy) / denom
    if t <= 0:
        # intersection is behind the ray origin
        return -1.0

    # compute intersection point P
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz

    # now test whether P lies inside the triangle hole v0,v1,v2
    v0 = (1.0, 1.0, 1.0)
    v1 = (-1.0, 0.0, 0.0)
    v2 = (0.0, -1.0, -1.0)

    # helper: vector subtraction
    def sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

    # helper: dot product
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    # we will do a simple barycentric‐coordinate test in the plane
    # using v0 as the origin of the triangle coordinate system:
    v0v1 = sub(v1, v0)
    v0v2 = sub(v2, v0)
    v0p  = sub((Px, Py, Pz), v0)

    dot00 = dot(v0v2, v0v2)
    dot01 = dot(v0v2, v0v1)
    dot11 = dot(v0v1, v0v1)
    dot02 = dot(v0v2, v0p)
    dot12 = dot(v0v1, v0p)

    # compute barycentric coordinates (u, v)
    denom_uv = dot00*dot11 - dot01*dot01
    if abs(denom_uv) < 1e-8:
        # degenerate triangle (should not happen in our problem)
        return t

    inv_denom = 1.0 / denom_uv
    u = ( dot11 * dot02 - dot01 * dot12) * inv_denom
    v = ( dot00 * dot12 - dot01 * dot02) * inv_denom

    # inside‐triangle test: u >= 0, v >= 0, u+v <= 1
    if u >= 0.0 and v >= 0.0 and (u + v) <= 1.0:
        # P lies in the triangular hole → no intersection
        return -1.0

    # otherwise P is on the plane but outside the hole → valid hit
    return t


#---------- 8.0 ----------#
import math

def ray_ellipsoid_intersection(O, D):
    """
    O: origin of the ray, tuple of 3 floats (Ox, Oy, Oz)
    D: direction of the ray, tuple of 3 floats (Dx, Dy, Dz)
    
    The ellipsoid is centered at C = (1,1,1) with semi-axes a=2, b=3, c=4:
       (x-1)^2/2^2 + (y-1)^2/3^2 + (z-1)^2/4^2 = 1.
    
    Returns the smallest t >= 0 such that O + t*D hits the ellipsoid,
    or -1 if there is no intersection.
    """

    # Ellipsoid center
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    # Semi-axes
    a2, b2, c2 = 2.0*2.0, 3.0*3.0, 4.0*4.0

    # Compute O - C
    ox, oy, oz = O[0] - Cx, O[1] - Cy, O[2] - Cz
    dx, dy, dz = D

    # Quadratic coefficients A t^2 + B t + C = 0
    A = (dx*dx)/a2 + (dy*dy)/b2 + (dz*dz)/c2
    B = 2.0 * (ox*dx/a2 + oy*dy/b2 + oz*dz/c2)
    C = (ox*ox)/a2 + (oy*oy)/b2 + (oz*oz)/c2 - 1.0

    # Discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0

    sqrt_disc = math.sqrt(disc)
    t0 = (-B - sqrt_disc) / (2.0 * A)
    t1 = (-B + sqrt_disc) / (2.0 * A)

    # We want the smallest non-negative root
    t_min = min(t0, t1)
    t_max = max(t0, t1)

    if t_min >= 0.0:
        return t_min
    elif t_max >= 0.0:
        return t_max
    else:
        return -1.0


#---------- 8.1 ----------#
import math

def ray_ellipsoid_intersect(O, D):
    """
    O: tuple of 3 floats, ray origin (Ox,Oy,Oz)
    D: tuple of 3 floats, ray direction  (Dx,Dy,Dz)
    Returns: float t >= 0 for the closest intersection, or -1 if none.
    """
    # Ellipsoid centre and semi‐axes
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    a, b, c    = 2.0, 3.0, 4.0

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Compute the coefficients A, B, C of the quadratic At^2 + Bt + C = 0
    # where we have ((Ox + t Dx − Cx)^2 / a^2) + ... = 1
    # Expand to standard form in t.
    A = (Dx*Dx) / (a*a) + (Dy*Dy) / (b*b) + (Dz*Dz) / (c*c)
    B = 2.0 * ( (Ox - Cx)*Dx / (a*a)
              + (Oy - Cy)*Dy / (b*b)
              + (Oz - Cz)*Dz / (c*c) )
    C = ((Ox - Cx)**2) / (a*a) \
      + ((Oy - Cy)**2) / (b*b) \
      + ((Oz - Cz)**2) / (c*c) \
      - 1.0

    # Compute discriminant
    disc = B*B - 4.0*A*C

    if disc < 0.0:
        # No real roots: no intersection
        return -1.0

    sqrt_disc = math.sqrt(disc)

    # Two possible solutions
    t0 = (-B - sqrt_disc) / (2.0*A)
    t1 = (-B + sqrt_disc) / (2.0*A)

    # We only want intersections "in front" of the ray origin, i.e. t >= 0
    ts = []
    if t0 >= 0.0:
        ts.append(t0)
    if t1 >= 0.0:
        ts.append(t1)

    if not ts:
        # Both intersections are behind the ray
        return -1.0

    # Return the nearest one
    return min(ts)

# Example usage:
if __name__ == "__main__":
    O = (0.0, 0.0, 0.0)
    D = (1.0, 1.0, 1.0)
    t = ray_ellipsoid_intersect(O, D)
    print("Intersection t =", t)


#---------- 8.2 ----------#
import math

def intersect_ray_ellipsoid(O, D):
    """
    O : tuple of 3 floats, ray origin (Ox, Oy, Oz)
    D : tuple of 3 floats, ray direction (Dx, Dy, Dz)
    Returns the smallest t >= 0 such that O + t*D lies on the ellipsoid
    centered at (1,1,1) with semi‐axes (2,3,4), or −1 if no intersection.
    """
    # Ellipsoid center and squared semi‐axes
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 2.0*2.0, 3.0*3.0, 4.0*4.0

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Move to ellipsoid's local coordinates: P(t) = O-C + t*D
    x0 = Ox - Cx
    y0 = Oy - Cy
    z0 = Oz - Cz

    # Build quadratic coefficients A t^2 + B t + C = 0
    A = (Dx*Dx)/a2 + (Dy*Dy)/b2 + (Dz*Dz)/c2
    B = 2.0*((x0*Dx)/a2 + (y0*Dy)/b2 + (z0*Dz)/c2)
    C = (x0*x0)/a2 + (y0*y0)/b2 + (z0*z0)/c2 - 1.0

    # If A is zero (degenerate), no real quadratic → no intersection
    if abs(A) < 1e-12:
        return -1.0

    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0   # no real roots → miss

    sqrt_disc = math.sqrt(disc)
    t1 = (-B - sqrt_disc) / (2.0*A)
    t2 = (-B + sqrt_disc) / (2.0*A)

    # We want the smallest non‐negative t
    t_min = float('inf')
    if t1 >= 0.0:
        t_min = min(t_min, t1)
    if t2 >= 0.0:
        t_min = min(t_min, t2)

    if t_min == float('inf'):
        return -1.0   # both intersections are "behind" the ray origin
    return t_min


#---------- 8.3 ----------#
import math

def intersect_ray_ellipsoid(O, D):
    """
    O: tuple of 3 floats, ray origin (Ox,Oy,Oz)
    D: tuple of 3 floats, ray direction (Dx,Dy,Dz)
    returns: float t ≥ 0 of first intersection, or -1 if none
    """
    # Ellipsoid center and squared semi‐axes
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 2.0**2, 3.0**2, 4.0**2

    # Translate ray origin into ellipsoid's local coords
    ox = O[0] - cx
    oy = O[1] - cy
    oz = O[2] - cz
    dx, dy, dz = D

    # Quadratic coefficients A t^2 + B t + C = 0
    A = (dx*dx)/a2 + (dy*dy)/b2 + (dz*dz)/c2
    B = 2.0 * (ox*dx/a2 + oy*dy/b2 + oz*dz/c2)
    C = (ox*ox)/a2 + (oy*oy)/b2 + (oz*oz)/c2 - 1.0

    # Solve discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0    # no real roots → no intersection

    # one or two real roots
    sqrt_disc = math.sqrt(disc)
    t1 = (-B - sqrt_disc) / (2.0 * A)
    t2 = (-B + sqrt_disc) / (2.0 * A)

    # we want the smallest non‐negative t
    tmin = float('inf')
    if t1 >= 0.0:
        tmin = t1
    if t2 >= 0.0 and t2 < tmin:
        tmin = t2

    return tmin if tmin != float('inf') else -1.0


#---------- 8.4 ----------#
import math

def intersect_ellipsoid(O, D):
    """
    O: ray origin, tuple of 3 floats (Ox, Oy, Oz)
    D: ray direction, tuple of 3 floats (Dx, Dy, Dz)
    returns: the smallest positive t so that O + t*D hits the ellipsoid
             centered at (1,1,1) with semi‐axes 2,3,4.  Returns -1 if none.
    """
    # ellipsoid center
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    # semi‐axes
    ax2, ay2, az2 = 2.0*2.0, 3.0*3.0, 4.0*4.0

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # shift ray origin into ellipsoid‐space
    x0 = Ox - Cx
    y0 = Oy - Cy
    z0 = Oz - Cz

    # build quadratic coefficients A t^2 + B t + C = 0
    A = (Dx*Dx)/ax2 + (Dy*Dy)/ay2 + (Dz*Dz)/az2
    B = 2.0 * ( x0*Dx/ax2 + y0*Dy/ay2 + z0*Dz/az2 )
    C = (x0*x0)/ax2 + (y0*y0)/ay2 + (z0*z0)/az2 - 1.0

    # if A is zero (degenerate), no intersection
    if abs(A) < 1e-12:
        return -1.0

    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0

    sqrt_disc = math.sqrt(disc)
    t1 = (-B - sqrt_disc) / (2.0*A)
    t2 = (-B + sqrt_disc) / (2.0*A)

    # we want the smallest positive t
    ts = []
    if t1 >= 0.0:
        ts.append(t1)
    if t2 >= 0.0:
        ts.append(t2)
    if not ts:
        return -1.0
    return min(ts)

# Example usage:
O = (0.0, 1.0, 1.0)
D = (1.0, 0.0, 0.0)
print(intersect_ellipsoid(O, D))   # should hit the ellipsoid at t≈?


#---------- 8.5 ----------#
You can get the intersection by plugging the parametric ray  
  R(t) = O + t D  
into the ellipsoid equation  
  (x–1)²/2² + (y–1)²/3² + (z–1)²/4² = 1  
and solving the resulting quadratic in t.  

Step by step:

 1. Shift the ray origin so that the ellipsoid is centered at the origin:  
    xc = O.x – 1  
    yc = O.y – 1  
    zc = O.z – 1  

 2. Write the ellipsoid in the form  
       x²/2² + y²/3² + z²/4² = 1.  

 3. Substitute x = xc + t·dx, y = yc + t·dy, z = zc + t·dz into that equation.  
    You get a quadratic  
       A t² + B t + C = 0  
    with  
       A = (dx*dx)/4  + (dy*dy)/9  + (dz*dz)/16  
       B = 2*(xc*dx/4 + yc*dy/9 + zc*dz/16)  
       C = (xc*xc)/4  + (yc*yc)/9  + (zc*zc)/16  – 1  

 4. Compute the discriminant  
       disc = B*B – 4*A*C.  
    If disc < 0, there is no intersection → return –1.  

 5. Otherwise compute the two roots  
       t0 = (–B – √disc)/(2*A)  
       t1 = (–B + √disc)/(2*A)  
    with t0 ≤ t1.  

 6. Since the ray starts outside and we only want intersections “in front” of the ray‐origin, pick the smallest positive root:  
       if t0 > 0: return t0  
       elif t1 > 0: return t1  
       else: return –1  

Here is a self‐contained Python function:

    from math import sqrt

    def intersect_ray_ellipsoid(O, D):
        # ellipsoid center and semi‐axes
        cx, cy, cz = 1.0, 1.0, 1.0
        a, b, c = 2.0, 3.0, 4.0

        # shift ray origin to ellipsoid‐centered coordinates
        xc = O[0] - cx
        yc = O[1] - cy
        zc = O[2] - cz
        dx, dy, dz = D

        # coefficients of A t^2 + B t + C = 0
        A = (dx*dx)/(a*a) + (dy*dy)/(b*b) + (dz*dz)/(c*c)
        B = 2.0*(xc*dx/(a*a) + yc*dy/(b*b) + zc*dz/(c*c))
        C = (xc*xc)/(a*a) + (yc*yc)/(b*b) + (zc*zc)/(c*c) - 1.0

        # solve quadratic
        disc = B*B - 4.0*A*C
        if disc < 0.0:
            return -1.0

        sqrt_disc = sqrt(disc)
        # two possible intersections
        t0 = (-B - sqrt_disc) / (2.0*A)
        t1 = (-B + sqrt_disc) / (2.0*A)

        # return the smallest positive t, or -1 if none
        if t0 > 0.0:
            return t0
        if t1 > 0.0:
            return t1
        return -1.0

Explanation:

- We form a quadratic by inserting the ray equation into the ellipsoid’s implicit formula.
- The coefficients A, B, C come straight out of expanding  
     (xc + t dx)²/a² + (yc + t dy)²/b² + (zc + t dz)²/c² = 1  
- A negative discriminant means no real roots (no intersection).  
- Otherwise the two roots t0, t1 are the “enter” and “exit” distances along the ray.  
- Since the ray starts outside, you pick the smaller positive root (if any).  

That single function will give you the distance to the nearest intersection or –1 if there is none.

#---------- 8.6 ----------#
import math

def intersect_ray_ellipsoid(O, D):
    """
    O: origin of ray, tuple of 3 floats (Ox, Oy, Oz)
    D: direction of ray, tuple of 3 floats (Dx, Dy, Dz)
    returns: t >= 0 to the first intersection, or -1.0 if none
    """
    # ellipsoid parameters
    cx, cy, cz = 1.0, 1.0, 1.0   # center
    a, b, c = 2.0, 3.0, 4.0      # semi-axes lengths

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # build quadratic coefficients A t^2 + B t + C = 0
    A = (Dx*Dx)/(a*a) + (Dy*Dy)/(b*b) + (Dz*Dz)/(c*c)
    B = 2.0 * ( Dx*(Ox-cx)/(a*a)
               + Dy*(Oy-cy)/(b*b)
               + Dz*(Oz-cz)/(c*c) )
    C = ((Ox-cx)**2)/(a*a) + ((Oy-cy)**2)/(b*b) + ((Oz-cz)**2)/(c*c) - 1.0

    # compute discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0     # no real intersections

    sqrt_disc = math.sqrt(disc)

    # two solutions
    t1 = (-B - sqrt_disc) / (2.0*A)
    t2 = (-B + sqrt_disc) / (2.0*A)

    # we only want non-negative t (ray starts outside the ellipsoid)
    ts = []
    if t1 >= 0.0:
        ts.append(t1)
    if t2 >= 0.0:
        ts.append(t2)

    if not ts:
        return -1.0

    # return the nearest intersection
    return min(ts)


# Example usage:
O = (0.0, 0.0, 0.0)
D = (1.0, 0.5, 0.7)
t = intersect_ray_ellipsoid(O, D)
print("first intersection t =", t)


#---------- 8.7 ----------#
import math

def intersect_ellipsoid(O, D):
    """
    Compute the intersection between a ray R(t) = O + t D and the
    axis-aligned ellipsoid centered at (1,1,1) with semi-axes a=2, b=3, c=4.
    O and D are 3‐tuples (ox, oy, oz), (dx, dy, dz).
    Returns the smallest t >= 0 if an intersection exists, otherwise -1.
    """
    # Unpack
    ox, oy, oz = O
    dx, dy, dz = D

    # Ellipsoid center and semi-axes
    cx, cy, cz = 1.0, 1.0, 1.0
    a, b, c = 2.0, 3.0, 4.0

    # Shift ray origin to ellipsoid coordinate system
    x0 = ox - cx
    y0 = oy - cy
    z0 = oz - cz

    # Precompute reciprocals of squared axes
    ia2 = 1.0 / (a * a)
    ib2 = 1.0 / (b * b)
    ic2 = 1.0 / (c * c)

    # Quadratic coefficients A t^2 + B t + C = 0
    A = (dx*dx)*ia2 + (dy*dy)*ib2 + (dz*dz)*ic2
    B = 2.0*((x0*dx)*ia2 + (y0*dy)*ib2 + (z0*dz)*ic2)
    C = (x0*x0)*ia2 + (y0*y0)*ib2 + (z0*z0)*ic2 - 1.0

    # If A is zero, the ray direction is degenerate in this metric:
    EPS = 1e-9
    if abs(A) < EPS:
        # Linear: B t + C = 0  => t = -C / B
        if abs(B) < EPS:
            return -1.0  # No solution or infinitely many if C==0, but ray start is outside so no hit
        t = -C / B
        return t if t >= 0.0 else -1.0

    # Compute discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0  # no real roots => no intersection

    sqrt_disc = math.sqrt(disc)
    t0 = (-B - sqrt_disc) / (2.0*A)
    t1 = (-B + sqrt_disc) / (2.0*A)

    # We want the smallest non‐negative root
    t_near = None
    if t0 >= 0.0:
        t_near = t0
    elif t1 >= 0.0:
        t_near = t1

    return t_near if t_near is not None else -1.0


# Example usage:
if __name__ == "__main__":
    # Ray from (0,0,0) in direction (1,1,1)
    t = intersect_ellipsoid((0,0,0), (1,1,1))
    if t >= 0:
        print("Hit at t =", t, "point =", (0 + t, 0 + t, 0 + t))
    else:
        print("No intersection")


#---------- 8.8 ----------#
def intersect_ellipsoid(O, D):
    """
    Compute the intersection of the ray R(t) = O + t D
    with the axis-aligned ellipsoid centered at (1,1,1)
    with semi-axes lengths a=2, b=3, c=4.
    Returns the smallest positive t, or -1 if no intersection.
    """
    # Unpack
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Translate so ellipsoid center is at the origin
    x0 = Ox - 1.0
    y0 = Oy - 1.0
    z0 = Oz - 1.0

    # Squares of semi-axes
    a2 = 2.0 * 2.0
    b2 = 3.0 * 3.0
    c2 = 4.0 * 4.0

    # Quadratic coefficients A t^2 + B t + C = 0
    A = (Dx*Dx)/a2 + (Dy*Dy)/b2 + (Dz*Dz)/c2
    B = 2.0 * (Dx*x0/a2 + Dy*y0/b2 + Dz*z0/c2)
    C = (x0*x0)/a2 + (y0*y0)/b2 + (z0*z0)/c2 - 1.0

    # Discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0   # no real roots → no intersection

    # Compute the two roots
    sqrt_disc = disc**0.5
    t1 = (-B - sqrt_disc) / (2.0 * A)
    t2 = (-B + sqrt_disc) / (2.0 * A)

    # We only want intersections in front of the ray origin: t > 0
    ts = [t for t in (t1, t2) if t > 0.0]
    if not ts:
        return -1.0

    # Return the nearest one
    return min(ts)


#---------- 8.9 ----------#
import math

def intersect_ellipsoid(O, D):
    """
    Ray: R(t) = O + t*D, t >= 0
    Ellipsoid: centered at C = (1,1,1), semi-axes a=2, b=3, c=4
       (x-1)^2/2^2 + (y-1)^2/3^2 + (z-1)^2/4^2 = 1

    Returns the smallest t >= 0 at which the ray hits the ellipsoid,
    or -1.0 if there is no intersection.
    """
    # Ellipsoid parameters
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 2.0*2.0, 3.0*3.0, 4.0*4.0

    # Translate ray origin into ellipsoid‐centered coords
    Px = O[0] - Cx
    Py = O[1] - Cy
    Pz = O[2] - Cz
    Dx, Dy, Dz = D

    # Quadratic coefficients A t^2 + B t + C = 0
    A = (Dx*Dx)/a2 + (Dy*Dy)/b2 + (Dz*Dz)/c2
    B = 2.0*(Px*Dx/a2 + Py*Dy/b2 + Pz*Dz/c2)
    C = (Px*Px)/a2 + (Py*Py)/b2 + (Pz*Pz)/c2 - 1.0

    eps = 1e-8
    # Degenerate (linear) case
    if abs(A) < eps:
        if abs(B) < eps:
            return -1.0      # no solution or ray lies exactly on surface
        t = -C / B
        return t if t >= 0.0 else -1.0

    # Compute discriminant
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0

    sqrt_disc = math.sqrt(disc)
    t1 = (-B - sqrt_disc) / (2.0*A)
    t2 = (-B + sqrt_disc) / (2.0*A)

    # We want the smallest non-negative t
    t = float('inf')
    if t1 >= 0.0:
        t = t1
    if t2 >= 0.0 and t2 < t:
        t = t2

    return t if t != float('inf') else -1.0


#---------- 9.0 ----------#
import math

def intersect_cut_sphere(O, D, eps=1e-9):
    """
    O: tuple of 3 floats = ray origin (Ox, Oy, Oz)
    D: tuple of 3 floats = ray direction   (Dx, Dy, Dz)
    returns: smallest t>=0 so that O + t D hits the solid (sphere ∩ (x+y<=1))
             or the cut‐disk (x+y=1, x^2+y^2+z^2<=1),
             or –1 if no intersection.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) Sphere intersection: solve |O + t D|^2 = 1
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    ts = []
    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        # two roots
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        for t in (t1, t2):
            if t >= 0:
                # compute intersection point
                x = Ox + t*Dx
                y = Oy + t*Dy
                z = Oz + t*Dz
                # check that point lies on the kept side of the cut: x+y <= 1
                if x + y <= 1 + eps:
                    ts.append(t)

    # 2) Plane intersection (the cut‐disk): plane is x + y = 1
    denom = Dx + Dy
    if abs(denom) > eps:
        t_plane = (1.0 - (Ox + Oy)) / denom
        if t_plane >= 0:
            x = Ox + t_plane*Dx
            y = Oy + t_plane*Dy
            z = Oz + t_plane*Dz
            # by construction x+y == 1; just check it's inside the sphere
            if x*x + y*y + z*z <= 1.0 + eps:
                ts.append(t_plane)

    if not ts:
        return -1.0
    # return the nearest intersection
    return min(ts)


#---------- 9.1 ----------#
import math

def intersect_cut_sphere(O, D, eps=1e-9):
    """
    O, D : 3‐tuples of floats, ray R(t)=O + t D
    Returns the smallest t>=0 where R(t) hits the solid { x^2+y^2+z^2<=1 and x+y<=1 },
    or -1 if no intersection.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    hits = []

    # 1) Sphere intersection: |O + t D|^2 = 1
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1

    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        for t in (t1, t2):
            if t >= 0:
                px = Ox + t*Dx
                py = Oy + t*Dy
                pz = Oz + t*Dz
                # check we are on the kept side of the plane x+y<=1
                if px + py <= 1 + eps:
                    hits.append(t)

    # 2) Plane intersection: x+y=1
    denom = Dx + Dy
    if abs(denom) > eps:
        t = (1 - Ox - Oy) / denom
        if t >= 0:
            px = Ox + t*Dx
            py = Oy + t*Dy
            pz = Oz + t*Dz
            # check point lies inside the sphere
            if px*px + py*py + pz*pz <= 1 + eps:
                hits.append(t)

    if not hits:
        return -1.0
    return min(hits)


# Example usage:
if __name__ == "__main__":
    O = (2.0, 0.0, 0.0)
    D = (-1.0, 0.0, 0.0)   # shooting towards the origin
    print(intersect_cut_sphere(O, D))  # should hit somewhere around t=1.0


#---------- 9.2 ----------#
import math

def intersect_cut_sphere(O, D):
    """
    O : tuple of 3 floats, ray origin (Ox,Oy,Oz)
    D : tuple of 3 floats, ray direction (Dx,Dy,Dz)
    returns: smallest positive t for R(t)=O+t*D hitting the cut sphere,
             or -1.0 if no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    best_t = float('inf')

    # 1) Ray-sphere intersection (centered at 0,0,0, radius=1)
    #    Solve a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    disc = b*b - 4*a*c
    if disc >= 0.0:
        sqrt_disc = math.sqrt(disc)
        # two possible roots
        for t in [(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)]:
            if t > 0.0:
                # point on sphere
                Px = Ox + t*Dx
                Py = Oy + t*Dy
                Pz = Oz + t*Dz
                # only keep it if it's on the kept side x+y <= 1
                if Px + Py <= 1.0 + 1e-9:
                    if t < best_t:
                        best_t = t

    # 2) Ray-plane intersection with x+y=1
    #    plane normal n=(1,1,0), plane eqn: x+y-1=0
    denom = Dx + Dy
    if abs(denom) > 1e-9:
        t_plane = (1.0 - (Ox + Oy)) / denom
        if t_plane > 0.0 and t_plane < best_t:
            # compute plane hit point
            Px = Ox + t_plane*Dx
            Py = Oy + t_plane*Dy
            Pz = Oz + t_plane*Dz
            # must lie inside the circle of radius 1
            if Px*Px + Py*Py + Pz*Pz <= 1.0 + 1e-9:
                # also ensure the ray is actually entering the kept halfspace x+y<1
                # i.e. just before the plane the ray was in x+y>1 region
                # since (d/dt)(x+y) = Dx+Dy = denom,
                # we need denom < 0 so that x+y is decreasing through 1
                if denom < 0.0:
                    best_t = t_plane

    if best_t == float('inf'):
        return -1.0
    else:
        return best_t


#---------- 9.3 ----------#
import math

def intersect_ray_cut_sphere(O, D, ε=1e-9):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: the smallest t >= 0 at which R(t)=O+tD first hits
             either the remaining spherical cap or the circular
             cut, or -1 if there is no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    ts = []

    # 1) Intersections with the full sphere x^2+y^2+z^2 = 1
    #    Solve a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        for s in (-1, +1):
            t_s = (-b + s*sqrt_disc) / (2*a)
            if t_s >= 0:
                # check that the intersection point lies in the half‐space x+y >= 1
                x = Ox + t_s*Dx
                y = Oy + t_s*Dy
                if x + y >= 1 - ε:
                    ts.append(t_s)

    # 2) Intersection with the cutting plane x + y = 1
    #    Solve Ox + t Dx + Oy + t Dy = 1  =>  t = (1 - Ox - Oy)/(Dx + Dy)
    denom = Dx + Dy
    if abs(denom) > ε:
        t_plane = (1.0 - Ox - Oy) / denom
        if t_plane >= 0:
            # compute the point
            x = Ox + t_plane*Dx
            y = Oy + t_plane*Dy
            z = Oz + t_plane*Dz
            # check that it lies within the circle x^2 + y^2 + z^2 <= 1
            if x*x + y*y + z*z <= 1 + ε:
                ts.append(t_plane)

    if not ts:
        return -1.0

    return min(ts)


#---------- 9.4 ----------#
import math

def intersect_cut_sphere(O, D):
    """
    O, D: 3‐tuples or lists of floats (ox,oy,oz), (dx,dy,dz)
    Sphere: center (0,0,0), radius=1
    Plane: x + y = 1, we remove the smaller cap (the side containing the sphere center),
           so the remaining solid is { x^2+y^2+z^2 <= 1  AND  x+y-1 >= 0 }.
    Returns: the smallest t>0 at which the ray hits the boundary of that solid,
             or -1 if there is no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    eps = 1e-8

    # 1) Intersect the infinite sphere x^2+y^2+z^2=1
    #    Solve (O + t D)·(O + t D) = 1
    #    -> (D·D) t^2 + 2 (O·D) t + (O·O - 1) = 0
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0

    disc = b*b - 4*a*c
    candidates = []

    if disc >= 0 and a > eps:
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        # t1 <= t2

        for t in (t1, t2):
            if t <= eps:
                continue
            # point on sphere
            px = ox + t*dx
            py = oy + t*dy
            pz = oz + t*dz
            # check it lies on the kept side of the plane: x+y-1 >= 0
            if px + py - 1.0 >= -eps:
                candidates.append(t)

    # 2) Intersect the infinite plane x + y = 1
    #    Plane normal (1,1,0), eq: (x+y-1)=0
    #    Solve (O + tD)·n - 1 = 0  =>  (ox+oy-1) + t*(dx+dy) = 0
    denom = dx + dy
    if abs(denom) > eps:
        t_plane = (1.0 - (ox + oy)) / denom
        if t_plane > eps:
            # check the ray is entering the half‐space x+y-1 >= 0
            #   that is, D·n > 0 (so that for t>0 we go from x+y<1 side to x+y>1 side)
            if denom > 0:
                # check that (px,py,pz) lies within the circle x^2+y^2+z^2 <= 1
                px = ox + t_plane*dx
                py = oy + t_plane*dy
                pz = oz + t_plane*dz
                if px*px + py*py + pz*pz <= 1.0 + eps:
                    candidates.append(t_plane)

    if not candidates:
        return -1.0

    return min(candidates)

# Example usage:
if __name__ == "__main__":
    # shoot a ray from (2,2,0) toward the origin
    O = (2.0, 2.0, 0.0)
    D = (-1.0, -1.0, 0.0)
    # normalize D
    L = math.sqrt(D[0]*D[0] + D[1]*D[1] + D[2]*D[2])
    D = (D[0]/L, D[1]/L, D[2]/L)
    t = intersect_cut_sphere(O, D)
    print("t =", t)


#---------- 9.5 ----------#
import math

def intersect_cut_sphere(O, D, eps=1e-9):
    """
    Compute intersection of ray R(t)=O + t*D with the solid
    { (x,y,z): x^2+y^2+z^2 <= 1  and  x+y <= 1 }.
    Return the smallest t>0 at which the ray hits the object,
    or -1 if there is no intersection.
    O and D are 3‐tuples or lists of floats.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) Intersections with the full sphere x^2+y^2+z^2=1
    #    Solve: ||O + tD||^2 = 1  -->  a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2 * (Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1

    ts = []

    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        for t in (t1, t2):
            if t > eps:
                # check that the hit point is on the kept side of the plane: x+y <= 1
                x = Ox + Dx*t
                y = Oy + Dy*t
                if x + y <= 1 + eps:
                    ts.append(t)

    # 2) Intersection with the cutting plane x+y=1 (we must also be inside the sphere there)
    denom = Dx + Dy
    if abs(denom) > eps:
        t_plane = (1.0 - (Ox + Oy)) / denom
        if t_plane > eps:
            x = Ox + Dx*t_plane
            y = Oy + Dy*t_plane
            z = Oz + Dz*t_plane
            # check that (x,y,z) lies inside the sphere
            if x*x + y*y + z*z <= 1 + eps:
                # and on the plane it automatically satisfies x+y=1
                ts.append(t_plane)

    if not ts:
        return -1.0
    return min(ts)


#---------- 9.6 ----------#
import math

def intersect_cut_sphere(O, D):
    """
    O: tuple of 3 floats = ray origin (Ox, Oy, Oz)
    D: tuple of 3 floats = ray direction (Dx, Dy, Dz)
    returns: the smallest t>=0 so that O + t D lies on the surface of the
             cut sphere (sphere of radius 1 at origin, sliced by x+y=1,
             the x+y>1 cap removed), or -1 if no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # Quadratic coefficients for |O + t D|^2 = 1
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Dx*Ox + Dy*Oy + Dz*Oz)
    c = (Ox*Ox + Oy*Oy + Oz*Oz) - 1.0

    ts = []

    # 1) Solve sphere intersection
    disc = b*b - 4*a*c
    if disc >= 0.0:
        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        # check the two roots in ascending order
        for t in sorted((t0, t1)):
            if t < 0.0:
                continue
            # hit point
            px = Ox + t*Dx
            py = Oy + t*Dy
            pz = Oz + t*Dz
            # only count it if it's on the kept side of the cut
            # kept side is x+y <= 1
            if px + py <= 1.0 + 1e-9:
                ts.append(t)
                break  # we only need the first valid sphere hit

    # 2) Solve plane intersection x+y=1
    denom = Dx + Dy
    if abs(denom) > 1e-9:
        t_plane = (1.0 - Ox - Oy) / denom
        if t_plane >= 0.0:
            px = Ox + t_plane*Dx
            py = Oy + t_plane*Dy
            pz = Oz + t_plane*Dz
            # check it's inside the circular cut (radius = sqrt(1 - (dist of plane)^2) = 1/sqrt(2), but
            # more simply check x^2+y^2+z^2 <= 1)
            if (px*px + py*py + pz*pz) <= 1.0 + 1e-9:
                # this is the planar face of the cut, so it is a valid boundary
                ts.append(t_plane)

    if not ts:
        return -1.0

    return min(ts)


#---------- 9.7 ----------#
import math

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def mul(s, v):
    return (s*v[0], s*v[1], s*v[2])

def intersect_cut_sphere(O, D, eps=1e-6):
    """
    O, D: 3‐tuples of floats
    Returns: the smallest t >= 0 at which the ray O + t D hits
             the unit sphere x^2+y^2+z^2 <= 1 cut by the half‐space x+y<=1,
             or -1 if there is no intersection.
    """
    hits = []

    # ---- 1) Ray‐sphere intersection ----
    A = dot(D, D)
    B = 2.0 * dot(O, D)
    C = dot(O, O) - 1.0

    disc = B*B - 4*A*C
    if disc >= 0.0:
        sqrt_disc = math.sqrt(disc)
        t1 = (-B - sqrt_disc) / (2*A)
        t2 = (-B + sqrt_disc) / (2*A)
        for t in sorted((t1, t2)):
            if t >= eps:
                P = add(O, mul(t, D))
                # keep only if on the x+y<=1 side
                if P[0] + P[1] <= 1.0 + eps:
                    hits.append(t)

    # ---- 2) Ray‐plane intersection (x+y=1) ----
    denom = D[0] + D[1]
    if abs(denom) > eps:
        t_plane = (1.0 - (O[0] + O[1])) / denom
        if t_plane >= eps:
            Pp = add(O, mul(t_plane, D))
            # check if point lies within the unit‐radius disk
            if dot(Pp, Pp) <= 1.0 + eps:
                # also we know x+y==1 here, and this is a boundary face,
                # so no further half‐space test is needed
                hits.append(t_plane)

    if not hits:
        return -1.0

    return min(hits)


#---------- 9.8 ----------#
import math

def intersect_cut_sphere(O, D):
    """
    O: ray origin, 3‐tuple of floats (ox, oy, oz)
    D: ray direction, 3‐tuple of floats (dx, dy, dz)
    returns: smallest positive t at which R(t)=O+t*D hits the cut sphere,
             or -1 if no intersection.
    """

    ox, oy, oz = O
    dx, dy, dz = D
    eps = 1e-7

    candidates = []

    # 1) Sphere intersection (unit sphere at origin)
    #    Solve ||O + t D||^2 = 1
    a = dx*dx + dy*dy + dz*dz
    b = 2.0*(dx*ox + dy*oy + dz*oz)
    c = ox*ox + oy*oy + oz*oz - 1.0

    disc = b*b - 4*a*c
    if disc >= 0.0:
        sq = math.sqrt(disc)
        t1 = (-b - sq) / (2*a)
        t2 = (-b + sq) / (2*a)
        for t in (t1, t2):
            if t > eps:
                # test whether the hit point lies on the kept part of the sphere:
                # we only keep the portion with x+y >= 1
                x = ox + t*dx
                y = oy + t*dy
                z = oz + t*dz
                if x + y >= 1.0 - eps:
                    candidates.append(t)

    # 2) Plane intersection (the cutting face x+y=1), but only if we are coming
    #    from the removed‐cap side (x+y<1) and we hit inside the circular disk
    denom = dx + dy
    if abs(denom) > eps:
        # solve (ox + t*dx) + (oy + t*dy) = 1  =>  t = (1 - ox - oy)/(dx+dy)
        t_plane = (1.0 - (ox + oy)) / denom
        if t_plane > eps:
            x = ox + t_plane*dx
            y = oy + t_plane*dy
            z = oz + t_plane*dz
            # must lie inside the circle of intersection radius sqrt(1 - (1/√2)^2) = 1/√2
            # easier: just check x^2 + y^2 + z^2 <= 1
            if x*x + y*y + z*z <= 1.0 + eps:
                # also must start on the removed side x+y<1
                if ox + oy < 1.0 - eps:
                    candidates.append(t_plane)

    if not candidates:
        return -1.0

    return min(candidates)


# --- some simple tests ---

# Ray from far below in negative z, pointing up, should hit the sphere body (if x+y>=1)
print(intersect_cut_sphere((0.5,0.5,-5),(0,0,1)))   # hits the spherical cap, ~4.2929

# Ray from the removed region side, aiming at the cut‐face disk
print(intersect_cut_sphere((0.2,0.2,-1),(0,0,1)))   # should hit the planar disk at z=?


#---------- 9.9 ----------#
import math

def intersect_cut_sphere(O, D):
    """
    Intersect a ray R(t) = O + t*D with a unit sphere at the origin
    cut by the plane x+y=1 (the small cap where x+y>1 is removed).
    Returns the smallest t>0 of either:
      - intersection with the remaining part of the sphere surface, or
      - intersection with the flat disc at x+y=1
    If no intersection, returns -1.
    O, D: 3-tuples or lists of floats.
    """

    # small epsilon to guard against floating‐point wobbles
    EPS = 1e-7

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # -- 1) Sphere intersections ------------------------------------------------
    # Solve |O + t D|^2 = 1
    #   A t^2 + B t + C = 0
    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    C = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    ts = []

    disc = B*B - 4*A*C
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t0 = (-B - sqrt_disc) / (2*A)
        t1 = (-B + sqrt_disc) / (2*A)
        # for each root, check if it actually lies on the kept part of the sphere
        for t in (t0, t1):
            if t > EPS:
                # compute hit point
                px = Ox + t*Dx
                py = Oy + t*Dy
                pz = Oz + t*Dz
                # only keep it if x+y <= 1 (i.e. not in the removed cap x+y>1)
                if px + py <= 1.0 + EPS:
                    ts.append(t)

    # -- 2) Disc (cut) intersection ---------------------------------------------
    # Solve x+y = 1  =>  (Ox + t Dx) + (Oy + t Dy) = 1
    denom = Dx + Dy
    if abs(denom) > EPS:
        t_plane = (1.0 - Ox - Oy) / denom
        if t_plane > EPS:
            # point on the plane
            px = Ox + t_plane*Dx
            py = Oy + t_plane*Dy
            pz = Oz + t_plane*Dz
            # we only want the circular disk of radius 1 in that plane
            if px*px + py*py + pz*pz <= 1.0 + EPS:
                ts.append(t_plane)

    # -- pick the smallest positive t -------------------------------------------
    if not ts:
        return -1.0
    return min(ts)


#---------- 10.0 ----------#
def intersect(O, D):
    """
    Compute the first intersection t>=0 of the ray R(t)=O + t D
    with the solid “sphere A radius 1 at (0,0,0)” minus the
    overlapping cap (sphere A∩sphere B), where
      sphere B is radius 1 at (0.5, 0.5, 0).

    Returns the smallest t>=0 if there is an intersection, or -1 otherwise.
    """

    import math

    def ray_sphere(O, D, C, R):
        # Solve ‖O + t D − C‖^2 = R^2  ⇒  a t^2 + b t + c = 0
        # return list of real roots (may be empty, or one if tangent, or two)
        Ox, Oy, Oz = O
        Dx, Dy, Dz = D
        Cx, Cy, Cz = C
        # compute O−C
        ocx, ocy, ocz = Ox - Cx, Oy - Cy, Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2*(Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - R*R
        disc = b*b - 4*a*c
        if disc < 0:
            return []
        sq = math.sqrt(disc)
        t1 = (-b - sq)/(2*a)
        t2 = (-b + sq)/(2*a)
        if abs(t1 - t2) < 1e-9:
            return [t1]
        return sorted([t1, t2])

    # Sphere A
    tA = ray_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    if not tA:
        return -1.0

    # We only care about the portion of the ray inside A:
    tA_in = tA[0]
    tA_out = tA[-1]
    # If the entire A‐interval is behind the ray origin, no hit:
    if tA_out < 0:
        return -1.0

    # Sphere B (the “hole”)
    tB = ray_sphere(O, D, (0.5, 0.5, 0.0), 1.0)

    # Build up the list of candidate intervals in A minus the overlap with B∩A
    intervals = []

    if not tB:
        # No intersection with B ⇒ nothing to remove
        intervals.append((tA_in, tA_out))
    else:
        # We have an interval for B
        tB_in = tB[0]
        tB_out = tB[-1]
        # We only remove where the ray is inside B *and* inside A,
        # but since we do A∖B, we simply subtract [tB_in, tB_out] from [tA_in, tA_out].
        # Four cases of overlap:
        #  1) no overlap
        if tB_out <= tA_in or tB_in >= tA_out:
            intervals.append((tA_in, tA_out))
        #  2) B fully covers A
        elif tB_in <= tA_in and tB_out >= tA_out:
            # everything in A is removed
            intervals = []
        #  3) overlap on the left
        elif tB_in <= tA_in < tB_out < tA_out:
            intervals.append((tB_out, tA_out))
        #  4) overlap on the right
        elif tA_in < tB_in < tA_out <= tB_out:
            intervals.append((tA_in, tB_in))
        #  5) B interval strictly inside A interval
        else:  # tA_in < tB_in < tB_out < tA_out
            intervals.append((tA_in, tB_in))
            intervals.append((tB_out, tA_out))

    # From all remaining intervals pick the smallest t>=0
    best_t = None
    for (t0, t1) in intervals:
        # we want the entry point of each interval that lies in front of the origin
        if t1 < 0:
            continue
        # if t0<0<t1, we are inside the shape at t=0 → first hit is at t1?  No,
        # but the problem guarantees the ray starts *outside* the object,
        # so we need not worry about that case.
        if t0 >= 0:
            cand = t0
        else:
            cand = t1  # should not happen here, but for safety
        if cand >= 0 and (best_t is None or cand < best_t):
            best_t = cand

    return best_t if best_t is not None else -1.0


#---------- 10.1 ----------#
import math

def intersect_sphere(O, D, C, r):
    """
    Ray‐sphere intersection.
    O, D, C are 3‐tuples, r is radius.
    Solve |O + t D − C|^2 = r^2 for t.
    Returns (t0, t1) with t0 <= t1, or None if no real roots.
    """
    # compute quadratic coefficients
    # (D·D) t^2 + 2 D·(O−C) t + |O−C|^2 − r^2 = 0
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    Cx, Cy, Cz = C
    OCx, OCy, OCz = Ox-Cx, Oy-Cy, Oz-Cz

    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2*(Dx*OCx + Dy*OCy + Dz*OCz)
    Cq = OCx*OCx + OCy*OCy + OCz*OCz - r*r

    disc = B*B - 4*A*Cq
    if disc < 0:
        return None
    sqrt_disc = math.sqrt(disc)
    t0 = (-B - sqrt_disc) / (2*A)
    t1 = (-B + sqrt_disc) / (2*A)
    if t0 > t1:
        t0, t1 = t1, t0
    return (t0, t1)

def ray_sphere_with_hole(O, D):
    """
    Returns the distance t to the first positive intersection of
    the ray (O, D) with the unit sphere at (0,0,0) minus the
    overlapping region with the unit sphere at (0.5,0.5,0).
    If no intersection, returns -1.
    """
    # 1) Intersection with outer sphere S0
    t0 = intersect_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    if t0 is None:
        return -1.0
    t0_enter, t0_exit = t0

    # 2) Intersection with the hole‐sphere S1
    t1 = intersect_sphere(O, D, (0.5, 0.5, 0.0), 1.0)
    if t1 is None:
        # no hole intersection → the object is just [t0_enter, t0_exit]
        if t0_exit < 0:
            return -1.0
        return max(t0_enter, 0.0)

    t1_enter, t1_exit = t1

    # 3) Build the “object intervals” = [t0_enter, t0_exit] minus [t1_enter, t1_exit]
    intervals = []
    a, b = t0_enter, t0_exit
    c, d = t1_enter, t1_exit

    # If the hole‐interval [c,d] does not overlap [a,b], keep [a,b]
    if d <= a or c >= b:
        intervals.append((a, b))
    else:
        # we have some overlap, carve it out
        if c > a:
            intervals.append((a, min(c, b)))
        if d < b:
            intervals.append((max(d, a), b))

    # 4) From those intervals pick the smallest t ≥ 0
    best_t = None
    for ta, tb in intervals:
        if tb < 0:
            continue
        # entry point into that slab
        t_candidate = max(ta, 0.0)
        if best_t is None or t_candidate < best_t:
            best_t = t_candidate

    return best_t if best_t is not None else -1.0

# Example usage:
O = ( 2.0, 0.0, 0.0 )      # ray origin
D = (-1.0, 0.1, 0.0)      # ray direction (need not be normalized)
t = ray_sphere_with_hole(O, D)
if t > 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", P)
else:
    print("No intersection.")


#---------- 10.2 ----------#
import math

def ray_sphere_intervals(O, D, C, R):
    """
    Compute the intersection interval [t0, t1] (entry and exit) of
    ray R(t) = O + t D with the sphere centered at C with radius R.
    Returns (t0, t1) with t0 <= t1, or None if no real intersections.
    """
    # Shift ray origin into sphere's local coordinates
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    Cx, Cy, Cz = C

    # Quadratic coefficients for |(O + t D) - C|^2 = R^2
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2 * ( Dx*(Ox - Cx) + Dy*(Oy - Cy) + Dz*(Oz - Cz) )
    c = (Ox - Cx)**2 + (Oy - Cy)**2 + (Oz - Cz)**2 - R*R

    disc = b*b - 4*a*c
    if disc < 0:
        return None
    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 > t1:
        t0, t1 = t1, t0
    return (t0, t1)

def intersect(O, D):
    """
    Computes the intersection of ray O + t D with the CSG shape:
      big_sphere(center=(0,0,0),radius=1)
      minus
      hole_sphere(center=(0.5,0.5,0),radius=1)
    Returns the smallest t >= 0 where the ray first hits the shape,
    or -1 if there is no intersection.
    """
    # 1) Find intersections with the big sphere
    S_big = ray_sphere_intervals(O, D, (0.0,0.0,0.0), 1.0)
    if S_big is None:
        return -1.0
    t0_b, t1_b = S_big

    # If both intersections are behind the ray, no hit
    if t1_b < 0:
        return -1.0

    # Clip to t >= 0
    t0_b = max(t0_b, 0.0)

    # 2) Find intersections with the hole sphere
    S_hole = ray_sphere_intervals(O, D, (0.5,0.5,0.0), 1.0)
    if S_hole is None:
        # No hole => the shape is just the big sphere
        return t0_b

    t0_h, t1_h = S_hole
    # If hole interval is entirely behind or entirely after the big-sphere interval,
    # it does not affect the first entry.
    if t1_h < 0 or t0_h > t1_b:
        return t0_b

    # 3) Decide whether the first big-sphere entry point is inside the hole
    #    Evaluate the point at t0_b and see if it's inside the hole
    Px = O[0] + t0_b*D[0]
    Py = O[1] + t0_b*D[1]
    Pz = O[2] + t0_b*D[2]
    # squared distance to hole center
    dx = Px - 0.5
    dy = Py - 0.5
    dz = Pz - 0.0
    if dx*dx + dy*dy + dz*dz < 1.0:
        # The ray enters the big sphere through the hole region,
        # so the first *visible* intersection is when it exits the hole:
        # that's t1_h, but only if that's still inside the big sphere.
        if t1_h > t1_b:
            return -1.0
        if t1_h < 0:
            return -1.0
        return t1_h

    # Otherwise the big-sphere entry is outside the hole, so that's our hit
    return t0_b


#---------- 10.3 ----------#
import math

def intersect_sphere_with_hole(O, D):
    """
    O : tuple of 3 floats = ray origin
    D : tuple of 3 floats = ray direction (need not be normalized)
    returns: smallest t>0 such that O + t*D lies inside the unit sphere at (0,0,0)
             but outside the unit sphere at (0.5,0.5,0).
             If no such t exists, returns -1.
    """
    # list to hold all intersection distances t
    ts = []
    
    # helper to solve ray-sphere intersection and append roots to ts
    def add_sphere_intersections(center, radius):
        # shift ray origin to sphere coord
        Ox, Oy, Oz = O
        Dx, Dy, Dz = D
        Cx, Cy, Cz = center
        # compute quadratic coefficients for |O + t D - C|^2 = R^2
        ocx, ocy, ocz = Ox - Cx, Oy - Cy, Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2*(Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - radius*radius
        disc = b*b - 4*a*c
        if disc < 0:
            return
        sqrt_d = math.sqrt(disc)
        t0 = (-b - sqrt_d) / (2*a)
        t1 = (-b + sqrt_d) / (2*a)
        ts.append(t0)
        # if disc > 0 we have two distinct roots
        if disc > 0:
            ts.append(t1)
    
    # intersect with the big sphere
    add_sphere_intersections((0.0, 0.0, 0.0), 1.0)
    # intersect with the hole sphere
    add_sphere_intersections((0.5, 0.5, 0.0), 1.0)
    
    # sort all candidate ts
    ts.sort()
    
    # epsilon to avoid numerical self‐hits at t≈0
    EPS = 1e-7
    
    # test each candidate
    for t in ts:
        if t <= EPS:
            continue   # behind or too close to origin
        # compute point P = O + t D
        Px = O[0] + t*D[0]
        Py = O[1] + t*D[1]
        Pz = O[2] + t*D[2]
        # test inside big sphere?
        if Px*Px + Py*Py + Pz*Pz > 1.0 + EPS:
            continue
        # test outside hole sphere?
        dx = Px - 0.5
        dy = Py - 0.5
        dz = Pz - 0.0
        if dx*dx + dy*dy + dz*dz < 1.0 - EPS:
            continue
        # passed both tests: this is the first hit on the shaped object
        return t
    
    # no valid intersection found
    return -1.0


#---------- 10.4 ----------#
import math

def intersect_ray_holey_sphere(O, D):
    """
    O : tuple of (Ox, Oy, Oz)    ray origin
    D : tuple of (Dx, Dy, Dz)    ray direction (need not be unit length)
    returns: smallest positive t for which O + t D hits
             the sphere-with-spherical-hole, or -1 if no hit.
    """
    # sphere1: center C1=(0,0,0), radius R1=1
    # sphere2 (the "cutter"): center C2=(0.5,0.5,0), radius R2=1
    
    # helper: solve a t^2 + b t + c = 0
    def solve_quadratic(a, b, c):
        disc = b*b - 4*a*c
        if disc < 0:
            return None
        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        if t0 > t1:
            t0, t1 = t1, t0
        return (t0, t1)
    
    # ---------------------------------------
    # 1) intersect ray with big sphere at origin, R=1
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Dx*Ox + Dy*Oy + Dz*Oz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0
    sol1 = solve_quadratic(a, b, c)
    if sol1 is None:
        return -1.0
    t0, t1 = sol1
    # we assume ray starts outside => t0 > 0, t1 > t0
    # but just in case, require t1>0
    if t1 <= 0:
        return -1.0
    # our candidate big-sphere interval is [t_enter, t_exit]
    t_enter = t0 if t0 > 0 else 0.0
    t_exit  = t1
    
    # ---------------------------------------
    # 2) intersect ray with cutter sphere (center=(0.5,0.5,0), R=1)
    C2x, C2y, C2z = 0.5, 0.5, 0.0
    Ox2 = Ox - C2x
    Oy2 = Oy - C2y
    Oz2 = Oz - C2z
    b2 = 2*(Dx*Ox2 + Dy*Oy2 + Dz*Oz2)
    c2 = Ox2*Ox2 + Oy2*Oy2 + Oz2*Oz2 - 1.0
    sol2 = solve_quadratic(a, b2, c2)
    
    # if no intersections with cutter, the whole [t_enter, t_exit] survives
    if sol2 is None:
        return t_enter
    
    u0, u1 = sol2
    # now [u0,u1] is the interval along the ray where we are inside the cutter
    
    # ---------------------------------------
    # 3) compute the set‐difference
    #    big_sphere_interval minus cutter_interval
    #    leaves up to two intervals, but we only care about the very first point
    
    # If the cutter interval lies entirely after our exit, or entirely before our entry,
    # it doesn’t affect our first entry:
    if u0 >= t_exit or u1 <= t_enter:
        return t_enter
    
    # Otherwise the cutter actually chews a piece out of [t_enter, t_exit].
    # If the cutter starts after we already entered the big sphere,
    # then the portion [t_enter, u0) is valid solid.
    if u0 > t_enter:
        return t_enter
    
    # otherwise the cutter starts at or before our entry, so we must wait until u1
    # to exit the cutter.  If u1 < t_exit, that is the next solid point:
    if u1 < t_exit:
        return u1
    
    # if u1 >= t_exit then the cutter covers the entire [t_enter, t_exit],
    # so there is no remaining solid:
    return -1.0


#---------- 10.5 ----------#
def intersect(O, D):
    """
    O : tuple of 3 floats = ray origin
    D : tuple of 3 floats = ray direction
    returns the distance t>=0 to the first hit of
    (sphere@0,r=1) minus (sphere@(.5,.5,0),r=1), or -1.
    """

    def solve_sphere(O, D, C, R):
        # solve |O + t D - C|^2 = R^2
        # returns sorted (t0, t1), or None if no real roots
        ocx = O[0] - C[0]
        ocy = O[1] - C[1]
        ocz = O[2] - C[2]
        a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
        b = 2.0*(D[0]*ocx + D[1]*ocy + D[2]*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - R*R
        disc = b*b - 4*a*c
        if disc < 0.0:
            return None
        sd = disc**0.5
        t0 = (-b - sd) / (2*a)
        t1 = (-b + sd) / (2*a)
        if t0 > t1:
            t0, t1 = t1, t0
        return t0, t1

    # 1) intersection with the big sphere S1 @ (0,0,0), r=1
    S1 = solve_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    if S1 is None:
        return -1.0
    t0, t1 = S1
    # if both intersections behind the ray, give up
    if t1 < 0.0:
        return -1.0
    # we're only interested in t>=0
    t0 = max(t0, 0.0)

    # 2) intersection with the cutter sphere S2 @ (0.5,0.5,0), r=1
    S2 = solve_sphere(O, D, (0.5, 0.5, 0.0), 1.0)
    if S2 is None:
        # no cutting region, entire [t0,t1] is valid
        return t0

    s0, s1 = S2

    # 3) subtract [s0,s1] from [t0,t1].
    # if they don't overlap at all, we can keep t0
    if s1 < t0 or s0 > t1:
        return t0

    # if there's a piece before the cut
    if s0 > t0:
        return t0

    # otherwise maybe there's a piece after the cut
    if s1 < t1:
        # start of the "after" piece is max(s1,t0), but s1>=t0 here
        return s1

    # completely removed
    return -1.0


#---------- 10.6 ----------#
import math

def intersect_sphere_with_hole(O, D):
    """
    O: ray origin, 3‐tuple floats
    D: ray direction, 3‐tuple floats
    returns: the smallest t >= 0 so that O + t D ∈ (sphere0 minus sphere1),
             or -1 if no such intersection
    """

    # helper dot and subtract
    def dot(u,v):       return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def sub(u,v):       return (u[0]-v[0], u[1]-v[1], u[2]-v[2])

    # Ray vs sphere0 (center C0=(0,0,0), radius R0=1)
    C0 = (0.0, 0.0, 0.0)
    R0 = 1.0
    OC = sub(O, C0)
    A0 = dot(D, D)
    B0 = 2.0 * dot(D, OC)
    Cc0 = dot(OC, OC) - R0*R0
    disc0 = B0*B0 - 4.0*A0*Cc0
    if disc0 < 0.0:
        return -1.0    # no intersection with the big sphere

    sqrt0 = math.sqrt(disc0)
    t0 = (-B0 - sqrt0) / (2.0*A0)
    t1 = (-B0 + sqrt0) / (2.0*A0)
    # sort so t0 <= t1
    if t0 > t1:
        t0, t1 = t1, t0

    # if the whole hit is behind origin, no intersection
    if t1 < 0.0:
        return -1.0

    # clamp entry to ≥ 0
    a = max(t0, 0.0)
    b = t1

    # Ray vs sphere1 (the “hole”): center C1=(0.5,0.5,0), radius R1=1
    C1 = (0.5, 0.5, 0.0)
    R1 = 1.0
    OC1 = sub(O, C1)
    A1 = A0  # same D
    B1 = 2.0 * dot(D, OC1)
    Cc1 = dot(OC1, OC1) - R1*R1
    disc1 = B1*B1 - 4.0*A1*Cc1

    # if no intersection with sphere1, then there is no hole on this ray
    if disc1 < 0.0:
        # the only solid interval is [a,b]
        return a if a <= b else -1.0

    sqrt1 = math.sqrt(disc1)
    s0 = (-B1 - sqrt1) / (2.0*A1)
    s1 = (-B1 + sqrt1) / (2.0*A1)
    if s0 > s1:
        s0, s1 = s1, s0

    # now [s0, s1] is the interval inside the hole
    # we subtract [s0,s1] from [a,b]:
    # compute overlap:
    over_start = max(a, s0)
    over_end   = min(b, s1)

    if over_end <= over_start:
        # no overlap → entire [a,b] is solid
        return a

    # otherwise we may have up to two surviving sub‐intervals:
    candidates = []
    # part before hole
    if a < over_start:
        candidates.append(a)
    # part after hole
    if over_end < b:
        candidates.append(over_end)

    # pick the smallest t ≥ 0
    candidates = [t for t in candidates if t >= 0.0]
    if not candidates:
        return -1.0
    return min(candidates)


#---------- 10.7 ----------#
import math

def intersect_ray_hollow_sphere(O, D):
    """
    O : tuple of 3 floats  -- ray origin
    D : tuple of 3 floats  -- ray direction (need not be unit length)
    returns t>=0 of the first intersection, or -1.0 if none.
    """

    # Sphere 1: center C1, radius R1
    C1 = (0.0, 0.0, 0.0)
    R1 = 1.0
    # Sphere 2 (the "hole"): center C2, radius R2
    C2 = (0.5, 0.5, 0.0)
    R2 = 1.0

    def solve_quadratic_sphere(O, D, C, R):
        """
        Solve (O + t D − C)·(O + t D − C) = R^2
        Returns (t0, t1) sorted if intersections exist, else None.
        """
        ox, oy, oz = O
        dx, dy, dz = D
        cx, cy, cz = C

        # shift so sphere is at origin
        ux = ox - cx
        uy = oy - cy
        uz = oz - cz

        a = dx*dx + dy*dy + dz*dz
        b = 2.0*(dx*ux + dy*uy + dz*uz)
        c = ux*ux + uy*uy + uz*uz - R*R

        disc = b*b - 4.0*a*c
        if disc < 0.0:
            return None

        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2.0*a)
        t1 = (-b + sqrt_disc) / (2.0*a)
        if t0 > t1:
            t0, t1 = t1, t0
        return (t0, t1)

    # 1) find entry/exit with the big sphere
    big = solve_quadratic_sphere(O, D, C1, R1)
    if big is None:
        return -1.0        # no hit on big sphere at all

    t_big0, t_big1 = big
    # if the entire big‐sphere interval is behind us, no intersection
    if t_big1 < 0.0:
        return -1.0

    # clamp entry to t≥0
    t_enter_big = max(t_big0, 0.0)
    t_exit_big  = t_big1

    # 2) find entry/exit with the small (hole) sphere
    small = solve_quadratic_sphere(O, D, C2, R2)
    if small is None:
        # no hole intersection, so the first hit is just the big‐sphere entry
        return t_enter_big

    t_hole0, t_hole1 = small

    # If the hole‐interval [t_hole0, t_hole1] does not overlap [t_enter_big, t_exit_big],
    # then the hole doesn't carve anything out of our visible segment.
    if t_hole1 < t_enter_big or t_hole0 > t_exit_big:
        return t_enter_big

    # Now they do overlap.  We have two possible “visible” sub‐intervals:
    #   A = [t_enter_big,  min(t_exit_big, t_hole0)]
    #   B = [max(t_enter_big, t_hole1),  t_exit_big]
    #
    # We just want the smallest t ≥ 0 in (A ∪ B), if either is non‐empty.

    # Case A in front of the hole:
    if t_enter_big < t_hole0:
        # A is non‐empty, and its first point is t_enter_big
        return t_enter_big

    # Otherwise t_enter_big ≥ t_hole0, so our entry into the big sphere
    # is already inside the hole.  We must wait until we exit the hole:
    if t_hole1 <= t_exit_big:
        # after t_hole1 we are back in the material
        # as long as that exit is ≥ 0
        return t_hole1 if t_hole1 >= 0.0 else -1.0

    # If we get here it means the hole covers the rest of the big sphere segment
    # (i.e. t_hole1 > t_exit_big), so there's no material left
    return -1.0


#---------- 10.8 ----------#
import math

def intersect_holed_sphere(O, D):
    """
    Ray: R(t) = O + t D,  t >= 0
    Object = unit sphere at (0,0,0) radius=1
             with the solid cap overlapping
             the unit sphere at (0.5,0.5,0) removed.
    O, D: 3‐tuples of floats.
    Returns: smallest t >= 0 of intersection, or -1 if none.
    """

    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

    def solve_sphere(O, D, C, R):
        """
        Solve ‖O + tD - C‖² = R² for t.
        Returns sorted list [t0, t1], or None if no real roots.
        """
        # Move to sphere‐center coords
        Ocx = O[0] - C[0]
        Ocy = O[1] - C[1]
        Ocz = O[2] - C[2]
        a = dot(D, D)
        b = 2*(D[0]*Ocx + D[1]*Ocy + D[2]*Ocz)
        c = Ocx*Ocx + Ocy*Ocy + Ocz*Ocz - R*R
        disc = b*b - 4*a*c
        if disc < 0:
            return None
        # one or two solutions
        sqrt_d = math.sqrt(disc)
        t0 = (-b - sqrt_d) / (2*a)
        t1 = (-b + sqrt_d) / (2*a)
        if t0 > t1:
            t0, t1 = t1, t0
        return (t0, t1)

    # 1) Intersection with S1 (center 0,0,0 radius=1)
    s1 = solve_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    if s1 is None:
        return -1.0
    t1_0, t1_1 = s1

    # 2) Intersection with S2 (center .5,.5,0 radius=1)
    s2 = solve_sphere(O, D, (0.5, 0.5, 0.0), 1.0)
    if s2 is None:
        # no removal → the whole [t1_0, t1_1] stands
        rem = None
    else:
        rem = s2  # removal interval

    # Now build the list of remaining intervals = [t1_0, t1_1] \ rem
    intervals = []
    eps = 1e-9

    if rem is None:
        # no subtraction
        intervals.append((t1_0, t1_1))
    else:
        t2_0, t2_1 = rem
        # before the removed chunk:
        if t2_0 > t1_0 + eps:
            intervals.append((t1_0, min(t1_1, t2_0)))
        # after the removed chunk:
        if t2_1 < t1_1 - eps:
            intervals.append((max(t1_0, t2_1), t1_1))

    # 3) From these intervals pick smallest t >= 0
    t_hit = None
    for (ta, tb) in intervals:
        if tb < 0:
            continue  # entire sub‐interval behind the ray origin
        t_candidate = max(ta, 0.0)
        # check validity
        if t_candidate <= tb + eps:
            if t_hit is None or t_candidate < t_hit:
                t_hit = t_candidate

    return t_hit if t_hit is not None else -1.0


#---------- 10.9 ----------#
def ray_sphere_hole(O, D):
    """
    O: ray origin, 3‐tuple of floats
    D: ray direction, 3‐tuple of floats (need not be normalized)
    returns: the smallest t>0 at which R(t)=O+tD hits the solid
             (sphere centered at (0,0,0) radius=1 with the
             intersection-with-sphere(center=(0.5,0.5,0),r=1) removed),
             or -1 if no intersection.
    """
    # small epsilon
    EPS = 1e-9

    # helper: dot product
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    # helper: subtract
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])

    # solve quadratic a t^2 + b t + c = 0
    def solve_quad(a, b, c):
        disc = b*b - 4*a*c
        if disc < -EPS:
            return None
        if disc < 0:
            disc = 0.0
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        if t0 > t1:
            t0, t1 = t1, t0
        return (t0, t1)

    # 1) Intersection with big sphere A (center Oa=(0,0,0), r=1)
    #    ||O + t D||^2 = 1  =>  a t^2 + b t + c = 0
    a = dot(D, D)
    b = 2 * dot(D, O)
    c = dot(O, O) - 1.0
    rootsA = solve_quad(a, b, c)
    if rootsA is None:
        return -1.0
    tA0, tA1 = rootsA

    # we assume the ray starts outside the final shape.
    # if tA1 < 0, both intersections are behind us => no hit
    if tA1 < EPS:
        return -1.0
    # clamp entry to zero if it straddles
    if tA0 < 0:
        tA0 = 0.0

    # 2) Intersection with cutting sphere B (center C=(0.5,0.5,0), r=1)
    C = (0.5, 0.5, 0.0)
    OminusC = sub(O, C)
    bB = 2 * dot(D, OminusC)
    cB = dot(OminusC, OminusC) - 1.0
    rootsB = solve_quad(a, bB, cB)

    # now we have A‐interval [tA0,tA1] and maybe B‐interval [tB0,tB1].
    # we want A minus (A∩B).  That can produce up to two intervals:
    #   I1 = [tA0, tB0], I2 = [tB1, tA1], provided those are valid.

    # if B has no intersection, the entire A‐interval is solid:
    if rootsB is None:
        return tA0 if tA0 > EPS else -1.0  # tA0≥0 by above clamp

    tB0, tB1 = rootsB

    # If the B‐interval doesn't overlap A at all, same
    if tB1 < tA0 + EPS or tB0 > tA1 - EPS:
        # no overlap, A is untouched
        return tA0

    # Otherwise they do overlap, carve out the overlap
    # overlap = [max(tA0,tB0), min(tA1,tB1)]
    ov0 = max(tA0, tB0)
    ov1 = min(tA1, tB1)

    # candidate intervals in the final solid:
    #   seg1 = [tA0, ov0]
    #   seg2 = [ov1, tA1]
    t_candidates = []
    if ov0 > tA0 + EPS:
        t_candidates.append(tA0)
    if ov1 < tA1 - EPS:
        # if we start inside the removed chunk, the next entry is at ov1
        # but we must make sure ov1 >= 0
        if ov1 > EPS:
            t_candidates.append(ov1)

    if not t_candidates:
        return -1.0

    # return the smallest positive candidate
    tmin = min(t_candidates)
    return tmin if tmin > EPS else -1.0


#---------- 11.0 ----------#
import math

def intersect_sphere_with_hole(O, D):
    """
    O, D: tuples or lists of 3 floats
    returns the smallest t >= 0 at which the ray O + t D hits the sphere-with-hole,
    or -1 if there is no intersection.
    """
    ox, oy, oz = O
    dx, dy, dz = D

    eps = 1e-6
    candidates = []

    # 1) intersections with the unit sphere x^2 + y^2 + z^2 = 1
    #    but only keep those points where y^2+z^2 >= 0.7^2 (i.e. outside the drilled‐out cylinder)
    A = dx*dx + dy*dy + dz*dz
    B = 2*(ox*dx + oy*dy + oz*dz)
    C = ox*ox + oy*oy + oz*oz - 1.0

    disc = B*B - 4*A*C
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t_s1 = (-B - sqrt_disc) / (2*A)
        t_s2 = (-B + sqrt_disc) / (2*A)
        for t in (t_s1, t_s2):
            if t >= eps:
                # point on sphere
                y = oy + t*dy
                z = oz + t*dz
                if y*y + z*z >= 0.7*0.7 - eps:
                    candidates.append(t)

    # 2) intersections with the infinite cylinder y^2 + z^2 = 0.7^2
    #    but only keep those points where x^2+y^2+z^2 <= 1 (i.e. inside the sphere)
    A = dy*dy + dz*dz
    B = 2*(oy*dy + oz*dz)
    C = oy*oy + oz*oz - 0.7*0.7

    # if A is zero then the ray is parallel to the cylinder axis (no cyl. intersection)
    if A > eps:
        disc = B*B - 4*A*C
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            t_c1 = (-B - sqrt_disc) / (2*A)
            t_c2 = (-B + sqrt_disc) / (2*A)
            for t in (t_c1, t_c2):
                if t >= eps:
                    # point on cylinder
                    x = ox + t*dx
                    y = oy + t*dy
                    z = oz + t*dz
                    if x*x + y*y + z*z <= 1.0 + eps:
                        candidates.append(t)

    if not candidates:
        return -1.0

    return min(candidates)


#---------- 11.1 ----------#
import math

def intersect(O, D):
    """
    O : tuple of 3 floats = ray origin (Ox,Oy,Oz)
    D : tuple of 3 floats = ray direction (Dx,Dy,Dz)
    returns the smallest t >= 0 such that
      P = O + t D lies on the sphere-with-hole boundary,
      or -1 if there is no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    eps = 1e-9
    hits = []

    # 1) Intersect with the sphere x^2 + y^2 + z^2 = 1
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    disc = b*b - 4*a*c
    if disc >= 0:
        sdisc = math.sqrt(disc)
        t_sph = [(-b - sdisc)/(2*a), (-b + sdisc)/(2*a)]
        for t in t_sph:
            if t >= 0:
                # point on sphere
                y = Oy + t*Dy
                z = Oz + t*Dz
                # must lie outside the cylindrical hole y^2+z^2 >= 0.7^2
                if y*y + z*z >= 0.49 - eps:
                    hits.append(t)

    # 2) Intersect with the infinite cylinder y^2 + z^2 = 0.7^2
    #    (axis along x).  We only keep those points that lie
    #    inside the sphere x^2+y^2+z^2 <= 1.
    A = Dy*Dy + Dz*Dz
    B = 2*(Oy*Dy + Oz*Dz)
    C = Oy*Oy + Oz*Oz - 0.49

    if abs(A) > eps:
        disc2 = B*B - 4*A*C
        if disc2 >= 0:
            s2 = math.sqrt(disc2)
            t_cyl = [(-B - s2)/(2*A), (-B + s2)/(2*A)]
            for t in t_cyl:
                if t >= 0:
                    x = Ox + t*Dx
                    y = Oy + t*Dy
                    z = Oz + t*Dz
                    # must lie inside the sphere
                    if x*x + y*y + z*z <= 1.0 + eps:
                        hits.append(t)

    if not hits:
        return -1.0
    return min(hits)


#---------- 11.2 ----------#
import math

def intersect_ray_sphere_with_hole(O, D):
    """
    O: (Ox,Oy,Oz) ray origin
    D: (Dx,Dy,Dz) ray direction (need not be normalized)
    returns: smallest t>=0 where ray hits the sphere-with-hole,
             or -1.0 if no intersection.
    """

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersect with the unit sphere |O+tD|^2 = 1
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    disc_s = b*b - 4*a*c
    if disc_s < 0.0:
        return -1.0        # no sphere intersection

    sqrt_s = math.sqrt(disc_s)
    t_s0 = (-b - sqrt_s) / (2*a)
    t_s1 = (-b + sqrt_s) / (2*a)
    if t_s0 > t_s1:
        t_s0, t_s1 = t_s1, t_s0

    # clamp to t >= 0
    if t_s1 < 0.0:
        return -1.0
    if t_s0 < 0.0:
        t_s0 = 0.0

    # 2) intersect with the infinite cylinder y^2+z^2 = 0.7^2
    a_c = Dy*Dy + Dz*Dz
    has_cyl = False
    if a_c > 1e-8:
        b_c = 2*(Oy*Dy + Oz*Dz)
        c_c = Oy*Oy + Oz*Oz - 0.7*0.7
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0.0:
            sqrt_c = math.sqrt(disc_c)
            t_c0 = (-b_c - sqrt_c) / (2*a_c)
            t_c1 = (-b_c + sqrt_c) / (2*a_c)
            if t_c0 > t_c1:
                t_c0, t_c1 = t_c1, t_c0
            has_cyl = True

    # 3) subtract [t_c0,t_c1] from [t_s0,t_s1]
    intervals = []
    if not has_cyl or t_c1 <= t_s0 or t_c0 >= t_s1:
        # no overlap ⇒ whole sphere interval survives
        intervals.append((t_s0, t_s1))
    else:
        # piece before the hole
        if t_c0 > t_s0:
            intervals.append((t_s0, min(t_c0, t_s1)))
        # piece after the hole
        if t_c1 < t_s1:
            intervals.append((max(t_c1, t_s0), t_s1))

    # 4) pick the smallest entry t among the remaining intervals
    if not intervals:
        return -1.0

    t_min = min(iv[0] for iv in intervals)
    return t_min


#---------- 11.3 ----------#
import math

def intersect_sphere_with_cyl_hole(O, D):
    """
    O: ray origin, tuple (Ox, Oy, Oz)
    D: ray direction, tuple (Dx, Dy, Dz)
    returns: the smallest t >= 0 at which R(t)=O+t*D hits the sphere-with-hole,
             or -1 if no intersection.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    eps = 1e-9
    hits = []

    # --- 1) Intersect with the unit sphere x^2 + y^2 + z^2 = 1 ---
    # solve a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0

    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        for t in (t0, t1):
            if t >= 0:
                # point on sphere
                y = Oy + Dy*t
                z = Oz + Dz*t
                # only keep it if it's outside the hole: y^2+z^2 >= 0.7^2
                if y*y + z*z >= 0.7*0.7 - eps:
                    hits.append(t)

    # --- 2) Intersect with the infinite cylinder y^2 + z^2 = (0.7)^2 ---
    # (axis along x), but only keep those points that lie inside the sphere
    a_c = Dy*Dy + Dz*Dz
    b_c = 2*(Oy*Dy + Oz*Dz)
    c_c = Oy*Oy + Oz*Oz - 0.7*0.7

    # if a_c is zero, the ray is parallel to the cylinder axis => no side‐wall hit
    if abs(a_c) > eps:
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0:
            sqrt_dc = math.sqrt(disc_c)
            tc0 = (-b_c - sqrt_dc) / (2*a_c)
            tc1 = (-b_c + sqrt_dc) / (2*a_c)
            for t in (tc0, tc1):
                if t >= 0:
                    # point on cylinder
                    x = Ox + Dx*t
                    y = Oy + Dy*t
                    z = Oz + Dz*t
                    # only keep it if it's inside the sphere: x^2+y^2+z^2 <= 1
                    if x*x + y*y + z*z <= 1.0 + eps:
                        hits.append(t)

    if not hits:
        return -1.0

    return min(hits)


#---------- 11.4 ----------#
import math

def intersect_sphere_with_hole(O, D,
                              R_sph=1.0,
                              R_cyl=0.7,
                              eps=1e-7):
    """
    O, D: 3‐tuples or lists of floats, ray origin and direction
    returns: smallest t>0 at which the ray hits
             the solid (sphere minus cylinder),
             or -1 if no hit.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) Sphere intersection
    #    Solve (D·D) t^2 + 2(O·D) t + (O·O - R_sph^2) = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - R_sph*R_sph

    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc)/(2*a)
        t2 = (-b + sqrt_disc)/(2*a)
        # try them in ascending order
        for t in sorted((t1, t2)):
            if t <= eps:
                continue
            # point on sphere
            Px = Ox + t*Dx
            Py = Oy + t*Dy
            Pz = Oz + t*Dz
            # is it outside the cylinder (i.e. in the kept part of the sphere)?
            if (Py*Py + Pz*Pz) > (R_cyl*R_cyl + eps):
                return t

    # 2) Cylinder‐side intersection (inner wall of the hole)
    #    Solve (Oy+t Dy)^2 + (Oz+t Dz)^2 = R_cyl^2
    a = Dy*Dy + Dz*Dz
    b = 2*(Oy*Dy + Oz*Dz)
    c = Oy*Oy + Oz*Oz - R_cyl*R_cyl

    cyl_ts = []
    if abs(a) > eps:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = math.sqrt(disc)
            tc1 = (-b - sqrt_disc)/(2*a)
            tc2 = (-b + sqrt_disc)/(2*a)
            cyl_ts = sorted((tc1, tc2))
    # else ray is parallel to cyl axis: no new side hits

    # Among the positive t's, pick the first one that lies inside the sphere
    best_tc = None
    for t in cyl_ts:
        if t <= eps:
            continue
        # point on cylinder
        Px = Ox + t*Dx
        Py = Oy + t*Dy
        Pz = Oz + t*Dz
        # only valid if also inside the sphere
        if Px*Px + Py*Py + Pz*Pz <= R_sph*R_sph + eps:
            best_tc = t
            break

    if best_tc is not None:
        return best_tc

    # no intersection
    return -1.0


#---------- 11.5 ----------#
import math

def intersect_sphere_with_hole(O, D):
    """
    O, D: 3‐tuples (Ox,Oy,Oz), (Dx,Dy,Dz)
    Ray: R(t) = O + t*D, t >= 0
    Object: sphere radius=1 minus cylinder radius=0.7 about x‐axis
    Returns: smallest t >= 0 of a valid intersection, or -1 if none.
    """
    Rs = 1.0    # sphere radius
    Rc = 0.7    # cylinder radius
    eps = 1e-9  # small tolerance
    
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    
    candidates = []
    
    # 1) Intersect with sphere: solve a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - Rs*Rs
    
    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        t_s1 = (-b - sqrt_disc)/(2*a)
        t_s2 = (-b + sqrt_disc)/(2*a)
        for t in (t_s1, t_s2):
            if t >= 0:
                # compute the hit point
                y = Oy + t*Dy
                z = Oz + t*Dz
                # accept this sphere‐hit only if it's outside the cylinder‐hole
                if y*y + z*z >= Rc*Rc - eps:
                    candidates.append(t)
    
    # 2) Intersect with cylinder y^2 + z^2 = Rc^2
    #    solve (Dy^2 + Dz^2) t^2 + 2(Oy Dy + Oz Dz)t + (Oy^2+Oz^2 - Rc^2) = 0
    a_c = Dy*Dy + Dz*Dz
    if a_c > eps:  # ray not parallel to cylinder axis in YZ‐plane
        b_c = 2*(Oy*Dy + Oz*Dz)
        c_c = Oy*Oy + Oz*Oz - Rc*Rc
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0:
            sqrt_c = math.sqrt(disc_c)
            t_c1 = (-b_c - sqrt_c)/(2*a_c)
            t_c2 = (-b_c + sqrt_c)/(2*a_c)
            for t in (t_c1, t_c2):
                if t >= 0:
                    # compute the hit point
                    x = Ox + t*Dx
                    y = Oy + t*Dy
                    z = Oz + t*Dz
                    # accept this cylinder‐hit only if it's inside the sphere
                    if x*x + y*y + z*z <= Rs*Rs + eps:
                        candidates.append(t)
    
    if not candidates:
        return -1.0
    return min(candidates)


#---------- 11.6 ----------#
import math

def intersect_holey_sphere(O, D):
    """
    O, D are 3‐tuples or lists of floats.
    Returns the smallest t>=0 so that O+t*D hits the unit sphere 
    with a cylindrical hole of radius 0.7 along the x-axis.
    Returns -1 if there is no such intersection.
    """
    # sphere radius
    R_s = 1.0
    # cylinder radius
    R_c = 0.7

    # --- 1) intersect ray with sphere: |O + t D|^2 = R_s^2 ---
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2*(O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = (O[0]*O[0] + O[1]*O[1] + O[2]*O[2]) - R_s*R_s

    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0    # no intersection with sphere at all

    sqrt_disc = math.sqrt(disc)
    t_s1 = (-b - sqrt_disc) / (2*a)
    t_s2 = (-b + sqrt_disc) / (2*a)

    # we assume the ray starts outside, so t_s2 is the "exit" and must be ≥ 0
    if t_s2 < 0:
        return -1.0

    # clamp the "entry" to t>=0
    t0 = max(t_s1, 0.0)
    t1 = t_s2

    # --- 2) intersect ray with infinite cylinder y^2+z^2 = R_c^2 ---
    # solve (O.y + t D.y)^2 + (O.z + t D.z)^2 = R_c^2
    a2 = D[1]*D[1] + D[2]*D[2]
    b2 = 2*(O[1]*D[1] + O[2]*D[2])
    c2 = (O[1]*O[1] + O[2]*O[2]) - R_c*R_c

    # special case: ray parallel to the x-axis => a2 == 0
    if abs(a2) < 1e-9:
        # if c2 > 0 then O is outside the cylinder and never hits it
        if c2 > 0:
            return t0
        else:
            # if c2 <=0 then the ray lies inside the cylinder for all t
            # so the entire [t0,t1] is punched out
            return -1.0

    disc2 = b2*b2 - 4*a2*c2
    if disc2 < 0:
        # no intersection with cylinder at all
        return t0

    sqrt2 = math.sqrt(disc2)
    t_c1 = (-b2 - sqrt2) / (2*a2)
    t_c2 = (-b2 + sqrt2) / (2*a2)
    if t_c1 > t_c2:
        t_c1, t_c2 = t_c2, t_c1

    # --- 3) subtract [t_c1,t_c2] from [t0,t1] and pick the smallest remaining t≥0 ---
    # if the cylinder interval lies completely outside the sphere interval, just return t0
    if t_c2 < t0 or t_c1 > t1:
        return t0

    # otherwise there is some overlap
    # if the ray hits the sphere *before* it hits the cylinder hole, that's the first point
    if t0 < t_c1:
        return t0

    # else the first bit of sphere surface is punched out, check the second bit
    if t_c2 < t1:
        return t_c2

    # otherwise the entire [t0,t1] got removed
    return -1.0


#---------- 11.7 ----------#
import math

def intersect_sphere(O, D, R=1.0):
    # Solve ||O + t D||^2 = R^2
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    b = 2*(O[0]*D[0] + O[1]*D[1] + O[2]*D[2])
    c = O[0]*O[0] + O[1]*O[1] + O[2]*O[2] - R*R
    disc = b*b - 4*a*c
    if disc < 0:
        return None
    s = math.sqrt(disc)
    t1 = (-b - s)/(2*a)
    t2 = (-b + s)/(2*a)
    return (min(t1, t2), max(t1, t2))

def intersect_cylinder_xaxis(O, D, rc=0.7, eps=1e-9):
    # Cylinder: y^2 + z^2 = rc^2
    # Solve (Oy + t Dy)^2 + (Oz + t Dz)^2 = rc^2
    A = D[1]*D[1] + D[2]*D[2]
    B = 2*(O[1]*D[1] + O[2]*D[2])
    C = O[1]*O[1] + O[2]*O[2] - rc*rc

    if abs(A) < eps:
        # Ray parallel to cylinder axis
        if C > 0:
            # Always outside cylinder → no intersection
            return None
        else:
            # Always inside cylinder → interval is whole real line
            return (-math.inf, math.inf)

    disc = B*B - 4*A*C
    if disc < 0:
        # No real roots → ray never crosses cylinder surface
        # Then either always outside (C>0) or always inside (C<0)
        if C > 0:
            # always outside → no intersection interval
            return None
        else:
            # always inside
            return (-math.inf, math.inf)

    s = math.sqrt(disc)
    t1 = (-B - s)/(2*A)
    t2 = (-B + s)/(2*A)
    return (min(t1, t2), max(t1, t2))


def ray_object_intersection(O, D):
    """
    Returns the smallest t>=0 for which O + t D lies inside the sphere
    of radius 1 but outside the cylinder of radius 0.7 around the x‐axis.
    Returns -1 if there is no such intersection.
    """
    eps = 1e-9

    # 1) find where ray enters/exits sphere
    sph = intersect_sphere(O, D, R=1.0)
    if sph is None:
        return -1.0
    ts0, ts1 = sph
    # If the entire sphere‐interval is behind us, no hit
    if ts1 < 0:
        return -1.0

    # 2) find where ray enters/exits cylinder (y^2+z^2=0.7^2)
    cyl = intersect_cylinder_xaxis(O, D, rc=0.7, eps=eps)
    # cyl is either None (ray always outside cyl),
    #       or (-inf,inf) (always inside cyl),
    #       or (tc0, tc1)

    # 3) form the intervals along t where the ray is inside the sphere:
    sphere_interval = (ts0, ts1)

    # 4) subtract out the cylinder‐interior interval from the sphere interval
    valid_intervals = []
    if cyl is None:
        # never inside cylinder → entire sphere interval counts
        valid_intervals.append(sphere_interval)
    elif cyl[0] == -math.inf and cyl[1] == math.inf:
        # always inside cylinder → sphere minus cylinder = empty
        valid_intervals = []
    else:
        tc0, tc1 = cyl
        # if no overlap
        if tc1 <= ts0 or tc0 >= ts1:
            valid_intervals.append((ts0, ts1))
        else:
            # overlap exists: carve out [max(ts0,tc0), min(ts1,tc1)]
            if tc0 > ts0:
                valid_intervals.append((ts0, min(tc0, ts1)))
            if tc1 < ts1:
                valid_intervals.append((max(tc1, ts0), ts1))

    # 5) among those remaining intervals, pick the smallest t>=0
    t_hit = math.inf
    for (a, b) in valid_intervals:
        # if interval is entirely behind us, skip
        if b < 0:
            continue
        # otherwise the first point ≥0 is max(a, 0)
        t0 = max(a, 0.0)
        if t0 <= b + eps and t0 < t_hit:
            t_hit = t0

    if t_hit == math.inf:
        return -1.0
    return t_hit


# Example usage:
O = (-2.0, 0.0, 0.0)       # ray origin
D = (1.0, 0.1, 0.2)        # ray direction (need not be normalized)
t = ray_object_intersection(O, D)
if t > 0:
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    print("Hit at t =", t, "point =", P)
else:
    print("No hit")


#---------- 11.8 ----------#
import math

def solve_quadratic(a, b, c, eps=1e-9):
    """
    Solve a*t^2 + b*t + c = 0.
    Returns a list of real roots (0, 1 or 2 entries).
    """
    if abs(a) < eps:
        # linear case b t + c = 0
        if abs(b) < eps:
            return []
        return [ -c/b ]
    disc = b*b - 4*a*c
    if disc < 0:
        return []
    elif abs(disc) < eps:
        return [ -b / (2*a) ]
    else:
        sd = math.sqrt(disc)
        return [ (-b - sd)/(2*a), (-b + sd)/(2*a) ]


def intersect_sphere_with_cylindrical_hole(O, D):
    """
    O, D: 3–tuples of floats, ray origin and direction.
    Returns the smallest t>0 at which R(t)=O+tD
    hits the sphere-with-hole, or -1 if no hit.
    """
    R_sphere = 1.0
    R_cyl    = 0.7

    Ox, Oy, Oz = O
    Dx, Dy, Dz = D

    # 1) intersections with the unit sphere x^2+y^2+z^2 = 1
    a_s = Dx*Dx + Dy*Dy + Dz*Dz
    b_s = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c_s = Ox*Ox + Oy*Oy + Oz*Oz - R_sphere*R_sphere
    ts_sphere = solve_quadratic(a_s, b_s, c_s)

    # 2) intersections with the infinite cylinder y^2+z^2 = 0.7^2
    a_c = Dy*Dy + Dz*Dz
    b_c = 2*(Oy*Dy + Oz*Dz)
    c_c = Oy*Oy + Oz*Oz - R_cyl*R_cyl
    ts_cyl = solve_quadratic(a_c, b_c, c_c)

    candidates = []

    # Check each sphere‐intersection: accept it only if
    # point is outside the cylinder (y^2+z^2 >= 0.7^2).
    for t in ts_sphere:
        if t <= 0:
            continue
        y = Oy + t*Dy
        z = Oz + t*Dz
        if y*y + z*z >= R_cyl*R_cyl - 1e-9:
            candidates.append(t)

    # Check each cylinder‐intersection: accept it only if
    # point is inside the sphere (x^2+y^2+z^2 <= 1).
    for t in ts_cyl:
        if t <= 0:
            continue
        x = Ox + t*Dx
        y = Oy + t*Dy
        z = Oz + t*Dz
        if x*x + y*y + z*z <= R_sphere*R_sphere + 1e-9:
            candidates.append(t)

    if not candidates:
        return -1.0

    return min(candidates)


#---------- 11.9 ----------#
import math

def intersect_sphere_with_hole(O, D,
                              sphere_radius = 1.0,
                              cyl_radius   = 0.7):
    """
    O, D : tuple or list of 3 floats  -- ray origin and (normalized) direction
    sphere_radius : float             -- radius of the outer sphere
    cyl_radius    : float             -- radius of the removed cylinder (about x-axis)
    
    Returns the smallest t>0 so that O + t D lies on the boundary of
      (Sphere of radius sphere_radius) MINUS (infinite cylinder about x with radius cyl_radius),
    or -1 if no such t exists.
    """
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    
    # ---- 1) Solve sphere intersection: |O + t D|^2 = sphere_radius^2
    #    a t^2 + b t + c = 0
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - sphere_radius*sphere_radius
    
    disc_s = b*b - 4*a*c
    if disc_s < 0:
        return -1.0     # no real intersection with sphere at all
    
    # we allow grazing => disc_s == 0 is fine
    sqrt_disc_s = math.sqrt(max(0.0, disc_s))
    t_s1 = (-b - sqrt_disc_s)/(2*a)
    t_s2 = (-b + sqrt_disc_s)/(2*a)
    # we want the smaller positive root first
    ts = []
    if t_s1 > 1e-8: ts.append(t_s1)
    if t_s2 > 1e-8: ts.append(t_s2)
    
    # ---- 2) Solve cylinder intersection: (Oy + t Dy)^2 + (Oz + t Dz)^2 = cyl_radius^2
    A = Dy*Dy + Dz*Dz
    B = 2*(Oy*Dy + Oz*Dz)
    C = Oy*Oy + Oz*Oz - cyl_radius*cyl_radius
    
    tc = []
    if abs(A) > 1e-12:
        disc_c = B*B - 4*A*C
        if disc_c >= 0:
            sqrt_disc_c = math.sqrt(max(0.0, disc_c))
            t_c1 = (-B - sqrt_disc_c)/(2*A)
            t_c2 = (-B + sqrt_disc_c)/(2*A)
            if t_c1 > 1e-8: tc.append(t_c1)
            if t_c2 > 1e-8: tc.append(t_c2)
    else:
        # Ray parallel to the cylinder axis (x-axis).
        # Then either never intersects (if Oy^2+Oz^2 != r^2)
        # or grazes everywhere (infinitely many), but grazing we treat as no boundary.
        tc = []
    
    # ---- 3) Collect all candidate hits:
    #   (a) sphere‐surface hits that lie *outside* the cylinder (i.e. hole region removed)
    #   (b) cylinder‐surface hits that lie *inside* the sphere
    candidates = []
    
    # (a) sphere hits
    for t in ts:
        x = Ox + t*Dx
        y = Oy + t*Dy
        z = Oz + t*Dz
        # only if outside the cylinder
        if (y*y + z*z) >= cyl_radius*cyl_radius:
            candidates.append(t)
    
    # (b) cylinder hits
    for t in tc:
        x = Ox + t*Dx
        y = Oy + t*Dy
        z = Oz + t*Dz
        # only if inside the sphere
        if (x*x + y*y + z*z) <= sphere_radius*sphere_radius:
            candidates.append(t)
    
    if not candidates:
        return -1.0
    
    return min(candidates)


# ------ quick smoke test ------
if __name__ == "__main__":
    # Ray straight at the center of the sphere along -z, starts at z=+5
    O = (0.0, 0.0, 5.0)
    D = (0.0, 0.0,-1.0)
    print(intersect_sphere_with_hole(O,D))   # should hit the sphere at t=4.0 (z=1)

    # Now nudge into the hole: start slightly off-center so that y^2+z^2<0.7^2
    O = (0.0, 0.5, 5.0)
    D = (0.0, 0.0,-1.0)
    print(intersect_sphere_with_hole(O,D))
    # This should return the t where the ray hits the cylinder wall inside the sphere,
    # i.e. solve y^2+z^2=0.7^2: 0.5^2 + (5-t)^2 = 0.7^2  => (5-t)^2 = 0.49-0.25=0.24  => t≈5±√0.24
    # the smaller positive ≈5 - 0.4899 = 4.5101



