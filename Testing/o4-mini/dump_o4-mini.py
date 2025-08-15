#---------- 0.0 ----------#
def intersect_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if a == 0 or disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 > 0:
        return t0
    if t1 > 0:
        return t1
    return -1.0

#---------- 0.1 ----------#
def intersect_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2*(Dx*Ox + Dy*Oy + Dz*Oz)
    C = Ox*Ox + Oy*Oy + Oz*Oz - 1.0
    disc = B*B - 4*A*C
    if disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t1 = (-B - sqrt_disc) / (2*A)
    t2 = (-B + sqrt_disc) / (2*A)
    t_near, t_far = (t1, t2) if t1 < t2 else (t2, t1)
    if t_near >= 0:
        return t_near
    if t_far >= 0:
        return t_far
    return -1.0

#---------- 0.2 ----------#
def intersect_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2 * (ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 > 0:
        return t0
    if t1 > 0:
        return t1
    return -1.0

#---------- 0.3 ----------#
def intersect_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2 * (dx*ox + dy*oy + dz*oz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 >= 0 and t1 >= 0:
        return min(t0, t1)
    if t0 >= 0:
        return t0
    if t1 >= 0:
        return t1
    return -1


#---------- 0.4 ----------#
def intersect_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2 * (Dx*Ox + Dy*Oy + Dz*Oz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    ts = sorted((t0, t1))
    for t in ts:
        if t >= 0:
            return t
    return -1


#---------- 0.5 ----------#
def intersect_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 >= 0:
        return t0
    if t1 >= 0:
        return t1
    return -1

#---------- 0.6 ----------#
def intersect_sphere(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1
    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    t_min = min(t0, t1)
    if t_min > 0:
        return t_min
    t_max = max(t0, t1)
    if t_max > 0:
        return t_max
    return -1


#---------- 0.7 ----------#
def intersect_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    t = None
    if t1 >= 0:
        t = t1
    if t2 >= 0 and (t is None or t2 < t):
        t = t2
    return t if t is not None else -1.0


#---------- 0.8 ----------#
def intersect_ray_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t0 = (-b - sqrt_disc) / (2*a)
    t1 = (-b + sqrt_disc) / (2*a)
    if t0 > 0 and t1 > 0:
        return min(t0, t1)
    if t0 > 0:
        return t0
    if t1 > 0:
        return t1
    return -1.0

#---------- 0.9 ----------#
def intersect_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    t = min(t1, t2)
    if t >= 0:
        return t
    t = max(t1, t2)
    if t >= 0:
        return t
    return -1.0

#---------- 1.0 ----------#
def intersect_cube(O, D):
    tmin, tmax = float('-inf'), float('inf')
    for o, d in zip(O, D):
        if abs(d) < 1e-8:
            if o < 0 or o > 1:
                return -1.0
        else:
            t1 = (0 - o) / d
            t2 = (1 - o) / d
            t_near, t_far = min(t1, t2), max(t1, t2)
            tmin = max(tmin, t_near)
            tmax = min(tmax, t_far)
            if tmin > tmax:
                return -1.0
    if tmax < 0:
        return -1.0
    return tmin if tmin >= 0 else tmax

#---------- 1.1 ----------#
def intersect_cube(O, D):
    tmin = -float('inf')
    tmax = float('inf')
    for o, d in zip(O, D):
        if d != 0.0:
            t1 = (0.0 - o) / d
            t2 = (1.0 - o) / d
            t_near = min(t1, t2)
            t_far = max(t1, t2)
            tmin = max(tmin, t_near)
            tmax = min(tmax, t_far)
            if tmin > tmax:
                return -1.0
        else:
            if o < 0.0 or o > 1.0:
                return -1.0
    return tmin if tmin >= 0.0 else -1.0

#---------- 1.2 ----------#
def intersect_cube(O, D):
    tmin = -float('inf')
    tmax = float('inf')
    for i in range(3):
        if D[i] != 0.0:
            t1 = (0.0 - O[i]) / D[i]
            t2 = (1.0 - O[i]) / D[i]
            tmin = max(tmin, min(t1, t2))
            tmax = min(tmax, max(t1, t2))
            if tmax < tmin:
                return -1.0
        else:
            if O[i] < 0.0 or O[i] > 1.0:
                return -1.0
    if tmax < 0.0:
        return -1.0
    return tmin if tmin >= 0.0 else tmax

#---------- 1.3 ----------#
def intersect_cube(O, D):
    tmin, tmax = float('-inf'), float('inf')
    bounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    for i in range(3):
        o_i, d_i = O[i], D[i]
        bmin, bmax = bounds[i]
        if abs(d_i) < 1e-8:
            if o_i < bmin or o_i > bmax:
                return -1.0
        else:
            t1 = (bmin - o_i) / d_i
            t2 = (bmax - o_i) / d_i
            t_near, t_far = min(t1, t2), max(t1, t2)
            tmin, tmax = max(tmin, t_near), min(tmax, t_far)
            if tmax < tmin:
                return -1.0
    return tmin if tmin >= 0.0 else -1.0


#---------- 1.4 ----------#
def intersect_cube(O, D):
    tmin = float('-inf')
    tmax = float('inf')
    for o, d in zip(O, D):
        if d == 0:
            if o < 0 or o > 1:
                return -1.0
        else:
            t1 = (0 - o) / d
            t2 = (1 - o) / d
            t_enter = min(t1, t2)
            t_exit = max(t1, t2)
            tmin = max(tmin, t_enter)
            tmax = min(tmax, t_exit)
            if tmin > tmax:
                return -1.0
    if tmax < 0:
        return -1.0
    return tmin if tmin >= 0 else tmax


#---------- 1.5 ----------#
def intersect_cube(O, D):
    inf = float('inf')
    tmin, tmax = -inf, inf
    for o, d in zip(O, D):
        if abs(d) < 1e-8:
            if o < 0 or o > 1:
                return -1.0
            continue
        t1 = (0 - o) / d
        t2 = (1 - o) / d
        t_near, t_far = min(t1, t2), max(t1, t2)
        tmin = max(tmin, t_near)
        tmax = min(tmax, t_far)
        if tmin > tmax:
            return -1.0
    if tmax < 0:
        return -1.0
    return tmin if tmin >= 0 else tmax

#---------- 1.6 ----------#
def intersect_ray_cube(O, D):
    import math
    tmin, tmax = -math.inf, math.inf
    for i in range(3):
        if abs(D[i]) < 1e-8:
            if O[i] < 0.0 or O[i] > 1.0:
                return -1.0
        else:
            t1 = (0.0 - O[i]) / D[i]
            t2 = (1.0 - O[i]) / D[i]
            tnear, tfar = min(t1, t2), max(t1, t2)
            tmin, tmax = max(tmin, tnear), min(tmax, tfar)
            if tmin > tmax:
                return -1.0
    if tmax < 0:
        return -1.0
    return tmin if tmin >= 0 else tmax


#---------- 1.7 ----------#
def intersect_cube(O, D):
    tmin, tmax = float('-inf'), float('inf')
    bounds = [(0.0, 1.0)] * 3
    for o, d, (b0, b1) in zip(O, D, bounds):
        if abs(d) < 1e-8:
            if o < b0 or o > b1:
                return -1.0
            continue
        t1 = (b0 - o) / d
        t2 = (b1 - o) / d
        t_near, t_far = min(t1, t2), max(t1, t2)
        tmin = max(tmin, t_near)
        tmax = min(tmax, t_far)
        if tmin > tmax:
            return -1.0
    if tmax < 0.0:
        return -1.0
    return tmin if tmin >= 0.0 else tmax


#---------- 1.8 ----------#
def intersect_cube(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    t_min, t_max = -float('inf'), float('inf')
    for o, d, mn, mx in ((ox, dx, 0.0, 1.0), (oy, dy, 0.0, 1.0), (oz, dz, 0.0, 1.0)):
        if d != 0.0:
            t1 = (mn - o) / d
            t2 = (mx - o) / d
            t_near, t_far = min(t1, t2), max(t1, t2)
        else:
            if o < mn or o > mx:
                return -1.0
            t_near, t_far = -float('inf'), float('inf')
        t_min = max(t_min, t_near)
        t_max = min(t_max, t_far)
        if t_min > t_max:
            return -1.0
    return t_min if t_min >= 0.0 else -1.0


#---------- 1.9 ----------#
def intersect_cube(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    t_min, t_max = -float('inf'), float('inf')
    for Oi, Di in ((Ox, Dx), (Oy, Dy), (Oz, Dz)):
        if Di != 0:
            t1 = (0 - Oi) / Di
            t2 = (1 - Oi) / Di
            t_near, t_far = min(t1, t2), max(t1, t2)
            t_min = max(t_min, t_near)
            t_max = min(t_max, t_far)
            if t_max < t_min:
                return -1
        else:
            if Oi < 0 or Oi > 1:
                return -1
    return t_min if t_min >= 0 else -1

#---------- 2.0 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Intersection with plane z=0
    if Dz == 0.0:
        return -1.0
    t = -Oz / Dz
    if t < 0.0:
        return -1.0
    # Compute intersection point
    x = Ox + t * Dx
    y = Oy + t * Dy
    # Check if inside square [0,1]x[0,1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.1 ----------#
def intersect_ray_square(O, D):
    eps = 1e-6
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    if abs(Dz) < eps:
        return -1
    t = -Oz / Dz
    if t < 0:
        return -1
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    if -eps <= Px <= 1 + eps and -eps <= Py <= 1 + eps:
        return t
    return -1


#---------- 2.2 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Ray-plane intersection with z=0 plane
    if abs(Dz) < 1e-8:
        return -1.0
    t = -Oz / Dz
    if t < 0.0:
        return -1.0
    # Intersection point
    x = Ox + t * Dx
    y = Oy + t * Dy
    # Check if point is inside the square [0,1] x [0,1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.3 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    if Dz == 0.0:
        return -1.0
    t = -Oz / Dz
    if t < 0.0:
        return -1.0
    x = Ox + t * Dx
    y = Oy + t * Dy
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.4 ----------#
def intersect_square(O, D):
    if D[2] == 0.0:
        return -1.0
    t = -O[2] / D[2]
    if t < 0.0:
        return -1.0
    x = O[0] + t * D[0]
    y = O[1] + t * D[1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.5 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Ray is parallel to the square's plane
    if abs(Dz) < 1e-8:
        return -1.0
    # Compute intersection with z=0 plane
    t = -Oz / Dz
    if t < 0:
        return -1.0
    x = Ox + t * Dx
    y = Oy + t * Dy
    # Check if the intersection point lies within the square bounds
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.6 ----------#
def intersect_ray_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    eps = 1e-6
    if abs(Dz) < eps:
        return -1.0
    t = -Oz / Dz
    if t < 0:
        return -1.0
    x = Ox + t * Dx
    y = Oy + t * Dy
    if -eps <= x <= 1 + eps and -eps <= y <= 1 + eps:
        return t
    return -1.0


#---------- 2.7 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    if Dz == 0:
        return -1.0
    t = -Oz / Dz
    if t <= 0:
        return -1.0
    x = Ox + t * Dx
    y = Oy + t * Dy
    if 0 <= x <= 1 and 0 <= y <= 1:
        return t
    return -1.0


#---------- 2.8 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Ray parallel to square's plane
    if Dz == 0:
        return -1.0
    # Compute intersection t with plane z=0
    t = -Oz / Dz
    if t < 0:
        return -1.0
    # Compute the intersection point
    x = Ox + t * Dx
    y = Oy + t * Dy
    # Check if the point lies within the square [0,1] x [0,1]
    if 0.0 <= x <= 1.0 and 0.0 <= y <= 1.0:
        return t
    return -1.0


#---------- 2.9 ----------#
def intersect_square(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Plane z = 0 normal (0,0,1)
    if abs(Dz) < 1e-8:
        return -1.0
    t = -Oz / Dz
    if t < 0:
        return -1.0
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    if 0.0 <= Px <= 1.0 and 0.0 <= Py <= 1.0:
        return t
    return -1.0


#---------- 3.0 ----------#
def intersect_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []
    # side intersections
    a = Dx*Dx + Dz*Dz
    b = 2*(Ox*Dx + Oz*Dz)
    c = Ox*Ox + Oz*Oz - 1
    if a != 0:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = disc**0.5
            for t in ((-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)):
                if t > 0:
                    y = Oy + t*Dy
                    if 0 <= y <= 1:
                        ts.append(t)
    # cap intersections y=0 and y=1
    if Dy != 0:
        for y_plane in (0, 1):
            t = (y_plane - Oy) / Dy
            if t > 0:
                x = Ox + t*Dx
                z = Oz + t*Dz
                if x*x + z*z <= 1:
                    ts.append(t)
    return min(ts) if ts else -1.0

#---------- 3.1 ----------#
def intersect_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    t_vals = []
    # Intersection with cylindrical surface
    a = Dx*Dx + Dz*Dz
    if a != 0.0:
        b = 2*(Ox*Dx + Oz*Dz)
        c = Ox*Ox + Oz*Oz - 1
        disc = b*b - 4*a*c
        if disc >= 0.0:
            sqrt_disc = disc**0.5
            inv2a = 0.5 / a
            for sign in (-1, 1):
                t = (-b + sign * sqrt_disc) * inv2a
                if t >= 0.0:
                    y = Oy + t * Dy
                    if 0.0 <= y <= 1.0:
                        t_vals.append(t)
    # Intersection with caps at y=0 and y=1
    if Dy != 0.0:
        for y_plane in (0.0, 1.0):
            t = (y_plane - Oy) / Dy
            if t >= 0.0:
                x = Ox + t * Dx
                z = Oz + t * Dz
                if x*x + z*z <= 1.0:
                    t_vals.append(t)
    if not t_vals:
        return -1.0
    return min(t_vals)


#---------- 3.2 ----------#
def intersect_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []
    # side intersections
    a = Dx*Dx + Dz*Dz
    if abs(a) > 1e-8:
        b = 2*(Ox*Dx + Oz*Dz)
        c = Ox*Ox + Oz*Oz - 1
        disc = b*b - 4*a*c
        if disc >= 0:
            sd = disc**0.5
            for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
                if t > 0:
                    y = Oy + t*Dy
                    if 0 <= y <= 1:
                        ts.append(t)
    # cap intersections
    if abs(Dy) > 1e-8:
        for y_plane in (0.0, 1.0):
            t = (y_plane - Oy)/Dy
            if t > 0:
                x = Ox + t*Dx
                z = Oz + t*Dz
                if x*x + z*z <= 1:
                    ts.append(t)
    return min(ts) if ts else -1.0


#---------- 3.3 ----------#
def intersect_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []
    # lateral surface
    A = Dx*Dx + Dz*Dz
    B = 2*(Ox*Dx + Oz*Dz)
    C = Ox*Ox + Oz*Oz - 1
    if abs(A) > 1e-8:
        disc = B*B - 4*A*C
        if disc >= 0:
            sqrt_disc = disc**0.5
            for t in [(-B - sqrt_disc)/(2*A), (-B + sqrt_disc)/(2*A)]:
                if t >= 0:
                    y = Oy + Dy*t
                    if 0 <= y <= 1:
                        ts.append(t)
    # caps
    if abs(Dy) > 1e-8:
        for y_plane in (0.0, 1.0):
            t = (y_plane - Oy)/Dy
            if t >= 0:
                x = Ox + Dx*t
                z = Oz + Dz*t
                if x*x + z*z <= 1:
                    ts.append(t)
    return min(ts) if ts else -1

#---------- 3.4 ----------#
def intersect_cylinder(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    t_min = float('inf')
    a = dx*dx + dz*dz
    b = 2*(ox*dx + oz*dz)
    c = ox*ox + oz*oz - 1
    if a != 0:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = disc**0.5
            for t in [(-b - sqrt_disc) / (2*a), (-b + sqrt_disc) / (2*a)]:
                if t > 0:
                    y = oy + t*dy
                    if 0 <= y <= 1 and t < t_min:
                        t_min = t
    if dy != 0:
        for y_cap in (0, 1):
            t = (y_cap - oy) / dy
            if t > 0:
                x = ox + t*dx
                z = oz + t*dz
                if x*x + z*z <= 1 and t < t_min:
                    t_min = t
    return t_min if t_min != float('inf') else -1

#---------- 3.5 ----------#
def intersect_cylinder(O, D):
    import math
    x0, y0, z0 = O
    dx, dy, dz = D
    ts = []
    a = dx*dx + dz*dz
    if a != 0.0:
        b = 2*(x0*dx + z0*dz)
        c = x0*x0 + z0*z0 - 1
        disc = b*b - 4*a*c
        if disc >= 0.0:
            sd = math.sqrt(disc)
            for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
                if t >= 0.0:
                    y = y0 + dy*t
                    if 0.0 <= y <= 1.0:
                        ts.append(t)
    if dy != 0.0:
        for y_plane in (0.0, 1.0):
            t = (y_plane - y0)/dy
            if t >= 0.0:
                x = x0 + dx*t
                z = z0 + dz*t
                if x*x + z*z <= 1.0:
                    ts.append(t)
    return min(ts) if ts else -1.0

#---------- 3.6 ----------#
def intersect_cylinder(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []

    # Side intersections (infinite cylinder x^2 + z^2 = 1)
    a = Dx*Dx + Dz*Dz
    if a != 0.0:
        b = 2.0*(Ox*Dx + Oz*Dz)
        c = Ox*Ox + Oz*Oz - 1.0
        disc = b*b - 4.0*a*c
        if disc >= 0.0:
            sd = math.sqrt(disc)
            for t in [(-b - sd)/(2.0*a), (-b + sd)/(2.0*a)]:
                if t > 0.0:
                    y = Oy + t*Dy
                    if 0.0 <= y <= 1.0:
                        ts.append(t)

    # Cap intersections (y = 0 and y = 1)
    if Dy != 0.0:
        for y_plane in (0.0, 1.0):
            t = (y_plane - Oy) / Dy
            if t > 0.0:
                x = Ox + t*Dx
                z = Oz + t*Dz
                if x*x + z*z <= 1.0:
                    ts.append(t)

    if not ts:
        return -1.0
    return min(ts)


#---------- 3.7 ----------#
def intersect_cylinder(O, D):
    x0, y0, z0 = O
    dx, dy, dz = D
    ts = []
    a = dx*dx + dz*dz
    b = 2*(x0*dx + z0*dz)
    c = x0*x0 + z0*z0 - 1
    if abs(a) > 1e-6:
        disc = b*b - 4*a*c
        if disc >= 0:
            sd = disc**0.5
            t1 = (-b - sd) / (2*a)
            t2 = (-b + sd) / (2*a)
            for t in (t1, t2):
                if t > 0:
                    y = y0 + dy*t
                    if 0 <= y <= 1:
                        ts.append(t)
    if abs(dy) > 1e-6:
        for y_plane in (0, 1):
            t = (y_plane - y0) / dy
            if t > 0:
                x = x0 + dx*t
                z = z0 + dz*t
                if x*x + z*z <= 1:
                    ts.append(t)
    return min(ts) if ts else -1.0


#---------- 3.8 ----------#
def intersect_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    t_vals = []
    A = Dx*Dx + Dz*Dz
    if A != 0.0:
        B = 2*(Ox*Dx + Oz*Dz)
        C = Ox*Ox + Oz*Oz - 1.0
        disc = B*B - 4*A*C
        if disc >= 0.0:
            sd = disc**0.5
            for t in [(-B - sd)/(2*A), (-B + sd)/(2*A)]:
                if t > 0.0:
                    y = Oy + t*Dy
                    if 0.0 <= y <= 1.0:
                        t_vals.append(t)
    for y_plane in (0.0, 1.0):
        if Dy != 0.0:
            t = (y_plane - Oy)/Dy
            if t > 0.0:
                x = Ox + t*Dx
                z = Oz + t*Dz
                if x*x + z*z <= 1.0:
                    t_vals.append(t)
    return min(t_vals) if t_vals else -1.0

#---------- 3.9 ----------#
def intersect_cylinder(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dz*dz
    b = 2*(ox*dx + oz*dz)
    c = ox*ox + oz*oz - 1
    ts = []
    # side intersections
    if abs(a) > 1e-6:
        disc = b*b - 4*a*c
        if disc >= 0:
            sqrt_disc = disc**0.5
            for t in [(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)]:
                if t > 0:
                    y = oy + t*dy
                    if 0 <= y <= 1:
                        ts.append(t)
    # cap intersections
    if abs(dy) > 1e-6:
        for y_plane in (0.0, 1.0):
            t = (y_plane - oy)/dy
            if t > 0:
                x = ox + t*dx
                z = oz + t*dz
                if x*x + z*z <= 1:
                    ts.append(t)
    return min(ts) if ts else -1

#---------- 4.0 ----------#
def intersect(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    best_t = float('inf')
    for cx, cy, cz in ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)):
        ocx = ox - cx
        ocy = oy - cy
        ocz = oz - cz
        a = dx*dx + dy*dy + dz*dz
        b = 2.0*(dx*ocx + dy*ocy + dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - 1.0
        disc = b*b - 4*a*c
        if disc < 0:
            continue
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        if t0 >= 0 and t0 < best_t:
            best_t = t0
        if t1 >= 0 and t1 < best_t:
            best_t = t1
    return best_t if best_t != float('inf') else -1.0

#---------- 4.1 ----------#
def intersect_two_spheres(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    best_t = float('inf')
    for Cx, Cy, Cz in ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)):
        ocx = Ox - Cx
        ocy = Oy - Cy
        ocz = Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2 * (Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - 1.0
        disc = b*b - 4*a*c
        if disc < 0:
            continue
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        if t0 >= 0 and t0 < best_t:
            best_t = t0
        if t1 >= 0 and t1 < best_t:
            best_t = t1
    return best_t if best_t != float('inf') else -1.0

#---------- 4.2 ----------#
def intersect_two_spheres(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    best_t = float('inf')
    for Cx, Cy, Cz in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]:
        ocx = Ox - Cx
        ocy = Oy - Cy
        ocz = Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2*(Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - 1.0
        disc = b*b - 4*a*c
        if disc < 0:
            continue
        sqrt_disc = disc**0.5
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        for t in (t1, t2):
            if t >= 0 and t < best_t:
                best_t = t
    return best_t if best_t != float('inf') else -1.0

#---------- 4.3 ----------#
def intersect_two_spheres(O, D):
    def intersect_sphere(O, D, C, r):
        Ox, Oy, Oz = O
        Dx, Dy, Dz = D
        Cx, Cy, Cz = C
        ocx, ocy, ocz = Ox - Cx, Oy - Cy, Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2*(Dx*ocx + Dy*ocy + Dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - r*r
        disc = b*b - 4*a*c
        if disc < 0:
            return None
        sqrt_disc = disc**0.5
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        ts = [t for t in (t1, t2) if t >= 0]
        return min(ts) if ts else None

    t0 = intersect_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    t1 = intersect_sphere(O, D, (1.0, 0.0, 0.0), 1.0)
    ts = [t for t in (t0, t1) if t is not None]
    return min(ts) if ts else -1.0


#---------- 4.4 ----------#
def intersect(O, D):
    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
    t_vals = []
    for C in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]:
        oc = (O[0] - C[0], O[1] - C[1], O[2] - C[2])
        b = 2 * (D[0]*oc[0] + D[1]*oc[1] + D[2]*oc[2])
        c = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2] - 1.0
        disc = b*b - 4*a*c
        if disc < 0:
            continue
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        if t0 > 0:
            t_vals.append(t0)
        if t1 > 0:
            t_vals.append(t1)
    if not t_vals:
        return -1.0
    return min(t_vals)


#---------- 4.5 ----------#
def intersect_spheres(O, D):
    import math
    def sphere_intersect(C, r):
        ox, oy, oz = O
        dx, dy, dz = D
        cx, cy, cz = C
        ocx, ocy, ocz = ox-cx, oy-cy, oz-cz
        a = dx*dx + dy*dy + dz*dz
        b = 2*(dx*ocx + dy*ocy + dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - r*r
        disc = b*b - 4*a*c
        if disc < 0: return None
        sq = math.sqrt(disc)
        t1 = (-b - sq)/(2*a); t2 = (-b + sq)/(2*a)
        if t1 >= 0: return t1
        if t2 >= 0: return t2
        return None
    t1 = sphere_intersect((0.0,0.0,0.0), 1.0)
    t2 = sphere_intersect((1.0,0.0,0.0), 1.0)
    ts = [t for t in (t1, t2) if t is not None]
    return min(ts) if ts else -1.0

#---------- 4.6 ----------#
def intersect_two_spheres(O, D):
    import math
    def intersect_sphere(C, r):
        ocx, ocy, ocz = O[0]-C[0], O[1]-C[1], O[2]-C[2]
        a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
        b = 2*(D[0]*ocx + D[1]*ocy + D[2]*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - r*r
        disc = b*b - 4*a*c
        if disc < 0:
            return None
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        tmin, tmax = min(t1, t2), max(t1, t2)
        if tmin >= 0:
            return tmin
        if tmax >= 0:
            return tmax
        return None
    t0 = intersect_sphere((0.0, 0.0, 0.0), 1.0)
    t1 = intersect_sphere((1.0, 0.0, 0.0), 1.0)
    ts = [t for t in (t0, t1) if t is not None]
    return min(ts) if ts else -1.0

#---------- 4.7 ----------#
def intersect_two_spheres(O, D):
    import math
    eps = 1e-6
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    ts = []
    for cx, cy, cz in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]:
        ocx, ocy, ocz = ox - cx, oy - cy, oz - cz
        b = 2 * (dx*ocx + dy*ocy + dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - 1.0
        disc = b*b - 4*a*c
        if disc < 0:
            continue
        s = math.sqrt(disc)
        t0 = (-b - s) / (2*a)
        t1 = (-b + s) / (2*a)
        if t0 >= 0: ts.append(t0)
        if t1 >= 0: ts.append(t1)
    if not ts:
        return -1.0
    ts.sort()
    for t in ts:
        px, py, pz = ox + t*dx, oy + t*dy, oz + t*dz
        d1 = math.sqrt(px*px + py*py + pz*pz)
        d2 = math.sqrt((px-1.0)**2 + py*py + pz*pz)
        if abs(d1 - 1.0) < eps and d2 > 1.0 + eps:
            return t
        if abs(d2 - 1.0) < eps and d1 > 1.0 + eps:
            return t
    return -1.0


#---------- 4.8 ----------#
def intersect(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    best_t = float('inf')
    for cx, cy, cz in ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)):
        ocx = ox - cx
        ocy = oy - cy
        ocz = oz - cz
        a = dx*dx + dy*dy + dz*dz
        b = 2.0 * (dx*ocx + dy*ocy + dz*ocz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - 1.0
        disc = b*b - 4.0*a*c
        if disc >= 0.0:
            sq = disc**0.5
            t1 = (-b - sq) / (2.0*a)
            t2 = (-b + sq) / (2.0*a)
            if 0.0 < t1 < best_t: best_t = t1
            if 0.0 < t2 < best_t: best_t = t2
    return best_t if best_t != float('inf') else -1.0

#---------- 4.9 ----------#
def intersect_two_spheres(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    best_t = float('inf')
    for Cx, Cy, Cz in ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)):
        r = 1.0
        ocx = Ox - Cx; ocy = Oy - Cy; ocz = Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2.0 * (ocx*Dx + ocy*Dy + ocz*Dz)
        c = ocx*ocx + ocy*ocy + ocz*ocz - r*r
        disc = b*b - 4.0*a*c
        if disc < 0.0: continue
        sqrt_disc = math.sqrt(disc)
        t1 = (-b - sqrt_disc) / (2.0*a)
        t2 = (-b + sqrt_disc) / (2.0*a)
        if t1 >= 0.0 and t1 < best_t: best_t = t1
        if t2 >= 0.0 and t2 < best_t: best_t = t2
    return best_t if best_t != float('inf') else -1.0

#---------- 5.0 ----------#
def intersect_ray_disc(O, D):
    n = (1.0, 1.0, 1.0)
    denom = n[0]*D[0] + n[1]*D[1] + n[2]*D[2]
    if abs(denom) < 1e-6:
        return -1
    numer = n[0]*O[0] + n[1]*O[1] + n[2]*O[2]
    t = -numer / denom
    if t < 0:
        return -1
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] <= 1.0:
        return t
    return -1


#---------- 5.1 ----------#
def intersect_disc(O, D):
    # Plane normal (normalized)
    inv_len = 1.0 / (3**0.5)
    n = (inv_len, inv_len, inv_len)
    # Dot product helper
    dot = lambda a, b: a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    # Compute denominator
    denom = dot(D, n)
    if abs(denom) < 1e-8:
        return -1.0
    # Compute t for intersection with plane
    t = -dot(O, n) / denom
    if t < 0.0:
        return -1.0
    # Compute intersection point P = O + t*D
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    # Check if within disc radius (radius = 1)
    if dot(P, P) <= 1.0:
        return t
    return -1.0

#---------- 5.2 ----------#
def intersect(O, D):
    import math
    # Disc parameters
    cx, cy, cz = 0.0, 0.0, 0.0
    radius = 1.0
    # Normalize disc normal
    nx, ny, nz = 1.0, 1.0, 1.0
    inv_len = 1.0 / math.sqrt(nx*nx + ny*ny + nz*nz)
    nx, ny, nz = nx*inv_len, ny*inv_len, nz*inv_len
    # Ray-plane intersection
    denom = D[0]*nx + D[1]*ny + D[2]*nz
    if abs(denom) < 1e-6:
        return -1.0
    t = -((O[0]-cx)*nx + (O[1]-cy)*ny + (O[2]-cz)*nz) / denom
    if t < 0.0:
        return -1.0
    # Compute intersection point
    px = O[0] + t*D[0]
    py = O[1] + t*D[1]
    pz = O[2] + t*D[2]
    # Check if within disc radius
    dx, dy, dz = px-cx, py-cy, pz-cz
    if dx*dx + dy*dy + dz*dz <= radius*radius:
        return t
    return -1.0


#---------- 5.3 ----------#
def intersect_disc(O, D):
    nx, ny, nz = 1.0, 1.0, 1.0
    denom = nx*D[0] + ny*D[1] + nz*D[2]
    if abs(denom) < 1e-6:
        return -1
    t = -(nx*O[0] + ny*O[1] + nz*O[2]) / denom
    if t < 0:
        return -1
    Px = O[0] + t * D[0]
    Py = O[1] + t * D[1]
    Pz = O[2] + t * D[2]
    if Px*Px + Py*Py + Pz*Pz <= 1.0:
        return t
    return -1


#---------- 5.4 ----------#
def intersect_disc(O, D):
    eps = 1e-6
    Nx, Ny, Nz = 1.0, 1.0, 1.0
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    denom = Dx*Nx + Dy*Ny + Dz*Nz
    if abs(denom) < eps:
        return -1.0
    t = -(Ox*Nx + Oy*Ny + Oz*Nz) / denom
    if t < 0.0:
        return -1.0
    Px = Ox + t*Dx
    Py = Oy + t*Dy
    Pz = Oz + t*Dz
    if Px*Px + Py*Py + Pz*Pz <= 1.0:
        return t
    return -1.0


#---------- 5.5 ----------#
def intersect_disc(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    nx, ny, nz = 1.0, 1.0, 1.0
    denom = dx*nx + dy*ny + dz*nz
    if abs(denom) < 1e-6:
        return -1.0
    t = (-(ox*nx + oy*ny + oz*nz)) / denom
    if t < 0.0:
        return -1.0
    px = ox + t*dx
    py = oy + t*dy
    pz = oz + t*dz
    if px*px + py*py + pz*pz <= 1.0:
        return t
    return -1.0

#---------- 5.6 ----------#
def intersect_disc(O, D):
    nx, ny, nz = 1.0, 1.0, 1.0
    denom = nx*D[0] + ny*D[1] + nz*D[2]
    if abs(denom) < 1e-6:
        return -1.0
    t = - (nx*O[0] + ny*O[1] + nz*O[2]) / denom
    if t < 0.0:
        return -1.0
    px = O[0] + t*D[0]
    py = O[1] + t*D[1]
    pz = O[2] + t*D[2]
    if px*px + py*py + pz*pz <= 1.0:
        return t
    return -1.0


#---------- 5.7 ----------#
def intersect_disc(O, D):
    N = (1.0, 1.0, 1.0)
    denom = D[0]*N[0] + D[1]*N[1] + D[2]*N[2]
    if abs(denom) < 1e-6:
        return -1.0
    t = -(O[0]*N[0] + O[1]*N[1] + O[2]*N[2]) / denom
    if t < 0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] <= 1.0:
        return t
    return -1.0


#---------- 5.8 ----------#
def intersect_disc(O, D):
    N = (1.0, 1.0, 1.0)
    denom = D[0]*N[0] + D[1]*N[1] + D[2]*N[2]
    if abs(denom) < 1e-6:
        return -1.0
    num = -(O[0]*N[0] + O[1]*N[1] + O[2]*N[2])
    t = num / denom
    if t < 0.0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] <= 1.0:
        return t
    return -1.0

#---------- 5.9 ----------#
def intersect_disc(O, D):
    Nx, Ny, Nz = 1.0, 1.0, 1.0
    denom = Nx*D[0] + Ny*D[1] + Nz*D[2]
    if abs(denom) < 1e-6:
        return -1.0
    t = -(Nx*O[0] + Ny*O[1] + Nz*O[2]) / denom
    if t < 0.0:
        return -1.0
    x, y, z = O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2]
    if x*x + y*y + z*z <= 1.0:
        return t
    return -1.0

#---------- 6.0 ----------#
def intersect_ray_quad_with_circle_hole(O, D):
    p0 = (1.0, 1.0, 1.0)
    p1 = (-1.0, 1.0, -1.0)
    p2 = (-2.0, -1.0, -2.0)
    p3 = (2.0, -1.0, 2.0)
    def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b): return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    n = cross(sub(p1, p0), sub(p3, p0))
    denom = dot(D, n)
    if abs(denom) < 1e-6:
        return -1.0
    t = dot(sub(p0, O), n) / denom
    if t < 0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    verts = [p0, p1, p2, p3]
    for i in range(4):
        vi = verts[i]
        vj = verts[(i+1) % 4]
        if dot(n, cross(sub(vj, vi), sub(P, vi))) < 0:
            return -1.0
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] <= 1.0:
        return -1.0
    return t


#---------- 6.1 ----------#
def intersect_quad_with_hole(O, D):
    def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    def inside_tri(P, A, B, C):
        v0, v1, v2 = sub(C, A), sub(B, A), sub(P, A)
        d00, d01, d11 = dot(v0, v0), dot(v0, v1), dot(v1, v1)
        d02, d12 = dot(v0, v2), dot(v1, v2)
        denom = d00*d11 - d01*d01
        if denom == 0: return False
        inv = 1.0/denom
        u = (d11*d02 - d01*d12)*inv
        v = (d00*d12 - d01*d02)*inv
        return u >= 0 and v >= 0 and u+v <= 1
    v0 = (1, 1, 1); v1 = (-1, 1, -1)
    v2 = (-2, -1, -2); v3 = (2, -1, 2)
    N = cross(sub(v1, v0), sub(v3, v0))
    denom = dot(N, D)
    if abs(denom) < 1e-6: return -1.0
    t = dot(N, sub(v0, O)) / denom
    if t < 0: return -1.0
    P = (O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2])
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] <= 1.0: return -1.0
    if inside_tri(P, v0, v1, v2) or inside_tri(P, v0, v2, v3):
        return t
    return -1.0


#---------- 6.2 ----------#
def intersect_quad_with_hole(O, D):
    v0 = (1.0, 1.0, 1.0)
    v1 = (-1.0, 1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = (2.0, -1.0, 2.0)
    eps = 1e-6
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

    # Plane normal
    N = cross(sub(v1, v0), sub(v3, v0))
    denom = dot(N, D)
    if abs(denom) < eps:
        return -1.0

    t = dot(N, sub(v0, O)) / denom
    if t < 0.0:
        return -1.0

    # Intersection point
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])

    # Inside quad test (edge-cross-halfspace)
    for A, B in ((v0, v1), (v1, v2), (v2, v3), (v3, v0)):
        if dot(cross(sub(B, A), sub(P, A)), N) < -eps:
            return -1.0

    # Circular hole test (radius = 1 at origin)
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] < 1.0 - eps:
        return -1.0

    return t


#---------- 6.3 ----------#
def intersect_quad_hole(O, D):
    epsilon = 1e-6
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    Ax, Ay, Az = 1.0, 1.0, 1.0
    Bx, By, Bz = -1.0, 1.0, -1.0
    Cx, Cy, Cz = -2.0, -1.0, -2.0
    Dx2, Dy2, Dz2 = 2.0, -1.0, 2.0
    # plane normal
    ABx, ABy, ABz = Bx-Ax, By-Ay, Bz-Az
    ADx, ADy, ADz = Dx2-Ax, Dy2-Ay, Dz2-Az
    Nx = ABy*ADz - ABz*ADy
    Ny = ABz*ADx - ABx*ADz
    Nz = ABx*ADy - ABy*ADx
    denom = Nx*Dx + Ny*Dy + Nz*Dz
    if abs(denom) < epsilon:
        return -1.0
    t = (Nx*(Ax-Ox) + Ny*(Ay-Oy) + Nz*(Az-Oz)) / denom
    if t <= epsilon:
        return -1.0
    Px = Ox + t*Dx; Py = Oy + t*Dy; Pz = Oz + t*Dz
    # inside quad
    verts = [(Ax,Ay,Az),(Bx,By,Bz),(Cx,Cy,Cz),(Dx2,Dy2,Dz2)]
    for i in range(4):
        x1,y1,z1 = verts[i]
        x2,y2,z2 = verts[(i+1)%4]
        ex, ey, ez = x2-x1, y2-y1, z2-z1
        vx, vy, vz = Px-x1, Py-y1, Pz-z1
        cx, cy, cz = ey*vz - ez*vy, ez*vx - ex*vz, ex*vy - ey*vx
        if cx*Nx + cy*Ny + cz*Nz < -epsilon:
            return -1.0
    # hole check
    if Px*Px + Py*Py + Pz*Pz < 1.0 - epsilon:
        return -1.0
    return t


#---------- 6.4 ----------#
def intersect_quad_with_hole(O, D):
    epsilon = 1e-6
    v0 = (1.0, 1.0, 1.0)
    v1 = (-1.0, 1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = (2.0, -1.0, 2.0)
    def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b): return (a[1]*b[2] - a[2]*b[1],
                              a[2]*b[0] - a[0]*b[2],
                              a[0]*b[1] - a[1]*b[0])
    # plane normal
    E1 = sub(v1, v0)
    E3 = sub(v3, v0)
    N = cross(E1, E3)
    denom = dot(D, N)
    if abs(denom) < epsilon:
        return -1.0
    t = dot(sub(v0, O), N) / denom
    if t < 0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    # inside quad test
    verts = [v0, v1, v2, v3]
    for i in range(4):
        a = verts[i]
        b = verts[(i+1) % 4]
        if dot(N, cross(sub(b, a), sub(P, a))) < -epsilon:
            return -1.0
    # circular hole
    if dot(P, P) >= 1.0:
        return t
    return -1.0

#---------- 6.5 ----------#
def intersect_ray_quad_with_circle_hole(O, D):
    # Vector operations
    def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b): return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def scale(v, s): return (v[0]*s, v[1]*s, v[2]*s)

    # Quad vertices
    v0 = (1.0, 1.0, 1.0)
    v1 = (-1.0, 1.0, -1.0)
    v2 = (-2.0, -1.0, -2.0)
    v3 = (2.0, -1.0, 2.0)

    # Compute plane normal
    N = cross(sub(v1, v0), sub(v3, v0))
    denom = dot(D, N)
    if abs(denom) < 1e-6:
        return -1.0
    t = dot(sub(v0, O), N) / denom
    if t < 0.0:
        return -1.0

    # Intersection point
    P = add(O, scale(D, t))

    # Inside quad test (convex polygon)
    verts = (v0, v1, v2, v3)
    for i in range(4):
        vi = verts[i]
        vj = verts[(i+1) % 4]
        if dot(cross(sub(vj, vi), sub(P, vi)), N) < 0.0:
            return -1.0

    # Circular hole test (radius = 1)
    if dot(P, P) < 1.0:
        return -1.0

    return t


#---------- 6.6 ----------#
def intersect_quad_with_hole(O, D):
    eps = 1e-6
    def add(a, b):
        return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    def sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def cross(a, b):
        return (a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0])
    def mul(a, s):
        return (a[0]*s, a[1]*s, a[2]*s)
    def length(a):
        return (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])**0.5

    p0 = (1.0, 1.0, 1.0)
    p1 = (-1.0, 1.0, -1.0)
    p2 = (-2.0, -1.0, -2.0)
    p3 = (2.0, -1.0, 2.0)
    # plane normal
    n = cross(sub(p1, p0), sub(p3, p0))
    ln = length(n)
    if ln < eps: return -1.0
    n = mul(n, 1.0/ln)
    denom = dot(n, D)
    if abs(denom) < eps: return -1.0
    t = dot(n, sub(p0, O)) / denom
    if t < eps: return -1.0
    P = add(O, mul(D, t))
    # check inside convex quad
    verts = [p0, p1, p2, p3]
    for i in range(4):
        vi, vj = verts[i], verts[(i+1)%4]
        edge = sub(vj, vi)
        vp = sub(P, vi)
        if dot(cross(edge, vp), n) < -eps:
            return -1.0
    # check circular hole of radius 1 at plane center
    vproj = sub(P, mul(n, dot(P, n)))
    if length(vproj) < 1.0 - eps:
        return -1.0
    return t

#---------- 6.7 ----------#
def intersect_quad_with_circle_hole(O, D):
    def dot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def cross(u, v): return (u[1]*v[2] - u[2]*v[1],
                              u[2]*v[0] - u[0]*v[2],
                              u[0]*v[1] - u[1]*v[0])
    def sub(u, v): return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    A = (1.0, 1.0,  1.0)
    B = (-1.0, 1.0, -1.0)
    C = (-2.0,-1.0, -2.0)
    Dv = (2.0,-1.0,  2.0)
    # Plane normal
    n = cross(sub(B, A), sub(Dv, A))
    denom = dot(D, n)
    if abs(denom) < 1e-6:
        return -1.0
    t = dot(sub(A, O), n) / denom
    if t < 0.0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    # Hole test (radius 1 at origin)
    if P[0]*P[0] + P[1]*P[1] + P[2]*P[2] < 1.0:
        return -1.0
    # Quad containment test (CCW order)
    verts = [A, B, C, Dv]
    for i in range(4):
        Vi, Vj = verts[i], verts[(i+1)%4]
        if dot(n, cross(sub(Vj, Vi), sub(P, Vi))) < 0.0:
            return -1.0
    return t


#---------- 6.8 ----------#
def intersect_ray_quad_with_hole(O, D):
    Ax,Ay,Az = 1.0,1.0,1.0
    Bx,By,Bz = -1.0,1.0,-1.0
    Cx,Cy,Cz = -2.0,-1.0,-2.0
    Dx,Dy,Dz = 2.0,-1.0,2.0
    Ox,Oy,Oz = O
    DxR,DyR,DzR = D
    # compute plane normal
    ux,uy,uz = Bx-Ax, By-Ay, Bz-Az
    vx,vy,vz = Cx-Ax, Cy-Ay, Cz-Az
    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx
    denom = nx*DxR + ny*DyR + nz*DzR
    if abs(denom) < 1e-6:
        return -1.0
    t = (nx*(Ax-Ox) + ny*(Ay-Oy) + nz*(Az-Oz)) / denom
    if t < 0:
        return -1.0
    Px,Py,Pz = Ox + t*DxR, Oy + t*DyR, Oz + t*DzR
    # circular hole test
    if Px*Px + Py*Py + Pz*Pz < 1.0:
        return -1.0
    def point_in_triangle(Px,Py,Pz, x0,y0,z0, x1,y1,z1, x2,y2,z2):
        u0,u1,u2 = x1-x0, y1-y0, z1-z0
        v0,v1,v2 = x2-x0, y2-y0, z2-z0
        w0,w1,w2 = Px-x0, Py-y0, Pz-z0
        dot00 = u0*u0 + u1*u1 + u2*u2
        dot01 = u0*v0 + u1*v1 + u2*v2
        dot02 = u0*w0 + u1*w1 + u2*w2
        dot11 = v0*v0 + v1*v1 + v2*v2
        dot12 = v0*w0 + v1*w1 + v2*w2
        denom2 = dot00*dot11 - dot01*dot01
        if abs(denom2) < 1e-6:
            return False
        alpha = (dot11*dot02 - dot01*dot12) / denom2
        beta  = (dot00*dot12 - dot01*dot02) / denom2
        return alpha >= 0 and beta >= 0 and alpha + beta <= 1
    if point_in_triangle(Px,Py,Pz, Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz) or \
       point_in_triangle(Px,Py,Pz, Ax,Ay,Az, Cx,Cy,Cz, Dx,Dy,Dz):
        return t
    return -1.0


#---------- 6.9 ----------#
def intersect_quad_hole(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    p0 = (1.0, 1.0, 1.0)
    p1 = (-1.0, 1.0, -1.0)
    p2 = (-2.0, -1.0, -2.0)
    p3 = (2.0, -1.0, 2.0)
    # Plane normal n = (-1, 0, 1)
    denom = -Dx + Dz
    if abs(denom) < 1e-6:
        return -1.0
    # n �� (O - p0) = Oz - Ox
    na = (Oz - Ox)
    t = -na / denom
    if t < 0.0:
        return -1.0
    # Intersection point
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz
    # Inside quad test
    n0, n1, n2 = -1.0, 0.0, 1.0
    for Vi, Vj in ((p0, p1), (p1, p2), (p2, p3), (p3, p0)):
        Ex, Ey, Ez = Vj[0] - Vi[0], Vj[1] - Vi[1], Vj[2] - Vi[2]
        VPx, VPy, VPz = Px - Vi[0], Py - Vi[1], Pz - Vi[2]
        cx = Ey * VPz - Ez * VPy
        cy = Ez * VPx - Ex * VPz
        cz = Ex * VPy - Ey * VPx
        if n0 * cx + n1 * cy + n2 * cz < 0.0:
            return -1.0
    # Circular hole test (radius 1 at origin)
    if Px * Px + Py * Py + Pz * Pz < 1.0:
        return -1.0
    return t

#---------- 7.0 ----------#
def intersect_ray_plane_with_triangle_hole(O, D):
    # Plane y = z �� normal n and plane constant d
    n = (0.0, 1.0, -1.0)
    d = 0.0
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    denom = n[1]*Dy + n[2]*Dz  # n��D = 1*Dy + (-1)*Dz
    if abs(denom) < 1e-8:
        return -1.0
    t = -(n[0]*Ox + n[1]*Oy + n[2]*Oz + d) / denom
    if t < 0.0:
        return -1.0
    Px, Py, Pz = Ox + t*Dx, Oy + t*Dy, Oz + t*Dz
    # Triangle hole vertices
    A = (1.0, 1.0, 1.0)
    B = (-1.0, 0.0, 0.0)
    C = (0.0, -1.0, -1.0)
    # check if P inside triangle ABC (hole) using edge tests
    def edge_test(U, V, P):
        # (V-U) x (P-U) �� n
        Ex, Ey, Ez = V[0]-U[0], V[1]-U[1], V[2]-U[2]
        Vx, Vy, Vz = P[0]-U[0], P[1]-U[1], P[2]-U[2]
        cx = Ey*Vz - Ez*Vy
        cy = Ez*Vx - Ex*Vz
        cz = Ex*Vy - Ey*Vx
        return cx*n[0] + cy*n[1] + cz*n[2]
    if (edge_test(A, B, (Px,Py,Pz)) <= 1e-8 and
        edge_test(B, C, (Px,Py,Pz)) <= 1e-8 and
        edge_test(C, A, (Px,Py,Pz)) <= 1e-8):
        return -1.0
    return t

#---------- 7.1 ----------#
def intersect_plane_with_tri_hole(O, D):
    eps = 1e-8
    # Plane: y - z = 0, normal N = (0,1,-1)
    N = (0.0, 1.0, -1.0)
    denom = N[0]*D[0] + N[1]*D[1] + N[2]*D[2]
    if abs(denom) < eps:
        return -1.0
    t = - (N[0]*O[0] + N[1]*O[1] + N[2]*O[2]) / denom
    if t < 0.0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    # Triangle hole vertices A, B, C in CCW order
    A = (1.0, 1.0, 1.0)
    B = (-1.0, 0.0, 0.0)
    C = (0.0, -1.0, -1.0)
    # Barycentric coordinate test
    v0 = (B[0]-A[0], B[1]-A[1], B[2]-A[2])
    v1 = (C[0]-A[0], C[1]-A[1], C[2]-A[2])
    v2 = (P[0]-A[0], P[1]-A[1], P[2]-A[2])
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    denom2 = dot00*dot11 - dot01*dot01
    if abs(denom2) < eps:
        return t
    u = (dot11*dot02 - dot01*dot12) / denom2
    v = (dot00*dot12 - dot01*dot02) / denom2
    # Strictly inside hole => no intersection
    if u > eps and v > eps and u + v < 1.0 - eps:
        return -1.0
    return t


#---------- 7.2 ----------#
def intersect_plane_with_triangle_hole(O, D):
    eps = 1e-6
    denom = D[1] - D[2]
    if abs(denom) < eps:
        return -1.0
    t = (O[2] - O[1]) / denom
    if t < 0:
        return -1.0
    # Intersection point
    Px, Py, Pz = O[0] + D[0]*t, O[1] + D[1]*t, O[2] + D[2]*t
    # Triangle vertices
    Ax, Ay, Az = 1.0, 1.0, 1.0
    Bx, By, Bz = -1.0, 0.0, 0.0
    Cx, Cy, Cz = 0.0, -1.0, -1.0
    # Vectors for barycentric
    v0 = (Cx-Ax, Cy-Ay, Cz-Az)
    v1 = (Bx-Ax, By-Ay, Bz-Az)
    v2 = (Px-Ax, Py-Ay, Pz-Az)
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    invDenom = 1.0 / (dot00*dot11 - dot01*dot01)
    u = (dot11*dot02 - dot01*dot12) * invDenom
    v = (dot00*dot12 - dot01*dot02) * invDenom
    # If inside triangle (including boundary), it's a hole -> no intersection
    if u >= 0.0 and v >= 0.0 and u + v <= 1.0:
        return -1.0
    return t


#---------- 7.3 ----------#
def intersect_plane_with_triangle_hole(O, D):
    x0, y0, z0 = O
    dx, dy, dz = D
    denom = dy - dz
    if abs(denom) < 1e-8:
        return -1.0
    t = (z0 - y0) / denom
    if t < 0.0:
        return -1.0
    px = x0 + t * dx
    py = y0 + t * dy
    pz = z0 + t * dz
    ax, ay, az = 1.0, 1.0, 1.0
    bx, by, bz = -1.0, 0.0, 0.0
    cx, cy, cz = 0.0, -1.0, -1.0
    v0x, v0y, v0z = cx - ax, cy - ay, cz - az
    v1x, v1y, v1z = bx - ax, by - ay, bz - az
    v2x, v2y, v2z = px - ax, py - ay, pz - az
    dot00 = v0x*v0x + v0y*v0y + v0z*v0z
    dot01 = v0x*v1x + v0y*v1y + v0z*v1z
    dot02 = v0x*v2x + v0y*v2y + v0z*v2z
    dot11 = v1x*v1x + v1y*v1y + v1z*v1z
    dot12 = v1x*v2x + v1y*v2y + v1z*v2z
    invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
    if u > 0.0 and v > 0.0 and u + v < 1.0:
        return -1.0
    return t

#---------- 7.4 ----------#
def intersect_ray_plane_with_triangle_hole(O, D):
    # plane: y = z  ->  n��X = 0 with n = (0,1,-1)
    n = (0.0, 1.0, -1.0)
    denom = n[1]*D[1] + n[2]*D[2]
    if abs(denom) < 1e-8:
        return -1.0
    t = -(n[1]*O[1] + n[2]*O[2]) / denom
    if t < 0.0:
        return -1.0
    P = (O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2])
    # triangle hole vertices
    A = (1.0, 1.0, 1.0)
    B = (-1.0, 0.0, 0.0)
    C = (0.0, -1.0, -1.0)
    # barycentric test
    v0 = (B[0]-A[0], B[1]-A[1], B[2]-A[2])
    v1 = (C[0]-A[0], C[1]-A[1], C[2]-A[2])
    v2 = (P[0]-A[0], P[1]-A[1], P[2]-A[2])
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    d00 = dot(v0, v0)
    d01 = dot(v0, v1)
    d11 = dot(v1, v1)
    d20 = dot(v2, v0)
    d21 = dot(v2, v1)
    denom2 = d00*d11 - d01*d01
    if abs(denom2) < 1e-8:
        return t
    v = (d11*d20 - d01*d21) / denom2
    w = (d00*d21 - d01*d20) / denom2
    u = 1.0 - v - w
    # inside hole if strictly inside
    if u > 0.0 and v > 0.0 and w > 0.0:
        return -1.0
    return t

#---------- 7.5 ----------#
def intersect(O, D):
    # Plane y = z: normal n = (0,1,-1)
    denom = D[1] - D[2]
    if abs(denom) < 1e-8:
        return -1.0
    t = (O[2] - O[1]) / denom
    if t < 0:
        return -1.0
    # Intersection point
    Px, Py, Pz = O[0] + t*D[0], O[1] + t*D[1], O[2] + t*D[2]
    # Triangle hole vertices
    Ax, Ay, Az = 1.0, 1.0, 1.0
    Bx, By, Bz = -1.0, 0.0, 0.0
    Cx, Cy, Cz = 0.0, -1.0, -1.0
    # Vectors for barycentric
    v0x, v0y, v0z = Bx - Ax, By - Ay, Bz - Az
    v1x, v1y, v1z = Cx - Ax, Cy - Ay, Cz - Az
    v2x, v2y, v2z = Px - Ax, Py - Ay, Pz - Az
    # Dot products
    dot00 = v0x*v0x + v0y*v0y + v0z*v0z
    dot01 = v0x*v1x + v0y*v1y + v0z*v1z
    dot02 = v0x*v2x + v0y*v2y + v0z*v2z
    dot11 = v1x*v1x + v1y*v1y + v1z*v1z
    dot12 = v1x*v2x + v1y*v2y + v1z*v2z
    denomB = dot00 * dot11 - dot01 * dot01
    if abs(denomB) < 1e-8:
        return t
    invDen = 1.0 / denomB
    u = (dot11 * dot02 - dot01 * dot12) * invDen
    v = (dot00 * dot12 - dot01 * dot02) * invDen
    # Strict inside triangle => hole
    if u > 0 and v > 0 and (u + v) < 1:
        return -1.0
    return t

#---------- 7.6 ----------#
def intersect_plane_with_triangle_hole(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    denom = dy - dz
    if abs(denom) < 1e-6:
        return -1.0
    t = (oz - oy) / denom
    if t < 0.0:
        return -1.0
    px = ox + t * dx
    py = oy + t * dy
    pz = oz + t * dz
    p1 = (1.0, 1.0, 1.0)
    p2 = (-1.0, 0.0, 0.0)
    p3 = (0.0, -1.0, -1.0)
    v0 = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])
    v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v2 = (px   - p1[0], py   - p1[1], pz   - p1[2])
    dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]
    dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]
    dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2]
    dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
    dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    if u >= 0.0 and v >= 0.0 and u + v <= 1.0:
        return -1.0
    return t


#---------- 7.7 ----------#
def intersect_plane_with_triangle_hole(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    denom = Dy - Dz
    if abs(denom) < 1e-8:
        return -1.0
    t = (Oz - Oy) / denom
    if t < 0.0:
        return -1.0
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    # Project to 2D (x,y), since on plane y=z we have Py == Pz
    P2x, P2y = Px, Py
    tri = [(1.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    for i in range(3):
        ax, ay = tri[i]
        bx, by = tri[(i + 1) % 3]
        cross = (ax - bx) * (P2y - by) - (ay - by) * (P2x - bx)
        if cross < 0.0:
            return t
    return -1.0


#---------- 7.8 ----------#
def intersect_plane_triangle_hole(O, D):
    tol = 1e-6
    ox, oy, oz = O
    dx, dy, dz = D
    den = dy - dz
    if abs(den) < tol:
        return -1.0
    t = (oz - oy) / den
    if t < 0:
        return -1.0
    px = ox + t * dx
    py = oy + t * dy
    pz = oz + t * dz
    Ax, Ay, Az = 1.0, 1.0, 1.0
    Bx, By, Bz = -1.0, 0.0, 0.0
    Cx, Cy, Cz = 0.0, -1.0, -1.0
    nx, ny, nz = 0.0, -3.0, 3.0
    # edge AB
    ex0, ey0, ez0 = Bx - Ax, By - Ay, Bz - Az
    vx0, vy0, vz0 = px - Ax, py - Ay, pz - Az
    cx0 = ey0 * vz0 - ez0 * vy0
    cy0 = ez0 * vx0 - ex0 * vz0
    cz0 = ex0 * vy0 - ey0 * vx0
    dot0 = cx0 * nx + cy0 * ny + cz0 * nz
    # edge BC
    ex1, ey1, ez1 = Cx - Bx, Cy - By, Cz - Bz
    vx1, vy1, vz1 = px - Bx, py - By, pz - Bz
    cx1 = ey1 * vz1 - ez1 * vy1
    cy1 = ez1 * vx1 - ex1 * vz1
    cz1 = ex1 * vy1 - ey1 * vx1
    dot1 = cx1 * nx + cy1 * ny + cz1 * nz
    # edge CA
    ex2, ey2, ez2 = Ax - Cx, Ay - Cy, Az - Cz
    vx2, vy2, vz2 = px - Cx, py - Cy, pz - Cz
    cx2 = ey2 * vz2 - ez2 * vy2
    cy2 = ez2 * vx2 - ex2 * vz2
    cz2 = ex2 * vy2 - ey2 * vx2
    dot2 = cx2 * nx + cy2 * ny + cz2 * nz
    if dot0 > 0 and dot1 > 0 and dot2 > 0:
        return -1.0
    return t

#---------- 7.9 ----------#
def intersect_plane_triangle_hole(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    eps = 1e-8
    # intersect with plane y = z
    denom = Dy - Dz
    if abs(denom) < eps:
        return -1.0
    t = -(Oy - Oz) / denom
    if t < 0.0:
        return -1.0
    # intersection point
    Px = Ox + t * Dx
    Py = Oy + t * Dy
    Pz = Oz + t * Dz
    # triangle vertices
    Ax, Ay, Az = 1.0, 1.0, 1.0
    Bx, By, Bz = -1.0, 0.0, 0.0
    Cx, Cy, Cz = 0.0, -1.0, -1.0
    # barycentric coordinates
    v0x, v0y, v0z = Bx - Ax, By - Ay, Bz - Az
    v1x, v1y, v1z = Cx - Ax, Cy - Ay, Cz - Az
    v2x, v2y, v2z = Px - Ax, Py - Ay, Pz - Az
    dot00 = v0x*v0x + v0y*v0y + v0z*v0z
    dot01 = v0x*v1x + v0y*v1y + v0z*v1z
    dot11 = v1x*v1x + v1y*v1y + v1z*v1z
    dot02 = v0x*v2x + v0y*v2y + v0z*v2z
    dot12 = v1x*v2x + v1y*v2y + v1z*v2z
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    # if inside triangle hole (open region), no intersection
    if u > eps and v > eps and u + v < 1.0 - eps:
        return -1.0
    return t


#---------- 8.0 ----------#
def intersect_ray_ellipsoid(O, D):
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    a, b, c = 2.0, 3.0, 4.0
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    A = Dx*Dx/(a*a) + Dy*Dy/(b*b) + Dz*Dz/(c*c)
    B = 2*((Ox-Cx)*Dx/(a*a) + (Oy-Cy)*Dy/(b*b) + (Oz-Cz)*Dz/(c*c))
    C = (Ox-Cx)*(Ox-Cx)/(a*a) + (Oy-Cy)*(Oy-Cy)/(b*b) + (Oz-Cz)*(Oz-Cz)/(c*c) - 1.0
    disc = B*B - 4*A*C
    if disc < 0.0:
        return -1.0
    sqrt_disc = disc**0.5
    t1 = (-B - sqrt_disc) / (2*A)
    t2 = (-B + sqrt_disc) / (2*A)
    t = float('inf')
    if t1 > 0 and t1 < t: t = t1
    if t2 > 0 and t2 < t: t = t2
    return t if t < float('inf') else -1.0

#---------- 8.1 ----------#
def intersect_ellipsoid(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    cx, cy, cz = 1.0, 1.0, 1.0
    a, b, c = 2.0, 3.0, 4.0
    ox = (ox - cx) / a
    oy = (oy - cy) / b
    oz = (oz - cz) / c
    dx /= a
    dy /= b
    dz /= c
    A = dx*dx + dy*dy + dz*dz
    B = 2.0*(ox*dx + oy*dy + oz*dz)
    C = ox*ox + oy*oy + oz*oz - 1.0
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0
    sqrt_disc = disc**0.5
    t0 = (-B - sqrt_disc) / (2.0*A)
    t1 = (-B + sqrt_disc) / (2.0*A)
    if t0 > 0.0 and t1 > 0.0:
        return min(t0, t1)
    elif t0 > 0.0:
        return t0
    elif t1 > 0.0:
        return t1
    else:
        return -1.0

#---------- 8.2 ----------#
def intersect_ellipsoid(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    ox -= 1.0; oy -= 1.0; oz -= 1.0
    A = (dx/2.0)**2 + (dy/3.0)**2 + (dz/4.0)**2
    B = 2.0 * (dx*ox/4.0 + dy*oy/9.0 + dz*oz/16.0)
    C = (ox/2.0)**2 + (oy/3.0)**2 + (oz/4.0)**2 - 1.0
    disc = B*B - 4.0*A*C
    if disc < 0.0 or A == 0.0:
        return -1.0
    sqrt_disc = disc**0.5
    t1 = (-B - sqrt_disc) / (2.0*A)
    t2 = (-B + sqrt_disc) / (2.0*A)
    t0 = t1 if t1 < t2 else t2
    if t0 >= 0.0:
        return t0
    t1, t2 = t2, t1  # ensure t1 is the larger
    return t1 if t1 >= 0.0 else -1.0

#---------- 8.3 ----------#
def intersect_ray_ellipsoid(O, D):
    import math
    Cx, Cy, Cz = 1.0, 1.0, 1.0
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    Px = Ox - Cx
    Py = Oy - Cy
    Pz = Oz - Cz
    a = Dx*Dx/4.0 + Dy*Dy/9.0 + Dz*Dz/16.0
    b = 2.0*(Px*Dx/4.0 + Py*Dy/9.0 + Pz*Dz/16.0)
    c = Px*Px/4.0 + Py*Py/9.0 + Pz*Pz/16.0 - 1.0
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t0 = (-b - sqrt_disc) / (2.0*a)
    t1 = (-b + sqrt_disc) / (2.0*a)
    t_near = min(t0, t1)
    t_far = max(t0, t1)
    if t_near >= 0.0:
        return t_near
    if t_far >= 0.0:
        return t_far
    return -1.0


#---------- 8.4 ----------#
def intersect_ellipsoid(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 4.0, 9.0, 16.0
    ex, ey, ez = Ox - cx, Oy - cy, Oz - cz
    A = Dx*Dx/a2 + Dy*Dy/b2 + Dz*Dz/c2
    B = 2.0*(ex*Dx/a2 + ey*Dy/b2 + ez*Dz/c2)
    C = ex*ex/a2 + ey*ey/b2 + ez*ez/c2 - 1.0
    disc = B*B - 4.0*A*C
    if disc < 0.0:
        return -1.0
    sqrt_d = math.sqrt(disc)
    t0 = (-B - sqrt_d) / (2.0*A)
    t1 = (-B + sqrt_d) / (2.0*A)
    tmin, tmax = min(t0, t1), max(t0, t1)
    if tmin >= 0.0:
        return tmin
    if tmax >= 0.0:
        return tmax
    return -1.0


#---------- 8.5 ----------#
def intersect_ray_ellipsoid(O, D):
    from math import sqrt
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 4.0, 9.0, 16.0
    Lx, Ly, Lz = O[0] - cx, O[1] - cy, O[2] - cz
    Dx, Dy, Dz = D
    A = Dx*Dx/a2 + Dy*Dy/b2 + Dz*Dz/c2
    B = Lx*Dx/a2 + Ly*Dy/b2 + Lz*Dz/c2
    C = Lx*Lx/a2 + Ly*Ly/b2 + Lz*Lz/c2 - 1.0
    disc = B*B - A*C
    if disc < 0 or A == 0:
        return -1.0
    sqrt_disc = sqrt(disc)
    t0 = (-B - sqrt_disc)/A
    t1 = (-B + sqrt_disc)/A
    if t0 > 0:
        return t0
    if t1 > 0:
        return t1
    return -1.0

#---------- 8.6 ----------#
def intersect_ellipsoid(O, D):
    import math
    # Ellipsoid parameters
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 2.0**2, 3.0**2, 4.0**2
    # Ray origin and direction
    ox, oy, oz = O
    dx, dy, dz = D
    # Coefficients for quadratic At^2 + Bt + C = 0
    xo, yo, zo = ox - cx, oy - cy, oz - cz
    A = dx*dx/a2 + dy*dy/b2 + dz*dz/c2
    B = 2*(dx*xo/a2 + dy*yo/b2 + dz*zo/c2)
    C = xo*xo/a2 + yo*yo/b2 + zo*zo/c2 - 1.0
    # Solve
    if abs(A) < 1e-8:
        if abs(B) < 1e-8:
            return -1.0
        t = -C / B
        return t if t >= 0.0 else -1.0
    disc = B*B - 4*A*C
    if disc < 0.0:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t0 = (-B - sqrt_disc) / (2*A)
    t1 = (-B + sqrt_disc) / (2*A)
    ts = [t for t in (t0, t1) if t >= 0.0]
    return min(ts) if ts else -1.0


#---------- 8.7 ----------#
def intersect_ray_ellipsoid(O, D):
    import math
    # Ellipsoid center and squared semi-axes
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 2.0**2, 3.0**2, 4.0**2
    # Translate ray origin
    ox, oy, oz = O[0] - cx, O[1] - cy, O[2] - cz
    dx, dy, dz = D
    # Quadratic coefficients
    A = dx*dx/a2 + dy*dy/b2 + dz*dz/c2
    B = 2*(ox*dx/a2 + oy*dy/b2 + oz*dz/c2)
    C = ox*ox/a2 + oy*oy/b2 + oz*oz/c2 - 1.0
    # Discriminant
    disc = B*B - 4*A*C
    if disc < 0 or abs(A) < 1e-8:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t0 = (-B - sqrt_disc) / (2*A)
    t1 = (-B + sqrt_disc) / (2*A)
    # Select smallest non-negative t
    if t0 >= 0 and t1 >= 0:
        return min(t0, t1)
    if t0 >= 0:
        return t0
    if t1 >= 0:
        return t1
    return -1.0


#---------- 8.8 ----------#
def intersect_ray_ellipsoid(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    # Ellipsoid center and squared semi-axes
    cx, cy, cz = 1.0, 1.0, 1.0
    a2, b2, c2 = 4.0, 9.0, 16.0
    # Translate ray origin to ellipsoid coordinate
    xo, yo, zo = ox - cx, oy - cy, oz - cz
    # Quadratic coefficients
    A = dx*dx/a2 + dy*dy/b2 + dz*dz/c2
    B = 2*(dx*xo/a2 + dy*yo/b2 + dz*zo/c2)
    C = xo*xo/a2 + yo*yo/b2 + zo*zo/c2 - 1.0
    # Discriminant
    disc = B*B - 4*A*C
    if A == 0 or disc < 0:
        return -1.0
    sqrt_disc = disc**0.5
    t0 = (-B - sqrt_disc) / (2*A)
    t1 = (-B + sqrt_disc) / (2*A)
    # Choose smallest positive t
    t = min(t0, t1) if min(t0, t1) > 0 else max(t0, t1)
    return t if t > 0 else -1.0


#---------- 8.9 ----------#
def intersect_ray_ellipsoid(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ocx = Ox - 1.0
    ocy = Oy - 1.0
    ocz = Oz - 1.0
    a = Dx*Dx/4.0 + Dy*Dy/9.0 + Dz*Dz/16.0
    b = 2.0*(ocx*Dx/4.0 + ocy*Dy/9.0 + ocz*Dz/16.0)
    c = ocx*ocx/4.0 + ocy*ocy/9.0 + ocz*ocz/16.0 - 1.0
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2.0*a)
    t2 = (-b + sqrt_disc) / (2.0*a)
    t = float('inf')
    if t1 >= 0.0:
        t = t1
    if t2 >= 0.0 and t2 < t:
        t = t2
    return t if t != float('inf') else -1.0

#---------- 9.0 ----------#
def intersect_cut_sphere(O, D):
    import math
    x0, y0, z0 = O
    dx, dy, dz = D
    # intersection with cut disc
    t_disc = None
    denom = dx + dy
    if denom != 0:
        t = (1 - (x0 + y0)) / denom
        if t > 0 and denom < 0:
            px, py, pz = x0 + t*dx, y0 + t*dy, z0 + t*dz
            if px*px + py*py + pz*pz <= 1:
                t_disc = t
    # intersection with sphere
    a = dx*dx + dy*dy + dz*dz
    b = 2*(x0*dx + y0*dy + z0*dz)
    c = x0*x0 + y0*y0 + z0*z0 - 1
    disc = b*b - 4*a*c
    t_sph = None
    if disc >= 0:
        sd = math.sqrt(disc)
        for t in sorted([(-b - sd)/(2*a), (-b + sd)/(2*a)]):
            if t > 0:
                px, py = x0 + t*dx, y0 + t*dy
                if px + py <= 1:
                    t_sph = t
                    break
    ts = [t for t in (t_disc, t_sph) if t is not None]
    return min(ts) if ts else -1.0

#---------- 9.1 ----------#
def intersect_ray_cut_sphere(O, D):
    import math
    ox, oy, oz = O
    dx, dy, dz = D
    # Sphere intersection (centered at origin, radius=1)
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    ts = []
    if disc >= 0:
        sd = math.sqrt(disc)
        for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
            if t > 0:
                px, py, pz = ox + t*dx, oy + t*dy, oz + t*dz
                if px + py <= 1:
                    ts.append(t)
    # Plane intersection (x + y = 1)
    nx, ny, nz = 1, 1, 0
    denom = nx*dx + ny*dy + nz*dz
    if abs(denom) > 1e-6:
        t_plane = (1 - (nx*ox + ny*oy + nz*oz)) / denom
        if t_plane > 0:
            px, py, pz = ox + t_plane*dx, oy + t_plane*dy, oz + t_plane*dz
            if px*px + py*py + pz*pz <= 1:
                ts.append(t_plane)
    return min(ts) if ts else -1.0

#---------- 9.2 ----------#
def intersect_cut_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    # Sphere intersection (center=0, radius=1)
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    t_hit = float('inf')
    if disc >= 0:
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        for t in (t0, t1):
            if t >= 0:
                x = ox + t*dx
                y = oy + t*dy
                z = oz + t*dz
                if x + y <= 1:
                    t_hit = t
                    break
    # Plane intersection (x+y=1)
    denom = dx + dy
    if denom != 0:
        t_plane = (1 - ox - oy) / denom
        if t_plane >= 0:
            x = ox + t_plane*dx
            y = oy + t_plane*dy
            z = oz + t_plane*dz
            if x*x + y*y + z*z <= 1:
                t_hit = min(t_hit, t_plane)
    return t_hit if t_hit != float('inf') else -1.0

#---------- 9.3 ----------#
def intersect_cut_sphere(O, D):
    import math
    eps = 1e-6
    ox, oy, oz = O
    dx, dy, dz = D

    # Sphere intersection
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    ts = []
    if disc >= 0:
        sd = math.sqrt(disc)
        t1 = (-b - sd) / (2*a)
        t2 = (-b + sd) / (2*a)
        for t in sorted((t1, t2)):
            if t >= 0:
                px = ox + t*dx
                py = oy + t*dy
                if px + py <= 1 + eps:
                    ts.append(t)

    # Plane intersection (x+y=1)
    nx, ny, nz = 1.0, 1.0, 0.0
    denom = nx*dx + ny*dy + nz*dz
    if abs(denom) > eps:
        t_plane = (1 - (nx*ox + ny*oy + nz*oz)) / denom
        if t_plane >= 0:
            px = ox + t_plane*dx
            py = oy + t_plane*dy
            pz = oz + t_plane*dz
            if (px-0.5)**2 + (py-0.5)**2 + pz*pz <= 0.5 + eps:
                ts.append(t_plane)

    return min(ts) if ts else -1


#---------- 9.4 ----------#
def intersect_cut_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    # Sphere intersection (centered at 0, radius = 1)
    a = dx*dx + dy*dy + dz*dz
    b = ox*dx + oy*dy + oz*dz
    c = ox*ox + oy*oy + oz*oz - 1.0
    disc = b*b - a*c
    t_sph = None
    if disc >= 0.0:
        sqrt_disc = disc**0.5
        t1 = (-b - sqrt_disc) / a
        t2 = (-b + sqrt_disc) / a
        ts = [t for t in (t1, t2) if t > 0.0]
        if ts:
            t_candidate = min(ts)
            px, py, pz = ox + t_candidate*dx, oy + t_candidate*dy, oz + t_candidate*dz
            if px + py <= 1.0:
                t_sph = t_candidate
    # Plane-cut intersection (plane x+y=1, disc radius 1)
    denom = dx + dy
    t_pl = None
    if denom != 0.0:
        t_plane = (1.0 - ox - oy) / denom
        if t_plane > 0.0 and ox + oy > 1.0 and denom < 0.0:
            px, py, pz = ox + t_plane*dx, oy + t_plane*dy, oz + t_plane*dz
            if px*px + py*py + pz*pz <= 1.0:
                t_pl = t_plane
    # Choose nearest positive intersection
    ts = [t for t in (t_sph, t_pl) if t is not None]
    return min(ts) if ts else -1.0


#---------- 9.5 ----------#
def intersect_cut_sphere(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    # Sphere intersection (x^2+y^2+z^2=1)
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    ts = []
    if disc >= 0:
        sd = disc**0.5
        for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
            if t >= 0:
                x, y, z = ox + t*dx, oy + t*dy, oz + t*dz
                if x + y <= 1:
                    ts.append(t)
    # Plane intersection (x+y=1) with disc
    denom = dx + dy
    if abs(denom) > 1e-6:
        t = (1 - (ox + oy)) / denom
        if t >= 0:
            x, y, z = ox + t*dx, oy + t*dy, oz + t*dz
            if x*x + y*y + z*z <= 1:
                ts.append(t)
    return min(ts) if ts else -1.0

#---------- 9.6 ----------#
def intersect_cut_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    C = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = B*B - 4*A*C
    ts = []
    eps = 1e-6
    if disc >= 0:
        sqrt_disc = disc**0.5
        for t in [(-B - sqrt_disc)/(2*A), (-B + sqrt_disc)/(2*A)]:
            if t > eps:
                Px = Ox + t*Dx
                Py = Oy + t*Dy
                if Px + Py <= 1:
                    ts.append(t)
    denom = Dx + Dy
    if abs(denom) > eps:
        t_plane = (1 - (Ox + Oy)) / denom
        if t_plane > eps:
            Px = Ox + t_plane*Dx
            Py = Oy + t_plane*Dy
            Pz = Oz + t_plane*Dz
            if Px*Px + Py*Py + Pz*Pz <= 1:
                ts.append(t_plane)
    return min(ts) if ts else -1

#---------- 9.7 ----------#
def intersect(O, D):
    import math
    ox, oy, oz = O
    dx, dy, dz = D
    # Sphere intersection
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    ts = []
    if disc >= 0:
        sd = math.sqrt(disc)
        for t in [(-b - sd) / (2*a), (-b + sd) / (2*a)]:
            if t > 0:
                px, py = ox + t*dx, oy + t*dy
                if px + py <= 1:
                    ts.append(t)
    # Plane intersection (x + y = 1)
    denom = dx + dy
    if abs(denom) > 1e-6:
        t_p = (1 - ox - oy) / denom
        if t_p > 0:
            px, py, pz = ox + t_p*dx, oy + t_p*dy, oz + t_p*dz
            if px*px + py*py + pz*pz <= 1:
                ts.append(t_p)
    return min(ts) if ts else -1.0


#---------- 9.8 ----------#
def intersect_cut_sphere(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    eps = 1e-6
    # Sphere intersection (centered at origin, radius=1)
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    ts = []
    if disc >= 0:
        sqrt_disc = disc**0.5
        t0 = (-b - sqrt_disc) / (2*a)
        t1 = (-b + sqrt_disc) / (2*a)
        for t in (t0, t1):
            if t > eps:
                x, y = Ox + t*Dx, Oy + t*Dy
                if x + y <= 1 + eps:
                    ts.append(t)
    # Plane intersection (x + y = 1) disk check
    denom = Dx + Dy
    if abs(denom) > eps:
        t_plane = (1 - (Ox + Oy)) / denom
        if t_plane > eps:
            x, y, z = Ox + t_plane*Dx, Oy + t_plane*Dy, Oz + t_plane*Dz
            if x*x + y*y + z*z <= 1 + eps:
                ts.append(t_plane)
    if not ts:
        return -1.0
    return min(ts)


#---------- 9.9 ----------#
def intersect_cut_sphere(O, D):
    from math import sqrt
    eps = 1e-6
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []
    # plane intersection (cut disc)
    denom = Dx + Dy
    if abs(denom) > eps:
        t_plane = (1 - (Ox + Oy)) / denom
        if t_plane > eps:
            Px = Ox + t_plane * Dx; Py = Oy + t_plane * Dy; Pz = Oz + t_plane * Dz
            vx = Px - 0.5; vy = Py - 0.5; vz = Pz
            if vx*vx + vy*vy + vz*vz <= 0.5 + eps:
                # outward normal is (1,1,0)
                if Dx + Dy < 0:
                    ts.append(t_plane)
    # sphere intersection
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2 * (Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = sqrt(disc)
        inv2a = 1 / (2*a)
        for t in ((-b - sqrt_disc) * inv2a, (-b + sqrt_disc) * inv2a):
            if t > eps:
                Px = Ox + t*Dx; Py = Oy + t*Dy; Pz = Oz + t*Dz
                if Px + Py <= 1 + eps:
                    # sphere normal = P (radius=1)
                    if Dx*Px + Dy*Py + Dz*Pz < 0:
                        ts.append(t)
                        break
    return min(ts) if ts else -1.0

#---------- 10.0 ----------#
def intersect(O, D):
    from math import sqrt
    def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def sphere_hits(O, D, C, r):
        OC = sub(O, C)
        a = dot(D, D)
        b = 2 * dot(D, OC)
        c = dot(OC, OC) - r*r
        disc = b*b - 4*a*c
        if disc < 0: return None
        s = sqrt(disc)
        t0 = (-b - s) / (2*a)
        t1 = (-b + s) / (2*a)
        if t0 > t1: t0, t1 = t1, t0
        return t0, t1
    h1 = sphere_hits(O, D, (0.0,0.0,0.0), 1.0)
    if h1 is None: return -1
    t01, t11 = h1
    h2 = sphere_hits(O, D, (0.5,0.5,0.0), 1.0)
    intervals = []
    if h2 is None:
        intervals.append((t01, t11))
    else:
        t02, t12 = h2
        a, b = t01, t11
        c, d = t02, t12
        if c > b or d < a:
            intervals.append((a, b))
        else:
            if a < c: intervals.append((a, min(b, c)))
            if d < b: intervals.append((max(a, d), b))
    t_hit = None
    for start, end in intervals:
        if end < 0: continue
        t = start if start >= 0 else end
        if t_hit is None or t < t_hit: t_hit = t
    return t_hit if t_hit is not None else -1

#---------- 10.1 ----------#
def intersect_sphere_hole(O, D):
    import math
    def dot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    a = dot(D, D)
    oc = O
    b = 2 * dot(oc, D)
    c = dot(oc, oc) - 1.0
    disc = b*b - 4*a*c
    if disc < 0: return -1.0
    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    t_min, t_max = min(t1, t2), max(t1, t2)
    oc2 = (O[0]-0.5, O[1]-0.5, O[2])
    b2 = 2 * dot(oc2, D)
    c2 = dot(oc2, oc2) - 1.0
    disc2 = b2*b2 - 4*a*c2
    segments = []
    if disc2 < 0:
        segments.append((t_min, t_max))
    else:
        sqrt_disc2 = math.sqrt(disc2)
        u1 = (-b2 - sqrt_disc2) / (2*a)
        u2 = (-b2 + sqrt_disc2) / (2*a)
        u_min, u_max = min(u1, u2), max(u1, u2)
        if u_max <= t_min or u_min >= t_max:
            segments.append((t_min, t_max))
        else:
            if u_min > t_min:
                segments.append((t_min, min(u_min, t_max)))
            if u_max < t_max:
                segments.append((max(u_max, t_min), t_max))
    t_hit = None
    for s, e in segments:
        if e < 0: continue
        t_candidate = s if s > 0 else e
        if t_candidate <= 0: continue
        if t_hit is None or t_candidate < t_hit:
            t_hit = t_candidate
    return t_hit if t_hit is not None else -1.0

#---------- 10.2 ----------#
def intersect_sphere_hole(O, D):
    import math
    def dot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def sub(u, v): return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    # Sphere0: center (0,0,0), r=1
    a = dot(D, D)
    b = 2 * dot(O, D)
    c = dot(O, O) - 1
    disc = b*b - 4*a*c
    if disc < 0: return -1.0
    sd = math.sqrt(disc)
    t0a = (-b - sd) / (2*a)
    t0b = (-b + sd) / (2*a)
    t0_near, t0_far = (t0a, t0b) if t0a < t0b else (t0b, t0a)
    if t0_far < 0: return -1.0
    # Sphere1 (hole): center (0.5,0.5,0), r=1
    C1 = (0.5, 0.5, 0.0)
    OC1 = sub(O, C1)
    b1 = 2 * dot(OC1, D)
    c1 = dot(OC1, OC1) - 1
    disc1 = b1*b1 - 4*a*c1
    intervals = []
    if disc1 < 0:
        intervals = [(t0_near, t0_far)]
    else:
        sd1 = math.sqrt(disc1)
        t1a = (-b1 - sd1) / (2*a)
        t1b = (-b1 + sd1) / (2*a)
        t1_near, t1_far = (t1a, t1b) if t1a < t1b else (t1b, t1a)
        if t1_near > t0_near:
            intervals.append((t0_near, min(t1_near, t0_far)))
        if t1_far < t0_far:
            intervals.append((max(t1_far, t0_near), t0_far))
    t_hit = None
    for u, v in intervals:
        if v < 0: continue
        t = max(u, 0.0)
        if t <= v and (t_hit is None or t < t_hit):
            t_hit = t
    return t_hit if t_hit is not None else -1.0

#---------- 10.3 ----------#
def intersect(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    def solve(Cx, Cy, Cz, r):
        rx, ry, rz = Ox - Cx, Oy - Cy, Oz - Cz
        a = Dx*Dx + Dy*Dy + Dz*Dz
        b = 2*(Dx*rx + Dy*ry + Dz*rz)
        c = rx*rx + ry*ry + rz*rz - r*r
        disc = b*b - 4*a*c
        if disc < 0:
            return ()
        sd = disc**0.5
        t1 = (-b - sd) / (2*a)
        t2 = (-b + sd) / (2*a)
        return (t1, t2) if t1 <= t2 else (t2, t1)

    # Primary sphere (center=(0,0,0), r=1)
    I1 = solve(0.0, 0.0, 0.0, 1.0)
    if not I1:
        return -1.0
    t1, t2 = I1
    if t2 < 0.0:
        return -1.0
    a = t1 if t1 > 0.0 else 0.0
    b = t2

    # Subtracting sphere (center=(0.5,0.5,0), r=1)
    I2 = solve(0.5, 0.5, 0.0, 1.0)
    if not I2:
        return a
    u1, u2 = I2
    if u2 < 0.0:
        return a
    c = u1 if u1 > 0.0 else 0.0
    d = u2

    # Compute set difference of intervals [a,b] \ [c,d]
    if b <= c or a >= d:
        return a
    if a < c:
        return a
    if d <= b:
        return d
    return -1.0


#---------- 10.4 ----------#
def intersect(O, D):
    import math
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2 * (dx*ox + dy*oy + dz*oz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    # sphere hole
    ox2, oy2, oz2 = ox - 0.5, oy - 0.5, oz
    b2 = 2 * (dx*ox2 + dy*oy2 + dz*oz2)
    c2 = ox2*ox2 + oy2*oy2 + oz2*oz2 - 1
    disc2 = b2*b2 - 4*a*c2
    hole = False
    if disc2 >= 0:
        sqrt2 = math.sqrt(disc2)
        u1 = (-b2 - sqrt2) / (2*a)
        u2 = (-b2 + sqrt2) / (2*a)
        if u1 > u2:
            u1, u2 = u2, u1
        hole = True
    if t1 > t2:
        t1, t2 = t2, t1
    for t in (t1, t2):
        if t >= 0 and (not hole or t < u1 or t > u2):
            return t
    return -1.0

#---------- 10.5 ----------#
def intersect_sphere_hole(O, D):
    import math
    def intersect_sphere(O, D, C, r):
        px, py, pz = O[0]-C[0], O[1]-C[1], O[2]-C[2]
        a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2]
        b = 2*(D[0]*px + D[1]*py + D[2]*pz)
        c = px*px + py*py + pz*pz - r*r
        disc = b*b - 4*a*c
        if disc < 0: return None
        s = math.sqrt(disc)
        t0 = (-b - s) / (2*a)
        t1 = (-b + s) / (2*a)
        if t0 > t1: t0, t1 = t1, t0
        return (t0, t1)
    s1 = intersect_sphere(O, D, (0.0, 0.0, 0.0), 1.0)
    if s1 is None or s1[1] < 0: return -1.0
    t1_0, t1_1 = max(s1[0], 0.0), s1[1]
    s2 = intersect_sphere(O, D, (0.5, 0.5, 0.0), 1.0)
    if s2 is None or s2[1] < 0: return t1_0
    t2_0, t2_1 = max(s2[0], 0.0), s2[1]
    if t1_1 <= t2_0 or t1_0 >= t2_1: return t1_0
    if t1_0 < t2_0: return t1_0
    return t2_1 if t2_1 <= t1_1 else -1.0

#---------- 10.6 ----------#
def intersect_sphere_with_hole(O, D):
    import math
    ox, oy, oz = O
    dx, dy, dz = D
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0
    disc = b*b - 4*a*c
    if disc < 0:
        return -1.0
    sqrt_disc = math.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    t_near1, t_far1 = (t1, t2) if t1 <= t2 else (t2, t1)
    if t_far1 < 0:
        return -1.0
    # Sphere 2 (hole) centered at (0.5,0.5,0)
    ox2, oy2, oz2 = ox - 0.5, oy - 0.5, oz
    b2 = 2*(ox2*dx + oy2*dy + oz2*dz)
    c2 = ox2*ox2 + oy2*oy2 + oz2*oz2 - 1.0
    disc2 = b2*b2 - 4*a*c2
    if disc2 < 0:
        return t_near1 if t_near1 >= 0 else (t_far1 if t_far1 >= 0 else -1.0)
    sqrt_disc2 = math.sqrt(disc2)
    s1 = (-b2 - sqrt_disc2) / (2*a)
    s2 = (-b2 + sqrt_disc2) / (2*a)
    t_near2, t_far2 = (s1, s2) if s1 <= s2 else (s2, s1)
    # No overlap between sphere1 interval and hole interval
    if t_far2 <= t_near1 or t_near2 >= t_far1:
        return t_near1 if t_near1 >= 0 else (t_far1 if t_far1 >= 0 else -1.0)
    # Compute valid segments = [t_near1, t_far1] \ [t_near2, t_far2]
    candidates = []
    # Left segment
    left_start, left_end = t_near1, min(t_far1, t_near2)
    if left_end >= left_start and left_end >= 0:
        candidates.append(max(left_start, 0.0))
    # Right segment
    right_start, right_end = max(t_near1, t_far2), t_far1
    if right_end >= right_start and right_end >= 0:
        candidates.append(max(right_start, 0.0))
    return min(candidates) if candidates else -1.0

#---------- 10.7 ----------#
def intersect_sphere_with_hole(O, D):
    import math
    def dot(u, v):
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def sub(u, v):
        return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def intersect_sphere(C, r):
        OC = sub(O, C)
        A = dot(D, D)
        B = 2 * dot(D, OC)
        Cc = dot(OC, OC) - r*r
        disc = B*B - 4*A*Cc
        if disc < 0:
            return None
        s = math.sqrt(disc)
        t0 = (-B - s) / (2*A)
        t1 = (-B + s) / (2*A)
        return (t0, t1) if t0 <= t1 else (t1, t0)
    s1 = intersect_sphere((0.0,0.0,0.0), 1.0)
    if s1 is None:
        return -1.0
    a, b = s1
    s2 = intersect_sphere((0.5,0.5,0.0), 1.0)
    if s2 is None:
        return a if a >= 0 else -1.0
    c, d = s2
    o_start, o_end = max(a, c), min(b, d)
    if o_start >= o_end:
        return a if a >= 0 else -1.0
    if o_start <= a and o_end >= b:
        return -1.0
    t = o_end if o_start <= a else a
    return t if t >= 0 else -1.0

#---------- 10.8 ----------#
def ray_sphere_hole_intersection(O, D):
    def intersect_sphere(O, D, C, r):
        ox, oy, oz = O; dx, dy, dz = D; cx, cy, cz = C
        Lx, Ly, Lz = ox - cx, oy - cy, oz - cz
        a = dx*dx + dy*dy + dz*dz
        b = 2*(dx*Lx + dy*Ly + dz*Lz)
        c = Lx*Lx + Ly*Ly + Lz*Lz - r*r
        disc = b*b - 4*a*c
        if disc < 0: return ()
        sqrt_disc = disc**0.5
        t1 = (-b - sqrt_disc)/(2*a)
        t2 = (-b + sqrt_disc)/(2*a)
        if abs(disc) < 1e-9: return (t1,)
        if t1 < t2: return (t1, t2)
        return (t2, t1)
    C1, r1 = (0.0, 0.0, 0.0), 1.0
    C2, r2 = (0.5, 0.5, 0.0), 1.0
    eps = 1e-6
    candidates = []
    for t in intersect_sphere(O, D, C1, r1):
        if t < 0: continue
        px, py, pz = O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2]
        dx2, dy2, dz2 = px-C2[0], py-C2[1], pz-C2[2]
        if dx2*dx2 + dy2*dy2 + dz2*dz2 >= r2*r2 - eps:
            candidates.append(t)
    for t in intersect_sphere(O, D, C2, r2):
        if t < 0: continue
        px, py, pz = O[0]+t*D[0], O[1]+t*D[1], O[2]+t*D[2]
        dx1, dy1, dz1 = px-C1[0], py-C1[1], pz-C1[2]
        if dx1*dx1 + dy1*dy1 + dz1*dz1 <= r1*r1 + eps:
            candidates.append(t)
    return min(candidates) if candidates else -1.0

#---------- 10.9 ----------#
def ray_sphere_hole_intersection(O, D):
    from math import sqrt
    eps = 1e-6
    def dot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    def sub(u, v): return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
    def add(u, v): return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
    def mul(u, s): return (u[0]*s, u[1]*s, u[2]*s)
    def sphere_intersect(C):
        OC = sub(O, C)
        a = dot(D, D)
        b = 2 * dot(D, OC)
        c = dot(OC, OC) - 1.0
        disc = b*b - 4*a*c
        if disc < 0: return []
        if abs(disc) < eps: return [-b / (2*a)]
        sd = sqrt(disc)
        return [(-b - sd) / (2*a), (-b + sd) / (2*a)]
    ts = sphere_intersect((0.0, 0.0, 0.0)) + sphere_intersect((0.5, 0.5, 0.0))
    ts.sort()
    prev = None
    for t in ts:
        if t < 0: continue
        if prev is not None and abs(t - prev) < eps: continue
        prev = t
        P = add(O, mul(D, t))
        d1 = sqrt(dot(P, P))
        v2 = sub(P, (0.5, 0.5, 0.0))
        d2 = sqrt(dot(v2, v2))
        if abs(d1 - 1.0) < eps and d2 > 1.0 - eps: return t
        if abs(d2 - 1.0) < eps and d1 < 1.0 + eps: return t
    return -1.0


#---------- 11.0 ----------#
def intersect(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    hits = []
    # Sphere intersection
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc >= 0:
        sd = math.sqrt(disc)
        for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
            if t >= 0:
                y = Oy + t*Dy
                z = Oz + t*Dz
                if y*y + z*z >= 0.49:
                    hits.append(t)
    # Cylinder intersection (axis along x)
    ac = Dy*Dy + Dz*Dz
    if ac != 0:
        bc = 2*(Oy*Dy + Oz*Dz)
        cc = Oy*Oy + Oz*Oz - 0.49
        disc_c = bc*bc - 4*ac*cc
        if disc_c >= 0:
            sc = math.sqrt(disc_c)
            for t in [(-bc - sc)/(2*ac), (-bc + sc)/(2*ac)]:
                if t >= 0:
                    x = Ox + t*Dx
                    y = Oy + t*Dy
                    z = Oz + t*Dz
                    if x*x + y*y + z*z <= 1:
                        hits.append(t)
    return min(hits) if hits else -1.0

#---------- 11.1 ----------#
def intersect_sphere_cylinder(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Sphere intersection
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0
    disc = b*b - 4*a*c
    if disc < 0: return -1.0
    sd = disc**0.5
    t0 = (-b - sd)/(2*a)
    t1 = (-b + sd)/(2*a)
    if t0 > t1: t0, t1 = t1, t0
    # Cylinder intersection (along x-axis)
    rc2 = 0.7*0.7
    ac = Dy*Dy + Dz*Dz
    bc = 2*(Oy*Dy + Oz*Dz)
    cc = Oy*Oy + Oz*Oz - rc2
    if ac == 0:
        if cc < 0: return -1.0
        return t0 if t0 > 0 else (t1 if t1 > 0 else -1.0)
    disc_c = bc*bc - 4*ac*cc
    if disc_c <= 0:
        if cc <= 0: return -1.0
        return t0 if t0 > 0 else (t1 if t1 > 0 else -1.0)
    sc = disc_c**0.5
    c0 = (-bc - sc)/(2*ac)
    c1 = (-bc + sc)/(2*ac)
    if c0 > c1: c0, c1 = c1, c0
    # Determine intersections with sphere minus cylinder
    ts = []
    # before cylinder hole
    if t0 < c0:
        t_end = min(t1, c0)
        if t_end > t0:
            if t0 > 0: ts.append(t0)
            elif t_end > 0: ts.append(t_end)
    # after cylinder hole
    if t1 > c1:
        t_start = max(t0, c1)
        if t1 > t_start:
            if t_start > 0: ts.append(t_start)
            elif t1 > 0: ts.append(t1)
    return min(ts) if ts else -1.0

#---------- 11.2 ----------#
def intersect(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    # Sphere intersection (radius 1 at origin)
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1.0
    disc = b*b - 4*a*c
    events = []
    if disc >= 0.0:
        sd = math.sqrt(disc)
        t0 = (-b - sd) / (2*a)
        t1 = (-b + sd) / (2*a)
        events.append((t0, "sphere"))
        events.append((t1, "sphere"))
    # Cylinder intersection (axis along x, radius 0.7)
    rr = 0.7*0.7
    a2 = Dy*Dy + Dz*Dz
    if a2 != 0.0:
        b2 = 2*(Oy*Dy + Oz*Dz)
        c2 = Oy*Oy + Oz*Oz - rr
        disc2 = b2*b2 - 4*a2*c2
        if disc2 >= 0.0:
            sd2 = math.sqrt(disc2)
            tc0 = (-b2 - sd2) / (2*a2)
            tc1 = (-b2 + sd2) / (2*a2)
            events.append((tc0, "cylinder"))
            events.append((tc1, "cylinder"))
    # Filter and sort events
    events = [(t, typ) for t, typ in events if t >= 0.0]
    events.sort(key=lambda e: e[0])
    # Test each event for actual boundary intersection
    for t, typ in events:
        x = Ox + Dx*t
        y = Oy + Dy*t
        z = Oz + Dz*t
        if typ == "sphere":
            # sphere surface outside hole
            if y*y + z*z > rr:
                return t
        else:
            # cylinder surface inside sphere
            if x*x + y*y + z*z < 1.0:
                return t
    return -1.0

#---------- 11.3 ----------#
def intersect(O, D):
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    candidates = []
    # Sphere intersection (radius=1)
    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    C = Ox*Ox + Oy*Oy + Oz*Oz - 1.0
    disc = B*B - 4*A*C
    if A != 0 and disc >= 0:
        sqrt_disc = disc**0.5
        for t in ((-B - sqrt_disc)/(2*A), (-B + sqrt_disc)/(2*A)):
            if t >= 0:
                y = Oy + t*Dy
                z = Oz + t*Dz
                if y*y + z*z >= 0.49:
                    candidates.append(t)
    # Cylinder intersection (radius=0.7 along x-axis)
    a = Dy*Dy + Dz*Dz
    if a != 0:
        b = 2*(Oy*Dy + Oz*Dz)
        c = Oy*Oy + Oz*Oz - 0.49
        disc_c = b*b - 4*a*c
        if disc_c >= 0:
            sqrt_c = disc_c**0.5
            for t in ((-b - sqrt_c)/(2*a), (-b + sqrt_c)/(2*a)):
                if t >= 0:
                    x = Ox + t*Dx
                    y = Oy + t*Dy
                    z = Oz + t*Dz
                    if x*x + y*y + z*z <= 1.0:
                        candidates.append(t)
    return min(candidates) if candidates else -1.0


#---------- 11.4 ----------#
def intersect_ray_sphere_with_cylinder_hole(O, D):
    import math
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    t_min = float('inf')
    # Sphere intersection
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = math.sqrt(disc)
        ts = [(-b - sqrt_disc) / (2*a), (-b + sqrt_disc) / (2*a)]
        for t in sorted(ts):
            if t >= 0:
                y = Oy + t*Dy
                z = Oz + t*Dz
                if y*y + z*z >= 0.49:
                    t_min = t
                    break
    # Cylinder intersection (along x-axis, radius 0.7)
    a_c = Dy*Dy + Dz*Dz
    if a_c != 0:
        b_c = 2*(Oy*Dy + Oz*Dz)
        c_c = Oy*Oy + Oz*Oz - 0.49
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0:
            sqrt_disc_c = math.sqrt(disc_c)
            ts_c = [(-b_c - sqrt_disc_c) / (2*a_c), (-b_c + sqrt_disc_c) / (2*a_c)]
            for t in sorted(ts_c):
                if t >= 0:
                    x = Ox + t*Dx
                    y = Oy + t*Dy
                    z = Oz + t*Dz
                    if x*x + y*y + z*z <= 1:
                        t_min = min(t_min, t)
                        break
    return t_min if t_min != float('inf') else -1


#---------- 11.5 ----------#
def intersect(O, D):
    from math import sqrt
    eps = 1e-6
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    def is_inside(px, py, pz):
        return px*px + py*py + pz*pz <= 1 and py*py + pz*pz >= 0.49
    candidates = []
    # Sphere intersection
    A = Dx*Dx + Dy*Dy + Dz*Dz
    B = 2 * (Ox*Dx + Oy*Dy + Oz*Dz)
    C = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = B*B - 4*A*C
    if disc >= 0:
        sd = sqrt(disc)
        for t in sorted([(-B - sd)/(2*A), (-B + sd)/(2*A)]):
            if t > eps:
                y = Oy + t*Dy
                z = Oz + t*Dz
                if y*y + z*z >= 0.49:
                    candidates.append(t)
    # Cylinder intersection
    A = Dy*Dy + Dz*Dz
    B = 2 * (Oy*Dy + Oz*Dz)
    C = Oy*Oy + Oz*Oz - 0.49
    disc = B*B - 4*A*C
    if A != 0 and disc >= 0:
        sd = sqrt(disc)
        for t in sorted([(-B - sd)/(2*A), (-B + sd)/(2*A)]):
            if t > eps:
                x = Ox + t*Dx
                y = Oy + t*Dy
                z = Oz + t*Dz
                if x*x + y*y + z*z <= 1:
                    candidates.append(t)
    if not candidates:
        return -1.0
    candidates.sort()
    for t in candidates:
        px, py, pz = Ox + t*Dx, Oy + t*Dy, Oz + t*Dz
        p0x, p0y, p0z = Ox + (t-eps)*Dx, Oy + (t-eps)*Dy, Oz + (t-eps)*Dz
        if not is_inside(p0x, p0y, p0z) and is_inside(px, py, pz):
            return t
    return -1.0

#---------- 11.6 ----------#
def intersect_sphere_cylinder_hole(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    best_t = float('inf')
    # Sphere intersection (radius=1)
    a_s = dx*dx + dy*dy + dz*dz
    b_s = 2*(ox*dx + oy*dy + oz*dz)
    c_s = ox*ox + oy*oy + oz*oz - 1.0
    disc_s = b_s*b_s - 4*a_s*c_s
    if disc_s >= 0:
        sqrt_s = disc_s**0.5
        for t in sorted([(-b_s - sqrt_s)/(2*a_s), (-b_s + sqrt_s)/(2*a_s)]):
            if t > 0:
                px, py, pz = ox + t*dx, oy + t*dy, oz + t*dz
                if py*py + pz*pz >= 0.7*0.7:
                    best_t = t
                    break
    # Cylinder intersection (axis=x, radius=0.7)
    a_c = dy*dy + dz*dz
    if a_c > 1e-8:
        b_c = 2*(oy*dy + oz*dz)
        c_c = oy*oy + oz*oz - 0.7*0.7
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0:
            sqrt_c = disc_c**0.5
            for t in sorted([(-b_c - sqrt_c)/(2*a_c), (-b_c + sqrt_c)/(2*a_c)]):
                if 0 < t < best_t:
                    px, py, pz = ox + t*dx, oy + t*dy, oz + t*dz
                    if px*px + py*py + pz*pz <= 1.0:
                        best_t = t
                        break
    return best_t if best_t < float('inf') else -1.0


#---------- 11.7 ----------#
def intersect_sphere_cylinder_hole(O, D):
    import math
    ox, oy, oz = O
    dx, dy, dz = D
    eps = 1e-6
    ts = []
    # Sphere intersection
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1.0
    disc = b*b - 4*a*c
    if disc >= 0:
        sd = math.sqrt(disc)
        for t in [(-b - sd)/(2*a), (-b + sd)/(2*a)]:
            if t > eps:
                y = oy + dy*t; z = oz + dz*t
                if y*y + z*z >= 0.49 - eps:
                    ts.append(t)
    # Cylinder intersection
    a_c = dy*dy + dz*dz
    if a_c > eps:
        b_c = 2*(oy*dy + oz*dz)
        c_c = oy*oy + oz*oz - 0.49
        disc_c = b_c*b_c - 4*a_c*c_c
        if disc_c >= 0:
            sd_c = math.sqrt(disc_c)
            for t in [(-b_c - sd_c)/(2*a_c), (-b_c + sd_c)/(2*a_c)]:
                if t > eps:
                    x = ox + dx*t; y = oy + dy*t; z = oz + dz*t
                    if x*x + y*y + z*z <= 1.0 + eps:
                        ts.append(t)
    return min(ts) if ts else -1.0

#---------- 11.8 ----------#
def intersect_sphere_cyl_hole(O, D):
    import math
    eps = 1e-6
    Ox, Oy, Oz = O
    Dx, Dy, Dz = D
    ts = []
    # Sphere intersection
    a = Dx*Dx + Dy*Dy + Dz*Dz
    b = 2*(Ox*Dx + Oy*Dy + Oz*Dz)
    c = Ox*Ox + Oy*Oy + Oz*Oz - 1
    disc = b*b - 4*a*c
    if disc >= 0:
        sd = math.sqrt(disc)
        t1 = (-b - sd) / (2*a)
        t2 = (-b + sd) / (2*a)
        for t in sorted((t1, t2)):
            if t >= eps:
                py = Oy + Dy*t
                pz = Oz + Dz*t
                if py*py + pz*pz >= 0.49 - eps:
                    ts.append(t)
                    break
    # Cylinder intersection
    a_c = Dy*Dy + Dz*Dz
    b_c = 2*(Oy*Dy + Oz*Dz)
    c_c = Oy*Oy + Oz*Oz - 0.49
    disc_c = b_c*b_c - 4*a_c*c_c
    if a_c > eps and disc_c >= 0:
        sd_c = math.sqrt(disc_c)
        tc1 = (-b_c - sd_c) / (2*a_c)
        tc2 = (-b_c + sd_c) / (2*a_c)
        for t in sorted((tc1, tc2)):
            if t >= eps:
                px = Ox + Dx*t
                py = Oy + Dy*t
                pz = Oz + Dz*t
                if px*px + py*py + pz*pz <= 1 + eps:
                    ts.append(t)
                    break
    return min(ts) if ts else -1


#---------- 11.9 ----------#
def intersect_sphere_with_cylinder_hole(O, D):
    ox, oy, oz = O
    dx, dy, dz = D
    eps = 1e-6
    ts = []
    # Sphere intersection
    a = dx*dx + dy*dy + dz*dz
    b = 2*(ox*dx + oy*dy + oz*dz)
    c = ox*ox + oy*oy + oz*oz - 1
    disc = b*b - 4*a*c
    if disc >= 0:
        sqrt_disc = disc**0.5
        t1 = (-b - sqrt_disc)/(2*a)
        t2 = (-b + sqrt_disc)/(2*a)
        for t in sorted((t1, t2)):
            if t >= eps:
                py = oy + t*dy
                pz = oz + t*dz
                if py*py + pz*pz >= 0.7*0.7 - eps:
                    ts.append(t)
                    break
    # Cylinder intersection (infinite along x)
    A = dy*dy + dz*dz
    if A > eps:
        B = 2*(oy*dy + oz*dz)
        C = oy*oy + oz*oz - 0.7*0.7
        disc_c = B*B - 4*A*C
        if disc_c >= 0:
            sqrt_dc = disc_c**0.5
            tc1 = (-B - sqrt_dc)/(2*A)
            tc2 = (-B + sqrt_dc)/(2*A)
            for t in sorted((tc1, tc2)):
                if t >= eps:
                    px = ox + t*dx
                    py = oy + t*dy
                    pz = oz + t*dz
                    if px*px + py*py + pz*pz <= 1 + eps:
                        ts.append(t)
                        break
    if not ts:
        return -1.0
    return min(ts)

