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