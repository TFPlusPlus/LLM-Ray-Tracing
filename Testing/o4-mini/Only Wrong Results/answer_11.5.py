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
