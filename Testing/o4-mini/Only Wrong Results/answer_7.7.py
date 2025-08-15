def wrong_solution(O, D):
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
        if cross < 0.0: # Cross product in 2D??
            return t
    return -1.0

def model_solution(origin, direction):
    """
    The object is a plane with a triangular hole. The plane is y = z. The triangular hole has corners (1, 1, 1), (-1, 0, 0), and (0, -1, -1) in counterclockwise order. Light can pass through the triangular hole but not the plane.
    """
    Ox, Oy, Oz = origin
    Dx, Dy, Dz = direction
    if Dy == Dz:
        return -1
    t = (Oy - Oz) / (Dz - Dy)
    if t < 0:
        return -1
    Px, Py = Ox + t * Dx, Oy + t * Dy
    # Barycentric coordinates: p = p0 + (p1 - p0) * a + (p2 - p0) * b
    # p0 = (-1, 0, 0), p1 = (1, 1, 1), p2 = (0, -1, -1)
    a = (Px + Py + 1) / 3
    b = a - Py
    if a < 0 or a > 1 or b < 0 or b > 1 or a + b > 1:
        return t
    return -1
