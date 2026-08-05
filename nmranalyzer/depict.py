"""2D coordinates and drawing for a parsed SMILES structure.

Rings are laid out as regular polygons and fused ring systems are built up by
reflecting each new ring across the bond it shares with one already placed;
acyclic parts are then grown outwards with the usual 120 degree zig-zag.  The
result is not a publication-grade depiction, but it is good enough to look at
a molecule and point at which proton is which.
"""

from __future__ import annotations

import math

BOND = 1.0                  # bond length in layout units
HETERO_ONLY = True          # carbons are drawn as vertices, not letters


# ---------------------------------------------------------------------------
# Ring perception
# ---------------------------------------------------------------------------


def find_rings(mol):
    """Smallest ring through each ring bond, de-duplicated.

    Not a strict SSSR, but for ordinary fused aromatics it finds the same
    rings a chemist would draw.
    """
    if not mol.bonds:
        return []

    bridges = _bridges(mol)
    rings = []
    seen = set()
    for bond in mol.bonds:
        if (bond.a.index, bond.b.index) in bridges:
            continue
        cycle = _smallest_cycle_through(bond)
        if not cycle:
            continue
        key = frozenset(cycle)
        if key not in seen:
            seen.add(key)
            rings.append(cycle)
    rings.sort(key=len)

    # Keep only as many independent rings as the cycle rank allows.
    wanted = mol.ring_count()
    if len(rings) > wanted:
        kept = []
        covered = set()
        for ring in rings:
            edges = _ring_edges(ring)
            if not edges <= covered:
                kept.append(ring)
                covered |= edges
            if len(kept) == wanted:
                break
        rings = kept
    return rings


def _ring_edges(ring):
    edges = set()
    for i, atom in enumerate(ring):
        other = ring[(i + 1) % len(ring)]
        edges.add(frozenset((atom.index, other.index)))
    return edges


def _bridges(mol):
    """Bonds that are in no cycle, found by depth-first low-link numbering."""
    result = set()
    order = {}
    low = {}
    counter = [0]

    def visit(atom, parent_bond):
        order[atom.index] = low[atom.index] = counter[0]
        counter[0] += 1
        for bond in atom.bonds:
            if bond is parent_bond:
                continue
            other = bond.other(atom)
            if other.index not in order:
                visit(other, bond)
                low[atom.index] = min(low[atom.index], low[other.index])
                if low[other.index] > order[atom.index]:
                    result.add((atom.index, other.index))
                    result.add((other.index, atom.index))
            else:
                low[atom.index] = min(low[atom.index], order[other.index])

    import sys
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(limit, len(mol.atoms) * 4 + 100))
    try:
        for atom in mol.atoms:
            if atom.index not in order:
                visit(atom, None)
    finally:
        sys.setrecursionlimit(limit)
    return result


def _smallest_cycle_through(bond):
    """Shortest path from one end of the bond to the other, avoiding the bond."""
    start, goal = bond.a, bond.b
    previous = {start.index: None}
    queue = [start]
    while queue:
        current = queue.pop(0)
        for edge in current.bonds:
            if edge is bond:
                continue
            other = edge.other(current)
            if other.index in previous:
                continue
            previous[other.index] = current
            if other is goal:
                path = [goal]
                node = current
                while node is not None:
                    path.append(node)
                    node = previous[node.index]
                return path
            queue.append(other)
    return None


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def layout(mol):
    """Return ``{atom_index: (x, y)}`` for every atom."""
    rings = find_rings(mol)
    systems = _ring_systems(rings)

    system_of = {}
    for sid, system in enumerate(systems):
        for ring in system:
            for atom in ring:
                system_of[atom.index] = sid
    local = [_layout_system_local(system) for system in systems]

    coords = {}
    done_systems = set()

    if systems:
        for index, point in local[0].items():
            coords[index] = point
        done_systems.add(0)
    elif mol.atoms:
        coords[mol.atoms[0].index] = (0.0, 0.0)

    ring_atoms = set(system_of)
    changed = True
    guard = 0
    while changed and guard < len(mol.atoms) * 4 + 20:
        guard += 1
        changed = False
        for atom in mol.atoms:
            if atom.index not in coords:
                continue
            ax, ay = coords[atom.index]
            for neighbour in atom.neighbours():
                if neighbour.index in coords:
                    continue
                angle = _free_angle(atom, coords, ring_atoms)
                target = (ax + BOND * math.cos(angle), ay + BOND * math.sin(angle))
                sid = system_of.get(neighbour.index)
                if sid is not None and sid not in done_systems:
                    _paste_system(local[sid], neighbour.index, target, angle, coords)
                    done_systems.add(sid)
                else:
                    coords[neighbour.index] = target
                changed = True

        if not changed:
            # A disconnected fragment needs its own starting point.
            missing = [a for a in mol.atoms if a.index not in coords]
            if missing:
                xs = [p[0] for p in coords.values()] or [0.0]
                start = (max(xs) + 2.5 * BOND, 0.0)
                sid = system_of.get(missing[0].index)
                if sid is not None and sid not in done_systems:
                    _paste_system(local[sid], missing[0].index, start, 0.0, coords)
                    done_systems.add(sid)
                else:
                    coords[missing[0].index] = start
                changed = True

    if not coords:
        for i, atom in enumerate(mol.atoms):
            coords[atom.index] = (i * BOND, 0.0)
    _relax(mol, coords)
    return coords


def _ring_systems(rings):
    """Group rings that share atoms into fused systems."""
    systems = []
    for ring in rings:
        members = set(a.index for a in ring)
        merged = [ring]
        rest = []
        for system in systems:
            if any(members & set(a.index for a in other) for other in system):
                merged.extend(system)
            else:
                rest.append(system)
        rest.append(merged)
        systems = rest
    return systems


def _layout_system_local(system):
    """Lay a fused ring system out on its own, in arbitrary position."""
    coords = {}
    _place_ring_system(system, coords, 0.0)
    return coords


def _paste_system(local, anchor_index, target, angle, coords):
    """Drop a pre-laid-out ring system onto the page.

    The system is rotated so that it grows away from the bond that reached it,
    then translated so ``anchor_index`` lands on ``target``.
    """
    if anchor_index not in local:
        coords[anchor_index] = target
        return
    ax, ay = local[anchor_index]
    cx = sum(p[0] for p in local.values()) / len(local)
    cy = sum(p[1] for p in local.values()) / len(local)
    current = math.atan2(cy - ay, cx - ax)
    rotation = angle - current
    cos_r, sin_r = math.cos(rotation), math.sin(rotation)
    for index, (x, y) in local.items():
        dx, dy = x - ax, y - ay
        coords[index] = (target[0] + dx * cos_r - dy * sin_r,
                         target[1] + dx * sin_r + dy * cos_r)


def _place_ring_system(system, coords, origin_x):
    system = sorted(system, key=len)
    placed = set()

    first = system[0]
    n = len(first)
    radius = BOND / (2.0 * math.sin(math.pi / n)) if n > 2 else BOND
    for k, atom in enumerate(first):
        angle = 2.0 * math.pi * k / n + math.pi / 2.0
        coords[atom.index] = (origin_x + radius * math.cos(angle),
                              radius * math.sin(angle))
        placed.add(atom.index)

    remaining = list(system[1:])
    guard = 0
    while remaining and guard < 200:
        guard += 1
        progressed = False
        for ring in list(remaining):
            shared = [a for a in ring if a.index in placed]
            if len(shared) < 2:
                continue
            if _place_fused_ring(ring, shared, coords, placed):
                remaining.remove(ring)
                progressed = True
        if not progressed:
            break
    return placed


def _place_fused_ring(ring, shared, coords, placed):
    """Place a ring that shares a bond with something already on the page."""
    n = len(ring)
    order = [a.index for a in ring]

    # Find two shared atoms that are adjacent in this ring.
    anchor = None
    for i in range(n):
        a = ring[i]
        b = ring[(i + 1) % n]
        if a.index in placed and b.index in placed:
            anchor = (i, a, b)
            break
    if anchor is None:
        return False
    start, a, b = anchor

    ax, ay = coords[a.index]
    bx, by = coords[b.index]
    mx, my = (ax + bx) / 2.0, (ay + by) / 2.0

    radius = BOND / (2.0 * math.sin(math.pi / n)) if n > 2 else BOND
    half = math.hypot(bx - ax, by - ay) / 2.0
    height = math.sqrt(max(radius * radius - half * half, 1e-6))

    dx, dy = bx - ax, by - ay
    length = math.hypot(dx, dy) or 1.0
    nx, ny = -dy / length, dx / length

    # Put the new ring's centre away from the existing atoms.
    candidates = [(mx + nx * height, my + ny * height),
                  (mx - nx * height, my - ny * height)]
    others = [coords[i] for i in placed if i not in (a.index, b.index)]

    def crowding(point):
        if not others:
            return 0.0
        return sum(1.0 / (1e-3 + math.hypot(point[0] - p[0], point[1] - p[1]))
                   for p in others)

    cx, cy = min(candidates, key=crowding)

    start_angle = math.atan2(ay - cy, ax - cx)
    end_angle = math.atan2(by - cy, bx - cx)
    step = 2.0 * math.pi / n
    # a and b are adjacent, so going round the ring from a to b is exactly one
    # step; whether that step is clockwise or not decides the walk direction.
    delta = (end_angle - start_angle) % (2.0 * math.pi)
    direction = 1.0 if abs(delta - step) < abs(delta - (2 * math.pi - step)) else -1.0

    for offset in range(1, n):
        atom = ring[(start + offset) % n]
        if atom.index in placed:
            continue
        angle = start_angle + direction * step * offset
        coords[atom.index] = (cx + radius * math.cos(angle),
                              cy + radius * math.sin(angle))
        placed.add(atom.index)
    return True


def _free_angle(atom, coords, ring_atoms):
    """Pick a direction that is as far as possible from existing neighbours."""
    taken = []
    ax, ay = coords[atom.index]
    for neighbour in atom.neighbours():
        if neighbour.index in coords:
            bx, by = coords[neighbour.index]
            taken.append(math.atan2(by - ay, bx - ax))
    if not taken:
        return math.radians(30.0)
    if len(taken) == 1:
        # Standard zig-zag; the sign alternates with the atom index so chains
        # do not fold back on themselves.
        turn = math.radians(120.0 if atom.index % 2 == 0 else -120.0)
        return taken[0] + turn

    best_angle = 0.0
    best_gap = -1.0
    for step in range(72):
        angle = 2.0 * math.pi * step / 72.0
        gap = min(abs(_wrap(angle - used)) for used in taken)
        if gap > best_gap:
            best_gap = gap
            best_angle = angle
    return best_angle


def _wrap(angle):
    while angle > math.pi:
        angle -= 2.0 * math.pi
    while angle < -math.pi:
        angle += 2.0 * math.pi
    return angle


def _relax(mol, coords, rounds=60):
    """Push overlapping atoms apart without disturbing bond lengths much."""
    ring_atoms = set()
    for ring in find_rings(mol):
        ring_atoms |= set(a.index for a in ring)

    indices = [a.index for a in mol.atoms]
    for _ in range(rounds):
        shift = {i: [0.0, 0.0] for i in indices}
        moved = False
        for i, first in enumerate(indices):
            x1, y1 = coords[first]
            for second in indices[i + 1:]:
                x2, y2 = coords[second]
                dx, dy = x2 - x1, y2 - y1
                distance = math.hypot(dx, dy)
                if distance >= 0.75 * BOND or distance == 0.0:
                    continue
                push = (0.75 * BOND - distance) / 2.0
                ux, uy = dx / distance, dy / distance
                if first not in ring_atoms:
                    shift[first][0] -= ux * push
                    shift[first][1] -= uy * push
                if second not in ring_atoms:
                    shift[second][0] += ux * push
                    shift[second][1] += uy * push
                moved = True
        if not moved:
            break
        for index in indices:
            dx, dy = shift[index]
            if dx or dy:
                x, y = coords[index]
                coords[index] = (x + dx * 0.5, y + dy * 0.5)


# ---------------------------------------------------------------------------
# Drawing primitives shared by the canvas and the SVG writer
# ---------------------------------------------------------------------------


def bounds(coords):
    xs = [p[0] for p in coords.values()]
    ys = [p[1] for p in coords.values()]
    if not xs:
        return 0.0, 0.0, 1.0, 1.0
    return min(xs), min(ys), max(xs), max(ys)


def transform(coords, width, height, margin=28.0):
    """Map layout units onto a pixel box, y flipped, aspect preserved."""
    x0, y0, x1, y1 = bounds(coords)
    span_x = max(x1 - x0, 1e-6)
    span_y = max(y1 - y0, 1e-6)
    usable_w = max(width - 2 * margin, 10.0)
    usable_h = max(height - 2 * margin, 10.0)
    scale = min(usable_w / span_x, usable_h / span_y)
    scale = min(scale, 62.0)                 # never blow a small molecule up
    off_x = (width - span_x * scale) / 2.0
    off_y = (height - span_y * scale) / 2.0

    out = {}
    for index, (x, y) in coords.items():
        out[index] = (off_x + (x - x0) * scale,
                      height - (off_y + (y - y0) * scale))
    return out, scale


def bond_segments(mol, screen, scale, rings=None):
    """Line segments for every bond, doubling up for double/aromatic bonds.

    The second line of a double or aromatic bond has to sit on the *inside* of
    its ring; drawing it always on the same side of the bond vector makes fused
    aromatics look like buckled polygons.
    """
    if rings is None:
        rings = find_rings(mol)
    centres = []
    for ring in rings:
        members = set(a.index for a in ring)
        cx = sum(screen[a.index][0] for a in ring) / len(ring)
        cy = sum(screen[a.index][1] for a in ring) / len(ring)
        centres.append((members, cx, cy))

    gap = max(2.5, scale * 0.09)
    shrink = max(6.0, scale * 0.19)
    segments = []
    for bond in mol.bonds:
        x1, y1 = screen[bond.a.index]
        x2, y2 = screen[bond.b.index]
        x1, y1, x2, y2 = _trim(bond, x1, y1, x2, y2, shrink)

        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length

        # Point the offset at the centre of whichever ring holds this bond.
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        inside = None
        for members, cx, cy in centres:
            if bond.a.index in members and bond.b.index in members:
                inside = (cx - mx, cy - my)
                break
        if inside is not None and (nx * inside[0] + ny * inside[1]) < 0:
            nx, ny = -nx, -ny

        if bond.aromatic:
            segments.append((x1, y1, x2, y2))
            if inside is not None:
                inner = 0.7
                hx, hy = (x2 - x1) * inner / 2.0, (y2 - y1) * inner / 2.0
                segments.append((mx - hx + nx * gap, my - hy + ny * gap,
                                 mx + hx + nx * gap, my + hy + ny * gap))
        elif bond.order >= 2:
            if inside is not None:
                # Ring double bond: main line on the bond, partner inside.
                segments.append((x1, y1, x2, y2))
                inner = 0.72
                hx, hy = (x2 - x1) * inner / 2.0, (y2 - y1) * inner / 2.0
                segments.append((mx - hx + nx * gap, my - hy + ny * gap,
                                 mx + hx + nx * gap, my + hy + ny * gap))
            else:
                segments.append((x1 + nx * gap, y1 + ny * gap,
                                 x2 + nx * gap, y2 + ny * gap))
                segments.append((x1 - nx * gap, y1 - ny * gap,
                                 x2 - nx * gap, y2 - ny * gap))
            if bond.order >= 3:
                segments.append((x1, y1, x2, y2))
        else:
            segments.append((x1, y1, x2, y2))
    return segments


def _trim(bond, x1, y1, x2, y2, shrink):
    """Pull a bond back from any atom that will be drawn as a letter."""
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy) or 1.0
    ux, uy = dx / length, dy / length
    if visible_label(bond.a):
        x1 += ux * shrink
        y1 += uy * shrink
    if visible_label(bond.b):
        x2 -= ux * shrink
        y2 -= uy * shrink
    return x1, y1, x2, y2


def visible_label(atom):
    """Carbon is left as a bare vertex; everything else gets its symbol."""
    if atom.symbol == "C" and atom.charge == 0 and not atom.isotope:
        return None
    text = atom.symbol
    if atom.n_hydrogens == 1:
        text += "H"
    elif atom.n_hydrogens > 1:
        text += "H%d" % atom.n_hydrogens
    if atom.charge:
        sign = "+" if atom.charge > 0 else "-"
        text += sign if abs(atom.charge) == 1 else "%s%d" % (sign, abs(atom.charge))
    return text


def render_svg(mol, width=340, height=260, highlight=(), title=""):
    """Draw the structure as a standalone SVG fragment."""
    coords = layout(mol)
    screen, scale = transform(coords, width, height)
    font = max(9.0, min(15.0, scale * 0.34))

    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
             'viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif">'
             % (width, height, width, height),
             '<rect width="%d" height="%d" fill="white"/>' % (width, height)]

    for index in highlight:
        if index in screen:
            x, y = screen[index]
            parts.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#ffe9a8"/>'
                         % (x, y, font * 1.25))

    for x1, y1, x2, y2 in bond_segments(mol, screen, scale):
        parts.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                     'stroke="#222" stroke-width="1.4"/>' % (x1, y1, x2, y2))

    for atom in mol.atoms:
        text = visible_label(atom)
        if not text:
            continue
        x, y = screen[atom.index]
        parts.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="white"/>'
                     % (x, y, font * 0.78))
        parts.append('<text x="%.1f" y="%.1f" font-size="%.1f" fill="#111" '
                     'text-anchor="middle" dominant-baseline="central">%s</text>'
                     % (x, y, font, text))

    if title:
        parts.append('<text x="6" y="14" font-size="11" fill="#666">%s</text>'
                     % title)
    parts.append("</svg>")
    return "\n".join(parts)
