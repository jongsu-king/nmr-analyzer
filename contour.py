"""Contour extraction by marching squares, tuned for 2D NMR maps.

A 2D spectrum is almost all baseline, so the inner loop first brackets each
cell by the min and max of its four corners and only visits the contour levels
that actually cross it.  That is what keeps a pure-Python contour of a
1024 x 1024 matrix interactive.
"""

from __future__ import annotations

import bisect


def levels(base, count=10, factor=1.5, negative=False):
    """A geometric ladder of contour levels starting at ``base``."""
    out = []
    value = float(base)
    for _ in range(max(1, count)):
        out.append(value)
        value *= factor
    return [-v for v in out] if negative else out


def downsample(data, rows, cols, max_rows, max_cols):
    """Block-maximum reduction, which keeps peaks rather than averaging them out.

    Returns ``(grid, row_step, col_step)`` so callers can map grid indices back
    to data indices.
    """
    row_step = max(1, rows // max_rows) if max_rows else 1
    col_step = max(1, cols // max_cols) if max_cols else 1
    if row_step == 1 and col_step == 1:
        return data, 1, 1

    grid = []
    for r0 in range(0, rows, row_step):
        block_rows = data[r0:r0 + row_step]
        out_row = []
        for c0 in range(0, cols, col_step):
            best = None
            for row in block_rows:
                chunk = row[c0:c0 + col_step]
                if not chunk:
                    continue
                local = max(chunk)
                if best is None or local > best:
                    best = local
            out_row.append(best if best is not None else 0.0)
        grid.append(out_row)
    return grid, row_step, col_step


def _interp(level, v0, v1):
    """Where between two corner values the contour crosses, as a fraction."""
    span = v1 - v0
    if span == 0.0:
        return 0.5
    t = (level - v0) / span
    return 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)


def segments(grid, level_list):
    """Contour line segments for every level.

    Returns ``{level: [((c0, r0), (c1, r1)), ...]}`` in fractional grid
    coordinates, where ``c`` runs along the row and ``r`` down the rows.
    """
    result = {level: [] for level in level_list}
    if not grid or len(grid) < 2 or len(grid[0]) < 2:
        return result

    ordered = sorted(level_list)
    rows = len(grid)
    cols = len(grid[0])

    for r in range(rows - 1):
        top = grid[r]
        bottom = grid[r + 1]
        for c in range(cols - 1):
            a = top[c]          # top-left
            b = top[c + 1]      # top-right
            d = bottom[c + 1]   # bottom-right
            e = bottom[c]       # bottom-left

            lo = a if a < b else b
            if d < lo:
                lo = d
            if e < lo:
                lo = e
            hi = a if a > b else b
            if d > hi:
                hi = d
            if e > hi:
                hi = e

            # Only the levels that actually cross this cell.
            start = bisect.bisect_left(ordered, lo)
            stop = bisect.bisect_right(ordered, hi)
            if start >= stop:
                continue

            for level in ordered[start:stop]:
                code = 0
                if a >= level:
                    code |= 8
                if b >= level:
                    code |= 4
                if d >= level:
                    code |= 2
                if e >= level:
                    code |= 1
                if code == 0 or code == 15:
                    continue

                # Crossing points on each edge, in cell-local coordinates.
                top_pt = (c + _interp(level, a, b), r)
                right_pt = (c + 1, r + _interp(level, b, d))
                bottom_pt = (c + _interp(level, e, d), r + 1)
                left_pt = (c, r + _interp(level, a, e))

                bucket = result[level]
                if code in (1, 14):
                    bucket.append((left_pt, bottom_pt))
                elif code in (2, 13):
                    bucket.append((bottom_pt, right_pt))
                elif code in (3, 12):
                    bucket.append((left_pt, right_pt))
                elif code in (4, 11):
                    bucket.append((top_pt, right_pt))
                elif code == 5:
                    bucket.append((left_pt, top_pt))
                    bucket.append((bottom_pt, right_pt))
                elif code in (6, 9):
                    bucket.append((top_pt, bottom_pt))
                elif code in (7, 8):
                    bucket.append((left_pt, top_pt))
                elif code == 10:
                    bucket.append((left_pt, bottom_pt))
                    bucket.append((top_pt, right_pt))
    return result
