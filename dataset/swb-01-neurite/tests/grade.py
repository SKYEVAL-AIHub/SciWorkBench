"""Measure geometric coverage and branch structure in physical space."""

from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree


def read_swc(path):
    points, parents = {}, {}
    scale = 1
    for line in path.read_text().splitlines():
        if line.strip().lower() == "# units: voxel":
            scale = np.array([0.525, 0.525, 1.5])
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        node, kind, x, y, z, radius, parent = fields
        node, parent = int(node), int(parent)
        values = np.array([float(x), float(y), float(z), float(radius)])
        if node <= 0 or node in points or not np.isfinite(values).all() or values[3] < 0:
            raise ValueError("Invalid SWC record")
        int(kind)
        points[node], parents[node] = values[:3], parent
    children = defaultdict(list)
    for node, parent in parents.items():
        if parent != -1:
            if parent not in points or node == parent:
                raise ValueError("Invalid parent")
            children[parent].append(node)
    roots = [node for node, parent in parents.items() if parent == -1]
    seen = set()
    pending = roots.copy()
    while pending:
        node = pending.pop()
        seen.add(node)
        pending.extend(children[node])
    if len(seen) != len(points):
        raise ValueError("Cycle in SWC")
    points = {node: value * scale for node, value in points.items()}
    branches = []
    for node in points:
        if parents[node] != -1 and len(children[node]) == 1:
            continue
        for child in children[node]:
            chain = [node, child]
            while len(children[chain[-1]]) == 1:
                chain.append(children[chain[-1]][0])
            vertices = np.array([points[i] for i in chain])
            lengths = np.r_[0, np.cumsum(np.linalg.norm(np.diff(vertices, axis=0), axis=1))]
            if lengths[-1] <= 0 or lengths[-1] > 3000:
                raise ValueError("Invalid branch length")
            steps = np.linspace(0, lengths[-1], int(np.ceil(lengths[-1])) + 1)
            samples = np.column_stack([np.interp(steps, lengths, axis) for axis in vertices.T])
            branches.append((samples, lengths[-1]))
    if not branches:
        raise ValueError("Empty skeleton")
    terminals = [points[node] for node in points if not children[node]]
    junctions = [points[node] for node in points if len(children[node]) > 1]
    return branches, terminals, junctions


def distances(first, second):
    return cKDTree(second).query(first)[0], cKDTree(first).query(second)[0]


def f1(a, b):
    return 2 * a * b / (a + b) if a + b else 0


def landmark_f1(first, second, radius):
    if not first or not second:
        return float(not first and not second)
    pairs = []
    for i, candidates in enumerate(cKDTree(second).query_ball_point(first, radius)):
        pairs.extend((np.linalg.norm(first[i] - second[j]), i, j) for j in candidates)
    matched_first, matched_second = set(), set()
    for _, i, j in sorted(pairs):
        if i not in matched_first and j not in matched_second:
            matched_first.add(i)
            matched_second.add(j)
    return 2 * len(matched_first) / (len(first) + len(second))


def grade(output, reference):
    try:
        first, tips_a, forks_a = read_swc(output / "reconstruction.swc")
    except (OSError, ValueError):
        return 0.0
    second, tips_b, forks_b = read_swc(reference / "reconstruction.swc")
    a, b = distances(np.concatenate([p for p, _ in first]), np.concatenate([p for p, _ in second]))
    cov2 = f1(np.mean(a <= 2), np.mean(b <= 2))
    cov5 = f1(np.mean(a <= 5), np.mean(b <= 5))
    ssd = np.clip((8 - (a.mean() + b.mean()) / 2) / 6.5, 0, 1)
    candidates = []
    for i, (a_path, _) in enumerate(first):
        for j, (b_path, _) in enumerate(second):
            if not (np.linalg.norm(a_path[[0, -1]] - b_path[[0, -1]], axis=1) <= 5).all():
                continue
            path_a, path_b = distances(a_path, b_path)
            overlap = min(np.mean(path_a <= 5), np.mean(path_b <= 5))
            if overlap >= .70:
                candidates.append((overlap, i, j))
    used_a, used_b = set(), set()
    for _, i, j in sorted(candidates, reverse=True):
        if i not in used_a and j not in used_b:
            used_a.add(i)
            used_b.add(j)
    length_a, length_b = (sum(length for _, length in paths) for paths in (first, second))
    length_score = np.clip((.30 - abs(length_a / length_b - 1)) / .25, 0, 1)
    path_score = 2 * len(used_a) / (len(first) + len(second))
    matched_length = f1(sum(first[i][1] for i in used_a) / length_a,
                        sum(second[j][1] for j in used_b) / length_b)
    return float((12 + 12 * cov2 + 6 * cov5 + 5 * ssd + 24 * path_score +
                  13 * landmark_f1(tips_a, tips_b, 5) + 13 * landmark_f1(forks_a, forks_b, 3) +
                  5 * length_score + 10 * matched_length) / 100)


if __name__ == "__main__":
    print(grade(*(Path(arg) for arg in sys.argv[1:])))
