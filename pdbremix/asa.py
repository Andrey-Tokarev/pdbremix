# encoding: utf-8

__doc__ = """ 

Calculate the accessible-surface area of atoms.

Uses the simple Shrake-Rupley algorithm, that generates a
relatively uniform density of dots over every atoms and
eliminates those within the sphere of another atom. The remaining
dots is used to calculate the area.

Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to
Solvent of Protein Atoms. Lysozyme and Insulin." J Mol Biol. 79
(1973) 351- 371. """

import math

import v3
import pdbatoms
from collections import defaultdict


def generate_sphere_points(n):
  """
  Returns list of coordinates on a sphere using the Golden-
  Section Spiral algorithm.
  """
  points = []
  inc = math.pi * (3 - math.sqrt(5))
  offset = 2 / float(n)
  for k in range(int(n)):
    y = k * offset - 1 + (offset / 2)
    r = math.sqrt(1 - y*y)
    phi = k * inc
    points.append(v3.vector(math.cos(phi)*r, y, math.sin(phi)*r))
  return points


def find_neighbor_indices(atoms, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  indices = range(k)
  indices.extend(range(k+1, len(atoms)))
  for i in indices:
    atom_i = atoms[i]
    dist = v3.distance(atom_k.pos, atom_i.pos)
    if dist < radius + atom_i.radius:
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa(atoms, probe, n_sphere_point=960):
  """
  Returns the accessible-surface areas of the atoms, by rolling a
  ball with probe radius over the atoms with their radius
  defined.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = find_neighbor_indices(atoms, probe, i)
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = range(j_closest_neighbor, n_neighbor)
      cycled_indices.extend(range(j_closest_neighbor))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff = v3.distance(atom_j.pos, test_point)
        if diff*diff < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)

  return areas


def makeBoxes(A, dmax):
    '''
    returns dictionary which keys are indecies of boxes (regions)
    with dmax length side and 
    values are indicies of atoms belonging to these boxes
    '''
    B = defaultdict(list) # space divided into boxes
    for i in xrange(len(A)):
        atom = A[i]
        boxCoor = tuple(int(math.floor(x / dmax)) for x in atom.pos)
        B[boxCoor].append(i)
    return B


def addBond(A, a1, a2, conn, dmax):
    '''
    add neighboring atoms (if distance to it is less than dmax)
    to atoms a1 and a2 in adjacency list conn
    '''
    atom1 = A[a1]
    atom2 = A[a2]
    if v3.mag2(atom1.pos - atom2.pos) <= dmax * dmax:  # connected
        conn[a1].append(a2)
        conn[a2].append(a1)


def neighbAtoms(B, box):
    '''
    returns list of atoms from half of neighbouring boxes of the box
    another half is accounted when symmetric (opposite) boxes considered
    '''
    na = [] # list for neighboring atoms
    x, y, z = box # coordinates of the box
    # top layer consisting of 9 boxes
    if (x + 1, y + 1, z +1) in B: na.extend(B[(x + 1, y + 1, z +1)])
    if (x, y + 1, z +1) in B: na.extend(B[(x, y + 1, z +1)])
    if (x + 1, y, z +1) in B: na.extend(B[(x + 1, y, z +1)])
    if (x, y, z +1) in B: na.extend(B[(x, y, z +1)])
    if (x - 1, y + 1, z +1) in B: na.extend(B[(x - 1, y + 1, z +1)])
    if (x + 1, y - 1, z +1) in B: na.extend(B[(x + 1, y - 1, z +1)])
    if (x, y - 1, z +1) in B: na.extend(B[(x, y - 1, z +1)])
    if (x - 1, y, z +1) in B: na.extend(B[(x - 1, y, z +1)])
    if (x - 1, y - 1, z +1) in B: na.extend(B[(x - 1, y - 1, z +1)])
    # half of the middle layer excluding the box itself (4 boxes)
    if (x + 1, y + 1, z) in B: na.extend(B[(x + 1, y + 1, z)])
    if (x, y + 1, z) in B: na.extend(B[(x, y + 1, z)])
    if (x + 1, y, z) in B: na.extend(B[(x + 1, y, z)])
    if (x + 1, y - 1, z) in B: na.extend(B[(x + 1, y - 1, z)])
    return na


def adjList(A, dmax):
    '''
    returns adjacency list from coordinate file
    in O(len(A)) time
    '''
    B = makeBoxes(A, dmax) # put atoms into the boxes with dmax length side
    # now go on boxes and check connections inside 3x3 superboxes
    conn = [[] for i in xrange(len(A))] # list of bond lengths each atom implicated
    for box in B:
        lb = len(B[box])
        for i in range(lb):
            a1 = B[box][i]
            # check possible connections inside the box
            for j in range(i+1, lb):
                a2 = B[box][j]
                addBond(A, a1, a2, conn, dmax)
            # check connections with atoms from neighbouring boxes
            na = neighbAtoms(B, box) # list of such atoms
            for a2 in na:
                addBond(A, a1, a2, conn, dmax)
    return conn


def find_neighbor_indices_mod(atoms, indices, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  for i in indices:
    if i == k: continue
    atom_i = atoms[i]
    dist2 = v3.mag2(atom_k.pos - atom_i.pos) # ToAn
    if dist2 < (radius + atom_i.radius) ** 2: # ToAn
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa_opt(atoms, probe, n_sphere_point=960):
  """
  Optimized version
  Returns the accessible-surface areas of the atoms, by rolling a
  ball with probe radius over the atoms with their radius
  defined.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  neighbor_list = adjList(atoms, 2 * (probe + max(atoms, key=lambda p: p.radius).radius))
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = [neig for neig in neighbor_list[i]]
    neighbor_indices = find_neighbor_indices_mod(atoms, neighbor_indices, probe, i) # even further narrow diapazon
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = range(j_closest_neighbor, n_neighbor)
      cycled_indices.extend(range(j_closest_neighbor))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff2 = v3.mag2(atom_j.pos - test_point)
        if diff2 < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)

  return areas

