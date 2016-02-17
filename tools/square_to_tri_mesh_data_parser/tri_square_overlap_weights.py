sqmesh_min_coord = [359919.189 - 360600.0, 3972158.559 - 3973000.0]
sqmesh_step = 2.0

import h5py
import math

tmesh_data = h5py.File("visdump_surface_mesh_jaramillo_384.h5",'r')
tmesh_key = '6234'

ME_len = len(tmesh_data[tmesh_key]['Mesh']['MixedElements'])
ntris = ME_len // 4;

tricells_inodes = [[0 for x in range(3)] for x in range(ntris)]

for i in range(ME_len):
  if i % 4 == 0:
    if tmesh_data[tmesh_key]['Mesh']['MixedElements'][i] != 4:
      raise RuntimeError("Mesh should only contain triangular cells!")
  else:
    tricells_inodes[i // 4][i % 4 - 1] = tmesh_data[tmesh_key]['Mesh']['MixedElements'][i]

tnodes = tmesh_data[tmesh_key]['Mesh']['Nodes']

def cohensutherland(left, top, right, bottom, x1, y1, x2, y2):
  """Clips a line to a rectangular area.
    
    This implements the Cohen-Sutherland line clipping algorithm.  left,
    top, right and bottom denote the clipping area, into which the line
    defined by x1, y1 (start point) and x2, y2 (end point) will be
    clipped.
    
    If the line does not intersect with the rectangular clipping area,
    four None values will be returned as tuple. Otherwise a tuple of the
    clipped line points will be returned in the form (cx1, cy1, cx2, cy2).
    """
  LEFT_, RIGHT_, BOTTOM_, TOP_ = 1, 2, 4, 8
  k1 = k2 = 0
          
  def _getclip(xa, ya):
    p = 0
    if xa < left:
      p |= LEFT_
    elif xa > right:
      p |= RIGHT_
    if ya < bottom:
      p |= BOTTOM_
    elif ya > top:
      p |= TOP_
    return p
                                
  k1 = _getclip(x1, y1)
  k2 = _getclip(x2, y2)
  while (k1 | k2) != 0:
    if (k1 & k2) != 0:
      return None, None, None, None
    opt = k1
    if k1 == 0:
      opt = k2
    if opt & TOP_:
      x = x1 + (x2 - x1) * (1.0*(top - y1)) / (y2 - y1)
      y = top
    elif opt & BOTTOM_:
      x = x1 + (x2 - x1) * (1.0*(bottom - y1)) / (y2 - y1)
      y = bottom
    elif opt & RIGHT_:
      y = y1 + (y2 - y1) * (1.0*(right - x1)) / (x2 - x1)
      x = right
    elif opt & LEFT_:
      y = y1 + (y2 - y1) * (1.0*(left - x1)) / (x2 - x1)
      x = left

    if opt == k1:
      x1 = x
      y1 = y
      k1 = _getclip(x1, y1)
    else:
      x2 = x
      y2 = y
      k2 = _getclip(x2, y2)

  return x1, y1, x2, y2

def get_intersect_poly(sq_i, sq_j, itri):
  left = sqmesh_min_coord[0] + sq_j*sqmesh_step
  right = sqmesh_min_coord[0] + (sq_j + 1)*sqmesh_step
  bottom = sqmesh_min_coord[1] + sq_i*sqmesh_step
  top = sqmesh_min_coord[1] + (sq_i + 1)*sqmesh_step

  #Triangle's nodes clockwise ordering: first node is the one with the min x coordinate
  ifir = 0
  for inode in range(1, 3):
    if tnodes[int(tricells_inodes[itri][inode])][0] < tnodes[int(tricells_inodes[itri][ifir])][0]:
      ifir = inode
  isec = -1
  clkw_ang = -math.pi
  for inode in range(3):
    if (inode != ifir):
      cur_ang = math.atan2(tnodes[int(tricells_inodes[itri][inode])][1] - tnodes[int(tricells_inodes[itri][ifir])][1], tnodes[int(tricells_inodes[itri][inode])][0] - tnodes[int(tricells_inodes[itri][ifir])][0])
      if cur_ang > clkw_ang:
        clkw_ang = cur_ang
        isec = inode

  inodes_clkw = [ifir, isec, 0]
  for inode in range(3):
    if (inode != ifir) and (inode != isec):
      inodes_clkw[2] = inode
      break

  for inode in range(3):
    inodes_clkw[inode] = int(tricells_inodes[itri][inodes_clkw[inode]])

  nclipped = 0
  seg_pts = [[[0.0 for x in range(2)] for x in range(2)] for x in range(3)]
  for iseg in range(3):
    inode1 = inodes_clkw[iseg]
    inode2 = inodes_clkw[(iseg + 1)%3]
    x1, y1, x2, y2 = cohensutherland(left, top, right, bottom, tnodes[inode1][0], tnodes[inode1][1], tnodes[inode2][0], tnodes[inode2][1])
    if x1 != None:
      seg_pts[nclipped][0][0] = x1
      seg_pts[nclipped][0][1] = y1
      seg_pts[nclipped][1][0] = x2
      seg_pts[nclipped][1][1] = y2
      nclipped += 1

  if nclipped == 0:
    return [[]]

  poly_nodes = [[0.0 for x in range(2)] for x in range(7)]
  poly_nodes[0] = seg_pts[0][0]
  inext_seg = 0
  npolynodes = 1;
  sides_cmp = [left, top, right, bottom]
  sq_nodes = [[left, top], [right, top], [right, bottom], [left, bottom]]

  if seg_pts[0][0] == seg_pts[nclipped - 1][1]:
    isq_side_start = -1
    isq_side_stop = -1
  else:
    node_cmp = [poly_nodes[0][0], poly_nodes[0][1], poly_nodes[0][0], poly_nodes[0][1]]
    for iside in range(4):
      if node_cmp[iside] == sides_cmp[iside]:
        isq_side_start = iside
        break
    for iside in range(4):
      if node_cmp[iside] == sides_cmp[iside]:
        isq_side_stop = iside

  while(1):
    if (inext_seg != nclipped) and (poly_nodes[npolynodes - 1] == seg_pts[inext_seg][0]):
      if (isq_side_stop == -1) and (inext_seg == nclipped - 1):
        break
      else:
        poly_nodes[npolynodes] = seg_pts[inext_seg][1]
        npolynodes += 1
        inext_seg += 1
        continue

    node_cmp = [poly_nodes[npolynodes - 1][0], poly_nodes[npolynodes - 1][1], poly_nodes[npolynodes - 1][0], poly_nodes[npolynodes - 1][1]]
    icurside = -1
    if isq_side_start != -1:
      for i in range(4):
        iside = (isq_side_start + i)%4
        if node_cmp[iside] == sides_cmp[iside]:
          icurside = iside
    else:
      for iside in range(4):
        if node_cmp[iside] == sides_cmp[iside]:
          icurside = iside
    isq_side_start = icurside

    if icurside == isq_side_stop:
      if inext_seg < nclipped:
        raise RuntimeError("Completed the intersection polygon before tracing all the clipped segments!")
      break

    if inext_seg < nclipped:
      next_node_cmp = [seg_pts[inext_seg][0][0], seg_pts[inext_seg][0][1], seg_pts[inext_seg][0][0], seg_pts[inext_seg][0][1]]
      if next_node_cmp[icurside] == sides_cmp[icurside]:
        poly_nodes[npolynodes] = seg_pts[inext_seg][0]
        npolynodes += 1
        continue

    poly_nodes[npolynodes] = sq_nodes[icurside]
    npolynodes += 1
      
  poly_nodes = poly_nodes[0:npolynodes]
  return poly_nodes

def get_poly_area(polynodes):
  nnodes = len(polynodes)
  poly_area = 0.0
  for itri in range(nnodes - 2):
    inodes = [0, itri + 1, itri + 2]
    tri_nodes = [[0.0 for x in range(2)] for x in range(3)]
    for i in range(3):
      tri_nodes[i] = polynodes[inodes[i]]
    tri_area = 0.5*abs((tri_nodes[0][0] - tri_nodes[2][0])*(tri_nodes[1][1] - tri_nodes[0][1]) - (tri_nodes[0][0] - tri_nodes[1][0])*(tri_nodes[2][1] - tri_nodes[0][1]))
    poly_area += tri_area
  return poly_area

def is_inside_tri(itri, sq_i, sq_j):
  tri_nodes = [[0.0 for x in range(2)] for x in range(3)]
  for inode in range(3):
    tri_nodes[inode] = tnodes[int(tricells_inodes[itri][inode])]
  sq_nodes = [[0.0 for x in range(2)] for x in range(4)]
  for i in range(2):
    for j in range(2):
      sq_nodes[2*i + j] = [sqmesh_min_coord[0] + (sq_j + j)*sqmesh_step, sqmesh_min_coord[1] + (sq_i + i)*sqmesh_step]
  is_inside = True
  for inode in range(4):
    det = (tri_nodes[1][1] - tri_nodes[2][1])*(tri_nodes[0][0] - tri_nodes[2][0]) + (tri_nodes[2][0] - tri_nodes[1][0])*(tri_nodes[0][1] - tri_nodes[2][1])
    s = ((tri_nodes[1][1] - tri_nodes[2][1])*(sq_nodes[inode][0] - tri_nodes[2][0]) + (tri_nodes[2][0] - tri_nodes[1][0])*(sq_nodes[inode][1] - tri_nodes[2][1])) / det
    t = ((tri_nodes[2][1] - tri_nodes[0][1])*(sq_nodes[inode][0] - tri_nodes[2][0]) + (tri_nodes[0][0] - tri_nodes[2][0])*(sq_nodes[inode][1] - tri_nodes[2][1])) / det
    if (s < 0) or (t < 0) or (s + t > 1):
      is_inside = False
      break
  return is_inside

fid = open('area_weights.dat', 'w')
fid.write(repr(ntris) + '\n')
fid.write(repr(sqmesh_min_coord[0]) + ' ' + repr(sqmesh_min_coord[1]) + ' ' + repr(sqmesh_step) + '\n')

for itri in range(ntris):
  print('Processing triangle ' + repr(itri + 1) + '/' + repr(ntris) + '\r'),
  xmin = tnodes[int(tricells_inodes[itri][0])][0]
  xmax = tnodes[int(tricells_inodes[itri][0])][0]
  ymin = tnodes[int(tricells_inodes[itri][0])][1]
  ymax = tnodes[int(tricells_inodes[itri][0])][1]
  for inode in range(1, 3):
    if tnodes[int(tricells_inodes[itri][inode])][0] < xmin:
      xmin = tnodes[int(tricells_inodes[itri][inode])][0]
    elif tnodes[int(tricells_inodes[itri][inode])][0] > xmax:
      xmax = tnodes[int(tricells_inodes[itri][inode])][0]
    if tnodes[int(tricells_inodes[itri][inode])][1] < ymin:
      ymin = tnodes[int(tricells_inodes[itri][inode])][1]
    elif tnodes[int(tricells_inodes[itri][inode])][1] > ymax:
      ymax = tnodes[int(tricells_inodes[itri][inode])][1]
  imin = int(math.floor((ymin - sqmesh_min_coord[1]) / sqmesh_step))
  imax = int(math.floor((ymax - sqmesh_min_coord[1]) / sqmesh_step))
  jmin = int(math.floor((xmin - sqmesh_min_coord[0]) / sqmesh_step))
  jmax = int(math.floor((xmax - sqmesh_min_coord[0]) / sqmesh_step))

  fid.write(repr(itri) + '\t')
  for i in range(imin, imax + 1):
    for j in range(jmin, jmax + 1):
      area_weight = 0.0
      if is_inside_tri(itri, i, j):
        area_weight = 1.0
      else:
        polygon = get_intersect_poly(i, j, itri)
        if polygon != [[]]:
          cur_poly_area = get_poly_area(polygon)
          area_weight = cur_poly_area / pow(sqmesh_step, 2)
      if area_weight > 0.0:
        fid.write(repr(i) + '\t' + repr(j) + '\t')
        fid.write(repr(area_weight) + '\t')
  fid.write('\n')

fid.close
