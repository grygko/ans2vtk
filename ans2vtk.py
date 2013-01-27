#!/usr/bin/env python

import sys
import numpy as np


def writefile(filename,nodes, nodedict, elements, elemdict, cell_type,
              nodedata={}, vecnodedata={}):
    """
    Notes:
        Format is based on "VTK 4.2 file formats" document
        currently supporting VTK cell types numbers 23 and 25 (see p. 10)
    """
    nonr=len(nodes)
    elnr=len(elements)
    print 'Nbr. Nodes: ',nonr
    print 'Nbr. Elements: ', elnr
    f=open(filename,'w')
    f.writelines('# vtk DataFile Version 2.0\nAnsys 2 VTK\nASCII\nDATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %d float\n'%nonr)
    for i in xrange(nonr):
        #f.write('%f %f %f\n'%tuple(nodes[i]))
        f.write('%s %s %s\n'%tuple(nodes[i]))

    # CELLS n size. n - num of cells, s - size of the cell list size
    size_CELLS = sum([len(nnodes) for nnodes in elements]) + elnr
    f.write('\nCELLS %d %d'%(elnr, size_CELLS)) 
    for e in range(elnr):
        numPoints = len(elements[e]) # num of nodes defining the element
        f.write("\n{0:d}".format(numPoints))
        for node in elements[e]:
            f.write(" {0:d}".format(nodedict[node]))

    #    f.write('8 %d %d %d %d %d %d %d %d\n'%tuple([nodedict[x]-1 for x in elements[n]]))
        #f.write('8 %d %d %d %d %d %d %d %d\n'%tuple([nodedict[i] for i in elements[n]]))

    # CELL_TYPES n
    #cell_types = {4: '9', 8 : '23', 20: '25'} # nodes num : VTK cell type
    f.write('\nCELL_TYPES %d\n'%elnr)
    for e in xrange(elnr):
        f.write('{0}\n'.format(cell_type[e]))

    # POINT_DATA n; n - number of points
    if len(nodedata) == 0: return 0
    f.write('\nPOINT_DATA %d\n'%nonr) # number of points
    for dat in nodedata:        # ['Nbr', 'Sin', 'Random']
        # SCALARS dataName dataType numComp (p.5)
        # LOOKUP_TABLE tableName
        f.write('\nSCALARS {0} float 1\nLOOKUP_Table default\n'.format(dat))
        for k in xrange(nonr):
            f.write('%g \n'%nodedata[dat][k])

    #for idx,case in enumerate(nodedata):
    for dat in vecnodedata:         # ['VecF1', '...']
        # VECTORS dataName dataType
        f.write('VECTORS {0} float\n'.format(dat))
        for j in xrange(nonr):
            f.write('{0[0]} {0[1]} {0[2]} \n'.format(vecnodedata[dat][j]))
    f.close()

def read_dot_node(filename, node_set):
    """Return node data from .NODE file written by ANSYS NWRITE command

    Each line in the .NODE file contains node number followed by x,y,z coordinates
    """
    with open(filename, 'rb') as f:
        ncount = 0
        node_dict = {}         # new : old numbers
        node_coord = []        # x, y, z coordinates
        for line in f:
            # get rid of "\n", empty spaces and split into a list
            entries = line.splitlines()[0].split()
            nnum = entries[0]  # node number
            if nnum in node_set:
                node_dict[nnum] = ncount
                node_coord.append(entries[1:4])
                ncount += 1
    return node_dict, node_coord

def read_dot_elem(filename, midnodes=False):
    """Return element data from .ELEM file written by ANSYS EWRITE command 

    Args:
        midnodes (bool): read mid-side nodes

    Description:
    The data description of each record is (EWRITE help): 
    I, J, K, L, M, N, O, P, MAT, TYPE, REAL, SECNUM, ESYS, IEL  -- nnodes <= 8
    Q, R, S, T, U, V, W, X, Y, Z, A, B                          -- nnodes > 8

    Returns:
        8 1 2 3 4 5 6 7 8
    """
    vtk_cell_types = {'VTK_QUAD': '9', 
                     'VTK_HEXAHEDRON': '12', 
                     'VTK_QUADRATIC_QUAD': '23', 
                     'VTK_QUADRATIC_HEXAHEDRON': '25'}

    with open(filename, 'rb') as f:
        ecount = 0
        elem_dict = {}         # new : old numbers
        elem_connect = []      # elements connectivity
        cell_type = {}
        for line in f:
            # get rid of "\n", empty spaces and split into a list
            entries = line.splitlines()[0].split()

            # record or not mid-side nodes without knowing if elem is 2d or 3d
            if not midnodes: # XXX stupid to check in the loop
                try:
                    enum = int(entries[13]) # el.num is at pos 13
                    elem_dict[enum] = ecount 
                    elem_connect.append(entries[:4])
                    cell_type[ecount] = vtk_cell_types['VTK_QUAD']
                    mnop_entries = entries[4:8]          # save M,N,O,P for later
                    ecount += 1
                except IndexError:                       # nodes > 8
                    cell_type[ecount-1] = vtk_cell_types['VTK_HEXAHEDRON']
                    elem_connect[-1] += mnop_entries
            else:
                try:
                    enum = int(entries[13]) # el.num is at pos 13
                    elem_dict[enum] = ecount 
                    elem_connect.append(entries[:8])
                    cell_type[ecount] = vtk_cell_types['VTK_QUADRATIC_QUAD']
                    ecount += 1
                except IndexError:                          # nodes > 8
                    cell_type[ecount-1] = vtk_cell_types[
                            'VTK_QUADRATIC_HEXAHEDRON']
                    elem_connect[-1] += entries

    node_set = [] # set of nodes contained in elem_connect
    for e in elem_connect:
        node_set += e           # 1-dim list
    node_set = set(node_set)
    return elem_dict, elem_connect, cell_type, node_set

def augment_node_data(swap_nodedict, part_data, pkey='NODE', 
        aug_dict={'rest': 0}):
    """Return augmented sorted list containing data for all nodes in node_dict

    Args:
        swap_nodedict (dict): dictionary with `new : old` node numbers
        part_data (dict): dictionary with lists partial data (i.e., for some nodes)
        pkey: key from `part_data` identifiyng a list ...
        aug_dict (dict): contains `key : augval` pairs i.e., to augment data 
        list `part_data[key]` with `augval`. Key `rest` to augment all unspecified
        keys
    """
    res_dict = {}
    data_keys = part_data.keys()
    data_keys.remove(pkey)

    for k in data_keys:
        res_dict[k] = []

    new_nodes = sorted(swap_nodedict.keys()) # int
    part_nodes = part_data[pkey]             # float

    # convert to the same type 
    obj_type = type(swap_nodedict.values()[0]) # str
    if  type(part_nodes[0]) != obj_type:
        part_nodes = [eval("obj_type(int(item))") for item in part_nodes]

    for k in data_keys:
        for node in new_nodes:
            old_node = swap_nodedict[node]
            if old_node not in part_nodes:
                res_dict[k].append(aug_dict['rest']) # TODO
            else:
                idx = part_nodes.index(old_node)
                res_dict[k].append(part_data[k][idx])

    return res_dict

if __name__=="__main__":
    
    nodefile, elfile = ('data2d/submodsm6_v01_struct.node', 
                        'data2d/submodsm6_v01_struct.elem') # 8 node
    nodefile, elfile = ('data3d/submodsm6_v01_therm.node', 
                        'data3d/submodsm6_v01_therm.elem') # 20 node

    nodes = []
    elements = []
    nodedict = {}
    elemdict = {}
    nodedata = {}
    vecnodedata = {}

    #read_nodes(nodefile,nodes,nodedict)
    #read_elements(elfile,elements,elemdict)
    elemdict, elements, cell_type, node_set = read_dot_elem(elfile, midnodes=False)
    nodedict, nodes = read_dot_node(nodefile, node_set)

    # some arbitrary data
    from numpy.random import random
    nodedata['Random'] = random(len(nodes))
    nodedata['Sin'] = np.sin(np.linspace(0,2*np.pi,len(nodes)))
    backlist = dict((v,k) for k,v in nodedict.iteritems()) # swap keys/values
    nodedata['Nbr'] = [int(backlist[i]) for i in range(len(nodes))]

    # vecdata
    vecnodedata['VecF1'] = zip(random(len(nodes)),np.zeros(len(nodes)), random(len(nodes)))

    ###########
    # REAL DATA (available not for all nodes)
    ###########
    sys.path.append('/usr2/kravchen/Documents/codes/pylib/plotCSV/')
    import csv_utils
    with open('data3d/test123.dat', 'rb') as f:
        data_dict = csv_utils.dict_from_csv(f)
    print data_dict.keys()

    aug_nodedata = augment_node_data(backlist, data_dict,
            aug_dict={'rest':0})
    nodedata['eps_a'] = aug_nodedata['eps_a']
    vecnodedata['eps'] = zip(aug_nodedata['eps_ax'],
                             aug_nodedata['eps_ay'],
                             aug_nodedata['eps_az'])


    # #######
    # EXPLORE
    # #######
    print '\nnode coordinates'
    for n in nodes[:5]:
        print n

    print '\nnode numbers: old : new'
    keys = sorted(nodedict.keys())
    for k in keys[:5]:
        print k, ':', nodedict[k]

    print '\nelements connectivity: old node numbers'
    for e in elements[:5]:
        print e

    print '\nelement numbers: old : new'
    keys = sorted(elemdict.keys())
    for k in keys[:5]:
        print k, ':', elemdict[k]

    print '\nswapped node keys:values new:old'
    keys_back = sorted(backlist.keys())
    for k in keys_back[:5]:
        print k, ':', backlist[k]

    print '\nvector field data'
    for v in vecnodedata['eps'][:5]:
        print v
    
    writefile('geo.vtk',nodes, nodedict, elements, elemdict, cell_type, nodedata, vecnodedata)
    print 'finished'

