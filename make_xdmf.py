from __future__ import division
import h5py as h5
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom

import os
import parameters as params

ele_coords, conductivity, path, sbspt = params.default_run('anisotropic')
#filepath = os.path.join(path, 'traub_pots_v_t.h5')
filepath = 'traub_pots_v_t.h5'
f = h5.File(filepath, 'r')
num_ele, dims = f['ele_pos'].value.shape
f.close()
start_t, end_t  = 2750, 3750

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def write_xdmf(filename='output.xdmf', num_ele=num_ele, dims=dims, start_t=start_t, end_t=end_t):
    """Write an XDMF file to accompany the HDF5 data file.
    """
    # Create the basic structure
    time_dims = str(end_t - start_t)
    root = Element('Xdmf')
    root.set('version', '1.0')
    domain = SubElement(root, 'Domain')

    time_grid = SubElement(domain, 'Grid', {'GridType': 'Collection',
                                            'CollectionType': 'Temporal',
                                            'Name': 'Potentials'})

    topology = SubElement(time_grid, 'Topology',
                          {'TopologyType': "Polyvertex",
                           'Name': "Main Topology",
                           'NodesPerElement': str(num_ele)})


    geometry = SubElement(time_grid, 'Geometry',
                          {'GeometryType': "X_Y_Z",
                           'Name': "Main Geometry"})

    ele_pos = SubElement(geometry, 'DataItem',
                         {'Dimensions': str(num_ele),
                          'Format': "HDF"})
    ele_pos.text = filepath + ":/x"
    ele_pos_y = SubElement(geometry, 'DataItem',
                         {'Dimensions': str(num_ele),
                          'Format': "HDF"})
    ele_pos_y.text = filepath + ":/y"
    ele_pos_z = SubElement(geometry, 'DataItem',
                         {'Dimensions': str(num_ele),
                          'Format': "HDF"})
    ele_pos_z.text = filepath + ":/z"



    for step in range(start_t, end_t):
        grid = SubElement(time_grid, 'Grid',
                          {'GridType': 'Uniform', 'Name': 'Structured Grid'})
        topology = SubElement(grid, 'Topology',
                              {'Reference': '//Topology[@Name="Main Topology"]'})
        geometry = SubElement(grid, 'Geometry',
                              {"Reference": '//Geometry[@Name="Main Geometry"]'})
        time_pt = SubElement(grid, 'Time',
                             {"Value": str(step)})
        density = SubElement(grid, "Attribute",
                             {"Name": 'phi',
                              "AttributeType": "Scalar",
                              "Center": "Node"})
        act_vals = SubElement(density, "DataItem",
                              {"Dimensions": str(num_ele),
                               "Format": "HDF"})
        act_vals.text = filepath + ':/' + str(step)

    with open(filename, 'w') as f:
        f.write(prettify(root))

write_xdmf()
