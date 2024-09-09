#!/usr/bin/env python
# Module for objects

import os
import re
import numpy as np
import rdkit.Chem as rdchem; from rdkit.Chem import rdDetermineBonds, rdmolops

class file_reader():
    lines=None
    lines_no=None
    comtag=None
    # Init function checks if file exists reads in the line from file
    def __init__(self, inp_fi):
        fi=os.path.realpath(inp_fi)
        if not os.path.isfile(fi):
            raise Exception(f"Not an existing file: \"{fi}\"")
        
        with open(fi, 'r') as rd:
            self.lines=[l.strip('\n') for l in rd.readlines()]
        if len(self.lines)<1:
            raise Exception(f"Empty File: \"{fi}\"")
    @property
    def lines_no(self):
        if type(self.lines)==type(None):
            raise Exception(f"Lines are not initialized")
        if type(self.lines)!=list:
            raise Exception(f"Lines are not in format list but {type(self.lines)}")
        if len(self.lines)<1:
            raise Exception(f"Lines is a list with no entries")

        self._lines_no=[ (i,l) for i,l in enumerate(self.lines) ]
        return self._lines_no
    @property
    def lines_nocom(self):
        if type(self.lines)==type(None):
            raise Exception(f"Lines are not initialized")
        if type(self.lines)!=list:
            raise Exception(f"Lines are not in format list but {type(self.lines)}")
        if len(self.lines)<1:
            raise Exception(f"Lines is a list with no entries")

        if type(self.comtag)==type(None):
            raise Exception(f"Comment tag is not defined")
        self._lines_nocom=self.lines
        for comtag in self.comtag:
            self._lines_nocom=[ l.split(comtag)[0] for i,l in enumerate(self._lines_nocom) ]

        return self._lines_nocom
    @property
    def lines_num_nocom(self):
        return [ (i,l) for i,l in enumerate(self.lines_nocom) ]
    def get_line_no(self, match=None, fullmatch=False, single=False):
        # Check that type is a string
        if type(match)!=None:
            if type(match)!=str:
                raise Exception()
        else:
            raise Exception()
        if not fullmatch:
            results=[ v[0] for v in self.lines_num_nocom if re.match( match, v[1] ) ]
        else:
            results=[ v[0] for v in self.lines_num_nocom if re.fullmatch( match, v[1] ) ]
        if single:
            if len(results)!=1:
                raise Exception()
            return results[0]
        else:
            return results
    def get_line_item(self,l, tot_leng=None, get_item=None, data_type=None, sep=None):
        if type(l)==str:
            pass
        elif type(l)==int:
            if l>=len(self.lines) or l<0:
                raise Exception(f"Requested line number \'{l}\' is not in allowed range ([0,{len(self.lines)}]).")
            else:
                l=self.lines_nocom[l]
        else:
            raise Exception()

        if type(sep)==type(None):
            ar=l.split()
        else:
            raise Exception(f"Not implemented")

        if type(tot_leng)!=type(None):
            if type(tot_leng)!=int:
                raise Exception(f"Wrong type for option \'tot_leng\':\n {tot_leng} ({type(tot_leng)})")
            else:
                if len(ar)!=tot_leng:
                    raise Exception(f"Expected line of length {tot_leng} but got {len(ar)}:\n{l}, {ar}")

        return_single=False
        if type(get_item)==int:
            get_item=[get_item]
            return_single=True
        elif type(get_item)!=list:
            raise Exception(f"Type of number of item to get is wrong, should be integer:\n{get_item} ({type(get_item)})")

        if any( [x>len(ar) for x in get_item ] ):
            raise Exception(f"Number of get item on this line is {get_item} longer than elements on line:\n{l}")
        if any( [x<0 for x in get_item] ):
            raise Exception()
        else:
            vals=[ ar[x] for x in get_item ]

        if type(data_type)!=type(None):
            try:
                vals=[ data_type(x) for x in vals ]
            except Exception as e:
                raise Exception(f"Could not convert value {val} for type conversion {data_type}")
        if return_single:
            vals=vals[0]
        return vals
    def read_data(self,start=None, end=None, data_type=None, object_type=None, line_items=None):
        if start<0 or start>=len(self.lines):
            raise Exception()
        if end<start+1 or end>len(self.lines):
            raise Exception()

        data=[]
        the_lines=self.lines_nocom[ start:end ]
        for l in the_lines:
            ar=l.split()

            #Check length
            if type(line_items)!=type(None):
                if type(line_items)!=int:
                    raise Exception()
                else:
                    if len(ar)!=line_items:
                        raise Exception()
            
            if type(data_type)!=type(None):
                try:
                    ar=[ data_type(x) for x in ar ]
                except Exception as e:
                    raise Exception(f"Could not read {ar} with {data_type}:\n{e}")
            data.append(ar)
        if type(object_type)==type(None):
            return data
        elif object_type=='np':
            try:
                data=np.array(data)
            except Exception as e:
                raise Exception(f"Could not convert to numpy array: {e}")
            return data
        else:
            raise Exception()
        
class read_map_file(file_reader):
    comtag=['!']
    dic={}
    def __init__(self, inp_fi):
        file_reader.__init__(self,inp_fi)
        self.extract_data()
    def extract_data(self):
        # Get the number of points
        point_lino=self.get_line_no(match='POINTS', single=True)
        point_line=self.lines_nocom[point_lino]
        no_pts=self.get_line_item( point_line, tot_leng=2, get_item=1, data_type=int)
        
        # Get the data
        no_data_start=self.get_line_no(match='BEGIN DATA', single=True)
        no_data_end=self.get_line_no(match='END DATA', single=True)
        #-- Get points
        data=self.read_data(start=no_data_start+1, end=no_data_start+1+no_pts, data_type=float, object_type='np', line_items=4)
        points=data[:,0:3]
        values=data[:,3]
        #-- Get the triangles
        triang_lino=self.get_line_no(match='TRIANGLES', single=True)
        no_triang=self.get_line_item( triang_lino, tot_leng=2, get_item=1, data_type=int)
        triangles=self.read_data(start=no_data_start+1+no_pts, end=no_data_start+1+no_pts+no_triang, 
                data_type=float, object_type='np', line_items=3)
        # Go to python index convention
        triangles=triangles-1


        self.dic.update({'no_pts':no_pts,
            'points': points,
            'values': values,
            'triangles': triangles,
        })

class geometry():
    ifi=None
    coordinates=None
    atom_types=None
    multiplicity=None
    multiplicity_ion=None
    charge=None
    adjacency_matrix=None
    def __init__(self,coordinates=None, atom_types=None):
        if type(coordinates) != type(None):
            self.coordinates=coordinates
        if type(atom_types) != type(None) :
            self.atom_types=atom_types
    def connec(self):
        #from rdkit.Chem import AllChem,rdDetermineBonds
        if type(self.ifi) != type(None):
            raw_mol = rdchem.rdmolfiles.MolFromXYZFile(self.ifi)
            rdmol = rdchem.Mol(raw_mol)
        else:
            mol_block=f"{len(self.coordinates)}\n\n"
            for ty,coor in zip( self.atom_types, self.coordinates):
                bohr_to_angstrom=0.529177249
                coor=[f"{float(x)*bohr_to_angstrom}" for x in coor]
                mol_block+=f"{ty} {' '.join(coor)}\n"
            rdmol=rdchem.MolFromXYZBlock(mol_block)
            #mol_block = rdchem.MolToMolBlock(rdmol)
        rdDetermineBonds.DetermineBonds(rdmol, charge=0)
        self.adjacency_matrix=np.array(rdmolops.GetAdjacencyMatrix(rdmol, useBO=True))
        #m_io.out().warn(f"Implement check for adjacency matrix, ie total valence electrons, number of bonds per species")
        #self.adj_mat=rdchem.rdmolops.GetAdjacencyMatrix(rdmol, useBO=True)
    #@classmethod




class read_mom_file(file_reader):
    comtag=['!']
    dic={}
    atom_types=None
    atom_coordinates=None
    def __init__(self,inp_fi):
        file_reader.__init__(self,inp_fi)
        self.extract_data()
    def extract_data(self):
        empty_lines=self.get_line_no(match=r'[ ]*', fullmatch=True)
        empty_lines+=self.get_line_no(match='', fullmatch=True)
        empty_lines=sorted( list(set(empty_lines)) )
        atom_lines=self.get_line_no(match=r'[A-Z][a-z]?[0-9]* ')

        start_end=[]
        atom_entries=[]
    
        # Check if atom lines were recognized correctly
        for i,e in enumerate(empty_lines):
            if i<len(empty_lines)-1:
                # Skip if consecutive empty line
                if empty_lines[i+1] == e+1:
                    pass
                # In case the document does not end after this line
                elif e<len(self.lines)-1:
                    if not e+1 in atom_lines:
                        raise Exception(f"Line {e} is empty, Total lines {len(self.lines)}")
            if i==len(empty_lines)-1:
                if e!=len(self.lines)-1:
                    raise Exception

        # Get the intervals
        atom_entries=[]
        for a in atom_lines:
            check=False
            for e in empty_lines:
                if e>a:
                    atom_entries.append( self.lines_num_nocom[a:e-1])
                    check=True
                    break
            if not check:
                raise Exception()

        # Process the intervals
        atoms=[]
        for atom in atom_entries:
            for i,line_pl_num in enumerate(atom):
                num, line =line_pl_num
                if i==0: # extract the atom data
                    xyz=self.get_line_item(num, tot_leng=8, get_item=[1,2,3], data_type=float)
                    xyz=np.array(xyz)

                    # Get atom type
                    dummy=self.get_line_item(num, tot_leng=8, get_item=4)
                    if dummy!='Type':
                        raise Exception(f"Expected string \'Type\' at position 4:\n{line}")
                    else:
                        at_type=self.get_line_item(num, tot_leng=8, get_item=5)
                        if not re.match(r'[A-Z][a-z]?', at_type):
                            raise Exception()
            atoms.append({'atom_type':at_type, 'xyz':xyz})

        self.atom_types=[ atom['atom_type'] for atom in atoms ]
        self.atom_coordinates=[ atom['xyz'] for atom in atoms ]


class grid():
    input_file=None
    grid_points=None
    values=None
    triangles=None
    statistical_indicators=None
    def __init__(self,inp=None):
        if type(None)==type(inp): # Start with empty object
            pass
        elif str==type(inp): # Read a file
            # Read file according to extension
            self.input_file=inp
            if inp.endswith('.map'):
                info_dict=read_map_file(inp).dic
                self.grid_points=info_dict['points']
                self.values=info_dict['values']
                self.triangles=info_dict['triangles']
            else:
                raise Exception(f"Unkown extension of file \"{inp}\"")
    def gen_statistical_indicators(self):
        def get_rmse(vals):
            avr=np.average(vals)
            return np.sqrt( np.average ( (vals-avr)**2 ) )
        average = np.average(self.values)
        minimum     = np.min(self.values)
        maximum     = np.max(self.values)
        min_dev = minimum-average
        max_dev = maximum-average
        rmse    = get_rmse(self.values)
        pos_dev = [ x for x  in self.values if x>=average ]
        neg_dev = [ x for x in self.values if x<=average ]
        rmse_pos=get_rmse(pos_dev)
        rmse_neg=get_rmse(neg_dev)
        self.statistical_indicators={
            'average':average, 'minimum':minimum, 'maximum':maximum, 
            'min_dev':min_dev, 'max_dev':max_dev, 'rmse':rmse
        }

class multipoles():
    atom_types=None
    atom_coordinates=None
    def __init__(self,inp_fi):
        if type(inp_fi)!=str:
            raise Exception()
        else:
            if inp_fi.endswith('.mom'):
                multipoles_raw=read_mom_file(inp_fi)
                self.atom_types=multipoles_raw.atom_types
                self.atom_coordinates=multipoles_raw.atom_coordinates
            else:
                raise Exception
    def return_geom(self, connectivity=False):
        geom=geometry(atom_types=self.atom_types, coordinates=self.atom_coordinates)
        geom.connec()
        if connectivity:
            geom.connec()
        return geom
