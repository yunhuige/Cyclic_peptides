#!/usr/bin/env python


import os, sys, glob, commands

class Mol2File:

    def __init__(self, infile):
        """Create a mol2 file object. """

        fin = open(infile,'r')
        self.rawchunks = ['@'+s for s in fin.read().split('@')]
        self.rawchunks.pop(0) # this first chunk is blank
        fin.close()

        self.chunks = {}  # a dictionary with the @<TRIPOS> keyword as the key
 			  # and chunks of text (list of lines) as the values
        self.chunkkeys = [] 
        for rawchunk in self.rawchunks:
            fields = rawchunk.split('\n') 
            header = fields.pop(0)
            keyword = header.split('>')[1] 
            self.chunks[keyword] = fields
            self.chunkkeys.append(keyword)

        # parse ATOM chunks into indiviual Mol2Atom objects
        self.atoms = [] # list of ordered Mol2Atom objects
        self.extract_atoms()

        self.translate = {}  # an atom index translation table, to be used later



    def extract_atoms(self):
        """Extract Mol2Atom atom objects from the ATOM chunk"""

        # lines should be: atom_id atom_name x y z atom_type [subst_id  [subst_name [charge [status_bit]]]]
        sample = '      1 N1          3.2160    2.9200    8.9280 NT        1 LYS     -0.802484'
        
        for i in range(len(self.chunks['ATOM'])):
            line = self.chunks['ATOM'][i]
            if len(line) > 8:
                fields = line.split()
                atom_id = int(fields[0])
                self.atoms.append( Mol2Atom() )
                index =  len(self.atoms)-1
                self.atoms[index].set_atom_id(atom_id)
                self.atoms[index].set_atom_name(fields[1])
                self.atoms[index].set_xyz( float(fields[2]), float(fields[3]), float(fields[4]) )
                self.atoms[index].set_atom_type(fields[5])
                if len(fields) > 6:
                    self.atoms[index].set_subst_id( int(fields[6]) )
                if len(fields) > 7:
                    self.atoms[index].set_subst_name( fields[7] )
                if len(fields) > 8:
                    self.atoms[index].set_charge( float(fields[8]) )
                if len(fields) > 9:
                    self.atoms[index].set_status_bit( fields[9] )

    def rewrite_atom_chunk(self):
        """Rewrite a Mol2 ATOM chunk from the information in the Mol2Atom objects"""

        # lines should be: atom_id atom_name x y z atom_type [subst_id  [subst_name [charge [status_bit]]]]
        sample = """
      1 N1          3.2160    2.9200    8.9280 NT        1 LYS     -0.802484
     35 H05        -2.9550    1.4110   -6.9120 H1        1 PHE      0.000000 
        """


        self.chunks['ATOM'] = []   # append lines to this list

        for i in range(len(self.atoms)):
            line = ''
            line += '%7d '%self.atoms[i].atom_id
            line += '%-8s '%self.atoms[i].atom_name
            line += '%9.4f %9.4f %9.4f '%(self.atoms[i].x, self.atoms[i].y, self.atoms[i].z)
            line += '%-7s '%self.atoms[i].atom_type
            if self.atoms[i].subst_id != None:
                line += '%3d '%self.atoms[i].subst_id
            else:
                line += '    '
            if self.atoms[i].subst_name != None:
                line += '%-6s '%self.atoms[i].subst_name
            else:
                line += '      '
            if self.atoms[i].charge != None:
                line += '%10.6f '%self.atoms[i].charge
            else:
                line += '          '

            self.chunks['ATOM'].append( line )



    def renumber(self):
        """Renumber the atom and bond numbering"""

        self.renumber_atoms()  # should fill self.translate{} needed for bond renumbering
        self.renumber_bonds()  

    def renumber_atoms(self):
        """Renumber the atom chunk of lines and atom objects.
        Also build a translation table that"""

        ### change the ATOM numbering
        sample = '      1 N1          3.2160    2.9200    8.9280 NT        1 LYS     -0.802484'

        # build a translation table
        new_indices = []
        old_indices = []
        for i in range(len(self.atoms)):
            new_indices.append( i+1 )
            old_indices.append( self.atoms[i].atom_id )
            self.atoms[i].atom_id = i+1


        # reorder the ATOM chunks
        self.chunks['ATOM'] = []
        for i in range(len(new_indices)):
            self.chunks['ATOM'].append( repr(self.atoms[i]) )
            self.translate[old_indices[i]] = new_indices[i] 

        # print the reordered cunks (for sanity) 
        for i in  range(len(self.chunks['ATOM'])):
            print self.chunks['ATOM'][i]

        print 'translation table', self.translate



    def renumber_bonds(self):
        """change the bond numbering"""

        sample = '     1    1    2 1'   # bond_id origin_atom_id target_atom_id bond_type
        bondnum = 1 
        for i in range(len(self.chunks['BOND'])):
            fields = self.chunks['BOND'][i].split()
            if len(fields) == 4:
                origin_atom_id = int(fields[1])
                target_atom_id = int(fields[2])
                bond_type = fields[3] 
                self.chunks['BOND'][i] = '%-6d %-4d %-4d %s'%(bondnum, self.translate[origin_atom_id], self.translate[target_atom_id], bond_type)
                bondnum += 1

        # print out the newly renumbered bonds
        for i in  range(len(self.chunks['BOND'])):
             print self.chunks['BOND'][i]

    def net_charge(self):
        """Return the net charge on the molecule."""

        charge = 0.0
        for i in range(len(self.atoms)):
            charge += self.atoms[i].charge
        return charge          

    def get_atom_index_by_name(self, name):

        for i in range(len(self.atoms)):
            if self.atoms[i].atom_name == name:
                return i
        return None

    def get_atom_id_by_name(self, name):

        for i in range(len(self.atoms)):
            if self.atoms[i].atom_name == name:
                return self.atoms[i].atom_id
        return None


    def edit_charge(self, name, charge):
        """Change the charge for an atom specified by atom name."""

        index = self.get_atom_index_by_name(name)
        print 'Changing %-7s charge from %10.4f to %10.4f'%(self.atoms[index].atom_name, self.atoms[index].charge, charge)
        self.atoms[index].charge = charge

    def make_net_charge_integer(self, target_charge, exclude_names=[]):
        """Shift all the atomic charges to achieve a net integer charge over all (non-excluded) atoms."""

        included_indices = []
        excluded_indices = []
        for i in range(len(self.atoms)):
            if exclude_names.count(self.atoms[i].atom_name) == 0:
                included_indices.append( i )
            else:
                excluded_indices.append( i )

        print 'included_indices', included_indices
        print 'excluded_indices', excluded_indices

        total_included_charge =  0.0
        for i in included_indices:
            total_included_charge += self.atoms[i].charge
        print 'total_included_charge', total_included_charge

        # distrubute a correction across all included atoms
        distribute_charge = (target_charge - total_included_charge)/float(len(included_indices))
        for i in included_indices:
            self.atoms[i].charge += distribute_charge


    def write(self, outfile):
        """Write *.mol2 file."""

        # make sure any updates to the ATOM chunk get written!
        self.rewrite_atom_chunk()
 
        fout = open(outfile,'w')
        fout.write(repr(self))
        fout.close()

        print '...Done.'

        print '*** Make sure the @<TRIPOS>MOLECULE header is correct -- MUST CHANGE THIS BY HAND!! ****    '
    
    def __repr__(self):
        """Return a representation of the Mol2File object, in the case, the text content of the the file."""

        outstr = ''
        for key in self.chunkkeys:
            outstr += '@<TRIPOS>%s\n'%key
            for i in range(len(self.chunks[key])):
                if self.chunks[key][i].strip() != '':
                    outstr += self.chunks[key][i]+'\n' 
        return outstr



class Mol2Atom:

    def __init__(self):
        """Initialize a Mol2Atom object."""

        self.atom_id = None
        self.atom_name = None
        self.x = None
        self.y = None
        self.z = None
        self.atom_type = None
        self.subst_id = None
        self.subst_name = None
        self.charge = None
        self.status_bit = None

    def set_atom_id(self, value):
        self.atom_id = value

    def set_atom_name(self, value):
        self.atom_name = value

    def set_xyz(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def set_atom_type(self, value):
        self.atom_type = value

    def set_subst_id(self, value):
        self.subst_id = value

    def set_subst_name(self, value):
        self.subst_name = value

    def set_charge(self, value):
        self.charge = value

    def set_status_bit(self, value):
        self.status_bit = value

    def __repr__(self):
        """Return a representation of the atom."""

        # lines should be: atom_id atom_name x y z atom_type [subst_id  [subst_name [charge [status_bit]]]]
        sample = '      1 N1          3.2160    2.9200    8.9280 NT        1 LYS     -0.802484'

        outstr = '%-7d %-7s  %9.4f %9.4f %9.4f %-7s '%(self.atom_id,
                                                                   self.atom_name,
                                                                   self.x, self.y, self.z,
                                                                   self.atom_type)
        if self.subst_id != None:
            outstr += '%3d '%self.subst_id
        else:
            return outstr

        if self.subst_name != None:
            outstr += '%-6s '%self.subst_name
        else:
            return outstr

        if self.charge != None:
            outstr += '%10f'%self.charge
        else:
            return outstr
          
        return outstr


"""
new file
     24 C11        4.0340    9.1000    2.8420      CA   1 EDA      -0.063700
     41 H10         0.8390   11.5540    2.3610 HA         1 EDA       0.177400
     42 H11        -0.6500   -2.6670   -0.0860 H        1 EDA       0.350700
     31 HG3         1.0940    1.3400    0.3390 HC        1 EDA       0.049400
old file
     24 C11         4.0340    9.1000    2.8420 CA        1 EDA      -0.0637

"""

