# -*- coding: utf-8 -*-
"""Module for manupulation of PDB files.

by Cristian Escobar
Last update 4/9/2020
"""


class PdbData():
    """Class that contains PDB data and functions to modify it.

    Parameters
    ----------
    pdb_filename : String
        PDB file name.

    Returns
    -------
        None.

    """

    def __init__(self, pdb_filename):
        """Start the PdbData instance.

        Instance needs the name of a file containing PDB data as string.

        Parameters
        ----------
        pdb_filename : String
            PDB file name.

        Returns
        -------
        None.

        """
        # Class variables
        self.pdb_filename = pdb_filename

        self.pdb_dictionary = {}

        self.chainIDs = []

    def open_pdbfile(self):
        """Open PDB file.

        Method to open a PDB file and creates a dictionary containing
        information about all protein atoms. Also collects information
        about the PDB file.

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        # open file
        pdb_file = open(self.pdb_filename, 'r')
        lines = pdb_file.readlines()
        pdb_file.close()

        pdb_dict = {}   # Dictionary will contain all atom data
        # {atom_num:{atom_property1: '', ...atom_propertyn:''},...}

        # get data
        atom_counter = 0

        for line in lines:

            if line[0:4] == 'ATOM':  # only reads ATOM data

                atom_counter += 1

                atom_dict = {}

                atom_dict['record_name'] = line[0:6]
                atom_dict['serial'] = int(line[6:11])  # %5d
                atom_dict['name'] = line[12:16]
                atom_dict['altLoc'] = line[16:17]
                atom_dict['resName'] = line[17:20]
                atom_dict['chainID'] = line[21:22]
                atom_dict['resSeq'] = int(line[22:26])  # %4d
                atom_dict['iCode'] = line[26:27]
                atom_dict['x'] = float(line[30:38])  # %8.3f
                atom_dict['y'] = float(line[38:46])  # %8.3f
                atom_dict['z'] = float(line[46:54])  # %8.3f
                atom_dict['occupancy'] = float(line[54:60])  # %6.2f
                atom_dict['tempFactor'] = float(line[60:66])  # %6.2f
                atom_dict['segid'] = line[72:76]  # 4 letter long
                atom_dict['element'] = line[76:78]
                atom_dict['charge'] = line[78:80]

                pdb_dict[atom_counter] = atom_dict

                # add chainIDs to list
                if atom_dict['chainID'] not in self.chainIDs:

                    self.chainIDs.append(atom_dict['chainID'])

        self.pdb_dictionary = pdb_dict

    def change_chainID(self, chainID, new_chainID):
        """Modify the chainID value.

        Method modifies the current ChainID for the new value provided.

        Parameters
        ----------
        chainID : String
            ChainID of the chain to modify, eg: 'A'
        new_chainID : String
            New ChainID, eg: 'B'

        Returns
        -------
        None.
        """
        if chainID in self.chainIDs:  # run only if chainID exist

            for atom_num in self.pdb_dictionary:  # loop over all atoms

                chainid = self.pdb_dictionary[atom_num]['chainID']

                if chainid == chainID:

                    self.pdb_dictionary[atom_num]['chainID'] = new_chainID

            # modify chainIDs list
            self.update_chainIDs()

        else:
            print("ChainID not found in the PDB file")

    def change_partial_chainID(self, start_res, end_res, chainID, new_chainID):
        """Change the chainID of a segment of residues.

        Method that changes the ChainID for a section of residues starting
        in start_res and ending in end_res for the new ChainID.

        Parameters
        ----------
        start_res : Integer
            Residue number of the first residue to modify.
        end_res : Integer
            Residue number of the last residue to modify.
        chainID : String
            ChainID of the residues to modify
        new_chainID : String
            New ChainID.

        Returns
        -------
        None.

        """
        for atom_num in self.pdb_dictionary:

            residue_num = self.pdb_dictionary[atom_num]['resSeq']

            chainid = self.pdb_dictionary[atom_num]['chainID']

            condition1 = residue_num >= start_res

            condition2 = residue_num <= end_res

            condition3 = chainid == chainID

            if condition1 and condition2 and condition3:

                self.pdb_dictionary[atom_num]['chainID'] = new_chainID

        # modify chainIDs list
        self.update_chainIDs()

    def update_chainIDs(self):
        """Update the list of chainID.

        Returns
        -------
        None.

        """
        new_chainid_list = []

        for atom_num in self.pdb_dictionary:

            chainid = self.pdb_dictionary[atom_num]['chainID']

            if chainid not in new_chainid_list:

                new_chainid_list.append(chainid)

        self.set_chainIds(new_chainid_list)

    def change_Aa_num(self, chainID, modifier):
        """Change the amino acid number by adding the modifier.

        Method that changes amino acid number of a given chain by a modifier
        factor.

        Parameters
        ----------
        chainID : String
            ChainID of the chain to be modified, eg. 'A'
        modifier : Integer
            Number used to modify the residue number

        Returns
        -------
        None.

        """
        for atom_num in self.pdb_dictionary:

            chainid = self.pdb_dictionary[atom_num]['chainID']

            if chainid == chainID:

                self.pdb_dictionary[atom_num]['resSeq'] += modifier

    def reset_Aa_num(self, chainID):
        """Reset amino acid numbering in a chain to start from 1.

        Parameters
        ----------
        chainID : String
            ChainID of the chain to be modified, eg: 'A'.

        Returns
        -------
        None.

        """
        Aa_counter = 0  # new residue number

        current_Aa = 0  # Number of current residue

        for atom_num in self.pdb_dictionary:

            chainid = self.pdb_dictionary[atom_num]['chainID']

            aa_number = self.pdb_dictionary[atom_num]['resSeq']

            if chainid == chainID:

                if aa_number != current_Aa:

                    current_Aa = aa_number

                    Aa_counter += 1

                    self.pdb_dictionary[atom_num]['resSeq'] = Aa_counter

                else:

                    self.pdb_dictionary[atom_num]['resSeq'] = Aa_counter

    def remove_segID(self):
        """Remove the SEGID from each atom.

        Returns
        -------
        None.

        """
        for atom_num in self.pdb_dictionary:

            self.pdb_dictionary[atom_num]['segid'] = '    '

    def change_segID(self, new_segid):
        """Change the SEGID for a new value.

        Parameters
        ----------
        new_segid : String
            SEGID should be up to 4 characteres long.

        Returns
        -------
        None.

        """
        for atom_num in self.pdb_dictionary:

            self.pdb_dictionary[atom_num]['segid'] = new_segid

    def add_gap(self, chainID, start_res, gap_length):
        """Add a gap in the amino acid sequence.

        Method to add a gap in the amino acid sequence of a given chain.
        The first residue after the gap correponds to start_res. Thus, all
        residues after start_res (including it) will be modified by adding
        gap_length to the residue number. Effectively, the gap is introduced
        after start_res-1.
        Example: if sequence: 1,2,3,4,5 and start_res: 3 and gap_length: 10,
        new seqeunce will be: 1,2,13,14,15.

        Parameters
        ----------
        chainID : String
            ChainID of chain to modify.
        start_res : Integer
            Residue number of residue after which the will be added.
        gap_length : Integer
            Lenght of gap added.

        Returns
        -------
        None.

        """
        for atom_num in self.pdb_dictionary:

            chainid = self.pdb_dictionary[atom_num]['chainID']

            aa_number = self.pdb_dictionary[atom_num]['resSeq']

            if chainid == chainID:

                if aa_number >= start_res:

                    self.pdb_dictionary[atom_num]['resSeq'] += gap_length

    def append_pdb(self, pdb_instance):
        """Append PDB data to the current data.

        Method to add PDB data from another instance of PdbData to the end of
        current data. Atom numbers of the new data will be updated con continue
        with the current numbers. ChainID of data added will be preserved.

        Parameters
        ----------
        pdb_instance : PdbData
            Instance of PdbData that will be added to the end of current
            instance.

        Returns
        -------
        None.

        """
        # set to the last current plus 1
        atom_counter = 1 + len(self.pdb_dictionary)

        for atom_num in pdb_instance.pdb_dictionary:

            atom_dict = pdb_instance.pdb_dictionary[atom_num]

            atom_dict['serial'] = atom_counter  # update to new atom number

            self.pdb_dictionary[atom_counter] = atom_dict

            atom_counter += 1

        # Check chainID if they are the same
        new_chainIDs = pdb_instance.chainIDs

        for chainid in new_chainIDs:

            if chainid in self.chainIDs:

                print('Warning! New chainID already exist.')

        self.update_chainIDs()

    def get_pdb_dictionary(self):
        """Return pdb dictionary data.

        Returns
        -------
        Dictionary.

        """
        return self.pdb_dictionary

    def set_pdb_dictionary(self, new_dictionary):
        """Set pdb_dictionary to a new dictionary.

        Parameters
        ----------
        new_dictionary : Dictionary
            Dictionary containing PDB data.

        Returns
        -------
        None.

        """
        self.pdb_dictionary = new_dictionary

    def get_chainIds(self):
        """Return ChainID list.

        Returns
        -------
        None.

        """
        return self.chainIDs

    def set_chainIds(self, chainid_list):
        """Set new chainID values.

        Parameters
        ----------
        chainid_list : List
            List of new ChainID values.

        Returns
        -------
        None.

        """
        self.chainIDs = chainid_list

    def get_pdb_filename(self):
        """Return PDB file name.

        Returns
        -------
        PDB file name.

        """
        return self.pdb_filename

    def set_pdb_filename(self, file_name):
        """Set PDB file name.

        Parameters
        ----------
        file_name : String
            File name of PDB file.

        Returns
        -------
        None.

        """
        self.pdb_filename = file_name
        self.pdb_filename = file_name

    def join_chains(self):
        """Join all chains into a single chain with continous residue numbers.

        Returns
        -------
        None.

        """
        chainid_list = list(self.chainIDs)

        # Change ChainIDs to a single value
        for chainid in chainid_list:

            self.change_chainID(chainid, 'A')

        # Reset amino acid numbering
        self.reset_Aa_num('A')

    def print_pdb(self, file_name):
        """Print data into a new PDB file.

        Parameters
        ----------
        file_name : String
            PDB file output name should have .pdb extension.

        Returns
        -------
        None.

        """
        out_file = open(file_name, 'w')

        for atom_num in self.pdb_dictionary:

            atom_dict = self.pdb_dictionary[atom_num]

            line = atom_dict['record_name']
            line += '%5d' % (atom_dict['serial'])
            line += ' '
            line += atom_dict['name']
            line += atom_dict['altLoc']
            line += atom_dict['resName']
            line += ' '
            line += atom_dict['chainID']
            line += '%4d' % (atom_dict['resSeq'])
            line += atom_dict['iCode']
            line += '   '
            line += '%8.3f' % (atom_dict['x'])
            line += '%8.3f' % (atom_dict['y'])
            line += '%8.3f' % (atom_dict['z'])
            line += '%6.2f' % (atom_dict['occupancy'])
            line += '%6.2f' % (atom_dict['tempFactor'])
            line += '      '
            line += atom_dict['segid']
            line += atom_dict['element']
            line += atom_dict['charge']

            out_file.write(line + '\n')

        out_file.write('END')

        out_file.close()

        print('DONE!')

    def print_chain(self, chainid, file_name):
        """Print a single chain from the PDB file.

        Parameters
        ----------
        chainid : String
            ChainID of the chain to be printed.
        file_name : String
            PDB file output name should have .pdb extension.

        Returns
        -------
        None.

        """
        if chainid in self.chainIDs:

            out_file = open(file_name, 'w')

            for atom_num in self.pdb_dictionary:

                atom_dict = self.pdb_dictionary[atom_num]

                if atom_dict['chainID'] == chainid:

                    line = atom_dict['record_name']
                    line += '%5d' % (atom_dict['serial'])
                    line += ' '
                    line += atom_dict['name']
                    line += atom_dict['altLoc']
                    line += atom_dict['resName']
                    line += ' '
                    line += atom_dict['chainID']
                    line += '%4d' % (atom_dict['resSeq'])
                    line += atom_dict['iCode']
                    line += '   '
                    line += '%8.3f' % (atom_dict['x'])
                    line += '%8.3f' % (atom_dict['y'])
                    line += '%8.3f' % (atom_dict['z'])
                    line += '%6.2f' % (atom_dict['occupancy'])
                    line += '%6.2f' % (atom_dict['tempFactor'])
                    line += '      '
                    line += atom_dict['segid']
                    line += atom_dict['element']
                    line += atom_dict['charge']

                    out_file.write(line + '\n')

            out_file.write('END')

            out_file.close()

            print('DONE!')

        else:

            print('ChainID not found in the PDB file.')


def combine_pdb(pdb_instance1, start_res1, end_res1, chainid1,
                pdb_instance2, start_res2, end_res2, chainid2,):
    """Combine pdb data in a new PdbData instance.

    Function to combine partial data from two different instances of PdbData
    to a new PdbData instance. Data to combine is defined by start_res, end_res
    and the respective ChainID. Data from instance1 will be added first and
    the instance2.

    Parameters
    ----------
    pdb_instance1 : PdbData instance
        Instance of PdbData number 1.
    start_res1 : Integer
        First residue of data to be  used from instance 1.
    end_res1 : Integer
        Last residue of data to be  used from instance 1.
    chainid1 : String
        ChainID of selected data of instance 1.
    pdb_instance2 : PdbData instance
        Instance of PdbData number 1.
    start_res2 : Integer
        First residue of data to be  used from instance 2.
    end_res2 : Integer
        Last residue of data to be  used from instance 2.
    chainid2 : String
        ChainID of selected data of instance 2.

    Returns
    -------
    New PdbData instance.

    """
    new_pdb_data = PdbData('')
    new_pdb_dict = {}
    chainID_list = []
    atom_counter = 1
    # add new chainIds
    if chainid1 != chainid2:
        chainID_list = [chainid1, chainid2]

    else:
        chainID_list = [chainid1]
        print('Warning! ChainID are the same.')

    # append first pdb
    for atom_num in pdb_instance1.pdb_dictionary:

        atom_dict = pdb_instance1.pdb_dictionary[atom_num]

        resnum = atom_dict['resSeq']
        chainid = atom_dict['chainID']

        if resnum >= start_res1 and resnum <= end_res1 and chainid == chainid1:

            atom_dict['serial'] = atom_counter

            new_pdb_dict[atom_counter] = atom_dict

            atom_counter += 1

    # append second pdb
    for atom_num in pdb_instance2.pdb_dictionary:

        atom_dict = pdb_instance2.pdb_dictionary[atom_num]

        resnum = atom_dict['resSeq']
        chainid = atom_dict['chainID']

        if resnum >= start_res2 and resnum <= end_res2 and chainid == chainid2:

            atom_dict['serial'] = atom_counter

            new_pdb_dict[atom_counter] = atom_dict

            atom_counter += 1

    new_pdb_data.set_pdb_dictionary(new_pdb_dict)
    new_pdb_data.set_chainIds(chainID_list)

    return new_pdb_data


def print_models(pdb_file_list, file_name):
    """Print PDBs to a single file as models.

    Functions takes a list of PDB files and combines them into a single
    PDB file where each PDB data correspond to one model in the colection.

    Parameters
    ----------
    pdb_file_list : List
        List of PDB file names.
    file_name : String
        File name for the output file with .pdb extension.

    Returns
    -------
    None.

    """
    out_file = open(file_name, 'w')

    model_num = 1

    for pdb_file_name in pdb_file_list:

        # open pdb data
        pdb_data = PdbData(pdb_file_name)
        pdb_data.open_pdbfile()
        pdb_dict = pdb_data.get_pdb_dictionary()

        # first line that defines a model
        model_line = 'MODEL '
        model_line += '    '
        model_line += '%4d' % (model_num)

        out_file.write(model_line + '\n')

        # loop over all atoms in a model
        for atom_num in pdb_dict:

            atom_dict = pdb_dict[atom_num]

            line = atom_dict['record_name']
            line += '%5d' % (atom_dict['serial'])
            line += ' '
            line += atom_dict['name']
            line += atom_dict['altLoc']
            line += atom_dict['resName']
            line += ' '
            line += atom_dict['chainID']
            line += '%4d' % (atom_dict['resSeq'])
            line += atom_dict['iCode']
            line += '   '
            line += '%8.3f' % (atom_dict['x'])
            line += '%8.3f' % (atom_dict['y'])
            line += '%8.3f' % (atom_dict['z'])
            line += '%6.2f' % (atom_dict['occupancy'])
            line += '%6.2f' % (atom_dict['tempFactor'])
            line += '      '
            line += atom_dict['segid']
            line += atom_dict['element']
            line += atom_dict['charge']

            out_file.write(line + '\n')

        out_file.write('ENDMDL')

        model_num += 1

    out_file.write('END')

    out_file.close()

    print('DONE_model!')
