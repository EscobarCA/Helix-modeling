# -*- coding: utf-8 -*-
"""Module for manupulation of structures using Pymol module.

by Cristian Escobar
Last update 2/11/2021
cristian.escobar.b@gmail.com
"""

from pymol import cmd
import PDB


def get_Ca_distance(pymol_object, chainID1, aa_number1, chainID2, aa_number2):
    """Get distance between CA atoms between two residues.

    Method to measure the distance between Ca carbon of different residues
    defined by their chainID and residue number.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID1 : String
        ChainID of the first residue.
    aa_number1 : Integer
        Residue number of the first residue.
    chainID2 : String
        ChainID of the second residue.
    aa_number2 : Integer
        Residue number of the second residue.

    Returns
    -------
    distance : Float
        Value of distance between the two Ca carbons selected.


    """
    # Atom selection following Pymol format
    atom1 = pymol_object + '//' + chainID1 + '/' + str(aa_number1) + '/CA'
    atom2 = pymol_object + '//' + chainID2 + '/' + str(aa_number2) + '/CA'

    distance = cmd.get_distance(atom1, atom2)

    return distance


def evaluate_clash(atom1, atom2, threshold):
    """Evaluate if atoms are closer that trashold.

    Parameters
    ----------
    atom1 : String
        String selection for atom1.
    atom2 : String
        String selection for atom1.
    threshold : Float
        Minimum distance to report a clash.

    Returns
    -------
    bool
        If distance is less than threshold, it returns True.

    """
    # get distance between atoms
    distance = cmd.get_distance(atom1, atom2)

    if distance <= threshold:

        return True

    return False


def evaluate_clash_objects(pymol_object1, chainID1, Aa_list1,
                           pymol_object2, chainID2, Aa_list2, threshold):
    """Evaluate whether there is a clash between two objects.

    Function to evaluate if there is a clash between CA carbons of one object
    and a second one. If the distance between the pair is less than threshold
    in Angstroms, it prints the chain and residue number of residues involved
    in the clash.

    Parameters
    ----------
    pymol_object1 : String
        Name of Pymol object 1 in the session.
    chainID1 : String
        ChainID of in object 1 to be used.
    Aa_list1 : List
        List of residues in object 1 to be evaluated.
    pymol_object2 : String
        Name of Pymol 2 object in the session.
    chainID2 : String
        ChainID of in object 2 to be used.
    Aa_list2 : List
        List of residues in object 2 to be evaluated.
    threshold : Float
        Minimum distance to report a clash.

    Returns
    -------
    None.

    """
    # Check clash between all atom pairs
    for residue1 in Aa_list1:

        for residue2 in Aa_list2:

            # select Ca carbons
            atom1 = pymol_object1 + '//' + chainID1 + '/' + str(residue1) + '/CA'

            atom2 = pymol_object2 + '//' + chainID2 + '/' + str(residue2) + '/CA'

            # check for clash
            if evaluate_clash(atom1, atom2, threshold):

                print(chainID1, residue1, '   ', chainID2, residue2)

    print('Finished!!')


def get_phi_angle(pymol_object, chainID, aa_number):
    """Get Phi dihedral angle of a given residue.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_number : Integer
        Residue number.

    Returns
    -------
    angle : Float
        Value of Phi angle.

    """
    # Atom selection for each atom that participates in the angle
    # Ci-1
    atom1 = pymol_object + '//' + chainID + '/' + str(aa_number-1) + '/C'
    # Ni
    atom2 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/N'
    # CAi
    atom3 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/CA'
    # Ci
    atom4 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/C'

    angle = cmd.get_dihedral(atom1, atom2, atom3, atom4)

    return angle


def get_psi_angle(pymol_object, chainID, aa_number):
    """Get Psi dihedral angle of a given residue.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_number : Integer
        Residue number.

    Returns
    -------
    angle : Float
        Value of Psi angle.

    """
    # Atom selection for each atom that participates in the angle
    # Ni
    atom1 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/N'
    # CAi
    atom2 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/CA'
    # Ci
    atom3 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/C'
    # Ni+1
    atom4 = pymol_object + '//' + chainID + '/' + str(aa_number+1) + '/N'

    angle = cmd.get_dihedral(atom1, atom2, atom3, atom4)

    return angle


def get_dihedral_angles(pymol_object, chainID, aa_list):
    """Get Phi and Psi dihedral angles.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_list : List
        List of amino acid numbers (int) used to get angles.

    Returns
    -------
    angle_list : List
        List of phi and psi angles (in a tuple) for each residue in list.
        angle_list = [(phi, psi), ....]

    """
    aa_counter = 1  # Track Aa number

    angle_list = []

    # sum over all amino acids
    phi_total = 0
    psi_total = 0

    # first line
    line = '  '
    line += '   Aa  '
    line += '  Phi  '
    line += '  Psi  '

    print(line)

    # get amino acid angles
    for aa_number in aa_list:

        phi = get_phi_angle(pymol_object, chainID, aa_number)
        psi = get_psi_angle(pymol_object, chainID, aa_number)

        angle_list.append((phi, psi))

        phi_total += phi
        psi_total += psi

        # new line for each amino acid
        line = '%2d ' % (aa_counter)
        line += '%4d ' % (aa_number)
        line += '%3.3f, ' % (phi)
        line += '%3.3f ' % (psi)

        print(line)

        aa_counter += 1

    # last line with average
    line = 'Ave '
    line += '%3.3f ' % (phi_total / float(len(aa_list)))
    line += '%3.3f ' % (psi_total / float(len(aa_list)))

    print(line)

    return angle_list


def set_phi_angle(pymol_object, chainID, aa_number, new_angle):
    """Set Phi dihedral angle of a given residue.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_number : Integer
        Residue number.
    new_angle : Float
        New phi angle to be set.

    Returns
    -------
    None.

    """
    # Atom selection for each atom that participates in the angle
    # Ci-1
    atom1 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/C'
    # Ni
    atom2 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/CA'
    # CAi
    atom3 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/N'
    # Ci
    atom4 = pymol_object + '//' + chainID + '/' + str(aa_number-1) + '/C'

    cmd.set_dihedral(atom1, atom2, atom3, atom4, new_angle)


def set_psi_angle(pymol_object, chainID, aa_number, new_angle):
    """Set Psi dihedral angle of a given residue.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_number : Integer
        Residue number.
    new_angle : Float
        New psi angle to be set.

    Returns
    -------
    None.

    """
    # Atom selection for each atom that participates in the angle
    # Ni+1
    atom1 = pymol_object + '//' + chainID + '/' + str(aa_number+1) + '/N'
    # Ci
    atom2 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/C'
    # CAi
    atom3 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/CA'
    # Ni
    atom4 = pymol_object + '//' + chainID + '/' + str(aa_number) + '/N'

    cmd.set_dihedral(atom1, atom2, atom3, atom4, new_angle)


def set_phi_psi_angles(pymol_object, chainID, aa_number, phi=0.0, psi=0.0):
    """Change both phi and psi angles of a given residue.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    aa_number : Integer
        Residue number.
    phi : Float
        Value for new Phi angle. The default is 0.0.
    psi : Float
        Value for new Psi angle. The default is 0.0.

    Returns
    -------
    None.

    """
    set_phi_angle(pymol_object, chainID, aa_number, phi)

    set_psi_angle(pymol_object, chainID, aa_number, psi)


def set_dihedral_angles_update(pymol_object, chainID, Aa_list,
                               phi_psi_list, protect_selection,
                               moving_atom, target_atom, start_num=1):
    """Set several dihedral angles to a set of amino acids.

    Function that will change the dihedral angles for each amino acid in
    Aa_list for the ones in the phi_psi_list. Amino acids selected in
    protect_selection will not move. After changing the angles the new
    structure is saved.
    This method is an updated version of "set_dihedral_angles". Instead of
    changing the angles of the structures and then setting them back to
    the original values, this method creates a copy of the pymol object
    and then changes its angles. Thus, this does not modify the initial
    structure. 
    This method also returns a list with information regarding the amino acid
    modified, the angle used, the distance between moving atom and target atom
    and the pdb file created. This is useful to create a log of the
    modifications made.

    Parameters
    ----------
    pymol_object : String
        Name of Pymol object in the session.
    chainID : String
        ChainID of residue.
    Aa_list : List
        List of amino acid numbers (int).
    phi_psi_list : List
        List of tuples containing phi and psi angles.
    protect_selection : String
        Slection string from pymol_object which will remain static.
        Example: '//X/44-316/*'
    moving_atom : String
        Atom selection for atom 1 used for distance mesurement, ussualy from
        chain being modified.
    target_atom : String
        Atom selection for atom 2 used for distance mesurement.
    start_num : Int, optional
        Number used to start the pdb files counter. The default is 1.

    Returns
    -------
    pdb_files : List
        List of tuples containing information about the output pdb file,
        distance between target-moving atoms, aminoacid modified and new phi
        and psi angles: (pdb_name,distance, aa_number, (new_phi, new_psi)).

    """
    # track of PDB file numbers
    pdb_num = start_num

    # pdb_file_list = [] #PDB file name list as output

    pdb_files = []

    for aa_number in Aa_list:

        for dihedral_angles in phi_psi_list:

            # Create a copy of object to be modified
            cmd.copy('obj_copy', pymol_object)

            # Protect residues
            cmd.protect('obj_copy'+protect_selection)

            new_phi, new_psi = dihedral_angles

            set_phi_psi_angles('obj_copy', chainID, aa_number,
                               phi=new_phi, psi=new_psi)

            # get distance
            distance = cmd.get_distance('obj_copy' + moving_atom, target_atom)

            # save new structure with pdb_num as name

            pdb_name = '{}_{}.pdb'.format(pymol_object, pdb_num)

            cmd.save(pdb_name, 'obj_copy')

            pdb_files.append((pdb_name, distance, aa_number, (new_phi, new_psi)))

            pdb_num += 1

            # Remove the object copy
            cmd.delete('obj_copy')

    print('Finished!')

    return pdb_files


def prepare_log(model_list, out_file_name):
    """Create output log with models created.

    Method takes information from set_dihedral_angles_update to prepare a log
    with the modifications made to the initial structure. It also combines the
    output pdb files into a single pdb file with several states.

    Parameters
    ----------
    data_list : List
        Otput from set_dihedral_angles_update.
    out_file_name : String
        file name for output log and pdb states file.

    Returns
    -------
    None.

    """
    out_file = open(out_file_name+'.log', 'w')

    models_name = []

    out_file.write('{:<3} {:<20} Distance\n'.format('#', 'File name'))

    total_models = len(model_list)

    for model in range(total_models):

        model_name = model_list[model][0]
        distance = model_list[model][1]
        aa_number = model_list[model][2]
        angles = model_list[model][3]

        line = '{:<3}  {:<20} {:<8.2f} {} {}\n'.format(model + 1, model_name,
                                                       distance, aa_number,
                                                       angles)

        out_file.write(line)

        models_name.append(model_name)

    out_file.close()

    # prepare pdb file with each model as a state

    out_file = out_file_name+'.pdb'

    PDB.print_models(models_name, out_file)


# List of phi and psi angles tested
phi_psi_list = [(-116.46990203857422, -41.119171142578125),
                (-100.94406127929688, -14.551700592041016),
                (-89.06173706054688, -81.47402954101562),
                (-86.58316802978516, -41.41770553588867),
                (-80.1707534790039, -50.459205627441406),
                (-79.94315338134766, -31.915937423706055),
                (-76.61148071289062, -37.83049011230469),
                (-76.44195556640625, -28.4343318939209),
                (-73.86849212646484, -28.97428321838379),
                (-71.58903503417969, -27.29974937438965),
                (-71.18693542480469, -37.75306701660156),
                (-70.86012268066406, -37.606258392333984),
                (-70.50489807128906, -25.322702407836914),
                (-70.42333221435547, -11.162373542785645),
                (-69.79801940917969, -41.28443145751953),
                (-68.54573059082031, -33.29055404663086),
                (-67.78120422363281, -41.32462692260742),
                (-66.43730163574219, -47.29899978637695),
                (-66.08631896972656, -23.022178649902344),
                (-65.08526611328125, -30.872766494750977),
                (-64.95262145996094, -48.96770477294922),
                (-63.79671096801758, -37.10122299194336),
                (-62.65116882324219, -32.652061462402344),
                (-62.33879470825195, -47.58149337768555),
                (-61.576995849609375, -44.33405303955078),
                (-60.91769790649414, -49.872581481933594),
                (-60.769203186035156, -55.0628776550293),
                (-59.36294174194336, -44.03278732299805),
                (-58.62078857421875, -48.016719818115234),
                (-58.105690002441406, -42.14945983886719),
                (-57.34004592895508, -39.09020233154297),
                (-57.048221588134766, -47.26634216308594),
                (-56.90250015258789, -48.021705627441406),
                (-56.62630081176758, -22.010334014892578),
                (-55.257667541503906, -48.15983581542969),
                (-55.248226165771484, -41.14173126220703),
                (-53.33686828613281, -39.00072479248047),
                (-44.463218688964844, -54.57436752319336),
                (-39.56572341918945, -50.36243438720703),
                (-38.418636322021484, -54.35171890258789)]


# Example of code used to modify a helix on residues 16 to 29

if __name__ == '__main__':

    # load Pymol sesion with object to be modified
    cmd.load('Helix_example.pse')

    #Run test of angles
    pymol_object = 'U_helix'    # Object in the session to be modified
    chainID = 'U'
    Aa_list = range(29, 15, -1)   # List of aminoacids
    protect_selection = '//U/32-161/*'
    target_atom = '/5wda/D/D/1/CA'  # target Ca atom in residue 1 in ChainID D
    moving_atom = '//U/1/CA'

    model_list = set_dihedral_angles_update(pymol_object, chainID, Aa_list,
                                            phi_psi_list, protect_selection,
                                            moving_atom, target_atom,
                                            start_num=1)

    # prepare log and output pdb states
    prepare_log(model_list, 'U_29-16')
