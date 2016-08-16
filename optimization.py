import rdkit.Chem.AllChem as ac
import rdkit.Chem as chem

import os
import subprocess as sp
from multiprocessing import Pool
from functools import partial

# More imports at the bottom of script.

def optimize_all(func_data, population):
    """
    Apply optimization function to all population members in parallel.

    Individual optimization functions defined within this module should 
    change the `optimized` attribute of ``MacroModel`` instances to 
    ``True``. They should also include the line
    
        if macro_mol.optimized:
            return None
    
    at the start. This prevents optimizing an already optimized
    structure again.
    
    If this function should be used, rather than its serial counterpart
    ``optimize_all_serial``, the ``optimize_population`` method in the
    ``Population`` class must be told to use it.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.
        
    population : Population
        The ``Population`` instance who's members must be optimized.
        
    Modifies
    --------
    MacroMolecule
        This function optimizes the structures of all the 
        ``MacroMolecule`` instances held in `population`. This means
        that their pristine ``.mol`` files are modified to their 
        optimized structures. However, only the content of these files
        is changed. The value of the `prist_mol_file` attributes remain 
        the same. 
    
    Returns
    -------
    None : NoneType
    
    """
    
    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the 
    # module.    
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may require.
    p_func = partial(func, **func_data.params)
    
    # Apply the function to every member of the population, in parallel.
    with Pool() as pool:
        pool.map(p_func, population)
    

def optimize_all_serial(func_data, population):
    """
    Apply optimization function to all population members, serially.

    Individual optimization functions defined within this module should 
    change the `optimized` attribute of ``MacroModel`` instances to 
    ``True``. They should also include the line
    
        if macro_mol.optimized:
            return None
    
    at the start. This prevents optimizing an already optimized
    structure again.
    
    If this function should be used, rather than its parallel 
    counterpart ``optimize_all``, the ``optimize_population`` method in 
    the ``Population`` class must be told to use it.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.
        
    population : Population
        The ``Population`` instance who's members must be optimized.
        
    Modifies
    --------
    MacroMolecule
        This function optimizes the structures of all the 
        ``MacroMolecule`` instances held in `population`. This means
        that their pristine ``.mol`` files are modified to their 
        optimized structures. However, only the content of these files
        is changed. The value of the `prist_mol_file` attributes remain 
        the same. 
    
    Returns
    -------
    None : NoneType
    
    """

    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the 
    # module.    
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may require.
    p_func = partial(func, **func_data.params)
    
    # Apply the function to every member of the population.    
    for member in population:
        p_func(member)

def update_prist_attrs_from_mol2(macro_mol):
    """
    Replaces instance in `prist_mol` from the optimized ``.mol2`` file.
    
    This function uses a ``.mol2`` file's structure to form a new rdkit 
    molecule instance. This new rdkit molecule instance is placed in the 
    `prist_mol` attribute of `macro_mol`.
    
    The ``.mol2`` file should be in the same location that the ``.mol``
    file is. It is converted to a ``.mol`` file so that the ``.mol``
    file in `prist_mol_file` holds the optimized structure.
        
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macro_molecule who's `prist_mol` and `prist_mol_file` 
        attributes are to be updated. Note that the `prist_mol_file`
        attribute itself is not changed. Only the data in the file it
        points to.
        
    Modifies
    --------
    macro_mol.prist_mol
        A new rdkit instance is placed in this attribute. The rdkit
        instances holds the molecule described by the ``.mol2`` file.
        
    macro_mol.prist_mol_file's content
        The content in this ``.mol`` file is replaced with the structure
        of the optimized molecule held in a ``.mol2`` file.
        
    Returns
    -------
    None : NoneType
    
    """
    
    # Get the name of the ``.mol2`` file. It should be in the same
    # directory and have the same name as the ``.mol`` file. Only a
    # different extension.
    mol2 = macro_mol.prist_mol_file.replace('.mol', '.mol2')    
    
    # Update the `prist_mol` attribute.
    print(os.path.isfile(mol2))
    macro_mol.prist_mol = chem.MolFromMol2File(mol2, removeHs=False,
                                              sanitize=False)

#    try:
#        chem.SanitizeMol(macro_mol.prist_mol)
#    except Exception as ex:
#        MacroMolError(ex, macro_mol, 'Sanitizing after optimization.')
    
    # Update content in ``prist_mol_file``.
    macro_mol.write_mol_file('prist')
    
def rdkit_optimization(macro_mol):
    """
    Optimizes the structure of the pristine molecule using rdkit.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure should be optimized.
        
    Modifies
    --------
    macro_mol.prist_mol
        The rdkit molecule held in this attribute has it's structure
        changed as a result of the optimization. This means the
        ``Conformer`` instance held by the rdkit molecule is changed.
    
    macro_mol.prist_mol_file's content
        The content of the ``.mol`` file located at 
        `macro_mol.prist_mol_file`, is changed so that it holds the
        structure of the optimized rdkit molecule.
    
    macro_mol.optimized
        After a successful optimization, this attribute is set to 
        ``True``.
    
    Returns
    -------
    None : NoneType   
    
    """
    
    # If `macro_mol` is already optmized, return.
    if macro_mol.optimized:
        return None
        
    # Sanitize then optimize the rdkit molecule in `prist_mol`.
    chem.SanitizeMol(macro_mol.prist_mol)
    ac.MMFFOptimizeMolecule(macro_mol.prist_mol)
    
    # Update the content of the ``.mol`` file.
    chem.MolToMolFile(macro_mol.prist_mol, macro_mol.prist_mol_file,
                      includeStereo=True, kekulize=False,
                      forceV3000=True)
    
    macro_mol.optimized = True
    
def macromodel_opt(macro_mol, 
                 macromodel_path=r"C:\Program Files\Schrodinger2016-2"):
    """
    Optimizes the molecule using MacroModel.

    This function runs a restricted optimization. The structures of the
    building blocks are frozen and only the new bonds formed between
    building blocks during assembly are optimized.    
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.
        
    macromodel_path : str
        The full path of the ``Schrodinger`` suite within the user's 
        machine. For example, in a default Microsoft installation the 
        folder will probably be something like
        ``C:\Program Files\Schrodinger2016-2``.
    
    Modifies
    --------
    macro_mol.prist_mol
        The rdkit molecule held in this attribute is replaced by an 
        rdkit molecule with an optimized structure.
    
    macro_mol.prist_mol_file's content
        The content of the ``.mol`` file located at 
        `macro_mol.prist_mol_file`, is changed so that it holds the
        structure of the optimized rdkit molecule.
    
    macro_mol.optimized
        After a successful optimization, this attribute is set to 
        ``True``.
    
    Returns
    -------
    None : NoneType       
    
    """

    # If the molecule is already optimized, return.
    if macro_mol.optimized:
        return None
    
    # MacroModel requires a ``.mae`` file as input. This creates a 
    # ``.mae`` file holding the molding the pristine molecule.    
    mae_file = _create_mae_file(macro_mol, macromodel_path)        

    # generate the ``.com`` file for the MacroModel run.
    _generate_COM(macro_mol)
    
    # To run MacroModel a command is issued to to the console via
    # ``subprocess.run``. The command is the full path of the ``bmin``
    # program. ``bmin`` is located in the Schrodinger installation
    # folder. On Windows, to run the software the ``.exe`` extension
    # must be added to the command and the entire path must be enclosed
    # in quotes. The path of the ``.mae`` file to be optimized is then
    # added to the command. On Windows and Unix machines the command
    # should look something like:
    #   "C:\\Program Files\\Schrodinger2016-2\\bmin.exe" mae_file_path
    # and
    #   $SCHRODINGER/bmin mae_file_path
    # respectively. Where ``mae_file_path`` does not include the 
    # ``.mae`` extension.
    file_root = macro_mol.prist_mol_file.replace(".mol", "")
    opt_cmd = os.path.join(macromodel_path, "bmin")
    if os.name == 'nt':
        opt_cmd = '"' + opt_cmd + '.exe"' 

    # Add the -WAIT option to the optimization command. This prevents
    # this means the optimization must finish before the next command
    # can be given to the console.
    opt_cmd = opt_cmd + " -WAIT " + file_root 
    # Run the optimization.
    sp.run(opt_cmd, shell=True)
    
    # Get the ``.mae`` file output from the optimization and convert it
    # to a ``.mol2`` file.
    _convert_mae_to_mol2(macro_mol, macromodel_path)
    
    update_prist_attrs_from_mol2(macro_mol) 

    # This command ensures that programs opened as a result of the 
    # optimization close. If this is not done after a population is
    # optimized, it is often the case the folders containing it cannot 
    # be moved or deleted. This is because some SCHRODINGER programs are
    # still running in the background and accessing them. This closes
    # all such programs.
    close_cmd = os.path.join(macromodel_path, "utilities", "jserver")
    if os.name == 'nt':
        close_cmd = '"' + close_cmd + '.exe"'
    close_cmd += " -cleanall" 
    sp.call(close_cmd, shell=True) 
   
    macro_mol.optimized = True       
    
def _generate_COM(macro_mol):
    """
    Create a ``.com`` file for a MacroModel optimization.

    The created ``.com`` file fixes all bond parameters which were not
    added during assembly. This means all bond distances, bond angles
    and torsional angles are fixed, except for cases where it involves
    a bond added during assembly of the macromolecule.
    
    This fixing is implemented by creating a ``.com`` file with various
    ``FX`` commands written within its body.
    
    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.

    Modifies
    --------
    This function creates a new ``.com`` file holding the instructions
    for optimizing the pristine macromolecule using MacroModel.

    Returns
    -------
    None : NoneType    
    
    """
    
    # This is the body of the ``.com`` file. The line that begins and
    # ends with exclamation lines is replaced with the various commands
    # that fix bond distances and angles.
    main_string= (" MMOD       0      1      0      0     0.0000     "
    "0.0000     0.0000     0.0000\n"
" DEBG      55      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" FFLD      16      1      0      0     1.0000     0.0000     "
"0.0000     0.0000\n"
" BDCO       0      0      0      0    41.5692 99999.0000     "
"0.0000     0.0000\n"
" CRMS       0      0      0      0     0.0000     0.5000     "
"0.0000     0.0000\n"
" BGIN       0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" READ       0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
"!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!\n"
" CONV       2      0      0      0     0.0500     0.0000     "
"0.0000     0.0000\n"
" MINI       1      0   2500      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" END        0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" ")

    # Create a path for the ``.com`` file. It is the same as that of the
    # ``.mol`` file but with a ``.com`` extension. Get the path of the
    # ``.mae`` file and the output file in the same way. 
    com_file = macro_mol.prist_mol_file.replace(".mol", ".com")
    mae = macro_mol.prist_mol_file.replace(".mol", ".mae")
    output = macro_mol.prist_mol_file.replace(".mol", "-out.maegz")
    
    # This function adds all the lines which fix bond distances and 
    # angles into ``main_string``.
    main_string = _fix_params_in_com_file(macro_mol, main_string)
    
    # Writes the ``.com`` file.
    with open(com_file, "w") as com:
        # The first line hold the ``.mae`` file containing the molecule
        # to be optimized.
        com.write(str(mae + "\n"))
        # The second line holds the name of the output file of the 
        # optimization.
        com.write(str(output + "\n"))
        # Next is the body of the ``.com`` file, held in
        # ``main_string``.
        com.write(main_string)


def _create_mae_file(macro_mol, macromodel_path):
    """
    Creates the ``.mae`` file holding the molecule to be optimized.    
    
    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized. Its ``.mol`` file is
        converted to a ``.mae`` file. The original ``.mol`` file is also
        kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".

    Modifies
    --------
    This function creates a new ``.mae`` file from the ``.mol`` file in
    `macro_mol.prist_mol_file`. This new file is placed in the same
    folder as the ``.mol`` file and has the same name. Only the 
    extensions are different.

    Returns
    -------
    str
        The full path of the newly created ``.mae`` file.     
    
    """
    
    # Create the name of the new ``.mae`` file. It is the same as the
    # ``.mol`` file, including the same path. Only the extensions are
    # different.
    mae_file = macro_mol.prist_mol_file.replace('.mol', 
                                                     '.mae')
    
    # ``convrt_cmd`` is the command entered into the console for turning
    # a ``.mol`` file to ``.mae``. It consists of the path to the 
    # program ``structconvert`` followed by the name of ``.mol`` file.
    # The option ``-omae`` specifies that the output should be a 
    # ``.mae`` file. This option is followed by the name of the ``.mae``
    # file. On a Windows machine the path must be placed in quotes and
    # include the ``.exe`` extension. Overall on a Windows and Unix
    # machine the line should look something like:
    #   C:\\Program Files\\Schrodinger2016-2\\utilities\\...
    #  ...structconvert.exe" mol_file.mol -omae mol_file.mae
    # and
    #   $SCHRODINGER/utilities/structconvert mol_file.mol -omae ...
    # ...mol_file.mae
    # respectively.
    convrt_cmd = os.path.join(macromodel_path, 'utilities',
                              'structconvert')  
    # For Windows systems add the ``.exe`` extension and encapsulate
    # path in quotes.                              
    if os.name == 'nt':
        convrt_cmd = '"' + convrt_cmd + '.exe"'    
    convrt_cmd += (" " + macro_mol.prist_mol_file + 
                   " -omae " + mae_file)

    sp.call(convrt_cmd, shell=True)    
    return mae_file

def _convert_mae_to_mol2(macro_mol, macromodel_path):
    """
    Converts a ``.mae`` file to a ``.mol2`` file.

    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized. The ``.mae`` file holding its
        optimized structure is converted to a ``.mol2`` file. Both
        versions are kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".   
    
    Modifies
    --------    
    This function creates a new ``.mol2`` file from the optimized 
    ``.mae`` file. This new file is placed in the same folder as the 
    ``.mae`` file.
    
    Returns
    -------
    None : NoneType
    
    """
    
    # Replace extensions to get the names of the various files.
    mol2 = macro_mol.prist_mol_file.replace(".mol", ".mol2")
    # ``out`` is the full path of the optimized ``.mae`` file.
    out = macro_mol.prist_mol_file.replace(".mol", 
                                                "-out.maegz")
    
    # ``convrt_cmd`` is the command entered into the console for turning
    # a ``.mae`` file to ``.mol2``. It consists of the path to the 
    # program ``structconvert`` followed by the option ``-imae`` and 
    # then the full path of the optimized ``.mae`` file. The option 
    # ``-omol2`` specifies that the output should be a ``.mol2`` file. 
    # This option is followed by the name of the ``.mol2`` file. On a 
    # Windows machine the path must be placed in quotes and include the 
    # ``.exe`` extension. Overall on a Windows and Unix machine the line 
    # should look something like:
    #   C:\\Program Files\\Schrodinger2016-2\\utilities\\...
    #  ...structconvert.exe" -imae mol_file.mae -omol2 mol_file.mol2
    # and
    #   $SCHRODINGER/utilities/structconvert -imae mol_file.mae... 
    # ... -mol2 mol_file.mol2
    # respectively.    
    convrt_cmd = os.path.join(macromodel_path, 'utilities', 
                                                     'structconvert')
    # For Windows systems add the ``.exe`` extension and encapsulate
    # path in quotes.
    if os.name == 'nt':
        convrt_cmd = '"' + convrt_cmd + '.exe"'                
    convrt_cmd = convrt_cmd + " -imae " + out + " -omol2 " + mol2
   
   # Execute the file conversion.
    sp.run(convrt_cmd, shell=True)

#        OPCD 1234567123456712345671234567 FFFFF.FFFF FFFFF.FFFF FFFFF.FFFF FFFFF.FFFF
def _fix_params_in_com_file(macro_mol, main_string):
    """
    Adds lines to the ``.com`` body fixing bond distances and angles.
    
    For each bond distance, bond angle and torisional angle that does
    not involve a bond created during assembly a ``FX`` command is added
    to the string holding holding the body of the ``.com`` file.
    
    These lines replace the filler line in the main string.
    
    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.    
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.
        
    main_string : str
        The body of the ``.com`` file which is to have fix commands
        added.
        
    Returns
    -------
    str
        A string holding the body of the ``.com`` file with instructions
        to fix the various bond distances and angles as described in the
        docstring.
    
    """
    
    # Make a string to hold all of the ``FX`` lines.
    fix_block = ""      
    # Add lines that fix the bond distance.
    fix_block = _fix_distance_in_com_file(macro_mol, fix_block)  
    # Add lines that fix the bond angles.                          
    fix_block = _fix_bond_angle_in_com_file(macro_mol, fix_block)
    # Add lines that fix the torsional angles.
    fix_block = _fix_torsional_angle_in_com_file(macro_mol, fix_block)
    
    return main_string.replace(("!!!BLOCK_OF_FIXED_PARAMETERS_"
                                "COMES_HERE!!!\n"), fix_block)

def _fix_distance_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing bond distances to ``.com`` body string.
    
    Only bond distances which do not involve bonds created during
    assembly are fixed.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond distances added to it
        by this function.
    
    """
    
    # This line holds the format for a line fixing the bond distance
    # between two atoms. The first two ``{...}`` are replaced with the
    # ids of atoms to be fixed. The last ``{...}`` is replaced with the
    # bond distance. Note that in the ``.mae`` files the indices of
    # atoms start at 1 while in rdkit they start at 0. As far as I can 
    # tell this corresponds to a shift of one for each atom index, with
    # the ordering being the same.
    fix_distance = (" FXDI {0:>7}{1:>7}      0      0"
                    "   100.0000 {2:>10.4f}     0.0000     0.0000")
   
    # Go through all the bonds in the heavy rdkit molecule. If the bond
    # is not between heavy atoms get its distance. Add a fix line using
    # the bond distance and atomic indices to the ``fix_block``. If the
    # bond does invovle two heavy atoms go to the next bond. This is
    # because a bond between 2 heavy atoms was added during assembly and
    # should therefore not be fixed.   
    for bond in macro_mol.heavy_mol.GetBonds():
        atom1 = bond.GetBeginAtom() 
        atom2 = bond.GetEndAtom()
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums and
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        atom1_id = atom1.GetIdx() 
        atom2_id = atom2.GetIdx() 
        
        bond_len = macro_mol.prist_distance(atom1_id, atom2_id)
        
        # Make sure that the indices are increased by 1 in the ``.mae``
        # file from their rdkit value.
        fix_block += (fix_distance.format(atom1_id+1, atom2_id+1,
                                         bond_len) + "\n")

    return fix_block
    
def _fix_bond_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing bond angles to ``.com`` body string.
    
    Only bond angles which do not involve bonds created during
    assembly are fixed. The exception to this is for bond angles next
    to bonds added during assembly. For example, consider:
        
        A-B-C-D=E
        
    if the bond between D and E (``=``) represents the bond added during
    assembly, the bond angle A-B-C will be fixed but B-C-D will not be.
    This is an artifact of the implementation but is not expected to
    play a significant role as the vast majority of bond angles which
    should be fixed, will be. The bond angle C-D=E will also not be
    fixed as that is the purpose of this function.         

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond angles added to it by
        this function.
    
    """
    # This line holds the format for a line fixing the bond angles
    # between 3 atoms. The first 3 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the bond angle. Note 
    # that in the ``.mae`` files the indices of atoms start at 1 while 
    # in rdkit they start at 0. As far as I can tell this corresponds to 
    # a shift of 1 for each atom index, with the ordering being the 
    # same.    
    fix_ba = (" FXBA {0:>7}{1:>7}{2:>7}      0   100.0000 "
              "{3:>10.4f}     0.0000     0.0000")
    
    # Create a substructure consisting of 3 dummy atoms bonded with 3
    # dummy bonds. This substructure will match with any 3 atoms which
    # are bonded together with any combination of bonds. These 3 atoms
    # will therefore have a bond angle.          
    ba_mol = chem.MolFromSmarts('[*]~[*]~[*]')
    
    # Get the indices of all atoms which have a bond angle. ``ba_atoms``
    # is a tuple of tuples of the form ((1,2,3), (4,5,6), (7,8,9), ...).
    # Each inner tuple holds the indicies of the atoms which form a bond
    # angle.
    ba_atoms = macro_mol.heavy_mol.GetSubstructMatches(ba_mol)
    
    # Get the conformer holding the atomic positions.
    conf = macro_mol.heavy_mol.GetConformer()    
    
    # For each bond angle check if a heavy atom is involved in forming
    # it. If no, a line fixing the bond angle is added to ``fix_block``.
    # If any atom of the 3 is a heavy atom the bond angle is not fixed.
    # This means that there will be some bond angles which consist of 2
    # bonds not added during assembly which will not be fixed. However,
    # it is assumed that the effect of this will be minimal.
    for atom1_id, atom2_id, atom3_id in ba_atoms:
        atom1 = macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.heavy_mol.GetAtomWithIdx(atom3_id)
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom3.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        ba = ac.GetAngleDeg(conf, atom1_id, atom2_id, atom3_id)
        
        fix_block += (fix_ba.format(atom1_id+1, atom2_id+1, 
                                    atom3_id+1, ba) + "\n")
    
    return fix_block
    
def _fix_torsional_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing torsional bond angles to ``.com`` body string.
    
    Only torsional angles which do not involve bonds created during
    assembly are fixed. The exception to this is for torsional angles 
    next to bonds added during assembly. For example, consider:
        
        A-B-C-D-E=F
        
    if the bond between E and F (``=``) represents the bond added during
    assembly, the torsional angle A-B-C-D will be fixed but B-C-D-E will 
    not be. This is an artifact of the implementation but is not 
    expected to play a significant role as the vast majority of bond 
    angles which should be fixed, will be. The bond angle C-D-E=F will 
    also not be fixed as that is the purpose of this function.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond angles added to it by
        this function.
    
    """

    # This line holds the format for a line fixing the torsional angles
    # between 4 atoms. The first 4 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the torsional angle. 
    # Note that in the ``.mae`` files the indices of atoms start at 1 
    # while in rdkit they start at 0. As far as I can tell this 
    # corresponds to a shift of 1 for each atom index, with the ordering 
    # being the same.    
    fix_ta = (" FXTA {0:>7}{1:>7}{2:>7}{3:>7}   100.0000 "
                "{4:>10.4f}     0.0000     0.0000")        

    # Create a substructure consisting of 4 dummy atoms bonded with 4
    # dummy bonds. This substructure will match with any 4 atoms which
    # are bonded together with any combination of bonds. These 4 atoms
    # will therefore have a torsinal angle.              
    ta_mol = chem.MolFromSmarts('[*]~[*]~[*]~[*]')
    
    # Get the indices of all atoms which have a torsional angle. 
    # ``ta_atoms`` is a tuple of tuples of the form ((1,2,3,4), 
    # (4,5,6,7), ...). Each inner tuple holds the indicies of the atoms 
    # which form a torsional angle.
    ta_atoms = macro_mol.heavy_mol.GetSubstructMatches(ta_mol)
    # Get the conformer holding the atomic positions.
    conf = macro_mol.heavy_mol.GetConformer()
    
    # For each torsional angle check if a heavy atom is involved in 
    # forming it. If no, a line fixing the torsional angle is added to 
    # ``fix_block``. If any atom of the 4 is a heavy atom the bond angle 
    # is not fixed. This means that there will be some bond angles which 
    # consist of 3 bonds not added during assembly which will not be 
    # fixed. However, it is assumed that the effect of this will be 
    # minimal.    
    for atom1_id, atom2_id, atom3_id, atom4_id in ta_atoms:
        atom1 = macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.heavy_mol.GetAtomWithIdx(atom3_id)
        atom4 = macro_mol.heavy_mol.GetAtomWithIdx(atom4_id)
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom3.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom4.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        ta = ac.GetDihedralDeg(conf, atom1_id, atom2_id, 
                                     atom3_id, atom4_id)
        
        fix_block += (fix_ta.format(atom1_id+1, atom2_id+1, 
                                atom3_id+1, atom4_id+1, ta) + "\n")

    return fix_block

def do_not_optimize(macro_mol):
    """
    Skips the optimization step.
    
    This is very useful when debugging so you do not waste your time
    waiting for molecules to get optimized. Use this in the input file
    in place of an optimization function when necessary.
    
    """
    
    return None
    
from .classes import FGInfo
from .classes.exception import MacroMolError