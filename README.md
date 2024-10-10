# ShowFeats
Tools to visualize pharmacophoric features.

This is a modification of Greg Landrum's original code that control's PyMOL. This version instead takes an SDF or MOL file and outputs a PDB formatted file by taking advantage of the default element colors in PyMOL to indicate different types of pharmacophoric features. It works best using CPK spheres, and ``set sphere_scale=0.25``.

## Usage:
``ShowFeatsToPDB.py [optional args] <filenames>``

Note: you must use ``--writeFeatsAsPDB`` to output a PDB file.

``<filenames>`` can be one or more sd or mol formatted files.

if "-" is provided as a filename, data will be read from stdin (the console)

## Example 1:
``python ShowFeatsToPDB.py --writeFeatsAsPDB --fdef=$RDBASE/Data/BaseFeatures.fdef input.sdf > output_feats.pdb``

 When the ``output_feats.pdb`` file is loaded into PyMOL and colored by element, the features will have the colors listed in the table "Color Key".
 
 Note that the ``--fdef=$RDBASE/Data/BaseFeatures.fdef`` flag is optional, and relies on RDKit's own feature definitions. If RDKit was installed using ``conda``, this may be in ``$HOME/opt/anaconda3//share/RDKit/Data/BaseFeatures.fdef``.

## Example 2:
If you give executable permission to ``ShowFeatsToPDB.py``, (using ``chmod +x ShowFeatsToPDB.py``, you can use a simpler command line call:

``ShowFeatsToPDB.py --writeFeatsAsPDB my_molecule.sdf > my_molecules_features.pdb``

## Color Key

|Pharmacophoric Feature |Element |Color  |
|-----------------|--------|-----------|
|Donor            |   F    |light blue |
|Acceptor         |   B    |pink       |
|NegIonizable     |   O    |red        |
|PosIonizable     |   N    |blue       |
|ZnBinder         |   D    |grey       |
|Aromatic         |   P    |orange     |
|LumpedHydrophobe |   Au   |gold       |
|Hydrophobe       |   S    |yellow     |
Read more about pharmacophores here: https://en.wikipedia.org/wiki/Pharmacophore
