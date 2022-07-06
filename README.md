# ShowFeats
Tools to visualize pharmacophoric features.

This is a modification of Greg Landrum's original code that control's PyMOL. This version instead takes an SDF or MOL file and outputs a PDB formatted file by taking advantage of the default element colors in PyMOL to indicate different types of pharmacophoric features. It works best using CPK spheres, and ``set sphere_scale=0.25``.

## Usage:
``ShowFeatsToPDB.py [optional args] <filenames>``

``<filenames>`` can be one or more sd or mol formatted files.

if "-" is provided as a filename, data will be read from stdin (the console)

## Example:
``python ShowFeatsToPDB.py --writeFeatsAsPDB --fdef=$RDBASE/Data/BaseFeatures.fdef input.sdf > output_feats.pdb``

 When this ``output_feats.pdb`` file is loaded into PyMOL and colored by element, the following features will have the following colors:

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
