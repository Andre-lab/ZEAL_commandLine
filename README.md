
![ZEAL graphical abstract](https://github.com/Andre-lab/ZEAL_commandLine/blob/main/graphicalAbstract_web.png)

This manual describes how to use the command-line version of ZEAL in MATLAB, 
available for download from github link XXXX.  For details on the method itself, 
see the publication [publication](https://academic.oup.com/bioinformatics):




        F. Ljung and I. André, ZEAL: Protein structure alignment based on shape similarity, *Bioinformatics* XX (2020)


# Requirements


ZEAL depends on the following classes 



   -  **PDB.m**             Handles reading/fetching/parsing of PDB files 
   -  **molShape.m**    Generates a molecular shape function (the molecular surface by default) from PDB object data 
   -  **ZC.m**                 Handles computation of Zernike-Canterakis moments and shape reconstruction 
   -  **ChiCoeffs.m**     Handles loading of pre-computed (object-independent) Chi-coefficients from .mat-files

## Toolboxes


ZEAL uses functions from these toolboxes which have to be installed (included in [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html) if ZEAL run as standalone program)


   -  `gads_toolbox`                                 G*lobal Optimization Toolbox* 
   -  `image_toolbox`              *Image Processing Toolbox* 
   -  `optimization_toolbox`       *Optimization Toolbox*
   -  `statistics_toolbox`         *Statistics and Machine Learning Toolbox* 
   -  `symbolic_toolbox`           *Symbolic Math Toolbox* 

  
# Abbreviations

   -  **ZC** Zernike-Canterakis 
   -  **MS** Molecular surface 
   -  **SA/SAS** Solvent-accessible surface 
   -  **vdW** van der Waals 
   -  **PDB** Protein Data Bank 

  
# Input


ZEAL can be run in two modes depending on if one or two protein structures 
are given as input. Protein structrues can either be given as files 
(`'filename.pdb'` or ` `'`filepath/filename.pdb'` if not in same directory) 
or as 4-letter PDB ID codes (eg. '5MOK'). If a 5-letter PDB ID code is given 
(eg. '5MOKA'), then the last letter is assumed to be the chain ID ('A') 
that should be selected for analysis. By default, all chains of the 
structure (first model if several exist) is used, excluding hydrogen atoms 
and any heteroatoms. If any alternate atom locations are defined, all states  
are selected. See the section \hyperref{H_56E64AE2}{Changing 
default parameters} to change what structure data should be used in the analysis.  


  


**Single mode** If one structure is given, then ZEAL will compute the 
Zernike-Canterkis shape descriptor for the structure 
(defined in the object as the `fixed` structure)



```matlab:Code
% Single mode
shapeData_pdbFile = ZEAL('5mok.pdb'); % -> reads pdb file 
```


```text:Output
Running ZEAL in single mode
 Importing fixed structure: 5mok.pdb
 Computing shape function for fixed structure
 Computing ZC moments for fixed structure
```


```matlab:Code
shapeData = ZEAL('5mokA'); % -> downloads the structure from PDB 
```


```text:Output
Running ZEAL in single mode
 Importing fixed structure: 5mokA
 Computing shape function for fixed structure
 Computing ZC moments for fixed structure
```

  


**Align mode** Adding a second structure with parameter `'rot'` will start 
the shape alignment search in ZEAL where the second structure 
(defined as `rotating` in the object) is rotated 



```matlab:Code
% Align mode
% shapeAlignData_pdbFile = ZEAL('5mok.pdb', 'rot', '2ho1.pdb'); % -> reads PDB files and and performs the shape alignment 
shapeAlignData = ZEAL('5mokA', 'rot', '2ho1A'); % -> Downloads structures and performs the shape alignment
```


```text:Output
Running ZEAL in Align mode
 Importing fixed structure: 5mokA
 Importing rotating structure: 2ho1A
 Computing shape function for fixed structure
 Computing shape function for rotating structure
 Computing ZC moments for fixed structure
 Computing ZC moments for rotating structure

 /////////////////////////////////////////////////////////////////////////////

	Searching for best shape superposition

	Stopping after 300 function evaluations

 ----------------------------------------------------------------------------
 Current best score      Euler (zyz)         iteration       time (s) 
 ----------------------------------------------------------------------------
	 0.31 		 0.00 0.00 0.00 	 0 		  0.2
	 0.51 		 3.14 1.57 3.14 	 1 		  0.4
	 0.58 		 4.82 2.03 2.56 	 20 		  3.5
	 0.61 		 4.91 1.92 2.18 	 23 		  4.0
	 0.62 		 4.63 2.09 1.78 	 27 		  4.7
	 0.70 		 4.45 2.04 2.08 	 31 		  5.4
	 0.74 		 4.15 1.92 1.89 	 34 		  5.9
	 0.78 		 3.87 1.75 1.87 	 35 		  6.1
	 0.79 		 4.03 1.71 2.05 	 43 		  7.4
	 0.80 		 3.83 1.63 2.12 	 51 		  9.0
	 0.81 		 3.77 1.62 2.12 	 109 		 21.9
 ----------------------------------------------------------------------------

	Search completed after 59.9 s.

	Best score 0.81 found after 109 iterations (21.9 s)

	using Euler angles (zyz) [3.766 1.616 2.120]

 /////////////////////////////////////////////////////////////////////////////
```

# Output

## Shape descriptors

The shape descriptors are accessed from the property `fixed.ZC.Descriptors` (single / align mode) 
or `rotating.ZC.Descriptors` (align mode only). 

```matlab:Code
shapeData.fixed.ZC.Descriptors
shapeData.rotating.ZC.Descriptors
```

Alternatively, they can be pulled directly using the `getShapeDescriptors` method


```matlab:Code
shapeData.getShapeDescriptors % -> fixed by default
shapeData.getShapeDescriptors('fixed')
shapeData.getShapeDescriptors('rotating')

% or 

getShapeDescriptors(shapeData)
getShapeDescriptors(shapeData, 'fixed')
getShapeDescriptors(shapeData, 'rottaing')

```

## ZC moments 

The moments are accessed from the property `fixed.ZC.Moments` (single / align mode) 
or `rotating.ZC.Moments` (align mode only), which will retirn a Matlab structure  
where the field `Values` contains the complex-valued moments, `Indiceslist` 
contains the n,l,m indices (col 1-3) for each moment with real part (col 4)
 and image part (col 5) separately. 


```matlab:Code
shapeData.fixed.ZC.Moments
shapeData.Rotating.ZC.Moments
```


```text:Output
ans = 
    IndicesList: [1771x5 double]
     CellValues: [21x21x41 double]
         Values: [1771x1 double]

```

Alternatively, they can be pulled directly using the `getMoments` method
```matlab:Code
shapeData.getMoments % -> fixed by default
shapeData.getMoments('fixed')
shapeData.getMoments('rotating')

% or 

getMoments(shapeData)
getMoments(shapeData, 'fixed')
getMoments(shapeData, 'rottaing')

```



## Shape similarity and ZEAL score (align mode only)
 The ZEAL score is accessed with

```matlab:Code
shapeAlignData.Score
```

```text:Output
ans = 0.8153
```


and the shape similairy (=Euclidean distance between the shape descriptorsm) with


```matlab:Code
shapeAlignData.ZCDdistance
```


```text:Output
ans = 0.0253
```

## Structures

Aligned structrues can be saved to PDB files using the save2pdb method. 
By default, both the fixed and rotating structure are exported to the current 
directory (see below on how to change). The names of the new pdb files have 
the format `originalname_ZEAL.pdb `(this can't be changed)`. `Also, HETATM 
records are omitted but this can be changed (se below).



```matlab:Code
save2pdb(shapeAlignData)

% or 

shapeAlignData.save2pdb

```


```text:Output
outputting PDB in file 5mokA_ZEAL.pdb ...
 done! closing file...
outputting PDB in file 2ho1A_ZEAL.pdb ...
 done! closing file...
```



To export ***all records***, use



```matlab:Code(Display)
save2pdb(shapeAlignData, 'includeAll', true)

% or 

shapeAlignData.save2pdb('includeAll', true)

```



If the original pdb file contains multiple chains, and the alignment was done 
with respect to chain X of the fixed structure and chain Y of the rotating 
structure, the coordinates are transformed so that the centroid of the X 
chain and Y chain are placed at the origo, and coordinates of the rotating 
structure are rotated relative that center. 


To ***include HETATM records***, use 



```matlab:Code(Display)
save2pdb(shapeAlignData, 'includeHetatoms', true)

% or 

shapeAlignData.save2pdb('includeHetatoms', true)

```



To ***save to specific directory***, use 



```matlab:Code(Display)
save2pdb(shapeAlignData, 'folderPath', '/Users/yourUserName/Desktop')

% or 

shapeAlignData.save2pdb('folderPath', '/Users/yourUserName/Desktop')
```

  
# Optional input



The parameters listed below are set by default, but can be changed (described in \hyperref{H_56E64AE2}{Changing default parameters}). 

### Shape 
|parameter|type|default value|expected values|description|
|:--:|:--:|:--:|:--:|:--:|
|'Order' | 'integer' | 20 | '>0' | The maximum expansion order of ZC moments.|
|'ChiCoeffPath' | 'char' | [pwd '/chi_coefficients'] | 'folder path' | Path to folder with pre-computed chi coefficients for order N with name  `chiCoeffs_order_N.mat`. Use ZC.computeChiCoeffs to compute them.|
|'GridRes' | 'integer'|64|'>0' | The side length of the cubic grid.|
|'Shape' | 'char'|'MS'| 'MS'/ 'SAS'/ 'vdw'/ 'electron_density' | The type of molecular shape function.|
|'ProbeRadius' | 'double'|1.4000| '>=0' | The radius of the probe in Å used for the surface computations.|
|'SmearFactor' | 'double'|0.3000| '>0, <1' | 'Fraction of grid to smear out atoms over when generating electron density.|
|'ShellThickness' | 'integer'|2|' >0' | Thickness of surfaces (shells) in grid units.|

### Search
|parameter|type|default value|expected values|description|
|:--:|:--:|:--:|:--:|:--:|
|'FunEvals' | 'integer' | 300 | '>0'| Number of ZEAL score evaluations until search stops..|
|'AlignLater' | 'logical' | 0 | 'true/false'| If false then ZEAL will not automatically align upon object construction. Manually start alignment using `shapeAlign(ZEALobject)` or `ZEALobject.shapeAlign()` .|

### PDB/CIF data
|parameter|type|default value|expected values|description|
|:--:|:--:|:--:|:--:|:--:|
|'fix_includeHetatoms | 'logical' |0| 'true/false' | Flag to indicate if HETATM records should be included in fixed structure.|
|'rot_includeHetatoms | 'logical' |0| 'true/false' | "---" in rotating structure.|
|'fix_includeHatoms' | 'logical' |0| 'true/false' | Flag to indicate if Hydrogen atoms should be included (if exists) in fixed structure..|
|'rot_includeHatoms' | 'logical' |0| 'true/false' | "---" in rotating structure.|
|'fix_chainID'|'char' | 'all' | 'all' or 1 letter' | The chain ID that should be selected (''all'' = all chains) in fixed structure.|
|'rot_chainID'|'char' | 'all' | 'all' or 1 letter' | "---" in rotating structure.|
|'fix_altLocID' |  'integer' | 'all' | 'all' or 1 letter' | The ID of any atom altlocs that should be selected in fixed structure. Use ''all'' to include all altlocs. |
|'rot_altLocID' |  'integer' | 'all' | 'all' or 1 letter' | "---" in rotating structure. |
  


## Changing default parameters


All parameters above can be changed by giving a comma-seperated argument list to ZEAL, or as a Matlab structure with fields equal to names of parameters. 


### Example: Using argument list

```matlab:Code(Display)
ZEAL('5mok','fix_chainID','B') % selects chain B for shape analysis 
ZEAL('5mok','fix_chainID','B', 'Shape', 'SAS', 'GridRes', 100) % + sets the solvent accessible surface as the shape function computed on a 100x100x100 grid
```

### Example: Using input structure

```matlab:Code(Display)
inputStruct.fix_chainID = 'B';
inputStruct.Shape = 'SAS';
inputStruct.GridRes = 100;

ZEAL('5mok', inputStruct)
```

  

