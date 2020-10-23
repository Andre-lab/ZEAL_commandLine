# ZEAL


This manual describes how to use the command-line version of ZEAL in MATLAB, available for download from github link XXXX.  For details on the method itself, see the publication [publication](https://academic.oup.com/bioinformatics):




        F. Ljung and I. André, ZEAL: Protein structure alignment based on shape similarity, *Bioinformatics* XX (2020)


# Requirements


ZEAL depends on the following classes 



   -  **PDB.m **             Handles reading/fetching/parsing of PDB files 
   -  **molShape.m**    Generates a molecular shape function (the molecular surface by default) from PDB object data 
   -  **ZC.m**                 Handles computation of Zernike-Canterakis moments and shape reconstruction 

## Toolboxes


ZEAL uses functions from these toolboxes which have to be installed (included in [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html) if ZEAL run as standalone program)



   -  `bioinformatics_toolbox`           *Bioinformatics toolbox* 
   -  `gads_toolbox`                                 G*lobal Optimization Toolbox* 
   -  `image_toolbox              `*Image Processing Toolbox* 
   -  `optimization_toolbox       `*Optimization Toolbox**`   `* 
   -  `statistics_toolbox         `*Statistics and Machine Learning Toolbox* 
   -  `symbolic_toolbox           `*Symbolic Math Toolbox* 

  
# Abbreviations

   -  ZC             Zernike-Canterakis 
   -  MS            Molecular surface 
   -  SA/SAS     Solvent-accessible surface 
   -  vdw            van der Waals 
   -  PDB          Protein Data Bank 

  
# Input


ZEAL can be run in two modes depending on if one or two protein structures are given as input. Protein structrues can either be given as files (`'filename.pdb'` or ` `'`filepath/filename.pdb'` if not in same directory ) or as 4-letter PDB ID codes (eg. '5MOK'). If a 5-letter PDB ID code is given (eg. '5MOKA'), then the last letter is assumed to be the chain ID ('A') that should be selected for analysis. By default, the A chain of the structure (first model if several exist) is used, excluding hydrogen atoms and any heteroatoms. If any alternate atom locations are defined, the 'A' state is selected. See the section \hyperref{H_56E64AE2}{Changing default parameters} to change what structure data should be used in the analysis.  


  


**Single mode** If one structure is given, then ZEAL will compute the Zernike-Canterkis shape descriptor for the structure (defined in the object as the `fixed` structure)



```matlab:Code
% Single mode
shapeData_pdbFile = ZEAL('5mok.pdb'); % -> reads pdb file 
```


```text:Output
 Importing fixed structure: 5mok.pdb
 Computing shape function for fixed structure
 Computing ZC moments for fixed structure
```


```matlab:Code
shapeData = ZEAL('5mokA'); % -> downloads the structure from PDB 
```


```text:Output
 Importing fixed structure: 5mokA
 Computing shape function for fixed structure
 Computing ZC moments for fixed structure
```

  


**Align mode** Adding a second structure with parameter `'rot'` will start the shape alignment search in ZEAL where the second structure (defined as `rotating` in the object) is rotated 



```matlab:Code
% Align mode
% shapeAlignData_pdbFile = ZEAL('5mok.pdb', 'rot', '2ho1.pdb'); % -> reads PDB files and and performs the shape alignment 
shapeAlignData = ZEAL('5mokA', 'rot', '2ho1A'); % -> Downloads structures and performs the shape alignment
```


```text:Output
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
	 0.30 		 0.00 0.00 0.00 	 0 		  0.2
	 0.50 		 3.14 1.57 3.14 	 1 		  0.3
	 0.56 		 4.73 2.03 2.59 	 20 		  3.3
	 0.56 		 4.23 1.84 2.76 	 22 		  3.7
	 0.68 		 4.29 1.85 2.36 	 27 		  4.4
	 0.76 		 4.02 1.76 2.12 	 31 		  5.1
	 0.77 		 3.99 1.61 2.10 	 39 		  6.3
	 0.79 		 3.77 1.65 2.03 	 47 		  7.6
 ----------------------------------------------------------------------------

	Search completed after 52.1 s.

	Best score 0.79 found after 47 iterations (7.6 s)

	using Euler angles (zyz) [3.770 1.649 2.029]

 /////////////////////////////////////////////////////////////////////////////
```

# Output


**Single mode **The shape descriptors and ZC moments are accessed from the property `fixed.ZC.Descriptors` 



```matlab:Code
shapeData.fixed.ZC.Descriptors
```


```text:Output
ans = 121x1    
    0.0158
    0.0000
    0.0251
    0.0041
    0.0013
    0.0016
    0.0185
    0.0109
    0.0010
    0.0039

```


```matlab:Code
shapeData.fixed.ZC.Moments
```


```text:Output
ans = 
    IndicesList: [1771x5 double]
     CellValues: [21x21x41 double]
         Values: [1771x1 double]

```



where the field `Values` contains the complex-valued moments, `Indiceslist` contains the n,l,m indices (col 1-3) for each moment with real part (col 4) and image part (col 5) separately. 


  


**Align mode **The shape descriptors and ZC moments for the rotating structure can be access as above, but from the property `rotating`  instead. The ZEAL score is accessed with



```matlab:Code
shapeAlignData.Score
```


```text:Output
ans = 0.7925
```



and the Euclidean distance between the shape descriptors with



```matlab:Code
shapeAlignData.ZCDdistance
```


```text:Output
ans = 0.0152
```

  
  
# Optional input


The parameters listed below are set by default, but can be changed (described in \hyperref{H_56E64AE2}{Changing default settings}). 



```matlab:Code
settingsTable = getZEALSettingsTable()
```

| |parameter|type|default value|expected values|description|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'Order'|'integer'|20|'>0'|'The maximum expansi...|
|2|'ChiCoeffPath'|'char'|'[pwd '/chi_coeffici...|'folder path'|'Path to folder with...|
|3|'GridRes'|'integer'|64|'>0'|'The side length of ...|
|4|'Shape'|'char'|'MS'|''MS'/ 'SAS'/ 'vdw'/...|'The type of molecul...|
|5|'ProbeRadius'|'double'|1.4000|'>=0'|'The radius of the p...|
|6|'SmearFactor'|'double'|0.3000|'>0, <1'|'Fraction of grid to...|
|7|'ShellThickness'|'integer'|2|'>0'|'Thickness of surfac...|
|8|'FunEvals'|'integer'|300|'>0'|'Number of ZEAL scor...|
|9|'AlignLater'|'logical'|0|'true/false'|'If false then ZEAL ...|
|10|'fix_includeHetatoms...|'logical'|0|'true/false'|'Flag to indicate if...|
|11|'rot_includeHetatoms...|'logical'|0|'true/false'|'"---" in rotating s...|
|12|'fix_includeHatoms'|'logical'|0|'true/false'|'Flag to indicate if...|
|13|'rot_includeHatoms'|'logical'|0|'true/false'|'"---" in rotating s...|
|14|'fix_chainID'|'char'|'A'|''all' or 1 letter'|'The chain ID that s...|
|15|'rot_chainID'|'char'|'A'|''all' or 1 letter'|'"---" in rotating s...|
|16|'fix_altLocID'|'char'|'A'|''all' or 1 letter'|'The ID of any atom ...|
|17|'rot_altLocID'|'char'|'A'|''all' or 1 letter'|'"---" in rotating s...|
|18|'fix_modelNumber'|'integer'|1|''|'The model number th...|
|19|'rot_modelNumber'|'integer'|1|''|'"---" in rotating s...|
|20|'LogLevel'|'char'|'standard'|''none'/ 'basic'/ 's...|'Determines the leve...|

  
  
# Changing default parameters


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

