classdef ZEALtest < matlab.unittest.TestCase
    %ZEALTEST class-based unit tests for ZEAL
    % unit =  smallest piece of code that can be logically isolated in a system
    %         = function, subroutine, method or property 
    %
    % USAGE
    %
    % Run all tests
    % -----------------------
    % testCase = ZEALtest;
    % res = run(testCase)
    %
    % Run single test method 
    % -----------------------
    % testCase = ZEALtest;
    % res = run(testCase, 'testPDBread'
    
    
    % Available tests
    % name                   description
    %                        |-> option different from default behaviour
        
    % testPDB_default        Test reading PDB file
    % testPDB_chain          |-> specific chain 
    % testPDB_hetatoms       |-> include hetatoms  
    
        
    % testCIF_default        Test reading CIF files
    % testCIF_chain          |-> specific chain 
    % testCIF_hetatoms       |-> include hetatoms  
    
    % testDownload_default   Test downloading (cif) structures from https://www.rcsb.org/
    % testDownload_chain     |-> specific chain 
    % testDownload_hetatoms  |-> include hetatoms  
    
    % testAlignment_default     Test performing an alignment 
    % testAlignment_funevals    |-> change stopping criterium 
    % testAlignment_alignlater  |-> hold alignment 
    
   

    
    
    % Test everything
    methods (Test)
        
    end
    
    
    % Test importing data 
    methods (Test, TestTags = {'class','PDB', 'CIF', 'Download'})
        
        function testPDB_default(testCase)
           % read PDB file containing 
            PDBstruct = PDB('5mok.pdb'); 
            
           % assert default selection criteria expected for ZEAL
           testCase.assertFalse(PDBstruct.Selection.includeHatoms, 'HETATM should not be selected by default');
           
           testCase.assertFalse(PDBstruct.Selection.includeHetatoms, 'Hydrogen atoms should not be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.chainID,'all', 'All chains in file should be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.altLocID, 'A', 'Altloc ID should be ''A'' by default');
           
           % assert data contains expected number of atoms
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 6358, 'Number of ATOM records differ from expected');
            
        end
        
        function testPDB_chain(testCase)
            
             PDBstruct = PDB('5mok.pdb', 'chainID', 'B');
             
             nAtoms = length(PDBstruct.Data.X);
             
             
             
            
        end
        
    end
        
    

    
    
    % Test alignment 
    
    
    
    
    

end

