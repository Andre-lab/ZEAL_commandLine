classdef ZEALtest < matlab.unittest.TestCase
    %ZEALTEST class-based unit tests for ZEAL
    % unit =  smallest piece of code that can be logically isolated in a system
    %         = function, subroutine, method or property 
    %
    % USAGE
    %
    % Run all tests
    % ----------------------------------------------
    % testCase = ZEALtest;
    % res = run(testCase)
    %
    % Run single test method 
    % ----------------------------------------------
    % testCase = ZEALtest;
    % res = run(testCase, 'testPDBread'
    %
    % Run groups of test based on category
    % ----------------------------------------------
    % results = runtests('ZEALtest','Tag','PDB');
    % 
    
    
    % Available tests
    % ----------------------------------------------
    % name                   description
    %                        |-> option different from default behaviour
    % ----------------------------------------------
    %
    %                        |-> option different from default behaviour
        
    % testPDB_default        Test reading PDB file
    % testPDB_chain          |-> specific chain 
    % testPDB_hetatoms       |-> include hetatoms  
    
        
    % testCIF_default        Test reading CIF files
    % testCIF_chain          |-> specific chain 
    % testCIF_hetatoms       |-> include hetatoms  
    
    % testDownload_default      Test downloading (cif) structures from https://www.rcsb.org/
    % testDownload_chain        |-> specific chain 
    % testDownload_hetatoms     |-> include hetatoms  
    % testDownload_chain_auto   |-> specific chain assumed from 5th letter
    % testDownload_auto_Hatoms  |-> try download if file doesn't exist,
    %                               include hydrogen atoms 
    
    % testAlignment_default     Test performing an alignment 
    % testAlignment_funevals    |-> change stopping criterium 
    % testAlignment_alignlater  |-> hold alignment 
    
    
    

    % ----------------------------------------------
    
    % Test everything
    methods (Test)
        
    end
    
    
    % Test importing data 
    % --------------------------------------------------------------------
    
    % PDB
    methods (Test, TestTags = {'class','PDB', 'import'})
        
        function testPDB_default(testCase)
           % read PDB file containing 
            PDBstruct = PDB('sample_data/5mok.pdb'); 
            
           % assert default selection criteria expected for ZEAL
           testCase.assertFalse(PDBstruct.Selection.includeHatoms, 'HETATM should not be selected by default');
           
           testCase.assertFalse(PDBstruct.Selection.includeHetatoms, 'Hydrogen atoms should not be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.chainID,'all', 'All chains in file should be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.altLocID, 'all', 'All altlocs should be selected by default');
           
           % assert data contains expected number of atoms
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 6358, 'Number of ATOM records differ from expected');
            
        end
        
        function testPDB_chain(testCase)
            
             PDBstruct = PDB('sample_data/5mok.pdb', 'chainID', 'B');
             
             nAtoms = length(PDBstruct.Data.X);
             testCase.assertEqual(nAtoms, 1594, 'Number of ATOM records differ from expected');                         
            
        end
        
        function testPDB_hetatoms(testCase)
            
             PDBstruct = PDB('sample_data/5mok.pdb', 'includeHetatoms', true);
             
             nAtoms = length(PDBstruct.Data.X);
             testCase.assertEqual(nAtoms, 6358+815, 'Number of ATOM+HETATM records differ from expected');                         
            
        end
        
    end
    

    % CIF
    methods (Test, TestTags = {'class' 'CIF', 'import'})
        
        function testCIF_default(testCase)
           % read PDB file containing 
            PDBstruct = PDB('sample_data/5mok.cif'); 
            
           % assert default selection criteria expected for ZEAL
           testCase.assertFalse(PDBstruct.Selection.includeHatoms, 'HETATM should not be selected by default');
           
           testCase.assertFalse(PDBstruct.Selection.includeHetatoms, 'Hydrogen atoms should not be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.chainID,'all', 'All chains in file should be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.altLocID, 'all', 'All altlocs should be selected by default');
           
           % assert data contains expected number of atoms
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 6358, 'Number of ATOM records differ from expected');
            
        end
        
        function testCIF_chain(testCase)
            
             PDBstruct = PDB('sample_data/5mok.cif', 'chainID', 'B');
             
             nAtoms = length(PDBstruct.Data.X);
             testCase.assertEqual(nAtoms, 1594, 'Number of ATOM records differ from expected');                         
            
        end
        
        function testCIF_hetatoms(testCase)
            
             PDBstruct = PDB('sample_data/5mok.cif', 'includeHetatoms', true);
             
             nAtoms = length(PDBstruct.Data.X);
             testCase.assertEqual(nAtoms, 6358+815, 'Number of ATOM+HETATM records differ from expected');                         
            
        end
        
    end
        
    % Download
    methods (Test, TestTags = {'class', 'Download', 'import'})
        
        function testDownload_default(testCase)
           % read PDB file containing 
            PDBstruct = PDB('5mok'); 
            
           % assert default selection criteria expected for ZEAL
           testCase.assertFalse(PDBstruct.Selection.includeHatoms, 'HETATM should not be selected by default');
           
           testCase.assertFalse(PDBstruct.Selection.includeHetatoms, 'Hydrogen atoms should not be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.chainID,'all', 'All chains in file should be selected by default');
           
           testCase.assertEqual(PDBstruct.Selection.altLocID, 'all', 'All altlocs should be selected by default');
           
           % assert data contains expected number of atoms
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 6358, 'Number of ATOM records differ from expected');
            
        end
        
        function testDownload_chain(testCase)
            
             PDBstruct = PDB('sample_data/5mok', 'chainID', 'B');
             
             nAtoms = length(PDBstruct.Data.X);
             testCase.assertEqual(nAtoms, 1594, 'Number of ATOM records differ from expected');
             
        end
        
        function testDownload_chain_auto(testCase)
            
            PDBstruct = PDB('sample_data/5mokB');
            
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 1594, 'Number of ATOM records differ from expected');
            
        end
        
        function testDownload_auto_Hatoms(testCase)
            
            PDBstruct = PDB('sample_data/5PTI.pdb', 'includeHatoms', true);
            
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 909, 'Number of ATOM records differ from expected');
            
        end
        
        function testDownload_hetatoms(testCase)
            
            PDBstruct = PDB('sample_data/5mok', 'includeHetatoms', true);
            
            nAtoms = length(PDBstruct.Data.X);
            testCase.assertEqual(nAtoms, 6358+815, 'Number of ATOM+HETATM records differ from expected');
            
        end
        
    end

    
    
    % Test alignment 
    
    methods (Test, TestTags = {'class', 'Alignment'})
        
        
    end   
    
    

end

