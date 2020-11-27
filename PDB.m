classdef PDB < handle
    % PDB Class to handle protein structures stored in the PDB file format.
    %
    % PDB files can be downloaded from the PDB using PDB.fetch('pdbid') or
    % parsed using PDB.read('file.pdb').
    %
    % An instance of this class is done giving a filepath or 4-letter PDB code, in
    % which case the constructor will automatically download or read the
    % file. If a 5-letter PDB ID code is given, then the last letter is assumed to be the chain ID.
    %
    % The PDB data can be filtered using optional selection criteria;
    % defaults are specified in the Selection property section below.
    %
    % NB: The methods fetch and read are slight changes of the getpdb() and
    % pdbread() functions from Matlab's Bioinformatics Toolbox. The
    % methods are only in this class for convience - license for the toolbox is
    % thus required - Copyright notices are in the method descriptions of
    % PDB.read and PDB.fetch.
    %
    %
    % Filip (Persson) Ljung
    %
    % filip.persson@gmail.com
    %
    % Last update: 2020-10-20
    %
    % -------------------------------------------------------------
    
    properties
        Name % filename or 4-letter PDB id code or 4+1-letter code (PDB id code)+(chain id)
        Source % source meta data for file/downloaded
        
        Selection % Filter the PDB data with the following selection fields : type : (defaults)
        % selection.ChainID     :  string      :  ('A')
        % Selection.Hetatoms    :  true/false  :  (false = hetatom records omitted)
        % Selection.Hatoms      :  true/false  :  (false = hydrogen atoms omitted)
        % Selection.AltLoc      :  string      :  ('A' locations kept)
        
        Data      % structured array containing ATOM and/or HETATOM records
        
    end
    
    properties (Hidden)
        AllData  % contains all records 
        ShowLog
    end
    
    methods
        
        function obj = PDB(IDorFileName, varargin)
            % PDB Construct an instance of this class
            % Download or read pdb file given file/PDB id code
            %
            % INPUT
            % -------------------------------------------------------------
            % <IDorFileName>               :      string
            % filepath or a 4-letter PDBid code. If a 5-letter code is given,
            % then the last letter is assumed to be the chain id.
            %
            % OPTIONAL
            % 'name':value pairs  or as a struct as given below where 'name'=fieldname
            %
            % <selection>                  :       struct
            %
            % <selection.includeHetatoms>  :       true/false   (false)
            %   Flag to indicate if HETATOM records should be kept/omitted.
            %   If not defined then hetatoms are omitted.
            %
            % <selection.includeHatoms>    :       true/false   (false)
            %   Flag to indicate if Hydrogen atoms  should be kept/omitted
            %   If not defined then Hydrogen atoms omitted.
            %
            % <selection.chainID>          :       string       ('all')
            %   The id of the chain to keep. If not defined then all included.
            %
            % <selection.altLocID>         :       string       ('A') {Use 'all' to include all altlocs}
            %   The id of any alternate locations to keep. If not defined then all
            %   included
            %
            % <selection.ShowLog>           :   true/false (false)
            %
            % -------------------------------------------------------------
            
            % Parse input
            default_includeHetatoms = false;
            default_includeHatoms = false;
            default_chainID = 'all';
            default_altLocID = 'A';
            default_ShowLog = false;
            
            p = inputParser;
            
            % Set required
            addRequired(p, 'pdbID_or_FileName', @(x)validateattributes(x,{'char'}, {'nonempty'}, 'pdbID_or_FileName'));
            
            % Set optional
            addOptional(p, 'includeHetatoms', default_includeHetatoms);
            addOptional(p, 'includeHatoms', default_includeHatoms);
            addOptional(p, 'chainID', default_chainID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'chainID'));
            addOptional(p, 'altLocID', default_altLocID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'altLocID'));
            addOptional(p, 'ShowLog', default_ShowLog);
            
            parse(p, IDorFileName, varargin{:});
            
            % Setup properties
            obj.Name = p.Results.pdbID_or_FileName;
            
            obj.Selection.includeHetatoms = p.Results.includeHetatoms;
            obj.Selection.includeHatoms = p.Results.includeHatoms;
            obj.Selection.chainID = p.Results.chainID;
            obj.Selection.altLocID = p.Results.altLocID;
            
            obj.ShowLog = p.Results.ShowLog;
            
            if isempty(obj.Selection.chainID)
                obj.Selection.chainID = 'all';
            end
            
            % Import data
            if exist(obj.Name,'file')>0
                                
                [obj.Source.Path, obj.Source.Name, obj.Source.Ext] = fileparts(obj.Name);
                
                if obj.ShowLog
                    startTime = tic;
                    fprintf('\n\t Reading PDB file:\n%s', obj.Name);
                end
                
                % Check if file is PDB or CIF
                
                switch obj.Source.Ext
                    case '.pdb'
                        obj.AllData = PDB.readPDBfile(fullfile(obj.Source.Path, [obj.Source.Name obj.Source.Ext] ));
                    case '.cif'
                        obj.AllData = PDB.readCIFdata(fullfile(obj.Source.Path, [obj.Source.Name obj.Source.Ext] ));
                end
                               
                if obj.ShowLog
                    fprintf('\n\t Done. Execution time: %2.2e s', toc(startTime));
                end
                
                if obj.ShowLog
                    fprintf('\n\t Filtering PDB data in model %d', obj.Selection.modelNumber);
                    fprintf('\n\t\t Keeping chain with ID:\t %s', obj.Selection.chainID);
                    
                    if obj.Selection.includeHetatoms
                        fprintf('\n\t\t Keeping HETATOM records  ');
                    else
                        fprintf('\n\t\t Leaving out HETATOM records ');
                    end
                    
                    if obj.Selection.includeHatoms
                        fprintf('\n\t\t Keeping Hydrogen atoms (if exist)');
                    else
                        fprintf('\n\t\t Leaving out Hydrogen atoms (if exist)');
                    end
                    
                    fprintf('\n\t\t Keeping alt locations with ID:\t %s', obj.Selection.altLocID);
                    
                end
                
                obj.Data = PDB.parsePDBstruct(obj.AllData, obj.Selection);
                
            else % download
                try
                    [~,name,~] = fileparts(obj.Name);
                    pdbid = name(1:4);
                    
                    queryStr = sprintf('https://models.rcsb.org/v1/%s/atoms?encoding=cif&copy_all_categories=false', pdbid);
                    cifDataPulled = webread(queryStr);
                                        
                    obj.AllData = PDB.readCIFdata(cifDataPulled);
                    
                    if numel(name)>4 % assume last letter specifies chain
                        obj.Selection.chainID = name(5);
                    end
                    
                    obj.Data = PDB.parsePDBstruct(obj.AllData, obj.Selection);
                    
                    obj.Source.Path = queryStr;
                    obj.Source.Name = pdbid;
                    obj.Source.Ext = 'cif';
                    
                catch ME
                    
                    fprintf('%s', ME.message);
                    error('\n Could not download structure %s in the RCSB PDB.', pdbid);
                    
                end
                
            end
            
        end
        
        function reparsePDB(obj)
            % Parse the original Matlab PDB struct again.
            % Run if any selection property is changed to update the PDB selection
            % data in obj.Data
            
            obj.Data = PDB.parsePDBstruct(obj.AllData, obj.Selection);
                        
        end
        
    end
    
    methods (Static)
        
        function PDBdata = readCIFdata(cifSource)
           
            % check if file
            if exist(cifSource,'file')               
                textLines = splitlines( fileread('1stm.cif') );                
            else                                
                textLines = splitlines(cifSource);                
            end
                    
            % get _atom_site labels for coordinate data             
            firstHeader = find(contains(textLines, '_atom_site.group_PDB'));
            lastHeader = find(contains(textLines, '_atom_site.pdbx_PDB_model_num'));
         
            % get ATOM and HETATM records
            sel = startsWith(textLines, {'ATOM', 'HETATM'});            
            dataLines = textLines(sel);
            
            % Create structure with dynamics fieldnames corresponding to
            % _atom_site labels the cif file 
            pdbDataCif = struct;
            for n=firstHeader:lastHeader
               
                headerName = split(textLines{n}, '.');
                pdbDataCif.(strip(headerName{2})) = cell(1,1);
                
            end
            
            % pre allocate cell sizes in struct            
            pdbDataFieldnames = fieldnames(pdbDataCif);
            nLines = numel(dataLines);            
            
            for n = 1:numel(pdbDataFieldnames)
                    pdbDataCif.(pdbDataFieldnames{n}) = cell(nLines,1);                
            end
                        
            % parse all record lines and assign to struct
            for n = 1:numel(dataLines)
                
                n_line = split(dataLines{n});
                
                pdbDataCif.(pdbDataFieldnames{1}){n} =  n_line{1};                
                pdbDataCif.(pdbDataFieldnames{2}){n} =  n_line{2};                
                pdbDataCif.(pdbDataFieldnames{3}){n} =  n_line{3};                
                pdbDataCif.(pdbDataFieldnames{4}){n} =  n_line{4};                
                pdbDataCif.(pdbDataFieldnames{5}){n} =  n_line{5};                
                pdbDataCif.(pdbDataFieldnames{6}){n} =  n_line{6};                
                pdbDataCif.(pdbDataFieldnames{7}){n} =  n_line{7};                
                pdbDataCif.(pdbDataFieldnames{8}){n} =  n_line{8};                
                pdbDataCif.(pdbDataFieldnames{9}){n} =  n_line{9};                
                pdbDataCif.(pdbDataFieldnames{10}){n} =  n_line{10};                
                pdbDataCif.(pdbDataFieldnames{11}){n} =  n_line{11};
                pdbDataCif.(pdbDataFieldnames{12}){n} =  n_line{12};
                pdbDataCif.(pdbDataFieldnames{13}){n} =  n_line{13};
                pdbDataCif.(pdbDataFieldnames{14}){n} =  n_line{14};
                pdbDataCif.(pdbDataFieldnames{15}){n} =  n_line{15};
                pdbDataCif.(pdbDataFieldnames{16}){n} =  n_line{16};
                pdbDataCif.(pdbDataFieldnames{17}){n} =  n_line{17};
                pdbDataCif.(pdbDataFieldnames{18}){n} =  n_line{18};
                pdbDataCif.(pdbDataFieldnames{19}){n} =  n_line{19};
                pdbDataCif.(pdbDataFieldnames{20}){n} =  n_line{20};
                pdbDataCif.(pdbDataFieldnames{21}){n} =  n_line{21};
                
            end
            
            % Create new struct to conform to legacy code
            PDBdata = struct;
            PDBdata.recordName = pdbDataCif.('group_PDB');
            PDBdata.atomNum = cellfun(@str2double, pdbDataCif.('id'));
            PDBdata.atomName = pdbDataCif.('label_atom_id');
          
            PDBdata.altLoc = pdbDataCif.('label_alt_id');
            PDBdata.resName = pdbDataCif.('label_comp_id');
           
            PDBdata.chainID = pdbDataCif.('label_asym_id');
            PDBdata.resNum = cellfun(@str2double, pdbDataCif.('label_seq_id'));
            PDBdata.X = cellfun(@str2double, pdbDataCif.('Cartn_x'));
            PDBdata.Y = cellfun(@str2double, pdbDataCif.('Cartn_y'));
            PDBdata.Z = cellfun(@str2double, pdbDataCif.('Cartn_z'));
            
            PDBdata.occupancy = cellfun(@str2double, pdbDataCif.('occupancy'));
            PDBdata.betaFactor = cellfun(@str2double, pdbDataCif.('B_iso_or_equiv'));
            PDBdata.element = pdbDataCif.('type_symbol');
            PDBdata.charge = cellfun(@str2double, pdbDataCif.('pdbx_formal_charge'));
            
        end
        
        function PDBdata = parsePDBstruct(pdbStruct, varargin)
            % Create a filtered struct from PDB data in (hidden) AllData
            % property
            %
            % INPUT
            % -------------------------------------------------------------
            % <pdbStruct>                 :       struct
            % the output from PDB.readCIFdata / PDB.readFile
            %
            % OPTIONAL
            % <selection>                 :       struct
            %
            % <selection.includeHetatoms> :       true/false
            %   Flag to indicate if HETATOM records should be kept/omitted.
            %   If not defined then hetatoms are omitted.
            %
            % <selection.includeHatoms>   :       true/false
            %   Flag to indicate if Hydrogen atoms  should be kept/omitted
            %   If not defined then Hydrogen atoms omitted.
            %
            % <selection.chainID>      :       string
            %   The id of the chain to keep. If not defined then all included.
            %
            % <selection.altLocID>       :       string
            %   The id of any alternate locations to keep. If not defined then all
            %   included
            %
            %
            % OUTPUT
            % -------------------------------------------------------------
            % PDBdata    :   Structured array with fields
            %
            % PDBdata.recordName
            % PDBdata.atomNum
            % PDBdata.atomName
            % PDBdata.altLoc
            % PDBdata.resName
            %
            % PDBdata.chainID
            % PDBdata.resNum
            % PDBdata.X
            % PDBdata.Y
            % PDBdata.Z
            
            % PDBdata.occupancy
            % PDBdata.betaFactor
            % PDBdata.element
            % PDBdata.charge
            %
            % -------------------------------------------------------------
            
            if ~isstruct(pdbStruct)
                error('\n Input must be a structure as created by the PDB.readCIFdata or PDB.readFile methods');
            end
            
            
            % parse input and setup selection filter varaibles
            default_includeHetatoms = false;
            default_includeHatoms = false;
            default_chainID = 'all';
            default_altLocID = 'all';
           
            
            p = inputParser;
            
            addRequired(p, 'pdbStruct');
            
            addOptional(p, 'includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric','logical'}, {'nonempty'}, 'includeHetatoms'));
            addOptional(p, 'includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric','logical'}, {'nonempty'}, 'includeHatoms'));
            addOptional(p, 'chainID', default_chainID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'chainID'));
            addOptional(p, 'altLocID', default_altLocID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'altLocID'));
            
            parse(p, pdbStruct, varargin{:});
            
            includeHetatoms = p.Results.includeHetatoms;
            includeHatoms = p.Results.includeHatoms;
            chainIDsel = p.Results.chainID;
            altLocID = p.Results.altLocID;
            
            if strcmp(chainIDsel, 'all')
                chainOp = false;
            else
                chainOp = true;
            end
            
            if strcmp(altLocID, 'all')
                altLocOp = false;
            else
                altLocOp = true;
            end
            
            
            % filter based on selection options: hetatoms, H-atoms, chain
            % and altlocs
            if includeHetatoms
                keepListTF = contains(pdbStruct.recordName, {'ATOM','HETATM'});
            else
                keepListTF = strcmp(pdbStruct.recordName, 'ATOM');
            end
            
            if ~includeHatoms
                keepHatomsTF = (strcmp(pdbStruct.element, 'H')) == 0;
                keepListTF = (keepListTF + keepHatomsTF) == 2;
            end
            
            if chainOp
                keepChainTF = strcmp(pdbStruct.chainID, chainIDsel);
                keepListTF = (keepListTF + keepChainTF) == 2;
            end
            
            if altLocOp
                keepAltLocTF = contains(pdbStruct.altLoc, {altLocID,''});
                keepListTF = (keepListTF + keepAltLocTF) == 2;
            end
            
            PDBdata.recordName = pdbStruct.recordName(keepListTF);
            PDBdata.atomNum    = pdbStruct.atomNum(keepListTF);
            PDBdata.atomName   = pdbStruct.atomName(keepListTF);
            PDBdata.altLoc     = pdbStruct.altLoc(keepListTF);
            PDBdata.resName    = pdbStruct.resName(keepListTF);
            
            PDBdata.chainID    = pdbStruct.chainID(keepListTF);
            PDBdata.resNum     = pdbStruct.resNum(keepListTF);
            PDBdata.X          = pdbStruct.X(keepListTF);
            PDBdata.Y          = pdbStruct.Y(keepListTF);
            PDBdata.Z          = pdbStruct.Z(keepListTF);
            
            PDBdata.occupancy  = pdbStruct.occupancy(keepListTF);
            PDBdata.betaFactor = pdbStruct.betaFactor(keepListTF);
            PDBdata.element    = pdbStruct.element(keepListTF);
            PDBdata.charge     = pdbStruct.charge(keepListTF);
            
        end
        
        
        
        
        function saveToFile(pdbData)
            % PDB.saveToFile save data to file in the PDB file format
            % Adopted from Evan (2020): read and write PDB files using matlab
            % (https://www.mathworks.com/matlabcentral/fileexchange/42957-read-and-write-pdb-files-using-matlab),
            % MATLAB Central File Exchange. Retrieved October 19, 2020.
            %
            % REQUIRED INPUT
            % -------------------------------------------------------------
            % pdbData.X            X coordinate data
            % pdbData.Y            Y coordinate data
            % pdbData.Z            Z coordinate data
            % pdbData.outfile      output file name
            
            % OPTIONAL INPUT
            % -------------------------------------------------------------
            % pdbData value        meaning                           default value
            %
            % pdbData.recordName   output record name of atoms      "ATOM"
            % pdbData.atomNum      atom serial number                sequential number
            % pdbData.atomName     name of atoms                    "OW" (water oxygen)
            % pdbData.altLoc       alt. location indicator          " "
            % pdbData.resName      name of residue                  "SOL" (water)
            %
            % pdbData.chainID      protein chain identifier         "A"
            % pdbData.resNum       residue sequence number           sequential number
            % pdbData.occupancy    occupancy factor                 "1.00"
            % pdbData.betaFactor   beta factor, temperature         "0.00"
            % pdbData.element      element symbol                   "O" (oxygen)
            % pdbData.charge       atomic charge                    " "
            %
            % -------------------------------------------------------------
            
            % coordinate data is required! Checking XYZ input
            if ~isfield(pdbData, 'X') || ~isfield(pdbData, 'Y') || ~isfield(pdbData, 'Z')
                error('Field(s) for XYZ coordinates not found.');
            end
            X = pdbData.X;
            Y = pdbData.Y;
            Z = pdbData.Z;
            if length(X) ~= length(Y) || length(X) ~= length(Z)
                error('XYZ coordinates not of equal dimension.');
            end
            
            % review optional data inputs
            % in case optional data  not given, fill in blanks
            if ~isfield(pdbData, 'outfile')
                pdbData.outfile = 'mat2PDB.pdb';
            end
            if ~isfield(pdbData, 'recordName')
                pdbData.recordName = cell(1,length(X));
                pdbData.recordName(1:end) = {'ATOM'};
            end
            if ~isfield(pdbData, 'atomNum')
                pdbData.atomNum = 1:length(X);
            end
            if ~isfield(pdbData, 'atomName')
                pdbData.atomName = cell(1,length(X));
                pdbData.atomName(1:end) = {'OW'};
            end
            if ~isfield(pdbData, 'altLoc')
                pdbData.altLoc = cell(1,length(X));
                pdbData.altLoc(1:end) = {' '};
            end
            if ~isfield(pdbData, 'resName')
                pdbData.resName = cell(1,length(X));
                pdbData.resName(1:end) = {'SOL'};
            end
            if ~isfield(pdbData, 'chainID')
                pdbData.chainID = cell(1,length(X));
                pdbData.chainID(1:end) = {'A'};
            end
            if ~isfield(pdbData, 'resNum')
                pdbData.resNum = 1:length(X);
            end
            if ~isfield(pdbData, 'occupancy')
                pdbData.occupancy = ones(1,length(X));
            end
            if ~isfield(pdbData, 'betaFactor')
                pdbData.betaFactor = zeros(1, length(X));
            end
            if ~isfield(pdbData, 'element')
                pdbData.element = cell(1,length(X));
                pdbData.element(1:end) = {'O'};
            end
            if ~isfield(pdbData, 'charge')
                pdbData.charge = cell(1,length(X));
                pdbData.charge(1:end) = {' '};
            end
            
            outfile    = pdbData.outfile;
            recordName = pdbData.recordName;
            atomNum    = pdbData.atomNum;
            atomName   = pdbData.atomName;
            altLoc     = pdbData.altLoc;
            resName    = pdbData.resName;
            chainID    = pdbData.chainID;
            resNum     = abs(pdbData.resNum);
            occupancy  = pdbData.occupancy;
            betaFactor = pdbData.betaFactor;
            element    = pdbData.element;
            charge     = pdbData.charge;
            
            % remove faulty inputs
            if length(recordName) ~= length(X)
                warning('recordName input is not the correct length!\n\tignoring user input\n');
                recordName = cell(1,length(X));
                recordName(1:end) = {'ATOM'};
            end
            if length(atomNum) ~= length(X)
                warning('atom serial number input is not the correct length!\n\tignoring user input\n');
                atomNum = 1:length(X);
            end
            if length(atomName) ~= length(X)
                warning('atom name input is not the correct length!\n\tignoring user input\n');
                atomName = cell(1,length(X));
                atomName(1:end) = {'OW'};
            end
            if length(altLoc) ~= length(X)
                warning('alternate location input is not the correct length!\n\tignoring user input\n');
                altLoc = cell(1,length(X));
                altLoc(1:end) = {' '};
            end
            if length(resName) ~= length(X)
                warning('residue name input is not the correct length!\n\tignoring user input\n');
                resName = cell(1,length(X));
                resName(1:end) = {'SOL'};
            end
            if length(chainID) ~= length(X)
                warning('chain ID input is not the correct length!\n\tignoring user input\n');
                chainID = cell(1,length(X));
                chainID(1:end) = {'A'};
            end
            if length(resNum) ~= length(X)
                warning('residue number input is not the correct length!\n\tignoring user input\n');
                resNum = 1:length(X);
            end
            if length(occupancy) ~= length(X)
                warning('occupancy input is not the correct length!\n\tignoring user input\n');
                occupancy = ones(1,length(X));
            end
            if length(betaFactor) ~= length(X)
                warning('beta factor input is not the correct length!\n\tignoring user input\n');
                betaFactor = zeros(1, length(X));
            end
            if length(element) ~= length(X)
                warning('element symbol input is not the correct length!\n\tignoring user input\n');
                element = cell(1,length(X));
                element(1:end) = {'O'};
            end
            if length(charge) ~= length(X)
                warning('charge input is not the correct length!\n\tignoring user input\n');
                charge = cell(1,length(X));
                charge(1:end) = {' '};
            end
            
            % fix atomName spacing
            for n = 1:length(atomName)
                atomName(n) = {sprintf('%-3s',cell2mat(atomName(n)))};
            end
            
            
            % create PDB
            
            % open file
            fprintf('outputting PDB in file %s\n', outfile);
            fid = fopen(outfile, 'w');
            if fid == -1
                
                error('Author:Function:OpenFile', 'Cannot open file: %s', outfile);
            end
            
            % output data
            for n = 1:length(atomNum)
                
                % standard PDB output line
                fprintf( fid, '%-6s%5u%5s%1.1s%3s %1.1s%4u%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
                    cell2mat(recordName(n)), atomNum(n), cell2mat(atomName(n)), ...
                    cell2mat(altLoc(n)), cell2mat(resName(n)), cell2mat(chainID(n)), ...
                    resNum(n), X(n), Y(n), Z(n), occupancy(n), betaFactor(n), ...
                    cell2mat(element(n)), cell2mat(charge(n)));
                
                % output progress in terminal
                if ~mod(n,400)
                    fprintf('   %6.2f%%', 100*n / length(atomNum));
                    if ~mod(n, 4000)
                        fprintf('\n');
                    end
                end
                
            end
            fprintf( fid, 'END\n');
            
            % close file
            fprintf('   %6.2f%%\n    done! closing file...\n', 100);
            
            fclose(fid);
            
        end
        
        
        function [PDBdata] = readPDBfile(readFile)
            %  -- readFromPDB.m --
            % This program is the most speedy way to read a PDB file that I could come
            % up with. It's function is simple: give it a PDB file and out comes a
            % matlab-friendly data structure. In cumbersomely large PDB's (such as those that
            % include solvent), this can shave off a good amount of time relative to
            % many programs. Unfortunately there is no easy way to hasten the slowest
            % step, which is turning strings into doubles.
            %
            % The output format is as given in online documentation
            % (as of July 2012 when writing this program)
            % http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
            %
            % It outputs 14 pieces total of information about the PDB.
            %
            % -- mandatory information (11) --
            %
            % outfile    (the name of the PDB, this is the only input on the command line)
            %
            % recordName (the class or type of atom, such as ATOM, HETATM, SOL, etc)
            % atomNum    (serial number of the atom)
            % atomName   (elemental identification of the atom)
            % altLoc     (alt. location indicator)
            % resName    (name of the amino acid/residue)
            %
            % chainID    (protein chain identifier)
            % resNum     (index number of the amino acid)
            % X          (X position of atom)
            % Y          (Y position of atom)
            % Z          (Z position of atom)
            %                      
            %  Evan (2020). read and write PDB files using matlab (https://www.mathworks.com/matlabcentral/fileexchange/42957-read-and-write-pdb-files-using-matlab),
            % MATLAB Central File Exchange. Retrieved November 5, 2020.
            
            PDBdata.outfile = readFile;
            % initialize file
            FileID = fopen(readFile);
            rawText = fread(FileID,inf,'*char');
            % parse lines by end-of-lines
            splitLines = strread(rawText, '%s', 'delimiter', '\n');
            % initialize variables
            numLines = length(splitLines);
            recordName = cell(1,numLines);
            atomNum    = cell(1,numLines);
            atomName   = cell(1,numLines);
            altLoc     = cell(1,numLines);
            resName    = cell(1,numLines);
            chainID    = cell(1,numLines);
            resNum     = cell(1,numLines);
            X          = cell(1,numLines);
            Y          = cell(1,numLines);
            Z          = cell(1,numLines);
            comment    = cell(1,numLines);
            % read each line
            m = 1;
            for n = 1:numLines
                
                thisLine = cell2mat(splitLines(n));
                
                if length(thisLine) > 53 && sum(isstrprop(thisLine(23:53), 'alpha')) == 0
                    
                    recordName(m) = {thisLine(1:6)};
                    atomNum(m)    = {thisLine(7:11)};
                    atomName(m)   = {thisLine(13:16)};
                    altLoc(m)     = {thisLine(17)};
                    resName(m)    = {thisLine(18:20)};
                    
                    chainID(m)    = {thisLine(22)};
                    resNum(m)     = {thisLine(23:26)};
                    X(m)          = {thisLine(31:38)};
                    Y(m)          = {thisLine(39:46)};
                    Z(m)          = {thisLine(47:54)};
                    
                    comment(m)            = {thisLine(55:end)};
                    
                    m = m + 1;
                end
                
            end
            % trim exess
            keepData = logical(strcmp(recordName,'ATOM  ') + strcmp(recordName,'HETATM'));
            recordName = recordName(keepData);
            atomNum    = atomNum(keepData);
            atomName   = atomName(keepData);
            altLoc     = altLoc(keepData);
            resName    = resName(keepData);
            chainID    = chainID(keepData);
            resNum     = resNum(keepData);
            X          = X(keepData);
            Y          = Y(keepData);
            Z          = Z(keepData);
            comment    = comment(keepData);
            % parse out "comment" section
            occupancy  = cell(1, length(recordName));
            betaFactor = cell(1, length(recordName));
            element    = cell(1, length(recordName));
            charge     = cell(1, length(recordName));
            % fix spacing
            for n = 1:length(recordName)
                thisLine = sprintf('%-26s',cell2mat(comment(n)));
                occupancy(n)  = {thisLine(1:6)};
                betaFactor(n) = {thisLine(7:12)};
                element(n)    = {thisLine(13:24)};
                charge(n)     = {thisLine(25:26)};
            end
            % reformat data for convenience
            PDBdata.recordName = (strtrim(recordName))';
            PDBdata.atomNum    = (str2double(atomNum))';
            PDBdata.atomName   = (strtrim(atomName))';
            PDBdata.altLoc     = (altLoc)';
            PDBdata.resName    = (strtrim(resName))';
            PDBdata.chainID    = (chainID)';
            PDBdata.resNum     = (str2double(resNum))';
            PDBdata.X          = (str2double(X))';
            PDBdata.Y          = (str2double(Y))';
            PDBdata.Z          = (str2double(Z))';
            PDBdata.occupancy  = (str2double(occupancy))';
            PDBdata.betaFactor = (str2double(betaFactor))';
            PDBdata.element    = (strtrim(element))';
            PDBdata.charge     = (strtrim(charge))';
            % I commented these lines out, since they cause more problems than they
            % solve. They do clean up the output for certain situations.
            % if isnan(PDBdata.occupancy(1))
            %     PDBdata.occupancy = strtrim(PDBdata.occupancy);
            % end
            % if isnan(PDBdata.betaFactor(1))
            %     PDBdata.occupancy = strtrim(PDBdata.betaFactor);
            % end
            % close file
            fclose(FileID);
            
        end        
    end
end

