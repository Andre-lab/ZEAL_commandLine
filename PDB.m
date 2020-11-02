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
        Name % filepath or 4-letter PDB id code or 4+1-letter code (PDB id code)+(chain id)
        
        Selection % Filter the PDB data with the following selection fields : type : (defaults)
                  % selection.ChainID     :  string      :  ('A')
                  % Selection.Hetatoms    :  true/false  :  (false = hetatom records omitted)
                  % Selection.Hatoms      :  true/false  :  (false = hydrogen atoms omitted)      
                  % Selection.AltLoc      :  string      :  ('A' locations kept)
                  % Selection.ModelNumber :  integer     :  (1)
                  
        Data      % structured array containing ATOM and/or HETATOM records      
        
    end
    
    properties (Hidden)
        MatlabData  % contains all records as processed by the PDB.read and PDB.fetch methods      
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
            % <selection.modelNumber>          :       integer    (1)   
            %   The model number to select. If not defined then the first
            %   model is selected
            %
            % <selection.ShowLog>           :   true/false (false)
            %
            % -------------------------------------------------------------
            
            % Parse input
            default_includeHetatoms = false;
            default_includeHatoms = false;            
            default_chainID = 'all';
            default_altLocID = 'A';
            default_modelID = 1;
            default_ShowLog = false;
                        
            p = inputParser;
            
            % Set required 
            addRequired(p, 'pdbID_or_FileName', @(x)validateattributes(x,{'char'}, {'nonempty'}, 'pdbID_or_FileName'));
            
            % Set optional
            addOptional(p, 'includeHetatoms', default_includeHetatoms);
            addOptional(p, 'includeHatoms', default_includeHatoms);            
            addOptional(p, 'chainID', default_chainID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'chainID'));
            addOptional(p, 'altLocID', default_altLocID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'altLocID'));
            addOptional(p, 'modelNumber', default_modelID);
            addOptional(p, 'ShowLog', default_ShowLog);         
            
            parse(p, IDorFileName, varargin{:});
            
            % Setup properties
            obj.Name = p.Results.pdbID_or_FileName;
            
            obj.Selection.includeHetatoms = p.Results.includeHetatoms;
            obj.Selection.includeHatoms = p.Results.includeHatoms;
            obj.Selection.chainID = p.Results.chainID;
            obj.Selection.altLocID = p.Results.altLocID;
            obj.Selection.modelNumber = p.Results.modelNumber;
            
            obj.ShowLog = p.Results.ShowLog;
            
            if isempty(obj.Selection.chainID)
                obj.Selection.chainID = 'all';
            end
            
            % Import data
            if exist(obj.Name,'file')>0
                            
                if obj.ShowLog
                    startTime = tic;
                    fprintf('\n\t Reading PDB file:\n%s', obj.Name);
                end
                
                obj.MatlabData = PDB.read(obj.Name);
                
                if obj.ShowLog
                    fprintf('\n\t Done. Execution time: %2.2e s', toc(startTime));
                end
                
                % [path, name, extension] = fileparts(obj.Name);
                [~, name, ~] = fileparts(obj.Name);
                
                obj.Name = name;
                
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
                
                obj.Data = PDB.parseMatlabPDBstruct(obj.MatlabData, obj.Selection);
                
            else % download
                try
                    [~,name,~] = fileparts(obj.Name);
                    pdbid = name(1:4);
                    
                    obj.MatlabData = PDB.fetch(pdbid);
                    
                    if numel(name)>4 % assume last letter specifies chain
                        obj.Selection.chainID = name(5);
                    end
                    
                    obj.Data = PDB.parseMatlabPDBstruct(obj.MatlabData, obj.Selection);
                    
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
           
           obj.Data = PDB.parseMatlabPDBstruct(obj.MatlabData, obj.Selection);
            
            
        end
        
    end
    
    methods (Static)                
        
        function PDBdata = parseMatlabPDBstruct(pdbStruct, varargin)
            % Create a simpler (filtered) struct for PDB data from Matlabs PDB object
            %
            % INPUT
            % -------------------------------------------------------------
            % <pdbStruct>                 :       struct      
            % the output from ZC.fetch/ ZC.read            
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
            % <selection.modelNumber>           :       integer       
            %   The model number to select. If not defined then the first
            %   model is selected
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
                error('\n Input must be a structure as created by the PDB.read or PDB.fetch methods of this class');                
            end
            
            
            % parse input and setup selection filter varaibles
            default_includeHetatoms = false;
            default_includeHatoms = false;            
            default_chainID = 'all';
            default_altLocID = 'all';
            default_modelID = 1;
                        
            p = inputParser;                    
            
            addRequired(p, 'pdbStruct');
            
            addOptional(p, 'includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric','logical'}, {'nonempty'}, 'includeHetatoms')); 
            addOptional(p, 'includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric','logical'}, {'nonempty'}, 'includeHatoms'));          
            addOptional(p, 'chainID', default_chainID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'chainID'));
            addOptional(p, 'altLocID', default_altLocID, @(x)validateattributes(x,{'char'}, {'nonempty'}, 'altLocID'));
            addOptional(p, 'modelNumber',  default_modelID, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','positive'}, 'modelNumber'));
             
            parse(p, pdbStruct, varargin{:});
            
            includeHetatoms = p.Results.includeHetatoms;
            includeHatoms = p.Results.includeHatoms;
            chainIDsel = p.Results.chainID;
            altLocID = p.Results.altLocID;
            modelNumber = p.Results.modelNumber;
                        
            nAtoms = size(pdbStruct.Model(1).Atom, 2);
                  
            % Set selection switches
            if isfield(pdbStruct.Model(1), 'HeterogenAtom')
                nHetAtoms = size(pdbStruct.Model(1).HeterogenAtom, 2);
            else
                nHetAtoms = 0;
            end
            
            if includeHetatoms                
                numLines = nAtoms + nHetAtoms;
            else
                numLines = nAtoms;
            end
            
            
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
            
            
            % pre-allocate cell arrays
            recordName = cell(1,numLines);
            atomNum    = zeros(1,numLines);
            atomName   = cell(1,numLines);
            altLoc     = cell(1,numLines);
            resName    = cell(1,numLines);
            
            chainID    = cell(1,numLines);
            resNum     = zeros(1,numLines);
            X          = zeros(1,numLines);
            Y          = zeros(1,numLines);
            Z          = zeros(1,numLines);
            
            occupancy  = zeros(1, numLines);
            betaFactor = zeros(1, numLines);
            element    = cell(1, numLines);
            charge     = cell(1, numLines);
            
            % get atom records
            m=1;
            for n = 1:nAtoms
                
                recordName(m) = {'ATOM  '};
                atomNum(m)    = pdbStruct.Model(modelNumber).Atom(n).AtomSerNo;
                atomName(m)   = {pdbStruct.Model(modelNumber).Atom(n).AtomName};
                altLoc(m)     = {pdbStruct.Model(modelNumber).Atom(n).altLoc};
                resName(m)    = {pdbStruct.Model(modelNumber).Atom(n).resName};
                
                chainID(m)    = {pdbStruct.Model(modelNumber).Atom(n).chainID};
                resNum(m)     = pdbStruct.Model(modelNumber).Atom(n).resSeq;
                X(m)          = pdbStruct.Model(modelNumber).Atom(n).X;
                Y(m)          = pdbStruct.Model(modelNumber).Atom(n).Y;
                Z(m)          = pdbStruct.Model(modelNumber).Atom(n).Z;
                
                occupancy(m)  = pdbStruct.Model(modelNumber).Atom(n).occupancy;
                betaFactor(m) = pdbStruct.Model(modelNumber).Atom(n).tempFactor;
                element(m)    = {pdbStruct.Model(modelNumber).Atom(n).element};
                charge(m)    = {pdbStruct.Model(modelNumber).Atom(n).charge};
                
                m = m + 1;
            end
            
            % get hetatom records
            if includeHetatoms
                for n = 1:nHetAtoms
                    
                    recordName(m) = {'HETATM'};
                    atomNum(m)    = pdbStruct.Model(modelNumber).HeterogenAtom(n).AtomSerNo;
                    atomName(m)   = {pdbStruct.Model(modelNumber).HeterogenAtom(n).AtomName};
                    altLoc(m)     = {pdbStruct.Model(modelNumber).HeterogenAtom(n).altLoc};
                    resName(m)    = {pdbStruct.Model(modelNumber).HeterogenAtom(n).resName};
                    
                    chainID(m)    = {pdbStruct.Model(modelNumber).HeterogenAtom(n).chainID};
                    resNum(m)     = pdbStruct.Model(modelNumber).HeterogenAtom(n).resSeq;
                    X(m)          = pdbStruct.Model(modelNumber).HeterogenAtom(n).X;
                    Y(m)          = pdbStruct.Model(modelNumber).HeterogenAtom(n).Y;
                    Z(m)          = pdbStruct.Model(modelNumber).HeterogenAtom(n).Z;
                    
                    occupancy(m)  = pdbStruct.Model(modelNumber).HeterogenAtom(n).occupancy;
                    betaFactor(m) = pdbStruct.Model(modelNumber).HeterogenAtom(n).tempFactor;
                    element(m)    = {pdbStruct.Model(modelNumber).HeterogenAtom(n).element};
                    charge(m)    = {pdbStruct.Model(modelNumber).HeterogenAtom(n).charge};
                    
                    m = m + 1;
                    
                end
            end
            
            % filter based on selection options: hetatoms, H-atoms, chain
            % and altlocs
            keepListTF = strcmp(recordName, 'ATOM  ');
            
            if includeHetatoms
                keepHetatomTF = strcmp(recordName, 'HETATM');
                keepListTF = (keepListTF + keepHetatomTF) == 2;
            end
            
            if ~includeHatoms
                keepHatomsTF = (strcmp(element, 'H')) == 0;
                keepListTF = (keepListTF + keepHatomsTF) == 2;
            end
            
            if chainOp
                keepChainTF = strcmp(chainID, chainIDsel);
                keepListTF = (keepListTF + keepChainTF) == 2;
            end
            
            if altLocOp
                keepAltLocTF = contains(altLoc, {altLocID,''});
                keepListTF = (keepListTF + keepAltLocTF) == 2;
            end
                        
            PDBdata.recordName = recordName(keepListTF);
            PDBdata.atomNum    = atomNum(keepListTF);
            PDBdata.atomName   = atomName(keepListTF);
            PDBdata.altLoc     = altLoc(keepListTF);
            PDBdata.resName    = resName(keepListTF);
            
            PDBdata.chainID    = chainID(keepListTF);
            PDBdata.resNum     = resNum(keepListTF);
            PDBdata.X          = X(keepListTF);
            PDBdata.Y          = Y(keepListTF);
            PDBdata.Z          = Z(keepListTF);
            
            PDBdata.occupancy  = occupancy(keepListTF);
            PDBdata.betaFactor = betaFactor(keepListTF);
            PDBdata.element    = element(keepListTF);
            PDBdata.charge     = charge(keepListTF);
            
        end
        
        function PDBstruct=fetch(pdbID,varargin)
            % fetch retrieves sequence information from the Protein Data Bank.
            %
            %   PDBSTRUCT = GETPDB(PDBID) searches for the ID in the Protein Data Bank
            %   (PDB) database and returns a structure containing information for the
            %   protein.
            %
            %   PDBSTRUCT = GETPDB(...,'TOFILE',FILENAME) saves the data returned from
            %   the database in the file FILENAME.
            %
            %   PDBSTRUCT = GETPDB(...,'SEQUENCEONLY',true) returns just the protein
            %   sequence. If the PDB file contains only one sequence then this will be
            %   returned as a character array. If more than one sequence is found, then
            %   these will be returned in a cell array.
            %
            %   Example:
            %
            %       pdbstruct = getpdb('2DHB')
            %
            %   This retrieves the structure information for horse deoxyhemoglobin
            %   (PDB ID 2DHB).
            %
            %   See also GETEMBL, GETGENBANK, GETGENPEPT, MOLVIEWER, PDBDISTPLOT,
            %   PDBREAD.
            
            %   Reference:
            %   H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H.
            %   Weissig, I.N. Shindyalov, P.E. Bourne: The Protein Data Bank. Nucleic
            %   Acids Research, 28 pp. 235-242 (2000)
            
            %   Copyright 2002-2012 The MathWorks, Inc.
            %
            
            %   Mirror sites are no longer supported after Jan 1st 2006.
            %   PDBSTRUCT = GETPDB(...,'MIRRORSITE',MIRROR) allows you to choose a
            %   mirror site for the PDB database. The default site is the San Diego
            %   Supercomputer Center, http://www.rcsb.org/pdb. Set MIRROR to
            %   http://rutgers.rcsb.org/pdb to use the Rutgers University Site or
            %   http://nist.rcsb.org/pdb for the National Institute of Standards and
            %   Technology site. Follow this link for a full list of PDB mirror sites:
            %   http://www.rcsb.org/pdb/static.do?p=general_information/mirror_sites/in
            %   dex.html
                        
            if nargin > 0
                pdbID = convertStringsToChars(pdbID);
            end
            
            if ~usejava('jvm')
                error(message('bioinfo:getpdb:NeedJVM', mfilename));
            end
            
            tofile = false;
            seqonly = false;
            PDBsite = 'http://www.rcsb.org/pdb';
            
            if nargin > 1
                if rem(nargin,2) == 0
                    error(message('bioinfo:getpdb:IncorrectNumberOfArguments', mfilename));
                end
                okargs = {'tofile','mirrorsite','sequenceonly'};
                for j=1:2:nargin-2
                    pname = varargin{j};
                    pval = varargin{j+1};
                    k = find(strncmpi(pname, okargs,length(pname)));
                    if isempty(k)
                        error(message('bioinfo:getpdb:UnknownParameterName', pname));
                    elseif length(k)>1
                        error(message('bioinfo:getpdb:AmbiguousParameterName', pname));
                    else
                        switch(k)
                            case 1    % tofile
                                if ischar(pval)
                                    tofile = true;
                                    filename = pval;
                                end
                            case 2    % mirrorsite
                                warning(message('bioinfo:getpdb:MirrorSiteNotSupported'))
                                
                            case 3  % sequenceonly
                                seqonly = bioinfoprivate.opttf(pval);
                                if isempty(seqonly)
                                    error(message('bioinfo:getpdb:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                                end
                        end
                    end
                end
            end
            
            
            % error if ID isn't a string
            if ~ischar(pdbID)
                error(message('bioinfo:getpdb:NotString'))
            end
            
            % get sequence from pdb.fasta if SEQUENCEONLY is true, otherwise full pdb
            if seqonly == true
                searchurl = [PDBsite '/downloadFile.do?fileFormat=FASTA&compression=NO&structureId=' pdbID];
                [~, pdb] = fastaread(searchurl);
            else
                searchurl = [PDBsite '/downloadFile.do?fileFormat=pdb&compression=NO&structureId=' pdbID];
                
                % get the html file that is returned as a string
                s=urlread(searchurl);
                
                % replace the html version of &
                s=strrep(s,'&amp;','&');
                
                % Find first line of the actual data
                start = strfind(s,'HEADER');
                
                if isempty(start)
                    % search for text indicating that there weren't any files found
                    notfound=regexp(s,'The file you requested does not exist.','once');
                    
                    % string was found, meaning no results were found
                    if ~isempty(notfound),
                        error(message('bioinfo:getpdb:PDBIDNotFound', pdbID)) ;
                    end
                    error(message('bioinfo:getpdb:PDBIDAccessProblem', pdbID));
                end
                
                [~, endOfFile] = regexp(s,'\nEND\s.*\n');
                
                % shorten string, to search for uid info
                s=s(start:endOfFile);
                
                %make each line a separate row in string array
                pdbdata = char(strread(s,'%s','delimiter','\n','whitespace',''));
                
                %pass to PDBREAD to create structure
                pdb=pdbread(pdbdata);
                
                
            end
            
            if nargout
                PDBstruct = pdb;
                if ~seqonly
                    % add URL
                    PDBstruct.SearchURL = searchurl;
                end
            else
                if seqonly || ~usejava('desktop')
                    disp(pdb);
                else
                    disp(pdb);
                    disp([char(9) 'SearchURL: <a href="' searchurl '"> ' pdbID ' </a>']);
                end
                
            end
            
            %  write out file
            if tofile == true
                fid = fopen(filename,'wt') ;
                if fid == (-1)
                    error(message('bioinfo:getpdb:CouldNotOpenFile', filename));
                else
                    rows = size(pdbdata,1);
                    
                    for rcount=1:rows-1,
                        fprintf(fid,'%s\n',pdbdata(rcount,:));
                    end
                    fprintf(fid,'%s',pdbdata(rows,:));
                    fclose(fid);
                end
            end
        end
        
        
        function PDBstruct = read(pdbfile,varargin)
            % read reads a Protein Data Bank (PDB) coordinate entry file.
            %
            %   PDBSTRUCT = PDBREAD(PDBFILE) reads in a PDB coordinate entry file
            %   PBDFILE and creates a structure PDBSTRUCT containing fields
            %   corresponding to the PDB records. The coordinate information is stored
            %   in the Model field. If the PDB file contains multiple models, then the
            %   coordinate information of each model will be stored in an array of
            %   structures in the Model field. The sequence information is stored in
            %   the Sequence field as structure where 'ResidueNames' contains the
            %   three-letter codes for the sequence residues and 'Sequence'
            %   contains the single-letter codes for the sequence.
            %
            %   PDBFILE can also be a URL, column vector of strings, or a MATLAB character
            %   array containing the text of a PDB file.
            %
            %   PDBREAD(PDBFILE,'MODELNUM',M) reads only the coordinate information of
            %   the M-th model and all other sections of the PDBFILE. If M does not
            %   correspond to an existing model serial number in the PDBFILE, the whole
            %   file is read.
            %
            %   Based on PDB format Version 2.0.
            %
            %   Examples:
            %
            %       % Get the Green Fluorescent Protein PDB data and save it to a file.
            %       gfl = getpdb('1GFL','TOFILE','1gfl.pdb');
            %
            %       % Then read it into MATLAB
            %       gfl = pdbread('1gfl.pdb')
            %
            %       % Read the Green Fluorescent Protein data from PDB
            %       gfl = pdbread('http://www.rcsb.org/pdb/files/1gfl.pdb')
            %
            %   See also GENPEPTREAD, GETPDB, MOLVIEWER, PDBDISTPLOT, PDBWRITE.
            
            %   Reference:
            %   The PDB format contents guide version 2.3.
            %   http://www.wwpdb.org/documentation/format23/v2.3.html
            
            %   Copyright 2002-2016 The MathWorks, Inc.
            
            % Validate input data
            if nargin > 0
                pdbfile = convertStringsToChars(pdbfile);
            end
            
            if nargin > 1
                [varargin{:}] = convertStringsToChars(varargin{:});
            end
            
            if nargin < 1
                error(message('bioinfo:pdbread:NotEnoughInputs', mfilename));
            end
            
            % Check pdbfile is a string.
            if ~ischar(pdbfile)
                error(message('bioinfo:pdbread:InvalidInput'));
            end
            
            % Get PDB filename
            useTempFile = true;
            if isvector(pdbfile) && exist(pdbfile,'file')
                filename = pdbfile;
                useTempFile = false;
            elseif ~isempty(pdbfile) && ~isempty(strfind(pdbfile(1,:),'HEADER')) && ~isempty(strfind(pdbfile(end,:),'END'))
                filename = savetempfile(pdbfile);
            elseif (strfind(pdbfile, '://'))
                if (~usejava('jvm'))
                    error(message('bioinfo:pdbread:NoJava'))
                end
                % Must be a URL
                pdbfile = urlread(pdbfile);
                % Clean up any &amp s
                pdbfile = strrep(pdbfile,'&amp;','&');
                filename = savetempfile(pdbfile);
            else
                error(message('bioinfo:pdbread:BadPDBFile'))
            end
            
            modelNum = 0;
            if  nargin > 1
                if rem(nargin,2) ~= 1
                    error(message('bioinfo:pdbread:IncorrectNumberOfArguments', mfilename));
                end
                okargs = {'modelnum'};
                for j=1:2:nargin-1
                    pname = varargin{j};
                    pval = varargin{j+1};
                    k = find(strncmpi(pname,okargs,numel(pname)));
                    if isempty(k)
                        error(message('bioinfo:pdbread:UnknownParameterName', pname));
                    elseif length(k)>1
                        error(message('bioinfo:pdbread:AmbiguousParameterName', pname));
                    else
                        switch(k)
                            case 1  % model number
                                modelNum = pval;
                                if ~isnumeric(modelNum) || numel(modelNum)> 1 || isempty(modelNum) || modelNum < 0
                                    error(message('bioinfo:pdbread:BadModelNumber'))
                                end
                        end
                    end
                end
            end
            
            % Open file
            fid = fopen(filename,'r');
            
            if fid == -1
                if ~useTempFile
                    error(message('bioinfo:pdbread:CouldNotOpenFile', filename));
                else
                    error(message('bioinfo:pdbread:CouldNotOpenTempFile', filename));
                end
            end
            
            theLines = textscan(fid,'%s','delimiter','\n');
            theLines = theLines{1};
            fclose(fid);
            
            if useTempFile
                delete(filename);
            end
            
            % Remove the lines corresponding to other models when requested a
            % particular model
            if modelNum ~= 0
                modelStarts = find(strncmp(theLines,'MODEL ',6));
                modelEnds = find(strncmp(theLines,'ENDMDL ',7));
                modelMatch = ['^MODEL \s*' num2str(modelNum) '\s'];
                h = find(~cellfun(@isempty,regexp(theLines(modelStarts),modelMatch,'once')));
                if isempty(h)
                    warning(message('bioinfo:pdbread:NoModelEntry', modelNum));
                else
                    q = [1:modelStarts(1)-1 modelStarts(h):modelEnds(h) modelEnds(end)+1:numel(theLines)];
                    theLines = theLines(q);
                end
            end
            
            %---  Initialize structure and variables
            
            PDBstruct = [];
            
            % Initializations necessary for various record types mentioned in the comments
            Journal_Entry = 1; % counter for the number of the JRNL records in the PDB file
            
            % Array of flags for the subrecords.
            % JRNLSubRecords(1) = AUTH
            % JRNLSubRecords(2) = TITL
            % JRNLSubRecords(3) = EDIT
            % JRNLSubRecords(4) = REF
            % JRNLSubRecords(5) = PUBL
            % JRNLSubRecords(6) = REFN
            % JournalSubRecords = zeros(1,6);
            % OldNumKeywds = 0; % KEYDS
            % OldNumEntries = 0; % OBSLTE
            % OldNumAuthors = 0; % AUTHOR
            
            NumOfRevisionDate = 0; % REVDAT
            NumOfDBReferences = 0; % DBREF
            NumOfSequenceConflicts = 0; % SEQADV
            NumOfModifiedResidues = 0; % MODRES
            NumOfHeterogen = 0; % HET
            NumOfHeterogenName = 0; % HETNAM
            NumOfHelix = 0; % HELIX
            NumOfSheet = 0; % SHEET
            NumOfTurn = 0; %TURN
            NumOfSSBond = 0; % SSBOND
            NumOfLink = 0; % LINK
            NumOfHydrogenBond = 0; % HYDBND
            NumOfSaltBridge =0; % SLTBRG
            NumOfCISPeptides = 0; % CISPEP
            NumOfTranslationVector = 0; % TVECT
            NumOfConnectivity = 0; % CONECT
            NumOfAtomSD = 0; % SIGATM
            NumOfAnisotropicTemp = 0; % ANISOU
            NumOfAnisotropicTempSD = 0; % SIGUIJ
            NumOfAtom = 0; % Atom
            NumOfHeterogenSynonym = 0; % HETSYN
            NumOfFormula = 0; % FORMUL
            NumOfHeterogenAtom = 0; % HETATM
            NumOfModel = 0; % MODEL
            NumOfTerminal = 0; % TER
            NumResChain = 0;
            JournalAuth = 0;
            ModelSerNum = 1;
            PrevRes = '';
            PrevSiteName = '';
            PrevHetIDHeterogenSynonym = '';
            PrevHetIDFormula = '';
            PrevHetIDHeterogenName = '';
            PrevRemark = 0;
            PrevFtnote = 0;
            
            % SITE record
            TmpStruct = struct('ResName',{},'ChainID',{''},'ResSeqNo',{},'InsCode',' ');
            TmpAtomStruct = allocateAtoms;
            TmpTerStruct = allocateTerminal;
            TmpHetStruct = allocateAtoms;
            TmpAtomSDStruct = allocateAtomSD;
            TmpAnisoTempStruct = allocateAnisoTemp;
            TmpAnisoTempSDStruct = allocateAnisoTempSD;
            
            % Flags used for every model
            modelWithAtom = false;
            modelWithHetAtom = false;
            modelWithTerminal = false;
            modelWithAtomSD = false;
            modelWithAnisoTemp = false;
            modelWithAnisoTempSD = false;
            
            numLines = numel(theLines);
            lineCounter = 1;
            blank = blanks(80);
            
            while lineCounter <= numLines
                tline = theLines{lineCounter};
                lineCounter = lineCounter+1;
                if ~ischar(tline)
                    break; % For end of file recognition
                end
                len = length(tline);
                if len == 0 % Omit the empty lines to avoid error of invalid matrix index.
                    continue
                end
                
                % RCSB web site format requires each line to have 80 characters. This
                % avoids exceeding the matrix dimension for lines with less than 80 characters.
                tline = [tline blank(len+1:80)];
                Record_name = upper(tline(1:6));
                
                % Remove trailing blanks, assuming that the record name will be left
                % aligned (as mentioned in the RCSB file format doc).
                Record_name = deblank(Record_name);
                
                % Take care of ORIGX1,ORIGX2,ORIGX3 as well as SCALE and MTRIX
                if strncmp(Record_name,'ORIGX',5) || strncmp(Record_name,'SCALE',5) || strncmp(Record_name,'MTRIX',5)
                    Record_name = Record_name(1:5);
                end
                
                switch Record_name
                    case 'HEADER' %Single/Mandatory
                        PDBstruct.Header = struct('classification',{deblank(tline(11:50))},...
                            'depDate',{tline(51:59)},...
                            'idCode',{tline(63:66)});
                        
                    case 'OBSLTE' % Single Continued/Optional : mandatory in withdrawn entries
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Obsolete = struct('repDate',{tline(12:20)},...
                                'idCode',{tline(22:25)},...
                                'rIdCode',{strtrim(tline(32:70))});
                        else
                            PDBstruct.Obsolete.rIdCode = strvcat(PDBstruct.Obsolete.rIdCode,...
                                strtrim(tline(32:70))); %#ok
                        end
                        
                    case 'TITLE'%Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Title = strtrim(tline(11:70));
                        else
                            PDBstruct.Title = strvcat(PDBstruct.Title,strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'CAVEAT' %Single Continued/Optional
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Caveat =  struct('idCode',{tline(12:15)},'comment',{strtrim(tline(20:70))});
                        else
                            PDBstruct.Caveat.comment = strvcat(PDBstruct.Caveat.comment,strtrim(tline(20:70))); %#ok
                        end
                        
                    case 'COMPND' %Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Compound = strtrim(tline(11:70));
                        else
                            PDBstruct.Compound = strvcat(PDBstruct.Compound,strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'SOURCE' %Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Source = strtrim(tline(11:70));
                        else
                            PDBstruct.Source = strvcat(PDBstruct.Source,strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'KEYWDS' %Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Keywords = strtrim(tline(11:70));
                        else
                            PDBstruct.Keywords = strvcat(PDBstruct.Keywords,strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'EXPDTA' %Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.ExperimentData = strtrim(tline(11:70));
                        else
                            PDBstruct.ExperimentData = strvcat(PDBstruct.ExperimentData,...
                                strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'AUTHOR' %Single Continued/Mandatory
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Authors = strtrim(tline(11:70));
                        else
                            PDBstruct.Authors = strvcat(PDBstruct.Authors,strtrim(tline(11:70))); %#ok
                        end
                        
                    case 'REVDAT' %Multiple Continued/Mandatory
                        % %             CurRevDate = strtrim(tline(14:22));
                        % %             if(isempty(CurRevDate))
                        % %                 PDB_struct.RevisionDate(NumOfRevisionDate).modType =...
                        % %                     [PDB_struct.RevisionDate(NumOfRevisionDate).modType;...
                        % %                     str2double(tline(32))];
                        % %                 PDB_struct.RevisionDate(NumOfRevisionDate).record =...
                        % %                     strvcat(PDB_struct.RevisionDate(NumOfRevisionDate).record,...
                        % %                     strtrim(tline(40:66)));%#ok
                        % %             else
                        NumOfRevisionDate = NumOfRevisionDate+1;
                        PDBstruct.RevisionDate(NumOfRevisionDate) =...
                            struct('modNum',{str2double(tline(8:10))},...
                            'modDate',{tline(14:22)},...
                            'modId',{tline(24:28)},...
                            'modType',{str2double(tline(32))},...
                            'record',{strtrim(tline(40:66))});
                        % %             end
                        
                    case 'SPRSDE' %Single Continued/Optional
                        if strcmp(tline(9:10),'  ') %check for continuation
                            PDBstruct.Superseded = struct('Supersededdate',{tline(12:20)},...
                                'idCode',{tline(22:25)},...
                                'sIdCode',{strtrim(tline(32:70))});
                        else
                            PDBstruct.Superseded.sIdCode =...
                                strvcat(PDBstruct.Superseded.sIdCode,strtrim(tline(32:70))); %#ok
                        end
                    case 'JRNL' %Other/Optional : This record has following sub-records: AUTH,TITL,EDIT,REF,PUBL,REFN
                        SubRecord = tline(13:16);
                        SubRecord = deblank(SubRecord); % Remove the trailing blanks. Needed for REF
                        
                        switch SubRecord
                            case 'AUTH'
                                if ~JournalAuth
                                    JournalAuth = 1;
                                    PDBstruct.Journal(Journal_Entry) = struct('Author',{''},...
                                        'Title',{''},...
                                        'Editor',{''},...
                                        'Reference',{''},...
                                        'Publisher',{''},...
                                        'CitationReference',{''});
                                    PDBstruct.Journal(Journal_Entry).Author =...
                                        strvcat(PDBstruct.Journal(Journal_Entry).Author,strtrim(tline(20:70))); %#ok
                                else
                                    PDBstruct.Journal(Journal_Entry).Author =...
                                        strvcat(PDBstruct.Journal(Journal_Entry).Author,strtrim(tline(20:70))); %#ok
                                end
                                
                            case 'TITL'
                                PDBstruct.Journal(Journal_Entry).Title =...
                                    strvcat(PDBstruct.Journal(Journal_Entry).Title,strtrim(tline(20:70))); %#ok
                                
                            case 'EDIT'
                                PDBstruct.Journal(Journal_Entry).Editor = ...
                                    strvcat(PDBstruct.Journal(Journal_Entry).Editor,strtrim(tline(20:70))); %#ok
                                
                            case 'REF'
                                PDBstruct.Journal(Journal_Entry).Reference = ...
                                    strvcat(PDBstruct.Journal(Journal_Entry).Reference,strtrim(tline(20:70))); %#ok
                                
                            case 'PUBL'
                                PDBstruct.Journal(Journal_Entry).Publisher = ...
                                    strvcat(PDBstruct.Journal(Journal_Entry).Publisher,strtrim(tline(20:70))); %#ok
                                
                            case 'REFN'
                                PDBstruct.Journal(Journal_Entry).CitationReference = ...
                                    strvcat(PDBstruct.Journal(Journal_Entry).CitationReference,tline(20:70)); %#ok
                                % REFN is the last subrecord and it is a single line record
                                Journal_Entry = Journal_Entry+1;
                                JournalAuth = 0;
                                
                            otherwise
                                %disp('Invalid subrecord type');
                        end
                        
                        % Some of the Remark records are mandatory and some are optional
                    case 'REMARK'
                        RemarkNo = str2double(tline(7:10));
                        
                        %Other/Optional
                        if RemarkNo == 1
                            CurRemark = RemarkNo;
                            if CurRemark ~= PrevRemark
                                PDBstruct.Remark1 = '';
                                PrevRemark = CurRemark;
                            end
                            if strcmp(tline(12:20),'REFERENCE')
                                PDBstruct.Remark1.NumJournals = str2double(tline(22:70));
                                PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals) = ...
                                    struct('Author',{''},...
                                    'Title',{''},...
                                    'Editor',{''},...
                                    'Reference',{''},...
                                    'Publisher',{''},...
                                    'CitationReference',{''});
                            else
                                SubRecord = tline(13:16);
                                SubRecord = deblank(SubRecord); % Remove the trailing blanks. Needed for REF
                                
                                switch SubRecord
                                    
                                    case 'AUTH'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Author = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Author,...
                                            strtrim(tline(20:70))); %#ok
                                        
                                    case 'TITL'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Title = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Title,...
                                            strtrim(tline(20:70)));  %#ok
                                        
                                    case 'EDIT'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Editor = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Editor,...
                                            strtrim(tline(20:70)));  %#ok
                                        
                                    case 'REF'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Reference = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Reference,...
                                            strtrim(tline(20:70)));  %#ok
                                        
                                    case 'PUBL'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Publisher = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).Publisher,...
                                            strtrim(tline(20:70)));  %#ok
                                        
                                    case 'REFN'
                                        PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).CitationReference = ...
                                            strvcat(PDBstruct.Remark1.JournalEntry(PDBstruct.Remark1.NumJournals).CitationReference,...
                                            strtrim(tline(20:70)));  %#ok
                                        
                                    otherwise
                                        %disp('Invalid subrecord type');
                                end
                            end
                        elseif RemarkNo == 2 % Other/Mandatory
                            CurRemark = RemarkNo;
                            if CurRemark ~= PrevRemark
                                PDBstruct.Remark2 = '';
                                PrevRemark = CurRemark;
                            end
                            if strcmp(tline(12:22),'RESOLUTION.')
                                
                                if strcmp(tline(29:38),'ANGSTROMS.')
                                    PDBstruct.Remark2.Resolution = str2double(tline(23:27));
                                else
                                    PDBstruct.Remark2.Detail = strtrim(tline(12:70));
                                end
                                
                            else
                                if isequal(tline(12:22),blanks(11))
                                    PDBstruct.Remark2.Detail = '';
                                else
                                    PDBstruct.Remark2.Detail = ...
                                        strvcat(PDBstruct.Remark2.Detail, strtrim(tline(12:70)));  %#ok
                                end
                            end
                            PDBstruct.Remark3.Refinement = '';
                            
                            % Other/Mandatory:  Following code assumes a single
                            % occurrence of the Remark3 record type
                        elseif RemarkNo == 3
                            CurRemark = RemarkNo;
                            if CurRemark ~= PrevRemark
                                PDBstruct.Remark3.Refinement = '';
                                PrevRemark = CurRemark;
                            else
                                PDBstruct.Remark3.Refinement = ...
                                    strvcat(PDBstruct.Remark3.Refinement,tline(12:70)); %#ok
                            end
                            
                            % Other/Optional
                        else
                            CurRemark = RemarkNo;
                            tmpRemark = sprintf('%d',CurRemark);
                            
                            if CurRemark ~= PrevRemark
                                PDBstruct.(['Remark' tmpRemark]) = tline(12:70);
                                PrevRemark = CurRemark;
                            else
                                PDBstruct.(['Remark' tmpRemark])  = ...
                                    strvcat(PDBstruct.(['Remark' tmpRemark]) ,tline(12:70)); %#ok
                            end
                            
                            
                        end
                        
                    case 'DBREF'%Multiple/Optional
                        NumOfDBReferences = NumOfDBReferences+1;
                        PDBstruct.DBReferences(NumOfDBReferences) = ...
                            struct('idCode',{tline(8:11)},...
                            'chainID',{tline(13)},...
                            'seqBegin',{str2double(tline(15:18))},...
                            'insertBegin',{tline(19)},...
                            'seqEnd',{str2double(tline(21:24))},...
                            'insertEnd',{tline(25)},...
                            'database',{tline(27:32)},...
                            'dbAccession',{tline(34:41)},...
                            'dbIdCode',{tline(43:54)},...
                            'dbseqBegin',{str2double(tline(56:60))},...
                            'idbnsBeg',{tline(61)},...
                            'dbseqEnd',{str2double(tline(63:67))},...
                            'dbinsEnd',{tline(68)});
                        
                    case 'SEQADV'%Multiple/Optional
                        NumOfSequenceConflicts = NumOfSequenceConflicts+1;
                        PDBstruct.SequenceConflicts(NumOfSequenceConflicts) = ...
                            struct('idCode',{tline(8:11)},...
                            'resName',{tline(13:15)},...
                            'chainID',{tline(17)},...
                            'seqNum',{str2double(tline(19:22))},...
                            'iCode',{tline(23)},...
                            'database',{tline(25:28)},...
                            'dbIdCode',{tline(30:38)},...
                            'dbRes',{tline(40:42)},...
                            'dbSeq',{str2double(tline(44:48))},...
                            'conflict',{strtrim(tline(50:70))});
                        
                    case 'SEQRES' %Multiple/Optional
                        if isspace(tline(12))
                            CurRes = sprintf('%d',str2double(tline(14:17)));
                        else
                            CurRes = tline(12);
                        end
                        
                        if ~isequal(CurRes,PrevRes)
                            NumResChain = NumResChain + 1;
                            PDBstruct.Sequence(NumResChain).NumOfResidues = str2double(tline(14:17));
                            PDBstruct.Sequence(NumResChain).ChainID = tline(12);
                            PDBstruct.Sequence(NumResChain).ResidueNames = tline(20:70);
                            PDBstruct.Sequence(NumResChain).Sequence = strtrim(GetSequence(tline(20:70)));
                            PrevRes = CurRes;
                        else
                            PDBstruct.Sequence(NumResChain).ResidueNames = ...
                                strcat(PDBstruct.Sequence(NumResChain).ResidueNames,[' ',tline(20:70)]);
                            PDBstruct.Sequence(NumResChain).Sequence = ...
                                strcat(PDBstruct.Sequence(NumResChain).Sequence,GetSequence(tline(20:70)));
                        end
                        
                    case 'MODRES'%Multiple/Optional
                        NumOfModifiedResidues = NumOfModifiedResidues+1;
                        PDBstruct.ModifiedResidues(NumOfModifiedResidues) = ...
                            struct('idCode',{tline(8:11)},...
                            'resName',{tline(13:15)},...
                            'chainID',{tline(17)},...
                            'seqNum',{str2double(tline(19:22))},...
                            'iCode',{tline(23)},...
                            'stdRes',{tline(25:27)},...
                            'comment',{strtrim(tline(30:70))});
                        
                    case 'HET'%Multiple/Optional
                        NumOfHeterogen = NumOfHeterogen +1;
                        PDBstruct.Heterogen(NumOfHeterogen) = ...
                            struct('hetID',{tline(8:10)},...
                            'ChainID',{tline(13)},...
                            'seqNum',{str2double(tline(14:17))},...
                            'iCode',{tline(18)},...
                            'numHetAtoms',{str2double(tline(21:25))},...
                            'text',{strtrim(tline(31:70))});
                        
                    case 'HETNAM'%Multiple Continued/Optional
                        CurHetIDHeterogenName = tline(12:14);
                        
                        if ~strcmp(CurHetIDHeterogenName,PrevHetIDHeterogenName)
                            NumOfHeterogenName = NumOfHeterogenName + 1;
                            PDBstruct.HeterogenName(NumOfHeterogenName).hetID = CurHetIDHeterogenName;
                            PDBstruct.HeterogenName(NumOfHeterogenName).ChemName = strtrim(tline(16:70));
                            PrevHetIDHeterogenName = CurHetIDHeterogenName;
                        else
                            PDBstruct.HeterogenName(NumOfHeterogenName).ChemName = ...
                                strvcat(PDBstruct.HeterogenName(NumOfHeterogenName).ChemName,...
                                strtrim(tline(16:70))); %#ok
                        end
                        
                        %Multiple/Optional
                    case 'HETSYN'
                        CurHetIDHeterogenSynonym = tline(12:14);
                        
                        if ~strcmp(CurHetIDHeterogenSynonym,PrevHetIDHeterogenSynonym)
                            NumOfHeterogenSynonym = NumOfHeterogenSynonym+1;
                            PDBstruct.HeterogenSynonym(NumOfHeterogenSynonym).hetID = CurHetIDHeterogenSynonym;
                            PDBstruct.HeterogenSynonym(NumOfHeterogenSynonym).hetSynonyms = strtrim(tline(16:70));
                            PrevHetIDHeterogenSynonym = CurHetIDHeterogenSynonym;
                        else
                            PDBstruct.HeterogenSynonym(NumOfHeterogenSynonym).hetSynonyms = ...
                                strvcat(PDBstruct.HeterogenSynonym(NumOfHeterogenSynonym).hetSynonyms,...
                                strtrim(tline(16:70))); %#ok
                        end
                        
                        %Multiple Continued/Optional
                    case 'FORMUL'
                        CurHetIDFormula = tline(13:15);
                        
                        if ~strcmp(CurHetIDFormula,PrevHetIDFormula)
                            NumOfFormula = NumOfFormula+1;
                            PDBstruct.Formula(NumOfFormula).CompNo = str2double(tline(9:10));
                            PDBstruct.Formula(NumOfFormula).hetID = tline(13:15);
                            PDBstruct.Formula(NumOfFormula).ChemForm = strtrim(tline(19:70));
                            PrevHetIDFormula = CurHetIDFormula;
                        else
                            PDBstruct.Formula(NumOfFormula).ChemForm = ...
                                strvcat(PDBstruct.Formula(NumOfFormula).ChemForm,strtrim(tline(19:70))); %#ok
                        end
                        
                        %Multiple/Optional
                    case 'HELIX'
                        NumOfHelix = NumOfHelix+1;
                        PDBstruct.Helix(NumOfHelix) = ...
                            struct('serNum',{str2double(tline(8:10))},...
                            'helixID',{tline(12:14)},...
                            'initResName',{tline(16:18)},...
                            'initChainID',{tline(20)},...
                            'initSeqNum',{str2double(tline(22:25))},...
                            'initICode',{tline(26)},...
                            'endResName',{tline(28:30)},...
                            'endChainID',{tline(32)},...
                            'endSeqNum',{str2double(tline(34:37))},...
                            'endICode',{tline(38)},...
                            'helixClass',{str2double(tline(39:40))},...
                            'comment',{tline(41:70)},...
                            'length',{str2double(tline(72:76))});
                        
                    case 'SHEET'%Multiple/Optional
                        NumOfSheet = NumOfSheet+1;
                        PDBstruct.Sheet(NumOfSheet) = ...
                            struct('strand',{str2double(tline(8:10))},...
                            'sheetID',{tline(12:14)},...
                            'numStrands',{str2double(tline(15:16))},...
                            'initResName',{tline(18:20)},...
                            'initChainID',{tline(22)},...
                            'initSeqNum',{str2double(tline(23:26))},...
                            'initICode',{tline(27)},...
                            'endResName',{tline(29:31)},...
                            'endChainID',{tline(33)},...
                            'endSeqNum',{str2double(tline(34:37))},...
                            'endICode',{tline(38)},...
                            'sense',{str2double(tline(39:40))},...
                            'curAtom',{tline(42:45)},...
                            'curResName',{tline(46:48)},...
                            'curChainId',{tline(50)},...
                            'curResSeq',{str2double(tline(51:54))},...
                            'curICode',{tline(55)},...
                            'prevAtom',{tline(57:60)},...
                            'prevResName',{tline(61:63)},...
                            'prevChainId',{tline(65)},...
                            'prevResSeq',{str2double(tline(66:69))},...
                            'prevICode',{tline(70)});
                        
                    case 'TURN' %Multiple/Optional
                        NumOfTurn = NumOfTurn+1;
                        PDBstruct.Turn(NumOfTurn) = ...
                            struct('seq',{str2double(tline(8:10))},...
                            'turnId',{tline(12:14)},...
                            'initResName',{tline(16:18)},...
                            'initChainId',{tline(20)},...
                            'initSeqNum',{str2double(tline(21:24))},...
                            'initICode',{tline(25)},...
                            'endResName',{tline(27:29)},...
                            'endChainId',{tline(31)},...
                            'endSeqNum',{str2double(tline(32:35))},...
                            'endICode',{tline(36)},...
                            'comment',{tline(41:70)});
                        
                        %Multiple/Optional
                    case 'SSBOND'
                        NumOfSSBond = NumOfSSBond+1;
                        PDBstruct.SSBond(NumOfSSBond) = ...
                            struct('serNum',{str2double(tline(8:10))},...
                            'resName1',{tline(12:14)},...
                            'chainID1',{tline(16)},...
                            'seqNum1',{str2double(tline(18:21))},...
                            'icode1',{tline(22)},...
                            'resName2',{tline(26:28)},...
                            'chainID2',{tline(30)},...
                            'seqNum2',{str2double(tline(32:35))},...
                            'icode2',{tline(36)},...
                            'sym1',{tline(60:65)},...
                            'sym2',{tline(67:72)});
                        
                        %Multiple/Optional
                    case 'LINK'
                        NumOfLink = NumOfLink+1;
                        PDBstruct.Link(NumOfLink) = ...
                            struct('remove1',{tline(13:16)},...
                            'altLoc1',{tline(17)},...
                            'resName1',{tline(18:20)},...
                            'chainID1',{tline(22)},...
                            'resSeq1',{str2double(tline(23:26))},...
                            'iCode1',{tline(27)},...
                            'AtomName2',{tline(43:46)},...
                            'altLoc2',{tline(47)},...
                            'resName2',{tline(48:50)},...
                            'chainID2',{tline(52)},...
                            'resSeq2',{str2double(tline(53:56))},...
                            'iCode2',{tline(57)},...
                            'sym1',{tline(60:65)},...
                            'sym2',{tline(67:72)});
                        
                        %Multiple/Optional
                    case 'HYDBND'
                        NumOfHydrogenBond = NumOfHydrogenBond+1;
                        PDBstruct.HydrogenBond(NumOfHydrogenBond) = ...
                            struct('AtomName1',{tline(13:16)},...
                            'altLoc1',{tline(17)},...
                            'resName1',{tline(18:20)},...
                            'Chain1',{tline(22)},...
                            'resSeq1',{str2double(tline(23:27))},...
                            'ICode1',{tline(28)},...
                            'nameH',{tline(30:33)},...
                            'altLocH',{tline(34)},...
                            'ChainH',{tline(36)},...
                            'resSeqH',{str2double(tline(37:41))},...
                            'iCodeH',{tline(42)},...
                            'name2',{tline(44:47)},...
                            'altLoc2',{tline(48)},...
                            'resName2',{tline(49:51)},...
                            'chainID2',{tline(53)},...
                            'resSeq2',{str2double(tline(54:58))},...
                            'iCode2',{tline(59)},...
                            'sym1',{tline(60:65)},...
                            'sym2',{tline(67:72)});
                        
                        %Multiple/Optional
                    case 'SLTBRG'
                        NumOfSaltBridge = NumOfSaltBridge+1;
                        PDBstruct.SaltBridge(NumOfSaltBridge) = ...
                            struct('AtomName1',{tline(13:16)},...
                            'altLoc1',{tline(17)},...
                            'resName1',{tline(18:20)},...
                            'chainID1',{tline(22)},...
                            'resSeq1',{str2double(tline(23:26))},...
                            'iCode1',{tline(27)},...
                            'AtomName2',{tline(43:46)},...
                            'altLoc2',{tline(47)},...
                            'resName2',{tline(48:50)},...
                            'chainID2',{tline(52)},...
                            'resSeq2',{str2double(tline(53:56))},...
                            'iCode2',{tline(57)},...
                            'sym1',{tline(60:65)},...
                            'sym2',{tline(67:72)});
                        
                        %Multiple/Optional
                    case 'CISPEP'
                        NumOfCISPeptides = NumOfCISPeptides+1;
                        PDBstruct.CISPeptides(NumOfCISPeptides) = ...
                            struct('serNum',{str2double(tline(8:10))},...
                            'ResName1',{tline(12:14)},...
                            'chainID1',{tline(16)},...
                            'seqNum1',{str2double(tline(18:21))},...
                            'icode1',{tline(22)},...
                            'ResName2',{tline(26:28)},...
                            'chainID2',{tline(30)},...
                            'seqNum2',{str2double(tline(32:35))},...
                            'icode2',{tline(36)},...
                            'modNum',{str2double(tline(44:46))},...
                            'measure',{str2double(tline(54:59))});
                        
                        %Multiple/Optional
                    case 'SITE'
                        CurSiteName = tline(12:14);
                        if ~isfield(PDBstruct,'Site') || isempty(PDBstruct.Site)
                            ResDetail = struct('ResName',{},...
                                'ChainID',{''},...
                                'ResSeqNo',{},...
                                'InsCode',{''});
                            PDBstruct.Site.SiteDetail = struct('SeqNo',{0},...
                                'SiteName',{''},...
                                'NoOfRes',{0},...
                                'ResDet',{ResDetail});
                            PDBstruct.Site.SiteNumber = 0;
                        end
                        if ~strcmp(CurSiteName,PrevSiteName)
                            ResNos = 0;
                            PDBstruct.Site.SiteNumber = PDBstruct.Site.SiteNumber+1;
                            PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).SeqNo = ...
                                str2double(tline(8:10));
                            PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).SiteName = ...
                                strtrim(tline(12:14));
                            PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).NumRes = ...
                                str2double(tline(16:17));
                            [PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).ResDet, ResNos] = ...
                                GetResidueStruct(TmpStruct,tline(19:61),ResNos);
                            PrevSiteName = CurSiteName;
                        else
                            PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).SeqNo = ...
                                [PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).SeqNo;str2double(tline(8:10))];
                            [PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).ResDet, ResNos] = ...
                                GetResidueStruct(PDBstruct.Site.SiteDetail(PDBstruct.Site.SiteNumber).ResDet,...
                                tline(19:61),ResNos);
                        end
                        
                        %Single/Mandatory
                    case 'CRYST1'
                        % Fields in this record: Header(record name), a,b,c (all 3 in
                        % Angstrom),alpha,beta,gamma(all 3 in degrees),sGroup,z
                        PDBstruct.Cryst1=struct('a',{str2double(tline(7:15))},...
                            'b',{str2double(tline(16:24))},...
                            'c',{str2double(tline(25:33))},...
                            'alpha',{str2double(tline(34:40))},...
                            'beta',{str2double(tline(41:47))},...
                            'gamma',{str2double(tline(48:54))},...
                            'sGroup',{tline(56:66)},...
                            'z',{str2double(tline(67:70))});
                        
                        %Single/Mandatory
                    case 'ORIGX'
                        %Fields in this record: Header(record name), O[n][1](O11),
                        %O[n][2](O12), O[n][3](O13), T[n](T1)
                        OrigNum = str2double(tline(6));
                        PDBstruct.OriginX(OrigNum).On1 = str2double(tline(11:20));
                        PDBstruct.OriginX(OrigNum).On2 = str2double(tline(21:30));
                        PDBstruct.OriginX(OrigNum).On3 = str2double(tline(31:40));
                        PDBstruct.OriginX(OrigNum).Tn = str2double(tline(46:55));
                        
                        %Single/Mandatory
                    case 'SCALE'
                        ScaleNum = str2double(tline(6));
                        PDBstruct.Scale(ScaleNum).Sn1 = str2double(tline(11:20));
                        PDBstruct.Scale(ScaleNum).Sn2 = str2double(tline(21:30));
                        PDBstruct.Scale(ScaleNum).Sn3 = str2double(tline(31:40));
                        PDBstruct.Scale(ScaleNum).Un = str2double(tline(46:55));
                        
                        %Single/Optional: Mandatory if the complete unit must be
                        %generated from the given coordinates using
                        %non-crystallographic symmetry
                    case 'MTRIX'
                        Matrix_num = str2double(tline(6));
                        PDBstruct.Matrix(Matrix_num).SerNo = str2double(tline(8:10));
                        PDBstruct.Matrix(Matrix_num).Mn1 = str2double(tline(11:20));
                        PDBstruct.Matrix(Matrix_num).Mn2 = str2double(tline(21:30));
                        PDBstruct.Matrix(Matrix_num).Mn3 = str2double(tline(31:40));
                        PDBstruct.Matrix(Matrix_num).Vn = str2double(tline(46:55));
                        PDBstruct.Matrix(Matrix_num).iGiven = str2double(tline(60));
                        
                        %Multiple/Optional
                    case 'TVECT'
                        NumOfTranslationVector = NumOfTranslationVector+1;
                        PDBstruct.TranslationVector(NumOfTranslationVector) = ...
                            struct('SerNo',{str2double(tline(8:10))},...
                            't1',{str2double(tline(11:20))},...
                            't2',{str2double(tline(21:30))},...
                            't3',{str2double(tline(31:40))},...
                            'text',{tline(41:70)});
                        
                        % Group/Optional
                    case 'MODEL'
                        NumOfModel = NumOfModel+1;
                        NumOfAtom = 0;
                        NumOfTerminal = 0;
                        NumOfHeterogenAtom = 0;
                        NumOfAtomSD = 0;
                        NumOfAnisotropicTemp = 0;
                        NumOfAnisotropicTempSD = 0;
                        ModelSerNum = str2double(tline(11:14));
                        
                        %Multiple/Optional
                    case 'ATOM'
                        % If there only one model
                        if NumOfModel == 0
                            NumOfModel = ModelSerNum;
                            PDBstruct.Model = [];
                        end
                        modelWithAtom = true;
                        NumOfAtom = NumOfAtom+1;
                        TmpAtomStruct(NumOfAtom) = ...
                            struct('AtomSerNo',{str2int(tline(7:11))},...
                            'AtomName',{strtrim(tline(13:16))},...
                            'altLoc',{strtrim(tline(17))},...
                            'resName',{strtrim(tline(18:20))},...
                            'chainID',{tline(22)},...
                            'resSeq',{str2int(tline(23:26))},...
                            'iCode',{strtrim(tline(27))},...
                            'X',{str2float(tline(31:38))},...
                            'Y',{str2float(tline(39:46))},...
                            'Z',{str2float(tline(47:54))},...
                            'occupancy',{str2int(tline(55:60))},...
                            'tempFactor',{str2float(tline(61:66))},...
                            'segID',{tline(73:76)},...
                            'element',{strtrim(tline(77:78))},...
                            'charge',{tline(79:80)},...
                            'AtomNameStruct',struct('chemSymbol',{strtrim(tline(13:14))},...
                            'remoteInd',{strtrim(tline(15))},...
                            'branch',{strtrim(tline(16))})); %#ok<*AGROW>
                        
                        %Multiple/Optional
                    case 'SIGATM'
                        modelWithAtomSD = true;
                        NumOfAtomSD = NumOfAtomSD+1;
                        TmpAtomSDStruct(NumOfAtomSD) = ...
                            struct('AtomSerNo',{str2double(tline(7:11))},...
                            'AtomName',{tline(13:16)},...
                            'altLoc',{tline(17)},...
                            'resName',{tline(18:20)},...
                            'chainID',{tline(22)},...
                            'resSeq',{str2double(tline(23:26))},...
                            'iCode',{tline(27)},...
                            'sigX',{str2double(tline(31:38))},...
                            'sigY',{str2double(tline(39:46))},...
                            'sigZ',{str2double(tline(47:54))},...
                            'sigOcc',{str2double(tline(55:60))},...
                            'sigTemp',{str2double(tline(61:66))},...
                            'segID',{tline(73:76)},...
                            'element',{tline(77:78)},...
                            'charge',{tline(79:80)},...
                            'AtomNameStruct',struct('chemSymbol',{strtrim(tline(13:14))},...
                            'remoteInd',{strtrim(tline(15))},...
                            'branch',{strtrim(tline(16))}));
                        
                        %Multiple/Optional
                    case 'ANISOU'
                        modelWithAnisoTemp = true;
                        NumOfAnisotropicTemp = NumOfAnisotropicTemp+1;
                        TmpAnisoTempStruct(NumOfAnisotropicTemp) = ...
                            struct('AtomSerNo',{str2double(tline(7:11))},...
                            'AtomName',{tline(13:16)},...
                            'altLoc',{tline(17)},...
                            'resName',{tline(18:20)},...
                            'chainID',{tline(22)},...
                            'resSeq',{str2double(tline(23:26))},...
                            'iCode',{tline(27)},...
                            'U00',{str2double(tline(29:35))},...
                            'U11',{str2double(tline(36:42))},...
                            'U22',{str2double(tline(43:49))},...
                            'U01',{str2double(tline(50:56))},...
                            'U02',{str2double(tline(57:63))},...
                            'U12',{str2double(tline(64:70))},...
                            'segID',{tline(73:76)},...
                            'element',{tline(77:78)},...
                            'charge',{tline(79:80)},...
                            'AtomNameStruct',struct('chemSymbol',{strtrim(tline(13:14))},...
                            'remoteInd',{strtrim(tline(15))},...
                            'branch',{strtrim(tline(16))}));
                        
                        %Multiple/Optional
                    case 'SIGUIJ'
                        modelWithAnisoTempSD = true;
                        NumOfAnisotropicTempSD = NumOfAnisotropicTempSD+1;
                        TmpAnisoTempSDStruct(NumOfAnisotropicTempSD) = ...
                            struct('AtomSerNo',{str2double(tline(7:11))},...
                            'AtomName',{tline(13:16)},...
                            'altLoc',{tline(17)},...
                            'resName',{tline(18:20)},...
                            'chainID',{tline(22)},...
                            'resSeq',{str2double(tline(23:26))},...
                            'iCode',{tline(27)},...
                            'SIG11',{str2double(tline(29:35))},...
                            'SIG22',{str2double(tline(36:42))},...
                            'SIG33',{str2double(tline(43:49))},...
                            'SIG12',{str2double(tline(50:56))},...
                            'SIG13',{str2double(tline(57:63))},...
                            'SIG23',{str2double(tline(64:70))},...
                            'segID',{tline(73:76)},...
                            'element',{tline(77:78)},...
                            'charge',{tline(79:80)},...
                            'AtomNameStruct',struct('chemSymbol',{strtrim(tline(13:14))},...
                            'remoteInd',{strtrim(tline(15))},...
                            'branch',{strtrim(tline(16))}));
                        
                        % Group/Optional
                    case 'TER'
                        modelWithTerminal = true;
                        NumOfTerminal = NumOfTerminal + 1;
                        TmpTerStruct(NumOfTerminal) = ...
                            struct('SerialNo',{str2double(tline(7:11))},...
                            'resName',{strtrim(tline(18:20))},...
                            'chainID',{strtrim(tline(22))},...
                            'resSeq',{str2double(tline(23:26))},...
                            'iCode',{strtrim(tline(27))});
                        
                        %Multiple Continued/Optional
                    case 'HETATM'
                        modelWithHetAtom = true;
                        NumOfHeterogenAtom = NumOfHeterogenAtom+1;
                        TmpHetStruct(NumOfHeterogenAtom) = ...
                            struct('AtomSerNo',{str2int(tline(7:11))},...
                            'AtomName',{strtrim(tline(13:16))},...
                            'altLoc',{strtrim(tline(17))},...
                            'resName',{strtrim(tline(18:20))},...
                            'chainID',{tline(22)},...
                            'resSeq',{str2double(tline(23:26))},...
                            'iCode',{strtrim(tline(27))},...
                            'X',{str2double(tline(31:38))},...
                            'Y',{str2double(tline(39:46))},...
                            'Z',{str2double(tline(47:54))},...
                            'occupancy',{str2double(tline(55:60))},...
                            'tempFactor',{str2double(tline(61:66))},...
                            'segID',{tline(73:76)},...
                            'element',{strtrim(tline(77:78))},...
                            'charge',{tline(79:80)},...
                            'AtomNameStruct',struct('chemSymbol',{strtrim(tline(13:14))},...
                            'remoteInd',{strtrim(tline(15))},...
                            'branch',{strtrim(tline(16))}));
                        % Group/Optional
                    case 'ENDMDL'
                        % End of model
                        PDBstruct.Model(NumOfModel).MDLSerNo = ModelSerNum;
                        if modelWithAtom
                            PDBstruct.Model(NumOfModel).Atom = TmpAtomStruct;
                            modelWithAtom = false;
                            TmpAtomStruct = allocateAtoms;
                        end
                        if modelWithAtomSD
                            PDBstruct.Model(NumOfModel).AtomSD = TmpAtomSDStruct;
                            modelWithAtomSD = false;
                            TmpAtomSDStruct = allocateAtomSD;
                        end
                        if modelWithAnisoTemp
                            PDBstruct.Model(NumOfModel).AnisotropicTemp = TmpAnisoTempStruct;
                            modelWithAnisoTemp = false;
                            TmpAnisoTempStruct = allocateAnisoTemp;
                        end
                        if modelWithAnisoTempSD
                            PDBstruct.Model(NumOfModel).AnisotropicTempSD = TmpAnisoTempSDStruct;
                            modelWithAnisoTempSD = false;
                            TmpAnisoTempSDStruct = allocateAnisoTempSD;
                        end
                        if modelWithTerminal
                            PDBstruct.Model(NumOfModel).Terminal = TmpTerStruct;
                            modelWithTerminal = false;
                            TmpTerStruct = allocateTerminal;
                        end
                        if modelWithHetAtom
                            PDBstruct.Model(NumOfModel).HeterogenAtom = TmpHetStruct;
                            modelWithHetAtom = false;
                            TmpHetStruct = allocateAtoms;
                        end
                        %Multiple/Optional
                    case 'CONECT'
                        NumOfConnectivity = NumOfConnectivity+1;
                        temp_a = str2double(tline(7:11));
                        temp_b = GetAtomList(tline(12:31));
                        
                        PDBstruct.Connectivity (NumOfConnectivity) = ...
                            struct('AtomSerNo',{temp_a},...
                            'BondAtomList',{temp_b},...
                            'HydrogenAtomList',[],...
                            'SaltBridgeAtom',[]);
                        
                        %Single/Mandatory
                    case 'MASTER'
                        PDBstruct.Master = struct('numREMARK',{str2double(tline(11:15))},...
                            'numHET',{str2double(tline(21:25))},...
                            'numHelix',{str2double(tline(26:30))},...
                            'numSheet',{str2double(tline(31:35))},...
                            'numTurn',{str2double(tline(36:40))},...
                            'numSite',{str2double(tline(41:45))},...
                            'numXform',{str2double(tline(46:50))},...
                            'numCoord',{str2double(tline(51:55))},...
                            'numTer',{str2double(tline(56:60))},...
                            'numConect',{str2double(tline(61:65))},...
                            'numSeq',{str2double(tline(66:70))});
                        
                        %Single/Mandatory
                    case 'END'
                        % Found end of file
                        
                        %Multiple/Optional
                    case 'FTNOTE'
                        FtnoteNo = str2double(tline(7:10));
                        CurFtnote = FtnoteNo;
                        tmpFtnote = sprintf('%d',CurFtnote);
                        
                        if CurFtnote ~= PrevFtnote
                            PDBstruct.(['Footnote' tmpFtnote]) = tline(12:70);
                            PrevFtnote = CurFtnote;
                        else
                            PDBstruct.(['Footnote' tmpFtnote])  = ...
                                strvcat(PDBstruct.(['Footnote' tmpFtnote]) ,tline(12:70)); %#ok
                        end
                        
                    otherwise
                        %disp('The file contains invalid record type');
                        
                end % for the SWITCH statement
                
            end % for the WHILE loop
            
            % If one of the following flags is still true then no 'ENDMDL' (End of
            % Model) line  has been found, this means that the file only contains one
            % model and we still need to assign the temporal structures to the output
            % structure:
            if modelWithAtom
                PDBstruct.Model.Atom = TmpAtomStruct;
            end
            if modelWithAtomSD
                PDBstruct.Model.AtomSD = TmpAtomSDStruct;
            end
            if modelWithAnisoTemp
                PDBstruct.Model.AnisotropicTemp = TmpAnisoTempStruct;
            end
            if modelWithAnisoTempSD
                PDBstruct.Model.AnisotropicTempSD = TmpAnisoTempSDStruct;
            end
            if modelWithTerminal
                PDBstruct.Model.Terminal = TmpTerStruct;
            end
            if modelWithHetAtom
                PDBstruct.Model.HeterogenAtom = TmpHetStruct;
            end
            
            function OutList = GetAtomList(InString)
                OutList = [];
                if ~isempty(InString)
                    try
                        OutList = str2num(reshape(InString, 5, [])');
                    catch allExceptions %#ok<NASGU>
                        OutList = [];
                    end
                end
            end
            
            function OutAcid = GetSequence(InAcid)
                % Residues in pdbfiles are arranged every four columns, there are three
                % cases:
                % 1. Residues are three letter code aminoacids:
                %    Ex: SEQRES 1 A 21 GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU
                % 2. Residues are deoxyribonucleotides:
                %    Ex: SEQRES 1 A 8   DA  DA  DC  DC  DG  DG  DT  DT
                % 3. Residues are ribonucleotides:
                %    Ex: SEQRES 1 X 39   U   C   C   C   C   C   G   U   G   C   C   C   A
                
                if all(isspace(InAcid(2:4:end)))
                    OutAcid = strtrim(InAcid(3:4:end));
                    return
                elseif all(isspace(InAcid(1:4:end)))
                    OutAcid = strtrim(InAcid(3:4:end));
                    return
                end
                
                OutAcid = strrep(upper(InAcid),'ALA','a');
                OutAcid = strrep(OutAcid,'ARG','r');
                OutAcid = strrep(OutAcid,'ASN','n');
                OutAcid = strrep(OutAcid,'ASP','d');
                OutAcid = strrep(OutAcid,'ASX','b');
                OutAcid = strrep(OutAcid,'CYS','c');
                OutAcid = strrep(OutAcid,'GLN','q');
                OutAcid = strrep(OutAcid,'GLU','e');
                OutAcid = strrep(OutAcid,'GLX','z');
                OutAcid = strrep(OutAcid,'GLY','g');
                OutAcid = strrep(OutAcid,'HIS','h');
                OutAcid = strrep(OutAcid,'ILE','i');
                OutAcid = strrep(OutAcid,'LEU','l');
                OutAcid = strrep(OutAcid,'LYS','k');
                OutAcid = strrep(OutAcid,'MET','m');
                OutAcid = strrep(OutAcid,'PHE','f');
                OutAcid = strrep(OutAcid,'PRO','p');
                OutAcid = strrep(OutAcid,'SER','s');
                OutAcid = strrep(OutAcid,'THR','t');
                OutAcid = strrep(OutAcid,'TRP','w');
                OutAcid = strrep(OutAcid,'TYR','y');
                OutAcid = strrep(OutAcid,'VAL','v');
                OutAcid = strrep(OutAcid,'UNK','x');
                OutAcid = regexprep(OutAcid,'[A-Z][A-Z][A-Z]','?');
                OutAcid = upper(OutAcid(~isspace(OutAcid)));
            end
            %---------------------------------------------------------------%
            function [OutStruct,OutNum] = GetResidueStruct(TmpStruct,InString,InNum)
                
                a=1; b=10;
                sz = size(InString);
                while b <= sz(2)
                    test_str = InString(a:b);
                    InNum = InNum + 1;
                    TmpStruct(InNum).ResName = test_str(1:3);
                    TmpStruct(InNum).ChainID = test_str(5);
                    TmpStruct(InNum).ResSeqNo = str2double(test_str(6:9));
                    TmpStruct(InNum).InsCode = test_str(10);
                    a=a+11;
                    b=b+11;
                end
                OutNum = InNum;
                OutStruct = TmpStruct;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function filename = savetempfile(pdbtext)
                
                filename =  [tempname '.spt'];
                fid=fopen(filename,'wb');
                
                rows = size(pdbtext,1);
                
                for rcount = 1:rows-1
                    fprintf(fid,'%s\n',pdbtext(rcount,:));
                end
                
                fprintf(fid,'%s',pdbtext(rows,:));
                fclose(fid);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function val = str2int(str)
                val = sscanf(str,'%d');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function val = str2float(str)
                val = sscanf(str,'%e');
            end
            %------------------------------------------------------------------%
            function Atom = allocateAtoms
                % Initialize Model.Atom, HeterogeneAtom structures
                % % Atom = struct('GlobalSerNo',{0},...
                Atom = struct('AtomSerNo',{0},...
                    'AtomName',{''},...
                    'altLoc',{''},...
                    'resName',{''},...
                    'chainID',{''},...
                    'resSeq',{0},...
                    'iCode',{''},...
                    'X',{0},...
                    'Y',{0},...
                    'Z',{0},...
                    'occupancy',{0},...
                    'tempFactor',{0},...
                    'segID',{''},...
                    'element',{''},...
                    'charge',{''},...
                    'AtomNameStruct',struct('chemSymbol',{''},...
                    'remoteInd',{''},...
                    'branch',{''}));
            end
            %-----------------------------------------------------------------%
            function Terminal = allocateTerminal
                % Initialize Model.Terminal structure
                Terminal = struct('SerialNo',{},...
                    'resName',{},...
                    'chainID',{},...
                    'resSeq',{0},...
                    'iCode',{});
            end
            %-------------------------------------------------------------------%
            function AtomSD = allocateAtomSD
                % Initialize Model.AtomSD
                AtomSD = struct('AtomSerNo',{},...
                    'AtomName',{},...
                    'altLoc',{},...
                    'resName',{},...
                    'chainID',{},...
                    'resSeq',{0},...
                    'iCode',{},...
                    'sigX',{0},...
                    'sigY',{0},...
                    'sigZ',{0},...
                    'sigOcc',{0},...
                    'sigTemp',{0},...
                    'segID',{},...
                    'element',{},...
                    'charge',{},...
                    'AtomNameStruct',struct('chemSymbol',{},...
                    'remoteInd',{},...
                    'branch',{}));
            end
            %-------------------------------------------------------------------%
            function AnisoTemp = allocateAnisoTemp
                % Initialize Model.AnisoTemp
                AnisoTemp = struct('AtomSerNo',{0},...
                    'AtomName',{},...
                    'altLoc',{},...
                    'resName',{},...
                    'chainID',{},...
                    'resSeq',{0},...
                    'iCode',{},...
                    'U00',{0},...
                    'U11',{0},...
                    'U22',{0},...
                    'U01',{0},...
                    'U02',{0},...
                    'U12',{0},...
                    'segID',{},...
                    'element',{},...
                    'charge',{},...
                    'AtomNameStruct',struct('chemSymbol',{},...
                    'remoteInd',{},...
                    'branch',{}));
            end
            %-------------------------------------------------------------------%
            function AnisoTempSD = allocateAnisoTempSD
                % Initialize Model.AnisoTempSD
                AnisoTempSD = struct('AtomSerNo',{0},...
                    'AtomName',{},...
                    'altLoc',{},...
                    'resName',{},...
                    'chainID',{},...
                    'resSeq',{0},...
                    'iCode',{},...
                    'SIG11',{0},...
                    'SIG22',{0},...
                    'SIG33',{0},...
                    'SIG12',{0},...
                    'SIG13',{0},...
                    'SIG23',{0},...
                    'segID',{},...
                    'element',{},...
                    'charge',{},...
                    'AtomNameStruct',struct('chemSymbol',{},...
                    'remoteInd',{},...
                    'branch',{}));
                
            end
            
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
        
        
        
    end
end

