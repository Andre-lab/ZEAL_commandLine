classdef ZEAL < handle
    %ZEAL Class to perform shape-based superposition of proteins.
    %   See description in the ZEAL method below for input details and
    %   usage or the manual at https://andre-lab.github.io/ZEAL/
    %
    %
    %   CITE AS
    %   ZEAL: F. Ljung and I. André, "Protein structure alignment based on
    %   shape similarity", Bioinformatics (2020)
    %
    %
    % Last modified 2021-01
    %
    % TODO: Add event listeners so that ZEAL responds automatically when
    % Settings property is chanegd -> recompute shapeFunction or ZC moments
    % or ZEALscore etc. Useful when we will build the GUI from this class
    %
    % TODO: Unit tests
    %
    % ---------------------------------------------------------------------
    
    properties
        fixed       % holder for data specific for fixed structure: PDB data, ZC moments etc
        rotating    % holder for data specific for rotating structure: PDB data, ZC moments etc
        Settings    % holder for various settings: PDB selection, ZC, shape function
        Search      % holder for data after alignment search
        ZCDdistance = [] % the Eucldiean distance between the shape descriptors of the fixed and rotating structure
        Score = []      % the ZEAL score
    end
    
    properties (Hidden)
        AlignMode % true/false
        % true if ZEAL should align (2 shapes defined), otherwise
        % compute shape descriptors (1 shape defined)
        
        Logging
        AlignLater % true/false (false by default)
        InputParserOnly % true/false ->
        % if true then ZEAL will not start aligning automatically
    end
    
    properties (Constant)
        ChiCoeffs = ChiCoeffs % A singleton  (=only one instance) object (ChiCoeffs class)
        % that handles the loading of pre-computed
        % chi-coefficients from .mat-files
        % The coefficients are visible (like a "global" variable) for all
        % instances of ZEAL, but they are only loaded
        % to memory once. See also the
        % reloadChiCoeffs method in this class.
    end
    
    methods
        function obj = ZEAL(fix, varargin)
            %ZEAL Construct an instance of this class
            %
            % INPUT
            % -------------------------------------------------------------
            % <fix>     :   string
            % filepath to pdb file or 4-letter PDB ID code which ZEAL will
            % try to download. If 4+1 -letter code, then ZEAL assumes last
            % letter is the chain ID that should be extracted
            %
            % OPTIONAL
            % 'name'-value pairs : type (default value)
            % (or given as a structure: the parser expands structures into separate inputs,
            % where each field name corresponds to an input parameter name)
            %
            % 'rot'             : string : same as for <fix>, but defines the
            %                              structure that will be rotated
            %
            % 'Order'           : integer (20)
            % 'GridRes'         : integer (64)
            % 'FunEvals'        : integer (300)
            %       Number of ZEAL score evaluations during search =
            %       the stopping criterium for search
            %
            % 'AlignLater'      : true/false (false)
            % If true then ZEAL will not start searching for best alignment
            % automatically. Manuall start the search using the method
            % shapeAlign(ZEALobject).
            %
            % 'InputParserOnly'       : true/false (false)
            % If true then only the option parsing is performed
            %
            % 'LogLevel'       : 'basic' / 'standard' / 'detailed' / 'none' ('standard')
            % Sets the level of information printed to the console.
            % 'detailed' not recommended unless option 'AlignLater' is set
            % to true;
            %
            % 'ShapeFunction'   : string ('MS') {choose between 'SAS','MS','vdw','electron_density'}
            % MS = molecular surface; SAS = solvent-accessible surface; vdw
            % = van der Waals; electron_density = Gaussian atom
            % representation
            %
            % 'ProbeRadius'     : integer (1.4)
            % 'SmearFactor'     : double (0.3) - for electron_density shape
            % 'ShellThickness'  : integer (1)
            %
            %
            % 'ChiCoeffPath '   : string ([pwd '/chi_coefficients/'])
            %       folder path to precomputed mat files ('chiCoeffs_order_N.mat') with the
            %       object-independent Chi coefficients computed for order N. These can be computed
            %       using ZC.computeChiCoeffs.
            %
            %
            % 'fix_includeHetatoms' :   true/false  (false)
            % 'fix_includeHatoms'   :   true/false  (false)
            % 'fix_chainID'         :   string      ('A')
            % 'fix_altLocID'        :   string      ('A')
            %
            % 'rot_includeHetatoms' :  true/false  (false)
            % 'rot_includeHatoms'   :  true/false  (false)
            % 'rot_chainID'         :  string      ('A')
            % 'rot_altLocID'        :  string       'A'
            %
            %
            % OUTPUT and USAGE
            % -------------------------------------------------------------
            %
            % 1) If only one structure is defined, then ZEAL will only compute the ZC moments
            % and shape descripors for that strucure.
            %
            % Ex:
            % ZCD_5MOK = ZEAL('5MOK');
            %
            % % get the Zernike-Canterakis shape descriptors
            % ZCD_5MOK.fixed.ShapeDescriptors
            %
            %
            % 2) If two structrues are given, then ZEAL will compute the shape similarity
            %    between the structures (the Euclidean distance between their shape descriptors (ZEAL.ZCDdistance)),
            %    and then search for the best superposition (unless 'AlignLater' = true )
            %
            % Ex:
            % aligned_shapeMatch = ZEAL('5MOKA', 'rot', '2HO1A');
            %
            % % get the Euclidean distance between the shape descriptors
            % aligned_shapeMatch.ZCDdistance
            %
            % % get the ZEAL score  for the superposition
            % aligned_shapeMatch.Score
            %
            %
            % ------------ CHANGING DEFAULT SEARCH SETTINGS --------------
            % Ex:
            %
            % - stopping-criteria for the search: change the maximum number of
            % function evaluations
            % ZEAL('5MOKA', 'rot', '2HO1A', 'FunEvals', 500);
            %
            %
            % - shape function to use: change the shape parameter
            % ZEAL('5MOKA', 'rot', '2HO1A', 'Shape', 'electron_density');
            %
            % -------------------------------------------------------------
            
            
            % --- parse input and set parameters ---
            
            % Set default values for parameters
            default_Order = 20;
            default_GridRes = 64;
            
            default_FunEvals = 300;
            
            expected_Shapes = {'SAS','MS','vdw','electron_density'};
            
            default_Shape = 'MS';
            default_ProbeRadius = 1.4;
            default_SmearFactor = 0.3;
            default_ShellThickness = 2;
            
            default_ChiCoeffPath = fullfile(pwd,'chi_coefficients');
            
            default_includeHetatoms = false;
            default_includeHatoms = false;
            default_chainID = 'all';
            default_altLocID = 'A';
            
            expected_ShowLog = {'basic', 'standard', 'detailed', 'none'};
            default_LogOption = 'standard';
            
            default_AlignLater = false;
            default_inputParserOnly = false;
            
            default_rot = ''; % = no rotating structure ->  zeal computes shape descriptors for fixed only
            
            
            % Setup input parser
            p = inputParser;
            
            addRequired(p, 'fix');
            addOptional(p, 'rot', default_rot);
            
            addOptional(p, 'Order', default_Order, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'Order'));
            addOptional(p, 'GridRes', default_GridRes, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', '>=', 32}, 'GridRes'));
            addOptional(p, 'FunEvals', default_FunEvals, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'Order'));
            
            addOptional(p, 'ShapeFunction', default_Shape, @(x) any(validatestring(x,expected_Shapes)));
            addOptional(p, 'ProbeRadius', default_ProbeRadius, @(x)validateattributes(x,{'numeric'}, {'nonempty', 'positive'}, 'ProbeRadius'));
            addOptional(p, 'SmearFactor', default_SmearFactor, @(x)validateattributes(x,{'numeric'}, {'nonempty', 'positive', '>',0, '<', 1}, 'SmearFactor'));
            addOptional(p, 'ShellThickness', default_ShellThickness, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'ShellThickness'));
            
            addOptional(p, 'fix_includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'fix_includeHetatoms'));
            addOptional(p, 'fix_includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'fix_includeHatoms'));
            addOptional(p, 'fix_chainID', default_chainID);
            addOptional(p, 'fix_altLocID', default_altLocID);
            
            addOptional(p, 'rot_includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'rot_includeHetatoms'));
            addOptional(p, 'rot_includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'rot_includeHatoms'));
            addOptional(p, 'rot_chainID', default_chainID);
            addOptional(p, 'rot_altLocID', default_altLocID);
            
            addOptional(p, 'ChiCoeffPath', default_ChiCoeffPath);
            
            addOptional(p, 'AlignLater', default_AlignLater);
            addOptional(p, 'InputParserOnly', default_inputParserOnly);
            
            addOptional(p, 'LogLevel', default_LogOption, @(x) any(validatestring(x,expected_ShowLog)));
            
            % Parse and set defaults
            parse(p, fix, varargin{:});
            
            
            
            % --- Assign parameters ---
            
            obj.Settings.Order = p.Results.Order;
            obj.Settings.GridRes = p.Results.GridRes;
            
            obj.Settings.ChiCoeffPath = p.Results.ChiCoeffPath;
            
            obj.Settings.molShape.FunctionType = p.Results.ShapeFunction;
            obj.Settings.molShape.GridRes = p.Results.GridRes;
            obj.Settings.molShape.ProbeRadius = p.Results.ProbeRadius;
            obj.Settings.molShape.SmearFactor = p.Results.SmearFactor;
            obj.Settings.molShape.ShellThickness = p.Results.ShellThickness;
            
            obj.AlignLater = p.Results.AlignLater;
            obj.InputParserOnly = p.Results.InputParserOnly;
            
            obj.Search.Performed = false;
            obj.Search.EulerAngles = [0 0 0];
            
            obj.Logging.Display = p.Results.LogLevel;
            obj.Settings.molShape.ShowLog = false;
            
            switch obj.Logging.Display
                case 'basic'
                    obj.Logging.Level = 1;
                case 'standard'
                    obj.Logging.Level = 2;
                case 'detailed'
                    obj.Logging.Level = 3;
                    obj.Settings.molShape.ShowLog = true;
                    obj.fixed.Selection.ShowLog = true;
                    obj.rotating.Selection.ShowLog = true;
                case 'none'
                    obj.Logging.Level = 0;
            end
            
            % Set surrogate optimization options
            obj.Settings.SurrOpt = optimoptions('surrogateopt');
            obj.Settings.SurrOpt.MaxFunctionEvaluations = double(p.Results.FunEvals);
            obj.Settings.SurrOpt.PlotFcn = [];
            obj.Settings.SurrOpt.Display = 'none';
            obj.Settings.SurrOpt.OutputFcn = @(a1,a2,a3)surrOptiOutputFcn(a1,a2,a3,obj);
            
            obj.Settings.EulerAngles.LowerBound = [0 0 0];
            obj.Settings.EulerAngles.UpperBound = [2*pi pi 2*pi];
            obj.Settings.EulerAngles.convention = 'zyz';
            
            % Fix
            obj.fixed.Name = p.Results.fix;
            obj.fixed.Selection.includeHetatoms = p.Results.fix_includeHetatoms;
            obj.fixed.Selection.includeHatoms = p.Results.fix_includeHatoms;
            obj.fixed.Selection.chainID = p.Results.fix_chainID;
            obj.fixed.Selection.altLocID = p.Results.fix_altLocID;
            
            if ~obj.InputParserOnly
                
                if isempty(p.Results.rot)
                    obj.AlignMode = false;
                    if obj.Logging.Level > 0
                        fprintf('\n\n Running ZEAL in single mode');
                    end
                else % rot structure given -> align mode activated
                    obj.AlignMode = true;
                    if obj.Logging.Level > 0
                        fprintf('\n\n Running ZEAL in Align mode');
                    end
                end
                
                
                % --- Import structures ---
                if obj.Logging.Level > 0
                    fprintf('\n - Importing fixed structure: %s', obj.fixed.Name);
                end
                obj.fixed.PDB = PDB(obj.fixed.Name, obj.fixed.Selection);
                
                obj.fixed.Rg = computeRadiusOfGyration(obj, obj.fixed.PDB.Data);
                
                if  obj.AlignMode
                    % Set PDB parameters for rotating if align mode activated
                    obj.rotating.Name = p.Results.rot;
                    obj.rotating.Selection.includeHetatoms = p.Results.rot_includeHetatoms;
                    obj.rotating.Selection.includeHatoms = p.Results.rot_includeHatoms;
                    obj.rotating.Selection.chainID = p.Results.rot_chainID;
                    obj.rotating.Selection.altLocID = p.Results.rot_altLocID;
                    
                    if obj.Logging.Level > 0
                        fprintf('\n - Importing rotating structure: %s', obj.rotating.Name);
                    end
                    
                    obj.rotating.PDB = PDB(obj.rotating.Name, obj.rotating.Selection);
                    
                    obj.rotating.Rg = computeRadiusOfGyration(obj, obj.rotating.PDB.Data);
                    
                end
                
                % --- Compute shape functions ---
                if obj.Logging.Level > 0
                    fprintf('\n - Computing shape function for fixed structure');
                end
                
                obj.fixed.molShape = molShape(obj.fixed.PDB.Data, obj.Settings.molShape);
                
                if obj.AlignMode
                    if obj.Logging.Level > 0
                        fprintf('\n - Computing shape function for rotating structure');
                    end
                    obj.rotating.molShape = molShape(obj.rotating.PDB.Data, obj.Settings.molShape);
                end
                
                
                % --- Computing ZC moments ---
                
                if obj.Logging.Level > 0
                    fprintf('\n - Computing ZC moments for fixed structure');
                end
                
                % Reload Chi Coeffs. if expansion order differs from order in
                % singelton property (ChiCoeffs)
                if ~isequal(obj.Settings.Order, obj.ChiCoeffs.Order)
                    obj.reloadChiCoeffs('Order',obj.Settings.Order);
                end
                
                obj.fixed.ZC = ZC(obj.fixed.molShape.FunctionData, obj.Settings.Order, obj.ChiCoeffs, 'ShowLog', obj.Settings.molShape.ShowLog);
                computeDescriptors(obj.fixed.ZC);
                
                if obj.AlignMode
                    if obj.Logging.Level > 0
                        fprintf('\n - Computing ZC moments for rotating structure');
                    end
                    
                    obj.rotating.ZC = ZC(obj.rotating.molShape.FunctionData, obj.Settings.Order, obj.ChiCoeffs, 'ShowLog', obj.Settings.molShape.ShowLog);
                    
                    computeDescriptors(obj.rotating.ZC);
                    
                    % Compute the Euclidean distance between shape descriptors
                    computeZCDdistance(obj);
                    
                    % Compute the initial ZEAL score
                    obj.Search.InitialScore = computeScore(obj);
                    
                    if ~obj.AlignLater
                        % Perform shape alignment
                        shapeAlign(obj);
                    end
                end
                
                fprintf('\n');
                
            end
            
        end
        
        function shapeAlign(obj)
            %shapeAlign Start searching for best shape alingment using
            % the surrogate optimization algorithm.
            
            [   obj.Search.SurrOptOut.EulerAngles, ...
                obj.Search.SurrOptOut.Score, ...
                obj.Search.SurrOptOut.exitflag, ...
                obj.Search.SurrOptOut.output, ...
                obj.Search.SurrOptOut.trials] = surrogateopt(@objFun, ...
                obj.Settings.EulerAngles.LowerBound, ...
                obj.Settings.EulerAngles.UpperBound, ...
                obj.Settings.SurrOpt);
            
            obj.Score = -1*obj.Search.SurrOptOut.Score;
            
            obj.Search.Performed = true;
            obj.Search.EulerAngles = obj.Search.SurrOptOut.EulerAngles;
            obj.Search.Score = -1*obj.Search.SurrOptOut.Score;
            
            function [d] = objFun(x)
                % objFun The objective function, i.e. the "expensive" ZEAL
                % score function
                
                % get rotation matrix from Euler angles
                R = ZEAL.euler2rotMat(x);
                
                % rotate orignal atom coordinates
                atomlistRot_rotated = obj.rotating.molShape.AtomData;
                atomlistRot_rotated(:,1:3) = atomlistRot_rotated(:,1:3) * R;
                
                % Compute new shape function
                updateShapeFunction(obj, atomlistRot_rotated)
                
                % Compute new ZC moments
                computeMoments(obj.rotating.ZC)
                
                % Compute new ZEALscore
                d = -1*computeScore(obj);
                
            end
            
            function updateShapeFunction(obj, atomList)
                % updateShapeFunction Computes new shape function from
                % transformed atom coordinates (atomList)
                
                switch obj.rotating.molShape.FunctionType
                    
                    case 'MS' % molecular surface
                        [obj.rotating.ZC.ShapeFunction, ~] = molShape.MSsurface(atomList, obj.Settings.molShape.GridRes, obj.Settings.molShape.ProbeRadius, obj.Settings.molShape.ShellThickness);
                        
                    case 'SAS' % solvent accessible surface
                        [obj.rotating.ZC.ShapeFunction, ~] = molShape.SAsurface(atomList, obj.Settings.molShape.GridRes, obj.Settings.molShape.ProbeRadius, obj.Settings.molShape.ShellThickness);
                        
                    case 'vdw' % van der Waals surface
                        probeRadius = 0;
                        [obj.rotating.ZC.ShapeFunction, ~] = molShape.SAsurface(atomList, obj.Settings.molShape.GridRes, probeRadius, obj.Settings.molShape.ShellThickness);
                        
                    otherwise % electron density
                        [obj.FunctionData, ~, ~] = molShape.electronDensity(atomLista, obj.Settings.molShape.GridRes, obj.Settings.molShape.SmearFactor);
                        
                end
                
            end
            
        end
        
        function  stop = surrOptiOutputFcn(x,optimvalues,state, obj)
            % surrOptiOutputFcn Output function that is called
            % during the surrogate optimization.
            
            stop = false;
            
            switch state
                case 'init' % surrogateopt is starting up
                    
                    if obj.Logging.Level > 0
                        fprintf('\n\n /////////////////////////////////////////////////////////////////////////////\n');
                        
                        fprintf('\n\tSearching for best shape superposition\n');
                        fprintf('\n\tStopping after %d function evaluations\n', obj.Settings.SurrOpt.MaxFunctionEvaluations);
                    end
                    
                    if obj.Logging.Level > 1
                        fprintf('\n ----------------------------------------------------------------------------');
                        fprintf('\n         Best score      Euler (zyz)         iteration       time (s) ');
                        fprintf('\n ----------------------------------------------------------------------------');
                    end
                    
                    obj.Search.History.Iteration = [];
                    obj.Search.History.Score = [];
                    obj.Search.History.EulerAngles = [];
                    
                    obj.Search.History.prev_val = 0;
                    obj.Search.History.startTime = tic;
                    
                case 'iter' % surrogateopt is running
                    
                    searchTime = toc(obj.Search.History.startTime);
                    
                    if optimvalues.fval < (obj.Search.History.prev_val*1.01)
                        
                        obj.Search.History.prev_val = optimvalues.fval;
                        
                        obj.Search.History.Iteration = [obj.Search.History.Iteration optimvalues.iteration];
                        obj.Search.History.EulerAngles = [obj.Search.History.EulerAngles; x];
                        obj.Search.History.Score = [obj.Search.History.Score -1*optimvalues.fval];
                        
                        obj.Search.History.BestScore = -1*optimvalues.fval;
                        obj.Search.History.BestTime = searchTime;
                        obj.Search.History.BestIteration = optimvalues.iteration;
                        obj.Search.History.BestEulerAngles = x;
                        
                        if obj.Logging.Level > 1
                            fprintf('\n\t %2.2f \t\t %2.2f %2.2f %2.2f \t %d \t\t %4.1f', obj.Search.History.BestScore, x(1), x(2), x(3), obj.Search.History.BestIteration, obj.Search.History.BestTime);
                        end
                    end
                    
                case 'done' % surrogateopt is finished
                    
                    searchTime = toc(obj.Search.History.startTime);
                    
                    obj.Search.History.Iteration = [obj.Search.History.Iteration optimvalues.iteration];
                    obj.Search.History.EulerAngles = [obj.Search.History.EulerAngles; x];
                    obj.Search.History.Score = [obj.Search.History.Score -1*optimvalues.fval];
                    
                    if obj.Logging.Level > 0
                        fprintf('\n ----------------------------------------------------------------------------');
                        
                        fprintf('\n\n\tSearch completed after %3.1f s.\n', searchTime);
                        fprintf('\n\tBest score %2.2f found after %d iterations (%3.1f s)\n', obj.Search.History.BestScore, obj.Search.History.BestIteration, obj.Search.History.BestTime);
                        fprintf('\n\tusing Euler angles (zyz) [%3.3f %3.3f %3.3f]', obj.Search.History.BestEulerAngles(1), obj.Search.History.BestEulerAngles(2), obj.Search.History.BestEulerAngles(3));
                        
                        fprintf('\n\n /////////////////////////////////////////////////////////////////////////////\n');
                    end
                    
                    stop = true;
            end
            
        end
        
        function reloadChiCoeffs(obj, varargin)
            % reloads chicoeffs using the master class ChiCoeff
            obj.ChiCoeffs.loadData(obj, varargin{:}); %all instances of ZEAL will see the newly loaded chi coeffs
        end
        
        
        function Score = computeScore(obj)
            % computeScore Compute the ZEAL score for the shape superposition.
            % The ZEAL score is the ZC moment correlation between the fixed and the
            % rotating structure, computed as the cosine of the angle
            % between the N-dimensional moment-vectors. Thus, cos(small
            % angle) corrsponds to high moment-correlation = the ZEAL
            % score. Because the moments are complex, the can angle can be defined
            % differently. Here we use the so called Euclidean angle, but a so called
            % Hermitian angle can also be defined.
            
            % ZEAL score = Cos(Euclidean angle)
            Score = real(dot(obj.rotating.ZC.Moments.Values, obj.fixed.ZC.Moments.Values)) / (norm(obj.rotating.ZC.Moments.Values)*norm(obj.fixed.ZC.Moments.Values));
            
            % Hermitian angle
            %         d = -1*abs(dot(ZC_mom_list, ZCref.mom_list)) / (norm(ZC_mom_list)*norm(ZCref.mom_list));
            
        end
        
        function computeZCDdistance(obj)
            % computeZCDdistance Computes the Euclidean distance between
            % shape descriptors of the fixed and rotating structure
            %
            % The distance is computed as the 2-norm of the difference
            % between the vector elements containing the shape-descriptor
            % values.
            
            if isempty(obj.fixed.ZC.Descriptors)
                computeDescriptors(obj.fixed.ZC);
            end
            
            if isempty(obj.rotating.ZC.Descriptors)
                computeDescriptors(obj.rotating.ZC);
            end
            
            obj.ZCDdistance = norm(obj.fixed.ZC.Descriptors -  obj.rotating.ZC.Descriptors, 2);
            
        end
        
        function saveSuccessFul = save2pdb(obj,varargin)
            % save2pdb Export strucure to pdb file
            % By default, if both structures exist in the object (fixed and rotated)
            % then those are exported to pdbfiles in the current folder.
            % The fixed strucure is exported if ZEAL is run in 'single
            % structure' mode (i.e. shape analysis of one structure = the fixed structure).
            
            % OPTIONAL INPUT
            %'name'-value pairs : type : default value {expected}
            %
            % 'structure'       :   char : 'all'  {'fixed','rotating','all'}
            %
            % 'includeHetatoms' :   logical : false {true/false}
            %
            % 'includeAll'      :   logical : false {true/false}
            %
            % 'folderPath'      :   char : current directory
            %
            % 'fixName'         :   char  : <source_name>_ZEAL.pdb
            %
            % 'rotName'         :   char  : <source_name>_ZEAL.pdb
            
            
            expectedStructures = {'fixed', 'rotating', 'all'};
            
            % set defaults
            defaultPath = pwd;
            
            defaultHetatoms = false;
            defaultAll = false;
            
            defaultFixName = [obj.fixed.PDB.Source.Name '_ZEAL.pdb'];
            
            if ~isempty(obj.rotating)
                defaultStructure = 'all';
                defaultRotName = [obj.rotating.PDB.Source.Name '_ZEAL.pdb'];
            else
                defaultStructure = 'fixed';
                defaultRotName = '';
            end
            
            % do option parsing
            p = inputParser;
            
            addOptional(p, 'structure', defaultStructure, @(x) any(validatestring(x,expectedStructures)));
            addOptional(p, 'includeHetatoms', defaultHetatoms);
            addOptional(p, 'includeAll', defaultAll);
            
            addOptional(p, 'fixName', defaultFixName);
            addOptional(p, 'rotName', defaultRotName);
            
            addOptional(p, 'folderPath', defaultPath);
            
            parse(p, varargin{:});
            
            selectedStructure = p.Results.structure;
            includeHetatoms = p.Results.includeHetatoms;
            includeAll = p.Results.includeAll;
            folderpath = p.Results.folderPath;
            
            fixName = p.Results.fixName;
            rotName = p.Results.rotName;
            
            saveSuccessFul = false;
            
            try
                
                if strcmp(selectedStructure,'fixed') ||  strcmp(selectedStructure,'all')
                    
                    selection = obj.fixed.PDB.Selection;
                    
                    if includeHetatoms || includeAll
                        selection.includeHetatoms = includeHetatoms;
                    end
                    
                    if includeAll
                        selection.chainID = 'all';
                        selection.altLocID = 'all';
                    end
                    
                    pdbData = PDB.parsePDBstruct(obj.fixed.PDB.AllData, selection);
                    
                    T = getTranslationMatrix(obj, 'fixed');
                    R = getRotationMatrix(obj, 'fixed');
                    
                    xyz = [pdbData.X pdbData.Y pdbData.Z];
                    
                    xyzRot = [ xyz ones(length(xyz),1)] * T * R;
                    
                    pdbData.X = xyzRot(:,1);
                    pdbData.Y = xyzRot(:,2);
                    pdbData.Z = xyzRot(:,3);
                    
                    filePath=fullfile(folderpath, fixName);
                    
                    exportPdb(obj, 'fixed', pdbData, filePath)
                    
                end
                
                
                if strcmp(selectedStructure,'rotating') ||  strcmp(selectedStructure,'all')
                    
                    selection = obj.rotating.PDB.Selection;
                    
                    if includeHetatoms || includeAll
                        selection.includeHetatoms = includeHetatoms;
                    end
                    
                    if includeAll
                        selection.chainID = 'all';
                        selection.altLocID = 'all';
                    end
                    
                    pdbData = PDB.parsePDBstruct(obj.rotating.PDB.AllData, selection);
                    
                    T = getTranslationMatrix(obj, 'rotating');
                    R = getRotationMatrix(obj, 'rotating');
                    
                    xyz = [pdbData.X pdbData.Y pdbData.Z];
                    
                    xyzRot = [ xyz ones(length(xyz),1)] * T * R;
                    
                    pdbData.X = xyzRot(:,1);
                    pdbData.Y = xyzRot(:,2);
                    pdbData.Z = xyzRot(:,3);
                    
                    filePath=fullfile(folderpath, rotName);
                    
                    exportPdb(obj, 'rotating', pdbData, filePath)
                    
                end
                
                saveSuccessFul = true;
                
            catch ME
                
                rethrow(ME)
                
            end
            
        end
        
        function T = getTranslationMatrix(obj, structure)
            
            switch structure
                
                case 'fixed'
                    XYZ = [obj.fixed.PDB.Data.X obj.fixed.PDB.Data.Y obj.fixed.PDB.Data.Z];
                case 'rotating'
                    XYZ = [obj.rotating.PDB.Data.X obj.rotating.PDB.Data.Y obj.rotating.PDB.Data.Z];
            end
            
            COM = mean(XYZ);
            
            T = eye(4);
            T(end, :) = [-COM 1];
            
        end
        
        
        function R = getRotationMatrix(obj, structure)
            
            R = eye(4);
            
            switch structure
                
                case 'fixed'
                    
                    R(end,end) = 0;
                    
                case 'rotating'
                    
                    R(1:3,1:3) = ZEAL.euler2rotMat(obj.Search.EulerAngles);
            end
            
        end
        
        function exportPdb(obj, structure, pdbData, fileSavePath)
            % Save pdbdata to pdbfile with ZEAL info in REMARK fields (if 2 structures loaded)
            %
            % INPUT
            %                   typr                 expected value
            % <structure>     :  char       :     {'fixed', 'rotating'}
            % Structure to save
            %
            % <pdbData>       : struct      :     structure from PDB class
            %
            % <fileSavePath>  : char
            % path to the file that we create and write data data
            
            
            % get file name from path
            [~, name, ~] = fileparts(fileSavePath);
            
            fprintf('outputting PDB in file %s.pdb ...', name);
            
            % create and open file
            fid = fopen(fileSavePath, 'w');
            if fid == -1
                error('Author:Function:OpenFile', 'Cannot open file: %s', name);
            end
            
            % Write
            if ~isempty(obj.fixed) && ~isempty(obj.rotating)
                
                % write settings used for the shape alignment to REMARK record
                fprintf(fid, '%-12s ZEAL Protein shape alignment %s\n', 'REMARK' , datestr(now));
                fprintf(fid, '%-12s \n', 'REMARK');
                
                if strcmp(structure,'fixed')
                    fprintf(fid, '%-12s ZEAL Fixed (this file): %s CHAIN %s\n', 'REMARK' , obj.fixed.Name, obj.fixed.PDB.Selection.chainID);
                    fprintf(fid, '%-12s ZEAL Rotating: %s CHAIN %s\n', 'REMARK' , obj.rotating.Name, obj.fixed.PDB.Selection.chainID);
                elseif strcmp(structure,'rotating')
                    fprintf(fid, '%-12s ZEAL Rotating (this file): %s CHAIN %s\n', 'REMARK' , obj.rotating.Name, obj.rotating.PDB.Selection.chainID);
                    fprintf(fid, '%-12s ZEAL Fixed: %s CHAIN %s\n', 'REMARK' , obj.fixed.Name, obj.fixed.PDB.Selection.chainID);
                end
                
                fprintf(fid, '%-12s ZEAL score: %5.5f \n', 'REMARK' , obj.Score);
                fprintf(fid, '%-12s \n', 'REMARK');
                fprintf(fid, '%-12s ZEAL SETTINGS\n', 'REMARK');
                fprintf(fid, '%-12s ZEAL Maximum order of Zernike-Canterakis moments: %d\n', 'REMARK', obj.Settings.Order);
                fprintf(fid, '%-12s ZEAL Grid resolution:\t%d\n', 'REMARK', obj.Settings.GridRes);
                fprintf(fid, '%-12s ZEAL Shape function:\t%s\n', 'REMARK', obj.Settings.molShape.FunctionType);
                
                if strcmp(obj.Settings.molShape.FunctionType, 'electron_density')
                    fprintf(fid, '%-12s ZEAL \t\t\t\t Smear factor: %2.2f\n', 'REMARK', obj.Settings.molShape.SmearFactor);
                else
                    fprintf(fid, '%-12s ZEAL \t\t\t\t Probe radius: %2.2f Å\n', 'REMARK', obj.Settings.molShape.ProbeRadius);
                    fprintf(fid, '%-12s ZEAL \t\t\t\t Shell thickness: %d\n', 'REMARK', obj.Settings.molShape.ShellThickness);
                end
                
                fprintf(fid, '%-12s \n', 'REMARK');
                
                fprintf(fid, '%-12s The 4x4 affine transformation matrix (Tmatrix)\n', 'REMARK');
                fprintf(fid, '%-12s to get new coordinates (this file) from original (orig) coordinates\n', 'REMARK');
                
                fprintf(fid, '%-12s \n', 'REMARK');
                
                T = getTranslationMatrix(obj, structure);
                R = getRotationMatrix(obj, structure);
                
                transMat = T*R;
                
                fprintf(fid, '%-12s Tmatrix = [\n', 'REMARK');
                
                fprintf(fid, '%-12s \t\t %5.5f %5.5f %5.5f %5.5f\n', 'REMARK', transMat(1,1), transMat(1,2), transMat(1,3), transMat(1,4));
                fprintf(fid, '%-12s \t\t %5.5f %5.5f %5.5f %5.5f\n', 'REMARK', transMat(2,1), transMat(2,2), transMat(2,3), transMat(2,4));
                fprintf(fid, '%-12s \t\t %5.5f %5.5f %5.5f %5.5f\n', 'REMARK', transMat(3,1), transMat(3,2), transMat(3,3), transMat(3,4));
                fprintf(fid, '%-12s \t\t %5.5f %5.5f %5.5f %5.5f\n', 'REMARK', transMat(4,1), transMat(4,2), transMat(4,3), transMat(4,4));
                fprintf(fid, '%-12s \t\t\t ] \n', 'REMARK');
                fprintf(fid, '%-12s \n', 'REMARK');
                
                fprintf(fid, '%-12s [ new_X new_Y new_Z ] = [ orig_X orig_Y orig_Z 1 ] * Tmatrix \n', 'REMARK');
                
                fprintf(fid, '%-12s \n', 'REMARK');
                fprintf(fid, '%-12s \n', 'REMARK');
                
                
                fprintf(fid, '%-12s Cite: \n', 'REMARK');
                fprintf(fid, '%-12s F. Ljung and I. André, \n', 'REMARK');
                fprintf(fid, '%-12s ZEAL: Structure alignment based on shape similarity, Bioinformatics (2020) \n', 'REMARK');
                fprintf(fid, '%-12s Bioinformatics (2020) \n', 'REMARK');
                
                fprintf(fid, '%-12s \n', 'REMARK');
                
            end
            
            % check integrity of some data fields
            nRecords = length(pdbData.atomNum);
            
            checkIfEmpty = @(x) isempty(x);
            
            % occupancy
            if isempty(pdbData.occupancy)
                pdbData.occupancy = zeros(nRecords,1);
            end
            
            % betaFactor
            if isempty(pdbData.betaFactor)
                pdbData.betaFactor = zeros(nRecords,1);
            end
            
            % element
            if nRecords == sum(cellfun(checkIfEmpty, pdbData.element))
                pdbData.element = cell(nRecords,1);
            end
            
            % charge
            if nRecords == sum(cellfun(checkIfEmpty, pdbData.charge))
                pdbData.charge = cell(nRecords,1);
            end
            
            ZEAL.writeModel(fid, pdbData)
            
            fprintf( fid, 'END\n');
            
            
            % close file
            fprintf('\n done! closing file...\n');
            fclose(fid);
            
        end
        
        function obj = mergeZEAL(fixObj, rotObj)
            % method to merge two ZEAL objects in "single mode" to a new
            % object in "align mode", i.e. a fixed and rotating structure
            % that can be aligned. The shape similarity (Euclidean ZCD
            % distance) is after the merge.
            %
            % Example:
            % Create two ZEAL objects in single mode:
            % a = ZEAL('1stmA');
            % b = ZEAL('2lisA')
            %
            % Merge them to a ZEAL object in align mode:
            % ab = mergeZEAL(a,b)
            %
            % Now you can do shape alignment of them using the shapeAlign
            % command:
            % shapeAlign(ab)
            
            obj = fixObj;
            obj.fixed = fixObj.fixed;
            obj.rotating = rotObj.fixed;
            
            computeZCDdistance(obj);
            
        end
        
        function ZCD = getShapeDescriptors(obj, varargin)
            
            if isempty(varargin)
                structure = 'fixed';
            else
                structure = varargin{1};
            end
            
            switch structure
                
                case 'fixed'
                    ZCD = obj.fixed.ZC.Descriptors;
                case 'rotating'
                    ZCD = obj.rotating.ZC.Descriptors;
            end
            
        end
        
        function ZCmoments = getMoments(obj, varargin)
            
            if isempty(varargin)
                structure = 'fixed';
            else
                structure = varargin{1};
            end
            
            switch structure
                
                case 'fixed'
                    ZCmoments = obj.fixed.ZC.Moments.Values;
                case 'rotating'
                    ZCmoments = obj.rotating.ZC.Moments.Values;
            end
            
        end
        
        function Rg = computeRadiusOfGyration(~,pdb_data)
            % Compute radius of gyration (Rg)
            % Input: the data field from the PDB class
            
            xyz = [pdb_data.X pdb_data.Y pdb_data.Z];
            natoms = size(xyz, 1);
            COM = mean(xyz);
            xyz_com = xyz-ones(natoms,1)*COM;
            
            Rg = sqrt(mean(sum(xyz_com.^2,2)));
        end
        
    end
    
    
    methods (Static)
        
        function Rmat = euler2rotMat(x)
            % euler2rotMat Get the rotation matrix corresponding to the
            % rotation parameterized by the Euler angles x and rotation
            % order (Euler convention) 'zyz'
            
            % Functions for Euler angles to rotation matrix
            Ry = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
            Rz = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
            
            Rmat = (Rz(x(3))*Ry(x(2))*Rz(x(1)))';
            
        end
        
        %         function writeModel(fid, model)
        %             % Write PDB structure data to file (id = fid). Structure format is that from PDB class.
        %             % Helper function to export2Pdb method in case we want to write
        %             % multiple models to same file
        %
        %             nAtoms = length(model.atomNum);
        %
        %             % output data
        %             try
        %                 for n = 1:nAtoms
        %
        %                     % fix atomName spacing
        %                     model.atomName(n) = {sprintf('%-3s',cell2mat(model.atomName(n)))};
        %
        %                     % standard PDB output line
        %                     fprintf( fid, '%-6s%5u%5s%1.1s%3s %1.1s%4i%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
        %                         cell2mat(model.recordName(n)), model.atomNum(n), cell2mat(model.atomName(n)), ...
        %                         cell2mat(model.altLoc(n)), cell2mat(model.resName(n)), cell2mat(model.chainID(n)), ...
        %                         model.resNum(n), model.X(n), model.Y(n), model.Z(n), model.occupancy(n), model.betaFactor(n), ...
        %                         cell2mat(model.element(n)), cell2mat(model.charge(n)));
        %                 end
        %
        %             catch
        %                 error('Failed to write PDB records.');
        %             end
        %
        %         end
        
        function writeModel(fid, pdbData)
            
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
            
            
            
        end
        
    end
    
end

