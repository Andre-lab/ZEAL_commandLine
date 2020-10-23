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
    % Last modified 2020-10-22
    %
    % TODO: Add event listeners so that ZEAL responds automatically when
    % Settings property is chanegd -> recompute shapeFunction or ZC moments
    % or ZEALscore etc. Useful when we will build the GUI from this class
    %
    % TODO: save2PDB method
    %
    % TODO: Unit tests
    %
    % ---------------------------------------------------------------------
    
    properties
        fixed       % holder for data specific for fixed structure: PDB data, ZC moments etc
        rotating    % holder for data specific for rotating structure: PDB data, ZC moments etc
        Settings    % holder for various settings: PDB selection, ZC, shape function 
        Search      % holder for data after alignment search 
        ZCDdistance % the Eucldiean distance between the shape descriptors of the fixed and rotating structure
        Score       % the ZEAL score
    end
    
    properties (Hidden)
        AlignMode % true/false
                    % true if ZEAL should align (2 shapes defined), otherwise
                    % compute shape descriptors (1 shape defined)
        ChiCoeffs   % holder for values and indices for chi coeffs
        Logging
        AlignLater % true/false (false by default)
                   % if true then ZEAL will not start aligning automatically                   
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
            % 'name'-value pairs (default value)
            % (or given as a structure: the parser expands structures into separate inputs,
            % where each field name corresponds to an input parameter name)
            %
            % 'rot'             : string : same as for <fix>, but defines the
            %                              structure that will be rotated
            %
            % 'Order'           : integer (20)
            % 'GridRes'         : integer (64)
            %
            % 'AlignLater'      : true/false (false) 
            % If true then ZEAL will not start searching for best alignment
            % automatically. Manuall start the search using the method
            % shapeAlign(ZEALobject).
            %
            % 'LogOption'       : 'basic' / 'standard' / 'detailed' / 'none' ('standard')
            % Sets the level of information printed to the console.
            % 'detailed' not recommended unless option 'AlignLater' is set
            % to true; 
            %
            % 'Shape'           : string ('MS') {choose between 'SAS','MS','vdw','electron_density'}
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
            % 'fix_modelNumber'     :   integer     1
            %
            % 'rot_includeHetatoms' :  true/false  (false)
            % 'rot_includeHatoms'   :  true/false  (false)
            % 'rot_chainID'         :  string      ('A')
            % 'rot_altLocID'        :  string       'A'
            % 'rot_modelNumber'     :  integer       1
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
            default_ShellThickness = 1;
            
            default_ChiCoeffPath = fullfile(pwd,'chi_coefficients');
            
            default_includeHetatoms = false;
            default_includeHatoms = false;
            default_chainID = 'A';
            default_altLocID = 'A';
            default_modelID = 1;
            
            expected_ShowLog = {'basic', 'standard', 'detailed', 'none'};
            default_LogOption = 'standard'; 
            
            default_AlignLater = false;
            
            default_rot = ''; % = no rotating structure ->  zeal computes shape descriptors for fixed only
            
            
            % Setup input parser
            p = inputParser;
            
            addRequired(p, 'fix');
            addOptional(p, 'rot', default_rot);
            
            addOptional(p, 'Order', default_Order, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'Order'));
            addOptional(p, 'GridRes', default_GridRes, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', '>=', 32}, 'GridRes'));
            addOptional(p, 'FunEvals', default_FunEvals, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'Order'));
            
            addOptional(p, 'Shape', default_Shape, @(x) any(validatestring(x,expected_Shapes)));
            addOptional(p, 'ProbeRadius', default_ProbeRadius, @(x)validateattributes(x,{'numeric'}, {'nonempty', 'positive'}, 'ProbeRadius'));
            addOptional(p, 'SmearFactor', default_SmearFactor, @(x)validateattributes(x,{'numeric'}, {'nonempty', 'positive', '>',0, '<', 1}, 'SmearFactor'));
            addOptional(p, 'ShellThickness', default_ShellThickness, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'ShellThickness'));
            
            addOptional(p, 'fix_includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'fix_includeHetatoms'));
            addOptional(p, 'fix_includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'fix_includeHatoms'));
            addOptional(p, 'fix_chainID', default_chainID);
            addOptional(p, 'fix_altLocID', default_altLocID);
            addOptional(p, 'fix_modelNumber', default_modelID, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'fix_modelNumber'));
            
            addOptional(p, 'rot_includeHetatoms', default_includeHetatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'rot_includeHetatoms'));
            addOptional(p, 'rot_includeHatoms', default_includeHatoms, @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'rot_includeHatoms'));
            addOptional(p, 'rot_chainID', default_chainID);
            addOptional(p, 'rot_altLocID', default_altLocID);
            addOptional(p, 'rot_modelNumber', default_modelID, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer', 'positive'}, 'rot_modelNumber'));
            
            addOptional(p, 'ChiCoeffPath', default_ChiCoeffPath);
            
            addOptional(p, 'AlignLater', default_AlignLater);
            addOptional(p, 'LogLevel', default_LogOption, @(x) any(validatestring(x,expected_ShowLog)));
            
            % Parse and set defaults
            parse(p, fix, varargin{:});
            
            if isempty(p.Results.rot) 
                obj.AlignMode = false;
                fprintf('\nRunning ZEAL in single mode');
            else % rot structure given -> align mode activated
                obj.AlignMode = true;
                fprintf('\nRunning ZEAL in Align mode');
            end
            
            % --- Assign parameters ---
            
            obj.Settings.Order = p.Results.Order;
            obj.Settings.GridRes = p.Results.GridRes;
            
            obj.Settings.ChiCoeffPath = p.Results.ChiCoeffPath;
            
            obj.Settings.molShape.Shape = p.Results.Shape;
            obj.Settings.molShape.GridRes = p.Results.GridRes;
            obj.Settings.molShape.ProbeRadius = p.Results.ProbeRadius;
            obj.Settings.molShape.SmearFactor = p.Results.SmearFactor;
            obj.Settings.molShape.ShellThickness = p.Results.ShellThickness;
            
            obj.AlignLater = p.Results.AlignLater;
            
            obj.Search.Performed = false;
            
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
            
            
            % Fix
            obj.fixed.Name = p.Results.fix;
            obj.fixed.Selection.includeHetatoms = p.Results.fix_includeHetatoms;
            obj.fixed.Selection.includeHatoms = p.Results.fix_includeHatoms;
            obj.fixed.Selection.chainID = p.Results.fix_chainID;
            obj.fixed.Selection.altLocID = p.Results.fix_altLocID;
            obj.fixed.Selection.modelNumber = p.Results.fix_modelNumber;                        
            
            % --- Import structures ---
            if obj.Logging.Level > 0
                fprintf('\n Importing fixed structure: %s', obj.fixed.Name);
            end
            obj.fixed.PDB = PDB(obj.fixed.Name, obj.fixed.Selection);
            
            
            if  obj.AlignMode
                % Set PDB parameters for rotating if align mode activated
                obj.rotating.Name = p.Results.rot;
                obj.rotating.Selection.includeHetatoms = p.Results.rot_includeHetatoms;
                obj.rotating.Selection.includeHatoms = p.Results.rot_includeHatoms;
                obj.rotating.Selection.chainID = p.Results.rot_chainID;
                obj.rotating.Selection.altLocID = p.Results.rot_altLocID;
                obj.rotating.Selection.modelNumber = p.Results.rot_modelNumber;
                
                if obj.Logging.Level > 0
                    fprintf('\n Importing rotating structure: %s', obj.rotating.Name);
                end
                
                obj.rotating.PDB = PDB(obj.rotating.Name, obj.rotating.Selection);
                
                % Set surrogate optimization options
                obj.Settings.SurrOpt = optimoptions('surrogateopt');
                obj.Settings.SurrOpt.MaxFunctionEvaluations = p.Results.FunEvals;
                obj.Settings.SurrOpt.PlotFcn = [];
                obj.Settings.SurrOpt.Display = 'none';
                obj.Settings.SurrOpt.OutputFcn = @(a1,a2,a3)surrOptiOutputFcn(a1,a2,a3,obj);
                
                obj.Settings.EulerAngles.LowerBound = [0 0 0];
                obj.Settings.EulerAngles.UpperBound = [2*pi pi 2*pi];
                obj.Settings.EulerAngles.convention = 'zyz';
                
            end
            
            % --- Compute shape functions ---
            if obj.Logging.Level > 0
                fprintf('\n Computing shape function for fixed structure');
            end
            
            obj.fixed.molShape = molShape(obj.fixed.PDB.Data, obj.Settings.molShape);
            
            if obj.AlignMode
                if obj.Logging.Level > 0
                    fprintf('\n Computing shape function for rotating structure');
                end
                obj.rotating.molShape = molShape(obj.rotating.PDB.Data, obj.Settings.molShape);
            end
            
            
            % --- Computing ZC moments ---
            loadChiCoeffs(obj);
            
            if obj.Logging.Level > 0
                fprintf('\n Computing ZC moments for fixed structure');
            end
            
            obj.fixed.ZC = ZC(obj.fixed.molShape.FunctionData, obj.Settings.Order, obj.ChiCoeffs, 'ShowLog', obj.Settings.molShape.ShowLog);
            computeDescriptors(obj.fixed.ZC);
            
            if obj.AlignMode
                if obj.Logging.Level > 0
                    fprintf('\n Computing ZC moments for rotating structure\n');
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
                atomlistRot_rotated(:,1:3) = atomlistRot_rotated(:,1:3) * R';
                
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
                        fprintf('\n Current best score      Euler (zyz)         iteration       time (s) ');
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
        
        function loadChiCoeffs(obj)
            % loadChiCoeffs Loads the object-independent chi coefficients
            % used for the ZC moment computation. The folder path is taken from 
            % obj.Settings.ChiCoeffPath and the files themself are assumed to
            % be mat-files with name 'chiCoeffs_order_X.mat', where X is the
            % ZC order.
            
            chiCoeffFilename = sprintf('chiCoeffs_order_%d.mat', obj.Settings.Order);
            
            chiCoeffDataPath = fullfile(obj.Settings.ChiCoeffPath, chiCoeffFilename);
            
            if obj.Logging.Level > 2
                fprintf('\n Loading Chi coefficients from file:\n\t %s', chiCoeffDataPath);
            end
            
            obj.ChiCoeffs=load(chiCoeffDataPath,'chi_coeff_cell', 'chi_nlm_rst_cell');
            
            obj.ChiCoeffs.Values = obj.ChiCoeffs.chi_coeff_cell;
            obj.ChiCoeffs.Indices = obj.ChiCoeffs.chi_nlm_rst_cell;
            obj.ChiCoeffs.Order = obj.Settings.Order;
                        
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
            
            obj.ZCDdistance = norm(obj.fixed.ZC.Descriptors -  obj.rotating.ZC.Descriptors, 2);
            
        end
               
        function save2pdb(obj,varargin)
            % save2pdb Export strucure to pdb file
            % By default all structures in the object (fixed and rotated)
            % are exported to pdbfiles in the current folder. 
            
            % OPTIONAL INPUT
            %
            %
            %
            %
            %
            %
            %
            
            
            
        end
        
    end
    
    methods (Static)
        
       function Rmat = euler2rotMat(x)
            % euler2rotMat Get the rotation matrix corresponding to the
            % rotation parameterized by the Euler angles x and rotation
            % order (Euler convention) 'zyz'
           
            % lambdas for Euler angles to rotation matrix
            Ry = @(theta) [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
            Rz = @(theta) [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
            
            Rmat = Rz(x(3))*Ry(x(2))*Rz(x(1));
          
        end 
                
    end
    
    
    
end

