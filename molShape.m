classdef molShape < handle
    % molShape Takes a PDB struct (created by pdb2mat.m) and computes the
    % molecular shape function (electron density, van der Waals, solvent-accessible surface and
    % the molecular (solvent excluded) surface). All surfaces are generated
    % in an NxNxN grid with grid side length N; this defines a shape function
    % which can be used as input to compute Zernike-Canterakis moments and
    % descriptors using the ZC class.
    %
    % Filip (Persson) Ljung
    %
    % filip.persson@gmail.com
    %
    % Last modified 2021-02-05
    %
    % TODO: Unit test 
    %
    % -------------------------------------------------------------
    
    properties
        FunctionType % molecular, solvent-accessible surface, vdw or electron density
        GridRes % side length of cubic grid 
        ProbeRadius % probe radius in Ångström
        SmearFactor % for electron density shape, this is the fraction of the grid that we smear out an atom over
        ShellThickness % the thickness in grid units for surfaces
        AtomData % Nx6 list with col X,Y,Z, positional varaince, element code, vdw radii
        ShowLog % flag to turn on/off standard output (logging to console)
        FunctionData % the voxelized (i.e. mapped to a grid) surface (molecular, solvent-accessible surface, vdw or electron density)
        PCAalign % flag to have atoms aligned using PCA: the greatest positional variation will be along the z-axis, then the y-axis, and the least variation along the x-axis
    end
    
    methods
        function obj = molShape(pdbStruct, varargin)   
            % molShape Construct an instance of this class -> compute shape
            % function from PDB data
            % 
            % INPUT
            % <pdbStruct>   :   struct 
            %  PDB.Data = the data property in the PDB class 
            %
            % OPTIONAL
            % 'name':value pairs (defaults)
            % 'FunctionType'           : string ('MS') {'SAS','MS','vdw','electron_density'}
            % 'GridRes'         : integer (64)
            % 'ProbeRadius'     : double (1.4)
            % 'SmearFactor'     : double (0.3)
            % 'ShellThickness'  : integer (2)
            % 'ShowLog'         : true/false (false)
            % 'defaultPCAalign' : true/false (false)
                % If true: atoms coordinates will be rotated so that the greatest positional variation will be along the z-axis, then the y-axis, and the least variation along the x-axis 
            
            expectedShapes = {'SAS','MS','vdw','electron_density'};
            defaultShape = 'MS';
            defaultGridSize = 64;
            defaultProbeRadius = 1.4;
            defaultSmearFactor = 0.3;
            defaultShellThickness = 2;
            defaultShowLog = false;
            defaultPCAalign = false;
            
            p = inputParser;
            
            % Set required
            addRequired(p, 'pdbStruct', @(x)validateattributes(x,{'struct'}, {'nonempty'}, 'pdbStruct'));
            
            % Set optional input
            addOptional(p, 'FunctionType', defaultShape, @(x) any(validatestring(x,expectedShapes)));
            addOptional(p, 'GridRes', defaultGridSize, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','positive'}, 'grid_size'));
            addOptional(p, 'ProbeRadius', defaultProbeRadius, @(x)validateattributes(x,{'numeric'}, {'nonempty','positive'}, 'probe_radius'));
            addOptional(p, 'SmearFactor', defaultSmearFactor,@(x)validateattributes(x,{'numeric'}, {'nonempty','nonnegative','<',1}, 'smear_factor'));
            addOptional(p, 'ShellThickness', defaultShellThickness, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','positive'}, 'shell_thickness'));
            addOptional(p, 'PCAalign', defaultPCAalign)
            addOptional(p, 'ShowLog', defaultShowLog);
            
            parse(p, pdbStruct, varargin{:});
            
            obj.PCAalign = p.Results.PCAalign;
            
            obj.GridRes = p.Results.GridRes;
            obj.FunctionType = p.Results.FunctionType;
            obj.ShowLog = p.Results.ShowLog;
            
            if strcmp(obj.FunctionType,'electron_density')
                obj.SmearFactor = p.Results.SmearFactor;
            else
                obj.SmearFactor = [];
                obj.ProbeRadius = p.Results.ProbeRadius;
                obj.ShellThickness = p.Results.ShellThickness;
            end
            
            % --- process PDB struct ---
            if obj.ShowLog
               startTime = tic;
               fprintf('\n\t Processing PDB data ... ');               
            end
            
            obj.AtomData = molShape.processPDBdata(pdbStruct);
            
            if obj.PCAalign 
                obj.AtomData = molShape.PCAalignAtoms(obj.AtomData);
            end
            
            if obj.ShowLog
               fprintf('Done. Execution time %2.2e s\n', toc(startTime));                    
            end
            
            % --- Compute the shape function ---
             if obj.ShowLog
                 startTime = tic;
                fprintf('\n\t Generating shape function of type "%s" on cubic grid with side length %d\n', obj.FunctionType, obj.GridRes);
                
                if contains(obj.FunctionType, {'MS', 'SAS', 'vdw'})
                    fprintf('\t and shell thickness %d using probe radius %2.2f Å', obj.ShellThickness, obj.ProbeRadius);
                else
                    fprintf('\t and smear factor %2.2f', obj.SmearFactor)
                end
             end
            
            computeShapeFunction(obj);
            
             if obj.ShowLog
               fprintf('\n\t Done. Execution time %2.2e s\n', toc(startTime));                
             end
            
        end
        
        function molSurfFV = getIsoSurface(obj, isovalue)
            % Generates the surface-mesh representation of the volume data
            % molSurfFV contains the faces and vertices of the isosurface
            % and can be passed directly to the PATCH command.
            
            gv = 1:Obj.GridRes;
            [xdata, ydata, zdata] = meshgrid(gv, gv, gv);
            molSurfFV = isosurface(ydata,xdata,zdata,round(obj.solidShape),isovalue);
            
        end
        
        function computeShapeFunction(obj)
                        
            switch obj.FunctionType
                
                case 'MS' % molecular surface
                    [obj.FunctionData, ~] = molShape.MSsurface(obj.AtomData, obj.GridRes, obj.ProbeRadius, obj.ShellThickness);
                    
                case 'SAS' % solvent accessible surface
                    [obj.FunctionData, ~] = molShape.SAsurface(obj.AtomData, obj.GridRes, obj.ProbeRadius, obj.ShellThickness);
                    
                case 'vdw' % van der Waals surface
                    obj.ProbeRadius = 0;
                    [obj.FunctionData, ~] = molShape.SAsurface(obj.AtomData, obj.GridRes, obj.ProbeRadius, obj.ShellThickness);                    
                    
                otherwise % electron density
                    [obj.FunctionData, ~, ~] = molShape.electronDensity(obj.AtomData, obj.GridRes, obj.SmearFactor);
                    
                    % alternatively use molShape.electronDensity_Morris which adds positional variance from B-factors 

            end      
            
        end
        
    end
    
    methods (Static)
        
        function atomList = PCAalignAtoms(atomList)
           
            % cooridnates are already translated so that their centroid is
            % placed at origo.
            XYZ = [atomList(:,1) atomList(:,2) atomList(:,3)];
                        
            % compute eigenvectors from covariance matrix 
            [V, D] = eig(cov(XYZ));
                        
            diagonalEigenvalues = diag(D);
                        
            % sort the eigenvectors based on size of eigenvalues
            [~, I] = sort(diagonalEigenvalues,'descend');
            V = V(:, I);
            
            
            % calculate the angles of the normal vector
            [alpha, beta] = unitVector2Angle(V(:,1));
            
            % align coordinates along:
            
            %  1) z direction
            [~, Ry, Rz] = rotMat(-alpha, pi-beta);
            XYZ = rotateAtoms(XYZ, Ry, Rz);
                     
            % 2) y-direction
            % calculate the angles of the normal vector
            [alpha, ~] = unitVector2Angle(V(:,2));
                     
            [~, Ry, Rz] = rotMat(pi/2 - alpha, 0);
            
            % rotate atoms
            XYZ = rotateAtoms(XYZ, Ry, Rz);
            
            atomList(:,1) = XYZ(:,1);
            atomList(:,2) = XYZ(:,2);
            atomList(:,3) = XYZ(:,3);
            
            
            function [alpha, beta] = unitVector2Angle(u)
                % compute rotational angle between the projected u on the xy plane and the x-axis
                alpha = atan2(u(2), u(1));
                % compute rotational angle between the u vector and the z-axis
                beta = atan2(sqrt(u(1)^2 + u(2)^2), u(3));
                
            end
            
            function [Rx, Ry, Rz] = rotMat(alpha, beta)
                
                Rx = [1 0 0; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)];
                Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
                Rz = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
                
            end
            
            function XYZ = rotateAtoms(XYZ, Ry, Rz)
               
                XYZ = (Ry*Rz*XYZ')';
              
            end
            
        end
        
        function [SAsolid, voxelResolution, scalingFactor] = createSAsolid(atomList, gridRes, probeRadius)
            % CREATE_SA_SOLID
            % Create a binary volumetric representation of the solvent-accessible (SA)
            % solid. Uses the coordinate extrema to find the minimal bounding cube
            % around the object
            
            % Step 1
            % Scale and translate all atoms to fit inside a bounding box with side
            % length L and with integer coordinates (voxels). The resolution of this
            % grid is thus L^3 voxels. After this scaling, the vdW radii r_i and
            % solvent probe radius r_p becomes sr_i and sr_p voxels.
            
            % Step 2
            % Create a solvent-accessible (SA) solid by assigning voxels a value of 1
            % if they are within sr_i+sr_p for each atom, or 0 otherwise.
            
            % INPUT
            % atom_list     :   Nx6 matrix generated from ZC.processPDBdata
                                % col 1-6: X,Y,Z,B-factor,
                      
            % grid_res      :   Side length of the cubic grid (in grid intervals)
            
            % probe_radius  :   Radius of the solvent-probe in Å
            
            % padding       :   padding to add as the fraction of L (0-1)
            
            
            % OUTPUT
            % SA_solid      :   NxNxN binary matrix where N = grid_res
            % Filled voxels = 1; Empty voxels = 0
            
            % TODO: Check input with parser
            %             p = inputParser;
            %
            %             addRequired(p, 'atomList', @(x)validateattributes(x,{'numeric'}, {'nonempty', 'ncols=6'}, 'atomList'));
            %
            %             addOptional(p, 'gridRes', @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','positive'}, 'gridRes'));
            
            
            % --- Step 1
            padding = 0.15;
            n_atoms = size(atomList,1);
            
            % compute centroid ("center of mass" (COM))
            COM = mean(atomList(:,1:3));
            
            % translate with COM at origo
            XYZ_com = atomList(:,1:3) - (ones(n_atoms,1) * COM);
            
            % fprintf('\n Using default scaling');
            % find longest edge of box containing all points
            Rmax = max( (max(XYZ_com)-min(XYZ_com))/2 ) + probeRadius + max(atomList(:,6));
            
            % resolution of voxelgrid in Å
            voxelResolution = Rmax / (gridRes/2) / (1-padding);
            
            scalingFactor = (1-padding) / (Rmax );
            
            
            % scaled probe radius
            sp_r = ceil(probeRadius / voxelResolution);
            sr_i = ceil(atomList(:,6) / voxelResolution);
            
            % scale and translate atom coordinates with center of mass placed at origo
            XYZ_comScaled =  XYZ_com * scalingFactor;
            
            % pre-allocate memory for 3D grid
            SAsolid = false(gridRes, gridRes, gridRes);
            
            % --- Step 2
            keep_on = true;
            
            while keep_on
                
                try
                    atomCrdGrid = round((XYZ_comScaled+1)*(gridRes/2));
                    
                    for i = 1:n_atoms
                        
                        xpos = atomCrdGrid(i,1);
                        ypos = atomCrdGrid(i,2);
                        zpos = atomCrdGrid(i,3);
                        
                        Rip_vox = (sp_r+sr_i(i));
                        Rip_vox2 = (Rip_vox^2);
                        
                        for x_fill = -Rip_vox:Rip_vox
                            for y_fill = -Rip_vox:Rip_vox
                                for z_fill = -Rip_vox:Rip_vox
                                    
                                    d = (x_fill)^2 + (y_fill)^2 + (z_fill)^2;
                                    
                                    if d <= Rip_vox2
                                        if (xpos+x_fill)>0 && (xpos+x_fill)<gridRes && ...
                                                (ypos+y_fill)>0 && (ypos+y_fill)<gridRes && ...
                                                (zpos+z_fill)>0 && (zpos+z_fill)<gridRes
                                            
                                            SAsolid(xpos+x_fill, ypos+y_fill, zpos+z_fill) = true;
                                            
                                        else
                                            error('grid boundary crossing')
                                        end
                                        
                                    end
                                    
                                end % z
                            end % y
                        end % x
                        
                    end
                    
                    keep_on = false;
                    
                catch
                    
                    keep_on = true;
                    %         fprintf('\n Could not fit SA solid on grid - rescaling again');
                    scalingFactor = 0.95 * scalingFactor;
                    
                    % rescale
                    XYZ_comScaled =  XYZ_com * scalingFactor;
                    
                end
                
            end
            
        end
        
        
        function [SAsurf] = buildBoundary(SA_solid)
            %BUILD_BOUNDARY Find boundary (perimeter) of the SA soild.
            % This is the solvent-accessible surface
            
            % Use a general technique that works for any dimensionality
            % and any connectivity.
            
            % INPUT
            % SA_solid  :   NxNxN binary matrix containing the voxelized representation
            % of the solvent-accessible solid (see create_SA_solid.m)
            
            % OUTPUT
            % SA_surf   :   NxNxN binary matrix containing the voxelized representation
            % of the solvent-accessible surface
            
            
            % define connectivity matrix
            conn = conndef(3, 'minimal');
            num_dims = max(3, ndims(conn));
            
            % add padding to SA_solid
            b = padarray(SA_solid,ones(1,3),0,'both');
            
            % erode the 3d image
            b_eroded = imerode(b,conn);
            p = b & ~b_eroded;
            idx = cell(1,3);
            
            for k = 1 : num_dims
                idx{k} = 2:(size(p,k) - 1);
            end
            
            SAsurf = p(idx{:});
            
        end
        
        function [surfaceShell, surfaceSolid] = MSsurface(atomList, gridRes, probeRadius, shellThickness)
            %EDTms Description
            
            [saSolid, voxRes, ~] = molShape.createSAsolid(atomList, gridRes, probeRadius);
            
            % Find boundary -> SAS
            saSurf = molShape.buildBoundary(saSolid);
            
            % Find interior
            M = logical(imfill(saSurf,'holes'));
            
            % Do Euclidean distance transform (EDT) -> Euclidean distance map (EDM)
            EDM = round(bwdist(saSurf,'euclidean'));
            
            % Change sign for interior voxels -> signed EDM (sEDM)
            sEDM=EDM;
            sEDM(M)=-sEDM(M);
            
            % scaled probe radius in voxel units
            sr_p = ceil(probeRadius/voxRes);
            
            % Create isosurface = molecular surface
            mol_surf = (sEDM == (-sr_p) );
            surfaceSolid = (sEDM <= (-sr_p) );
            
            % Create thicker shells
            surfaceShell = mol_surf;
            
            if shellThickness > 1
                for s = 1:shellThickness-1
                    surfaceShell = surfaceShell + (sEDM == (-sr_p - s) );
                end
            end
            
        end
        
        function [surfaceShell, surfaceSolid] = SAsurface(atomList, gridRes, probeRadius, shellThickness)
            
            [surfaceSolid, ~, ~ ] = molShape.createSAsolid(atomList, gridRes, probeRadius);
            
            % Find boundary -> SAS
            saSurf = molShape.buildBoundary(surfaceSolid);
            
            % dilate surface (i.e. thicken the surface)
            se = strel('cube', shellThickness);
            
            surfaceShell = imdilate(saSurf,se);
            
        end
        
        function atomList = processPDBdata(pdbData)
            
            % PARSE_PDB_STRUCT
            % Create list with atom coordinates together
            % with their positional variance (from B-factors) and vdW radii
            
            % INPUT
            % pdbData  :   structured array from pdbread()/pdb2mat (Matlab bioinformatics toolbox )
            
            % OUTPUT
            % atom_list :   Nx6 matrix where N are the number of heavy atoms in the
            % structrure
            % column 1 - X coordinate
            % column 2 - Y coordinate
            % column 3 - Z coordinate
            % column 4 - positional variance (Å)
            % column 5 - element code: 1=H; 2=C; 3=N; 4=O; 5=S
            % column 6 - vdW radii (Å)
            
            % References
            % Bondi, A. (1964). "Van der Waals Volumes and Radii". J. Phys. Chem. 68 (3): 441–451. doi:10.1021/j100785a001
            
            
            
            Xcrds = pdbData.X;
            Ycrds = pdbData.Y;
            Zcrds = pdbData.Z;
            
            atomPositionVariance = (3* [pdbData.betaFactor]) /  (8.0 * pi^2);
            atomSernumber = pdbData.atomNum;
            
            XYZ = [Xcrds Ycrds Zcrds];
            
            COM = mean(XYZ);
            
            n_atoms = size(Xcrds,1);
            
            elements = {'H','C','N','O','S','X'};
            elementCodes = [1 2 3 4 5 6];
            vdwRadii = [1.2 1.7 1.55 1.52 1.8];
            avgRadius = mean(vdwRadii(2:end));
            
            count = 0;
            Hcount = 0;
            
            atomList = zeros(n_atoms,7);
            hydrogenList = zeros(n_atoms,1);
            
            for i = 1:n_atoms
                
                element_i = pdbData.element{i};
                if isempty(element_i) % no element annotation, check resname instead to get element
                    element_i = pdbData.atomName{i}(1);
                    % Take of cases such as 1H..; 2H..
                    if ~isletter(element_i)
                        element_i = pdbData.atomName{i}(2);
                    end
                end
                
                % ONLY HEAVY ATOMS
                try
                    if (element_i=='C') || (element_i=='N')  || (element_i=='O') || (element_i=='S') || (element_i=='X')
                        
                        count = count + 1;
                        %
                        atomList(count,1) = Xcrds(i)-COM(1);
                        atomList(count,2) = Ycrds(i)-COM(2);
                        atomList(count,3) = Zcrds(i)-COM(3);
                        
                        atomList(count,4) = atomPositionVariance(i);
                        atomList(count,7) = atomSernumber(i);
                        
                        element_code = elementCodes(strcmp(element_i, elements));
                        
                        switch element_code
                            
                            case 2 % C
                                atomList(count,6) = vdwRadii(2);                                
                            case 3 % N
                                atomList(count,6) = vdwRadii(3);
                            case 4 % O
                                atomList(count,6) = vdwRadii(4);
                            case 5 % S
                                atomList(count,6) = vdwRadii(5);
                            case 6
                                atomList(count,6) = avgRadius;
                        end
                        
                        atomList(count,5) = element_code;
                        
                    elseif (element_i == 'H') || (element_i == 'D')
                        
                        Hcount = Hcount + 1;
                        hydrogenList(Hcount) = i;
                    end
                    
                catch
                    warning('Could not parse atom name.');
                end
            end
            
            atomList(count+1:end,:) = [];
            %hydrogenList(Hcount+1:end,:) = [];
            %varargout(1) = {hydrogen_list};
            
            
        end
        
        function [density, scalingFactor, voxelRes] = electronDensity_Morris(atomList, gridRes, smearFactor)
            % CREATE_ELECDENSITY Generate shape as an electron density of a molecule using
            % Gaussian-atom representation as implemented in the Python code by 
            % Grandison S., Roberts C., Morris R. J. (2009) in
            % The application of 3D zernike moments for the description of "model-free" molecular
            % structure, functional motion, and structural reliability
            % Journal of Computational Biology 16 (3) 487-500
            % DOI:10.1089/cmb.2008.0083
            
            % The Gaussians are weighted by the B-factor (positional variance)
            
            % INPUT
            % atomList: The Nx6 list of atom coordinates etc created by the
            %           ZC.processPDBdata method in this class.
            % gridRes:  The side length (resolution) of the cubic grid
            % smear_factor: the fraction of the grid that we smear out an atom over
            
            % Step 1
            padding = 0.6 * smearFactor;
            
            nAtoms = size(atomList,1);
            
            COM = mean(atomList(:,1:3));
            
            % translate with COM at origo
            XYZ_com = atomList(:,1:3) -  (ones(nAtoms,1) * COM);

            smearRange = round(smearFactor * gridRes/2);

            % find longest edge of box containing all points
            Rmax = max((max(XYZ_com)-min(XYZ_com))/2) + max(atomList(:,6));
            
            % resolution of voxelgrid in Å
            voxelRes = Rmax / (gridRes/2) / (1-padding);
            scalingFactor = (1-padding) / (Rmax );
            
            % scaled atoms and positional variance (B-factors)
            sr_i = ceil(atomList(:,6) / voxelRes); % scaled vdw
            %sr_B = ceil(atomList(:,4) / voxelRes); % scaled positional var
            
            % scale and translate atom coordinates with center of mass placed at origo
            XYZ_comScaled =  XYZ_com * scalingFactor;
            
            % pre-allocate memory for 3D grid
            density = zeros(gridRes, gridRes, gridRes);
            
            % Step 2
            keep_on = true;
            
            while keep_on
                
                try
                    
                    for i = 1:nAtoms
                        
                        %sigma = sqrt(sr_B(i)) + sr_i(i); % scaled [position variance from B-factor] + scaled [vdw radius]
                        
                        sigma = sr_i(i);
                        
                        xpos = round(0.5 * (XYZ_comScaled(i,1) + 1) * gridRes);
                        ypos = round(0.5 * (XYZ_comScaled(i,2) + 1) * gridRes);
                        zpos = round(0.5 * (XYZ_comScaled(i,3) + 1) * gridRes);
                        
                        for x_smear = -smearRange:smearRange
                            for y_smear = -smearRange:smearRange
                                for z_smear = -smearRange:smearRange
                                    
                                    d = (x_smear)^2 + (y_smear)^2 + (z_smear)^2;
                                    
                                    if d <= smearRange^2
                                        if (xpos+x_smear)>0 && (xpos+x_smear)<gridRes && ...
                                                (ypos+y_smear)>0 && (ypos+y_smear)<gridRes && ...
                                                (zpos+z_smear)>0 && (zpos+z_smear)<gridRes
                                            
                                            val = 8 * exp(-d^2 / (2.0 * (sigma)^2));
                                            
                                            val = val / (sigma * 2 * pi)^0.2; % normalize 
                                            
                                            density(xpos+x_smear, ypos+y_smear, zpos+z_smear) = density(xpos+x_smear, ypos+y_smear, zpos+z_smear) + val;
                                            
                                        else
                                            error('Boundary crossing, rescaling atoms on grid.');
                                        end
                                        
                                    end
                                    
                                end % z
                            end % y
                        end % x
                        
                        keep_on = false;
                        
                    end
                    
                catch % grid overflow, rescale
                    
                    keep_on = true;
                    
                    scalingFactor = 0.95* scalingFactor;
                    
                    XYZ_comScaled =  XYZ_com * scalingFactor;
                    
                end
                
            end
            
        end
        
        function [density, scalingFactor, voxelRes] = electronDensity(atomList, gridRes, smearFactor)
            %CREATE_ELECDENSITY Generate shape as an electron density of a molecule using
            % Gaussian-atom representation as described in 
            
            % J. A. Grant, M. A. Gallardo, and B. T. Pickup, 
            % “A fast method of molecular shape comparison: A simple application 
            % of a Gaussian description of molecular shape,”
            % Journal of Computational Chemistry. 1996.
            
            % INPUT
            % atomList: The Nx6 list of atom coordinates etc created by the
            %           ZC.processPDBdata method in this class.
            % gridRes:  The side length (resolution) of the cubic grid
            % smear_factor: the fraction of the grid that we smear out an atom over
            
            % Step 1
            padding = 0.6 * smearFactor;
            
            weight = 2*sqrt(2); 
            
            nAtoms = size(atomList,1);
            
            COM = mean(atomList(:,1:3));
            
            % translate with COM at origo
            XYZ_com = atomList(:,1:3) -  (ones(nAtoms,1) * COM);

            smearRange = round(smearFactor * gridRes/2);

            % find longest edge of box containing all points
            Rmax = max((max(XYZ_com)-min(XYZ_com))/2) + max(atomList(:,6));
            
            % resolution of voxelgrid in Å
            voxelRes = Rmax / (gridRes/2) / (1-padding);
            scalingFactor = (1-padding) / (Rmax );
            
            % scaled atoms and positional variance (B-factors)
            sr_i = ceil(atomList(:,6) / voxelRes); % scaled vdw
            %sr_B = ceil(atomList(:,4) / voxelRes); % scaled positional var
            
            % scale and translate atom coordinates with center of mass placed at origo
            XYZ_comScaled =  XYZ_com * scalingFactor;
            
            % pre-allocate memory for 3D grid
            density = zeros(gridRes, gridRes, gridRes);
            
            % Step 2
            keep_on = true;
            
            while keep_on
                
                try
                    
                    for i = 1:nAtoms
                        
                        %sigma = sqrt(sr_B(i)) + sr_i(i); % scaled [position variance from B-factor] + scaled [vdw radius]
                        
                        sigma_i = sr_i(i);
                        
                        alpha_i = pi*((3*pi) / (4*pi*sigma_i^3) )^(2/3);
                        
                        
                        xpos = round(0.5 * (XYZ_comScaled(i,1) + 1) * gridRes);
                        ypos = round(0.5 * (XYZ_comScaled(i,2) + 1) * gridRes);
                        zpos = round(0.5 * (XYZ_comScaled(i,3) + 1) * gridRes);
                        
                        for x_smear = -smearRange:smearRange
                            for y_smear = -smearRange:smearRange
                                for z_smear = -smearRange:smearRange
                                    
                                    d = (x_smear)^2 + (y_smear)^2 + (z_smear)^2;
                                    
                                    if d <= smearRange^2
                                        if (xpos+x_smear)>0 && (xpos+x_smear)<gridRes && ...
                                                (ypos+y_smear)>0 && (ypos+y_smear)<gridRes && ...
                                                (zpos+z_smear)>0 && (zpos+z_smear)<gridRes
                                            
                                            val = weight * exp(-d^2 * alpha_i);
                                            
                                            
                                            density(xpos+x_smear, ypos+y_smear, zpos+z_smear) = density(xpos+x_smear, ypos+y_smear, zpos+z_smear) + val;
                                            
                                        else
                                            error('Boundary crossing, rescaling atoms on grid.');
                                        end
                                        
                                    end
                                    
                                end % z
                            end % y
                        end % x
                        
                        keep_on = false;
                        
                    end
                    
                catch % grid overflow, rescale
                    
                    keep_on = true;
                    
                    scalingFactor = 0.95* scalingFactor;
                    
                    XYZ_comScaled =  XYZ_com * scalingFactor;
                    
                end
                
            end
            
        end
        
        
    end % methods (Static)
    
    
end

