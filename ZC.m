classdef ZC < handle
    % ZC Class with methods to compute Zernike-Canterakis moments and
    % associated shape descripors. An object is created, but the class
    % contains static methods that can be used without having an instance of
    % the object. See the ZC (constructor method) for input details.
    %
    % Many methods are based on the C++ library by Novotni & Klein
    % described in:
    % Novotni, M., & Klein, R. (2003).
    % 3D Zernike Descriptors for Content Based Shape Retrieval.
    % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
    % 216–225. https://doi.org/10.1145/781636.781639
    %
    %
    % Filip (Persson) Ljung
    %
    % filip.persson@gmail.com
    %
    % Last modified 2021-02-08
    %
    % IMPORTANT: The bug in the original C++ library that caused the shape
    % descriptors to be cumulative has been fixed in this library.
    %
    % TODO: Unit test
    %
    % -------------------------------------------------------------
    
    properties
        Order       % integer: max ZC expansion order
        Settings    % struct :
        Moments     % struct :
        % Moments.Values = complex-valued vector with ZC moments
        % Moments.CellValues = 3d cell for complex-valued ZC moments arranged  based on their labels,
        %                      i.e. ZCmoments3Dcell{n,l,m} = ZCmom_nlm
        % Moments.IndicesList = Nx5 col 1-3 are the n,l,m
        %                       indices and col 4-5 are the real and imaginary
        %                       components of the moments respectively.
        Descriptors % vector: the real-valued ZC shape descriptor
        ShapeFunction % 3D matrix: the cubic grid with voxel-values defining the shape function
        ChiCoeffs   % : object (ChiCoeffs class)
        %
        %      Values  : 3d cell array
        %      Indices : 3d cell arrayx
        %      Order   : integer
        GridRes     % integer: side length (in grid units) of the cubic grid containing the shape function
        
    end
    
    methods
        
        function obj = ZC(shapeFunction, order, chiCoeffs, varargin)
            % ZC Construct an instance of this class
            %
            % INPUT
            % -------------------------------------------------------------
            % <shapeFunction> : NxNxN double/logical
            %                 (cubic) grid with voxel-values defining the shape function
            % <order> : integer
            %           The max ZC expansion order
            %
            % <ChiCoeffs> : struct with fields (from ZC.computeChiCoeffs in this class)
            %
            % <ChiCoeffs.Values>  : 3d cell array
            % <ChiCoeffs.Indices> : 3d cell arrayx
            % <ChiCoeffs.Order>   : integer
            %
            % OPTIONAL
            % 'name'-value pairs (default)
            %
            % 'scaleOption' : 1<= integer <= 3 (1)
            %                 Option for object scaling during normalization
            %
            % 'showLog' : true/false (false)
            %
            % -------------------------------------------------------------
            
            defaultScaleOption = 1;
            defaultShowLogOp = false; % for debugging
            
            % ----- Option parsing and parameter setup -----
            p = inputParser;
            
            % set required
            addRequired(p, 'shapeFunction', @(x)validateattributes(x,{'numeric', 'logical'}, {'nonempty'}, 'shapeFunction'));
            addRequired(p, 'order', @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','positive'}, 'order'));
            addRequired(p, 'chiCoeffs');
            
            % set optional
            addOptional(p, 'scaleOption', defaultScaleOption, @(x)validateattributes(x,{'numeric'}, {'nonempty','integer','>=1','<=3'}, 'scaleOption'));
            addOptional(p, 'ShowLog', defaultShowLogOp);
            
            parse(p, shapeFunction, order, chiCoeffs, varargin{:});
            
            % ----- Store parameters to object -----
            obj.Order = p.Results.order;
            
            obj.Settings.scaleOption = p.Results.scaleOption;
            
            obj.Settings.ShowLog = p.Results.ShowLog;
            
            obj.ShapeFunction = p.Results.shapeFunction;
            
            obj.ChiCoeffs.Values = p.Results.chiCoeffs.Values;
            obj.ChiCoeffs.Indices = p.Results.chiCoeffs.Indices;
            obj.ChiCoeffs.Order = p.Results.chiCoeffs.Order;
            
            obj.GridRes = size(obj.ShapeFunction, 1);
            
            % Check that we have the right chi coeffcients to compute the
            % ZC moments (at least to the degree that their size/assigned order
            % matches the order).
            if ~isequal(obj.ChiCoeffs.Order, obj.Order)
                error('The order of the Chi coefficients is %d, but needs to be equal to the desired order=%d.', obj.ChiCoeffs.Order, obj.Order);
            end
            
        end
        
        function computeMoments(obj)
            
            if obj.Settings.ShowLog
                fprintf('\n ZC.computeMoments(obj): Perform ZC moment computation');
                startTime = tic;
            end
            
            
            % flatten 3d grid to 1d grid
            if obj.Settings.ShowLog
                fprintf('\n\t Flattening cubic grid to linear grid');
            end
            linearGrid = ZC.cubic2LinearGrid(obj.ShapeFunction, obj.GridRes);
            
            % Compute first order geometric moments
            if obj.Settings.ShowLog
                fprintf('\n\t Computing geometric moments up to first order (the center of gravity)');
            end
            [~, geoMom_COM] = ZC.computeGeoMoments(linearGrid, obj.GridRes , 0, 0, 0, 1, 1);
            
            % 0'th order moments -> normalization
            null_moment = geoMom_COM(1,1,1);
            
            % 1'st order moments -> center of gravity
            xCOM = geoMom_COM(2,1,1) / null_moment;
            yCOM = geoMom_COM(1,2,1) / null_moment;
            zCOM = geoMom_COM(1,1,2) / null_moment;
            
            % Scale to fit inside unit ball
            if obj.Settings.ShowLog
                fprintf('\n\t Computing scale factor to fit inside unit sphere using option=%d', obj.Settings.scaleOption);
            end
            
            [scale, ~] = ZC.computeScaling(obj.Settings.scaleOption, linearGrid, obj.GridRes, xCOM, yCOM, zCOM);
            
            % Normalize: cut off object function for values outide unit
            if obj.Settings.ShowLog
                fprintf('\n\t Normalizing shape function');
            end
            linearGridNorm = ZC.normalizeGrid(linearGrid, obj.GridRes , xCOM, yCOM, zCOM, scale);
            
            % Compute geometric moments
            if obj.Settings.ShowLog
                fprintf('\n\t Computing geometric moments up to order %d', obj.Order);
            end
            
            [~, geoMoments] = ZC.computeGeoMoments(linearGridNorm, obj.GridRes , xCOM, yCOM, zCOM, scale, obj.Order);
            
            % Compute ZC moments
            if obj.Settings.ShowLog
                fprintf('\n\t Computing ZC moments up to order %d\n', obj.Order);
            end
            
            [obj.Moments.IndicesList, obj.Moments.CellValues] = ZC.computeZCmoments(obj.Order, obj.ChiCoeffs.Values, obj.ChiCoeffs.Indices, geoMoments);
            
            obj.Moments.Values = complex(obj.Moments.IndicesList(:,4), obj.Moments.IndicesList(:,5));
            
            if obj.Settings.ShowLog
                stopTime = toc(startTime);
                fprintf('\n Done. Execution time %2.2e s', stopTime);
            end
            
        end
        
        
        function computeDescriptors(obj)
            % COMPUTE_ZCINVARIANTS
            % Computes the Zernike-Canterakis shape descriptors, which are the norms of vectors with
            % components of Z_nl^m with m being the running index.
            
            % The algorithm is based on the C++ code by Novotni & Klein described in
            %
            % % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            %
            % The bug in the original C++ library that caused the shape
            % descriptors to be cumulative has been fixed here.
            %
            % -------------------------------------------------------------
            
            if isempty(obj.Moments)
                computeMoments(obj);
            end
            
            if obj.Settings.ShowLog
                fprintf('\n\t Computing ZC shape descriptors from moments up to order %d\n', obj.Order);
            end
            
            nInvariants = ZC.numberOfInvariants(obj.Order);
            
            obj.Descriptors = zeros(nInvariants,1);
            
            inv_count = 0;
            
            for n = 0:obj.Order
                
                % sum_tmp = 0; % NB This is the bug in the original N&K
                % code that caused the invariants to be cumulatative. Keep
                % this information for legacy reasons.
                
                for l = mod(n,2):2:n
                    
                    sum_tmp = 0;
                    
                    for m=-l:l
                        
                        absM = abs(m);
                        
                        % The ZC_nlm moment
                        mom = obj.Moments.CellValues(n+1, l+1, absM+1);
                        
                        %conjugate if m negative
                        if m<0
                            mom = conj(mom);
                            % take care of sign for odd m
                            if mod(absM,2)
                                mom = -1*mom;
                            end
                        end
                        
                        sum_tmp = sum_tmp + norm(mom)^2;
                        % the C++ std:norm function gives the square of the L2
                        %(euclidian norm), which is the so called field norm
                        
                    end
                    
                    inv_count = inv_count + 1;
                    obj.Descriptors(inv_count,1) =  sqrt(sum_tmp);
                    
                end
            end
            
        end
        
        function computeDescriptors_NKbugVersion(obj)
            % COMPUTE_ZCINVARIANTS
            % Computes the Zernike-Canterakis shape descriptors, which are the norms of vectors with
            % components of Z_nl^m with m being the running index.
            
            % The algorithm is based on the C++ code by Novotni & Klein described in
            %
            % % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            %
            %
            % NB: There is a known bug in this code which
            %
            % -------------------------------------------------------------
            
            if isempty(obj.Moments)
                computeMoments(obj);
            end
            
            if obj.Settings.ShowLog
                fprintf('\n\t Computing ZC shape descriptors from moments up to order %d\n', obj.Order);
            end
            
            nInvariants = ZC.numberOfInvariants(obj.Order);
            
            obj.Descriptors = zeros(nInvariants,1);
            
            inv_count = 0;
            
            for n = 0:obj.Order
                
                sum_tmp = 0; % NB This is a bug that causes the invariants to be cumulatative
                
                for l = mod(n,2):2:n
                    
                    for m=-l:l
                        
                        absM = abs(m);
                        
                        % The ZC_nlm moment
                        mom = obj.Moments.CellValues(n+1, l+1, absM+1);
                        
                        %conjugate if m negative
                        if m<0
                            mom = conj(mom);
                            % take care of sign for odd m
                            if mod(absM,2)
                                mom = -1*mom;
                            end
                        end
                        
                        sum_tmp = sum_tmp + norm(mom)^2;
                        % the C++ std:norm function gives the square of the L2
                        %(euclidian norm), which is the so called field norm
                        
                    end
                    
                    inv_count = inv_count + 1;
                    obj.Descriptors(inv_count,1) =  sqrt(sum_tmp);
                    
                end
            end
            
        end
        
        function setShapeFunction(obj, cubicGrid)
            obj.cubicGrid = cubicGrid;
        end
        
        
    end
    
    
    methods (Static)
        
        
        function [linearGrid] = cubic2LinearGrid(cubicGrid, gridRes)
            % ZC.cubic2LinearGrid
            % Flatten 3D (cubic) voxel-grid to 1D vector.
            % xyz coordinates are mapped to the vector as
            % element_id = dim*((z-1)*dim+(x-1)) + (y-1)+1;
            % where dim is the side length in grid intervals of the cubic grid.
            
            % INPUT
            % -------------------------------------------------------------
            % meshgrid_3D : NxNxN matrix containing the voxelized object
            % N = grid_res
            
            % OUTPUT
            % -------------------------------------------------------------
            % linearGrid       : 1D vector
            
            
            % -------------------------------------------------------------
            
            linearGrid = zeros(gridRes^3,1);
            
            for x = 1:gridRes
                for y = 1:gridRes
                    for z = 1:gridRes
                        
                        if cubicGrid(x,y,z)>0
                            
                            % id = grid_res*((z-1)*grid_res+(y-1)) + (x-1)+1;
                            id = gridRes*((z-1)*gridRes+(x-1)) + (y-1)+1;
                            
                            linearGrid(id) = cubicGrid(x,y,z);
                            
                        end
                        
                    end
                end
            end
            
        end
        
        function [geo_mom_list, geo_moments] = computeGeoMoments(LinearGrid, gridRes , xCOM, yCOM, zCOM, scale, order)
            % ZC.computeGeoMoments
            % Computes geometric moments M_rst up to order n for a voxelized object
            % within a grid. The object, with center of mass at (x|y|z)COG, is scaled to fit within
            % the unit ball before M_rst is computed for each combination of indices,
            % such that r,s,t>0 and r+s+t<n. The algorithm is based on the C++ code by
            % Novotni & Klein described in
            %
            % % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            
            % INPUT
            % -------------------------------------------------------------
            % <LinearGrid> : 1D vector
            % The flatted voxelized object, where x, y and z coordinates are mapped to
            % the array elements as ((z-1) * dim + y-1) * dim + x-1+1
            
            % <dim> : integer
            % The side-length resolution of the cubic grid containing the voxelized
            % object.
            
            % <xCOG>, <yCOG>, <zCOG> : float
            % The center of gravity of the object in the x, y and z dimension
            % respectively
            
            % <scale> : float
            % The scaling factor to fit object into unit ball
            
            % <order> : integer
            % The maximum order for geometric moments.
            
            % OUTPUT
            % -------------------------------------------------------------
            % <geo_mom_list> : [n x 4] matrix, float
            % The geometric moments for each combination of indices r, s, t.
            % Column 1 = index r
            % Column 2 = index s
            % Column 3 = index t
            % Column 4 = geo metric moment
            
            % <geo_moments> : [(order+1) x (order+1) x (order +1)] matrix
            % Each element contains the geometric moments for a given combinatino of
            % indices r,s,t, i.g. M_{r=2,s=2,t=2} = geo_moments(r=2,s=2,t=2)
            
            % -------------------------------------------------------------
            
            % Scale object to fit into unit ball
            xDim = gridRes;
            yDim = gridRes;
            zDim = gridRes;
            
            dimVec = [xDim yDim zDim];
            
            min_v = zeros(3,1);
            
            min_v(1) = (-xCOM+1)*scale;
            min_v(2) = (-yCOM+1)*scale;
            min_v(3) = (-zCOM+1)*scale;
            
            samples = zeros(max(dimVec)+1,3);
            
            for r=1:3
                for s=0:dimVec(r)
                    samples(s+1,r)=min_v(r) + s*scale;
                end
            end
            
            % Compute geometric moments
            arrayDim = zDim;
            layerDim = yDim*zDim;
            
            diffArrayDim = zDim+1;
            diffLayerDim = (yDim+1)*zDim;
            diffGridDim = (xDim+1)*layerDim;
            
            diffGrid = zeros(diffGridDim,1);
            diffLayer = zeros(diffLayerDim,1);
            diffArray = zeros(diffArrayDim,1);
            
            layer = zeros(layerDim,1);
            array = zeros(arrayDim,1);
            
            geo_moments = zeros(order+1, order+1, order+1);
            geo_mom_list = zeros(ZC.numberOfMoments(order),4);
            
            % generate the diff version of the voxel grid in x direction
            iter = 1;
            diffIter = 1;
            
            for x = 0:layerDim-1
                % for x = 0:0
                
                diffGrid = ZC.ComputeDiffFunction(LinearGrid, diffGrid, iter, diffIter, xDim);
                
                iter = iter + xDim;
                diffIter = diffIter + xDim + 1;
                
            end
            
            count = 0;
            
            for r=0:order
                
                % diffGrid
                diffIter = 1;
                
                for p=0:layerDim-1
                    
                    sampleIter = 1;
                    [layer(p+1), diffGrid] = ZC.Multiply(diffGrid, samples(:,1), diffIter, sampleIter, xDim+1);
                    
                    diffIter = diffIter + xDim + 1;
                    
                end
                
                % layer
                iter = 1;
                % diffLayer
                diffIter = 1;
                
                for y=0:arrayDim-1
                    
                    diffLayer = ZC.ComputeDiffFunction(layer, diffLayer, iter, diffIter, yDim);
                    
                    iter = iter + yDim;
                    diffIter = diffIter + yDim + 1;
                    
                end
                
                for s=0:order-r
                    
                    % diffLayer
                    diffIter = 1;
                    
                    for p=0:arrayDim-1
                        
                        sampleIter = 1;
                        [array(p+1), diffLayer] = ZC.Multiply(diffLayer, samples(:,2), diffIter, sampleIter, yDim+1);
                        diffIter = diffIter + yDim + 1;
                    end
                    
                    % array
                    iter = 1;
                    % diffarray
                    diffIter = 1;
                    
                    diffArray = ZC.ComputeDiffFunction(array, diffArray, iter, diffIter, zDim);
                    
                    for t=0:order-r-s
                        
                        count = count+1;
                        sampleIter = 1;
                        [moment, diffArray] = ZC.Multiply(diffArray, samples(:,3), diffIter, sampleIter, zDim+1);
                        geo_moments(r+1,s+1,t+1) = moment / ( (1+r)*(1+s)*(1+t) );
                        geo_mom_list(count,:) = [r s t geo_moments(r+1,s+1,t+1)];
                        
                    end
                    
                end % j
            end %i
            
            geo_mom_list(count+1:end,:) = [];
            
        end
        
        function [B] = ComputeDiffFunction(A, B, Aiter, Biter, dim)
            
            B(Biter) = -1*A(Aiter);
            
            for i = 1:dim-1
                
                B(Biter+i) = A(Aiter+i-1)-A(Aiter+i);
                
            end
            
            B(Biter+dim) = A(Aiter+(dim-1));
            
        end
        
        
        function [sum_val, diffGrid] = Multiply(diffGrid, sample, diffIter, sampleIter, dim)
            
            sum_val = 0;
            
            for i=0:dim-1
                
                diffGrid(diffIter+i) = diffGrid(diffIter+i) * sample(sampleIter+i);
                sum_val = sum_val + diffGrid(diffIter+i);
                
            end
            
        end
        
        function [scaling_factor, max_extent] = computeScaling(scale_option, voxels, dim,  xCOG, yCOG, zCOG)
            % ZC.computeScaling
            % To compute Zernike-Canterakis moments for a 3d object, the object has to
            % be scaled so that it can fit inside the unit-ball (where the ZC functions live).
            
            % This function gives that scaling factor by two different methods. Reproduction accuracy
            % falls off close to the ball surface (due to discretization effects?) and both methods make sure
            % the scaling avoids this. Novotni and Klein used method 2 as per their code (not described in paper)
            % without any justification. Method 1 is used by Grandison, Roberts and Morris as per the paper,
            % but no justification for why the factor of 0.6 was chosen.
            
            
            % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            
            % Grandison, S., Roberts, C., & Morris, R. J. (2009).
            % The Application of 3D Zernike Moments for the Description of
            % “Model-Free” Molecular Structure, Functional Motion, and Structural Reliability.
            % Journal of Computational Biology, 16(3), 487–500. https://doi.org/10.1089/cmb.2008.0083
            
            
            % INPUT
            % -------------------------------------------------------------
            % <scale_option> : 1,2,3
            % 1  =  scaling so that the maximum distance between a filled voxel
            %       and the grid centroid is 70 % of the unit ball radius.
            
            % 2  =  (Used by Novotni & Klein) =  scaling by 2 x Radius of gyration
            %       (average distance between filled voxels and the grid
            %       centroid).
            
            % 3  =  scaling so that the maximum distance between a
            %       filled voxel and the grid centroid is 100 % of the unit ball radius.
            
            % <voxels> : 1D vector
            % The flatted voxelized object, where x, y and z coordinates are mapped to
            % the array elements as ((z-1) * dim + y-1) * dim + x-1+1
            
            % <dim> : integer
            % The side-length resolution of the cubic grid containing the voxelized
            % object.
            
            % <xCOG>, <yCOG>, <zCOG> : float
            % The center of gravity of the object in the x, y and z dimension
            % respectively
            
            % OUTPUT
            % -------------------------------------------------------------
            % <scaling_factor> : float
            % The scaling factor for fitting object voxels into the unit ball.
            % <Rmax> : float
            % Maximum distance from grid centroid
            % <max_extent> : float
            % For scaling option 1: The fraction of the unit ball radius where the largest distance from the
            % centroid is
            % <voxel_resolution> : float
            % The resolution (in Ångström) captured by a voxel
            
            % -------------------------------------------------------------
            
            switch scale_option
                
                case 1
                    Rmax = 0;
                    for x = 1:dim
                        for y = 1:dim
                            for z = 1:dim
                                
                                
                                id = ((z-1) * dim + y-1) * dim + x-1+1;
                                
                                if voxels(id) > 0
                                    
                                    mx = x - xCOG;
                                    my = y - yCOG;
                                    mz = z - zCOG;
                                    
                                    temp = mx^2+my^2+mz^2;
                                    
                                    if temp > Rmax
                                        Rmax = temp;
                                    end
                                    
                                end
                                
                            end
                        end
                    end
                    
                    max_extent = 0.7;
                    
                    scaling_factor = max_extent*(1 / sqrt(Rmax));
                    %         Rmax = sqrt(Rmax);
                    
                    
                case 2
                    temp=0;
                    nvox = 0;
                    for x = 1:dim
                        for y = 1:dim
                            for z = 1:dim
                                
                                id = ((z-1)* dim + y-1) * dim + x-1+1;
                                
                                if  voxels(id) > 0
                                    
                                    mx = x - xCOG;
                                    my = y - yCOG;
                                    mz = z - zCOG;
                                    
                                    temp = temp + mx^2+my^2+mz^2;
                                    nvox=nvox+1;
                                    
                                end
                            end
                        end
                    end
                    
                    Rmean = sqrt(temp/nvox);
                    
                    max_extent = 0.5;
                    scaling_factor = 1/(2* Rmean);
                    
                case 3
                    
                    Rmax = 0;
                    for x = 1:dim
                        for y = 1:dim
                            for z = 1:dim
                                
                                
                                id = ((z-1) * dim + y-1) * dim + x-1+1;
                                
                                if voxels(id) > 0
                                    
                                    mx = x - xCOG;
                                    my = y - yCOG;
                                    mz = z - zCOG;
                                    
                                    temp = mx^2+my^2+mz^2;
                                    
                                    if temp > Rmax
                                        Rmax = temp;
                                    end
                                    
                                end
                                
                            end
                        end
                    end
                    
                    max_extent = 1;
                    scaling_factor = max_extent*(1 / sqrt(max(Rmax)));
                    
            end
            
            % voxel_resolution = Rmax / ( (dim/2));
            
        end
        
        function [zernikeMoments_list, zernikeMoments] = computeZCmoments(order, chi_coeff_cell, chi_nlm_rst_cell, geo_moments)
            % ZC.computeZCmom
            % Computes the Zernike moments of order n from geometric moments
            % and the object-independent chi-coefficients of the same order.
            % This computation is data dependent and has to be performed for
            % each new object and/or transformation.
            
            % The algorithm is based on the C++ code by Novotni % Klein,
            % described in
            % % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            
            % INPUT
            % -------------------------------------------------------------
            % <order> : integer
            % The maximum expansion order
            
            % <chi_coeff_cell> : cell array
            % The chi coefficients from ZC.computeChiCoeffs (precomputed and
            % loaded from mat files typically (can be expensive to compute))
            
            % <chi_nlm_rst_cell> : cell array
            % The indices to label the chi coefficients from ZC.computeChiCoeffs
            
            % OUTPUT
            % -------------------------------------------------------------
            % <zernikeMoments_list> : NxM matrix
            
            % <zernikeMoments> :
            
            % -------------------------------------------------------------
            
            
            zernikeMoments = zeros(order+1, order+1, order+1);
            
            zernikeMoments_list = zeros(ZC.numberOfMoments(order),5);
            momcount = 0;
            
            pi_factor = (3/(4*pi));
            
            for n=0:order
                
                l0 = mod(n,2);
                
                for l = l0:2:n
                    
                    for m=-l:l
                        
                        M=abs(m);
                        
                        zm = complex(0);
                        
                        % get chi coeffs
                        chi_nlm_rst = chi_nlm_rst_cell{n+1,l+1,M+1};
                        chi_values = chi_coeff_cell{n+1,l+1,M+1};
                        nCoeffs = size(chi_nlm_rst,1);
                        
                        for i = 1:nCoeffs
                            
                            r = chi_nlm_rst(i,4)+1;
                            s = chi_nlm_rst(i,5)+1;
                            t = chi_nlm_rst(i,6)+1;
                            
                            %                 geoMoms_test(count2,:) = [n l m r-1 s-1 t-1 real(chi_values(i)) imag(chi_values(i)) geo_moments(r,s,t)];
                            
                            zm = zm + conj(chi_values(i)) * geo_moments(r,s,t);
                            
                        end
                        
                        zm = zm * pi_factor;
                        
                        if n == 0 && l == 0 && m == 0
                            nullMoment = real(zm);
                        end
                        
                        if m<0
                            
                            zernikeMoments(n+1,l+1,m+l+1) = (-1)^M * conj(zm);
                            
                        else
                            
                            zernikeMoments(n+1,l+1,m+l+1) = zm;
                            
                        end
                        
                        
                        momcount = momcount + 1;
                        
                        zernikeMoments_list(momcount,1) = n;
                        zernikeMoments_list(momcount,2) = l;
                        zernikeMoments_list(momcount,3) = m;
                        
                        zernikeMoments_list(momcount,4) = real(zernikeMoments(n+1,l+1,m+l+1));
                        zernikeMoments_list(momcount,5) = imag(zernikeMoments(n+1,l+1,m+l+1));
                        
                    end % m
                end % l
            end % n
            
            zernikeMoments_list(momcount+1:end,:) = [];
            
        end
        
        function [ChiCoeff, chiCoeffList] = computeChiCoeffs(order)
            % ZC.computeChiCoeffs
            % Computes the coefficients associated with geometric moments used in the computation
            % of Zernike-Canterakis (ZC) moments. These coefficients do not depend on the object
            % and should be pre-computed and stored for computational efficiency
            
            % The algorithm is based on the C++ code by Novotni % Klein described in
            %
            % % Novotni, M., & Klein, R. (2003).
            % 3D Zernike Descriptors for Content Based Shape Retrieval.
            % Proceedings of the Eighth ACM Symposium on Solid Modeling and Applications,
            % 216–225. https://doi.org/10.1145/781636.781639
            
            % INPUT
            % -------------------------------------------------------------
            % <order>    :   integer
            % order for the ZC moments
            
            % OUTPUT
            % -------------------------------------------------------------
            % ChiCoeff struct with field
            % Values : 3d cell array with dim (order+1)x(order+1)x(order+1)
            % Each cell contains coefficients for all [r,s,t] (indices) triplets associated with an [n,l,m] triplet
            % These are obtained using the syntax ChiCoeff.Values{n,l,m}
            
            % Indices  : 3d cell array with dim (order+1)x(order+1)x(order+1)
            % Each cell contains the [n,l,m] indices (col 1-3)
            % and the combinations of [r,s,t] indices associated with that [n,l,m] triplet.
            % These are obtained using the syntax ChiCoeff.Indices{n,l,m}
            
            % chi_coeff         : Nx8 matrix for the N combinations of [n,l,m] (col 1-3)
            % and [r,s,t] (col 4-6) triplets.
            % col 7 and 8 contains the real and imaginary values
            % associated with the triplet combinations.
            %
            % -------------------------------------------------------------
            
            
            count_coeff = 0;
            setStart = 1;
            set_count = 0;
            
            chi_coeff = zeros((order+1)^4,8);
            
            chi_coeff_cell = cell(order+1, order+1, order+1);
            chi_nlm_rst_cell = cell(order+1, order+1, order+1);
            
            chi_ind_map = zeros((order+1)^3,6);
            
            fprintf('\n Computing chi coefficients for order = %d', order);
            
            n_count_vec = zeros(order+1,1);
            
            for n = 0:order
                
                fprintf('\n Doing order %d', n);
                
                li=0;
                l0 = mod(n,2);
                
                n_count = 0;
                
                for l = l0:2:n % only even values of l since Zernike functions require that (n-l) is even
                    li=li+1;
                    
                    for m = 0:l
                        
                        c_set_count=0;
                        
                        cs = c_lm(l,m) ;
                        
                        w = cs/ 2^(m);
                        
                        k=round( (n-l)/2 );
                        
                        for nu = 0:k
                            
                            qs = q_klnu(k,l,nu);
                            w_Nu = w * qs;
                            
                            for alpha = 0:nu
                                
                                w_NuA = w_Nu * nchoosek(nu,alpha);
                                
                                for beta = 0:(nu-alpha)
                                    
                                    w_NuAB = w_NuA * nchoosek(nu-alpha, beta);
                                    
                                    for p = 0:m
                                        
                                        w_NuABP = w_NuAB * nchoosek(m,p);
                                        
                                        for mu = 0:floor((l-m)/2)
                                            
                                            w_NuABPMu = w_NuABP ...
                                                * nchoosek(l, mu) ...
                                                * nchoosek(l-mu, m+mu) ...
                                                / 2^(2 * mu);
                                            
                                            for q = 0:mu
                                                
                                                w_NuABPMuQ = w_NuABPMu * nchoosek(mu, q);
                                                
                                                % the sign
                                                if mod((m-p+mu),2)
                                                    w_NuABPMuQ = -1 * w_NuABPMuQ;
                                                end
                                                
                                                rest = mod(p,4);
                                                
                                                switch rest
                                                    case 0
                                                        c = complex(w_NuABPMuQ, 0);
                                                    case 1
                                                        c = complex(0, w_NuABPMuQ);
                                                    case 2
                                                        c = complex(-1 * w_NuABPMuQ, 0);
                                                    case 3
                                                        c = complex(0, -1 * w_NuABPMuQ);
                                                end
                                                
                                                t_i = l - m + 2 * (nu - alpha - beta - mu);
                                                s_i = 2 * (mu - q + beta) + m - p;
                                                r_i = 2 * q + p + 2 * alpha;
                                                
                                                c_set_count = c_set_count +1;
                                                
                                                count_coeff = count_coeff + 1;
                                                n_count = n_count + 1;
                                                
                                                chi_coeff(count_coeff,1) = n;
                                                chi_coeff(count_coeff,2) = l;
                                                chi_coeff(count_coeff,3) = m;
                                                
                                                chi_coeff(count_coeff,4) = r_i;
                                                chi_coeff(count_coeff,5) = s_i;
                                                chi_coeff(count_coeff,6) = t_i;
                                                
                                                chi_coeff(count_coeff,7) = real(c);
                                                chi_coeff(count_coeff,8) = imag(c);
                                                
                                            end %q
                                        end % mu
                                    end % p
                                end % beta
                            end % alpha
                        end % nu
                        
                        if c_set_count > 0
                            
                            set_count = set_count + 1;
                            
                            chi_ind_map(set_count,1) = n;
                            chi_ind_map(set_count,2) = l;
                            chi_ind_map(set_count,3) = m;
                            
                            chi_ind_map(set_count, 4) = setStart;
                            chi_ind_map(set_count, 5) = setStart + c_set_count - 1;
                            
                            setStart = chi_ind_map(set_count, 5) + 1;
                        end
                        
                        sel_int = chi_ind_map(set_count, 4): chi_ind_map(set_count, 5);
                        
                        chi_coeff_cell{n+1,l+1,m+1} = complex( chi_coeff(sel_int, 7), chi_coeff(sel_int, 8) );
                        chi_nlm_rst_cell{n+1,l+1,m+1} = chi_coeff(sel_int, 1:6);
                        
                    end % m
                end % l
                
                n_count_vec(n+1,1) = n_count;
                
                
            end % n
            
            chi_coeff(count_coeff+1:end,:) = [];
            chi_ind_map(set_count+1:end,:) = [];
            
            chi_ind_map(:, 6) = chi_ind_map(:, 5) - chi_ind_map(:, 4)+1;
            
            ChiCoeff.Values = chi_coeff_cell;
            ChiCoeff.Indices = chi_nlm_rst_cell;
            
            chiCoeffList = chi_coeff;
            
            
            
            % Help functions
            function c_lm_val = c_lm(l,m)
                % c_l^m = c_l^-m
                
                c_lm_val = sqrt( (2 * l + 1) * factorial(l + m) * factorial(l - m) ) / factorial(l);
                
                
            end
            
            function q_klnu_val = q_klnu(k,l,nu)
                
                % nominator of straight part
                nom = nchoosek(2*k,k) * (nchoosek(k, nu)) * (nchoosek(2 * (k + l + nu) + 1, 2 * k));
                
                if mod(k+nu,2)
                    nom = -1*nom;
                end
                
                % denominator of straight part
                den = 2^(2*k) * nchoosek(k+l+nu,k);
                
                % nominator of sqrt part
                n_sqrt = 2*l + 4*k + 3;
                
                % denominator of sqrt part
                d_sqrt = 3;
                
                q_klnu_val =  nom / den * sqrt(n_sqrt / d_sqrt);
                
            end
            
        end
        
        function [linearGrid] = normalizeGrid(linearGrid, dim, xCOG, yCOG, zCOG, scale)
            %NORMALIZE_GRID Cuts of the object function for values outside the unit
            %ball
            
            radius = 1/scale;
            sqrRadius = radius^2;
            
            for x = 1:dim
                for y = 1:dim
                    for z = 1:dim
                        
                        
                        id = ((z-1) * dim + y-1) * dim + x-1+1;
                        
                        if linearGrid(id) > 0
                            
                            mx = x - xCOG;
                            my = y - yCOG;
                            mz = z - zCOG;
                            
                            sqrLen = mx^2+my^2+mz^2;
                            
                            if sqrLen > sqrRadius
                                linearGrid(id) = 0;
                            end
                            
                        end
                        
                    end
                end
            end
            
        end
        
        function [nInvariants] = numberOfInvariants(order)
            
            nInvariants = floor( ((order+2)/2)^2 );
            
        end
        
        function [nMoments] = numberOfMoments(order)
            
            nMoments = (order+1)*(order+2)*(order+3)/6;
            
        end
        
        function [order] = expansionOrder(nMoments)
            
            p = [1 6 11 (6-nMoments*6)];
            
            solutions = roots(p);
            
            for n = 1:3               
                if isreal(solutions(n))                    
                    order = solutions(n);                    
                end                
            end
            
        end
        
        function [nlmList] = NLMlabels(order)
            
            nlmList = zeros(int64(ZC.numberOfMoments(order)),3);
            
            momCount = 0;
            
            for n = 0:order
                
                for l = mod(n,2):2:n
                    
                    for m=-l:l
                        
                        momCount = momCount + 1;
                        
                        nlmList(momCount,:) = [n l m];
                        
                    end
                    
                end
                
            end
            
        end
        
        function shapeRecon = reconstructShape(gridRes, ZCm, ZPgrid)
            % Reconstruct object from ZC moments and pre-computed
            % ZC-functions for the specified grid size
            
            % -------------------------------------------------------------
            
            fprintf('\n Reconstructing shape function from moments... ');
            shapeRecon = zeros(gridRes,gridRes,gridRes);
            %grid_recon2 = grid_recon;
            
            %             nmoms = size(ZCm,1);
            for x = 1:gridRes
                fprintf('\n Doing layer %d/%d', x, gridRes);
                for y = 1:gridRes
                    for z = 1:gridRes
                        
                        %                         zp_tmp = reshape(ZPgrid(x,y,z,:), [nmoms 1 1 1]);
                        zp_tmp = squeeze(ZPgrid(x,y,z,:));
                        shapeRecon(x,y,z) = real(sum(zp_tmp.*ZCm));
                        
                    end
                end
            end
            
            %grid_recon2 = permute(grid_recon,[2 1 3]);
            fprintf('done. \n ');
        end
        
        function mat2DX(DXdata, varargin)
            % this function outputs a DX file for viewing in VMD or PyMOL. It uses
            % six inputs to function properly (subvariables of "DXdata"). While only
            % one is necessary (the 3D matrix called "densityMatrix"), the others shape
            % and position the 3D matrix in space. A minimum point and a voxel length
            % are required to do this. The minimum point is recorded as the "minX"
            % "minY" and "minZ" variables.
            % As of July 2012, file output format is as defined on
            % http://www.poissonboltzmann.org/file-formats/mesh-and-data-formats/opendx-scalar-data
            %
            %
            % -- required inputs (1) --
            %
            % input value           meaning
            %
            % input.densityMatrix  rectangular-prism 3D matrix holding density values.
            %                       All XYZ values correspond to spatial location of
            %                       data.
            %
            % -- optional inputs (5): generates defaults when not user-specified --
            %
            % input value        meaning                           default value
            %
            % input.outfile      output file name                 "mat2dx.dx"
            % input.minX         origin (min) X coord (angstroms)  0.00
            % input.minY         origin (min) Y coord (angstroms)  0.00
            % input.minZ         origin (min) Z coord (angstroms)  0.00
            % input.voxelLength  length to side of voxel (angst.)  1.00
            %
            %
            % -- example usage: output a spherical gradient 10 angstroms in radius --
            %
            % [X Y Z] = meshgrid(-10:10);                        % make a 3D gradient
            % input.densityMatrix = ((X.^2 + Y.^2 + Z.^2).^0.5); % construct 3D space
            % isosurface(input.densityMatrix);                   % plot in matlab
            % mat2dx(input);                                     % output as DX file
            %
            % -- example usage: normalize DX output --
            %
            % DXdata = dx2mat('DXfile.dx');                      % interpret DX data
            % mat3D = DXdata.densityMatrix;                      % make name convenient
            % maxValue = max(max(max(mat3D)));                   % find max value
            % minValue = min(min(min(mat3D)));                   % find min value
            % valueRange = maxValue - minValue;                  % find range of data
            % DXdata.densityMatrix = (mat3D-minValue)/valueRange; % normalize and save
            % DXdata.outfile = 'DXfile_normaized.dx';            % rename output
            % mat2dx(DXdata);                                    % write to DX file
            
            p = inputParser;
            
            default_NormalizeOp = true;
            default_FilePath = fullfile(pwd,'ZC_shapereconstruxtion.dx');            % rename output
            
            addOptional(p, 'Normalize', default_NormalizeOp);
            addOptional(p, 'FilePath', default_FilePath);
            
            
            parse(p, varargin{:});
            
            normalizeOp = p.Results.Normalize;
            filePath = p.Results.FilePath;
            
            if normalizeOp
                
                mat3D = DXdata.densityMatrix;                      % make name convenient
                maxValue = max(max(max(mat3D)));                   % find max value
                minValue = min(min(min(mat3D)));                   % find min value
                valueRange = maxValue - minValue;                  % find range of data
                DXdata.densityMatrix = (mat3D-minValue)/valueRange; % normalize and save
                DXdata.outfile = filePath;
                
            end
            
            % review input data
            
            % check required input
            if ~isfield(DXdata, 'densityMatrix')
                fprintf(['this profram requires a 3D matrix "densityMatrix"\n' ...
                    'to function properly\n\texiting...']);
                return;
            end
            
            % check optional inputs
            if ~isfield(DXdata, 'minX')
                fprintf('no minX variable given, default is 0.00\n');
                DXdata.minX = 0;
            end
            if ~isfield(DXdata, 'minY')
                fprintf('no minY variable given, default is 0.00\n');
                DXdata.minY = 0;
            end
            if ~isfield(DXdata, 'minZ')
                fprintf('no minZ variable given, default is 0.00\n');
                DXdata.minZ = 0;
            end
            if ~isfield(DXdata, 'outfile')
                fprintf('no output file name given, default is mat2dx.dx\n');
                DXdata.outfile = 'mat2dx.dx';
            end
            if ~isfield(DXdata, 'voxelLength')
                fprintf('no voxel length given, default is 1.00 angstroms\n');
                DXdata.voxelLength = 1;
            end
            
            fprintf('outputting 3D data to %s\n', DXdata.outfile);
            
            % prepare variables for output function
            
            % relabel variables
            outfile = DXdata.outfile;
            minX = DXdata.minX;
            minY = DXdata.minY;
            minZ = DXdata.minZ;
            voxelLength = DXdata.voxelLength;
            densityMatrix = DXdata.densityMatrix;
            
            % convenience variables
            dimen = size(densityMatrix);
            totalElements = prod(dimen);
            overFlowVals = mod(totalElements,3);
            numRows = floor(totalElements / 3);
            
            % reshape data for fast output
            densityMatrix = reshape(permute(densityMatrix, [3 2 1]), [1 totalElements]);
            out3DMat = reshape(densityMatrix(1:end - overFlowVals), [3 numRows])';
            
            % output data into file
            
            % begin file creation
            FILE = fopen(outfile, 'w');
            
            % output header information
            fprintf( FILE, 'object 1 class gridpositions counts%8.0f%8.0f%8.0f\n', ...
                dimen(1), dimen(2), dimen(3));
            fprintf( FILE, 'origin%16g%16g%16g\n', minX, minY, minZ);
            fprintf( FILE, 'delta%16g 0 0\n', voxelLength);
            fprintf( FILE, 'delta 0 %16g 0\n', voxelLength);
            fprintf( FILE, 'delta 0 0 %16g\n', voxelLength);
            fprintf( FILE, 'object 2 class gridconnections counts%8.0f%8.0f%8.0f\n', ...
                dimen(1), dimen(2), dimen(3));
            fprintf( FILE, 'object 3 class array type double rank 0 items  %g data follows\n', ...
                totalElements);
            
            % output density information
            newLine = 0;
            for n = 1:numRows
                fprintf(FILE,'%16E%16E%16E\n',out3DMat(n,:));
                % output status
                if ~mod(n, floor(numRows/20))
                    fprintf('   %6.0f%%',n/numRows*100);
                    newLine = newLine + 1;
                    if ~mod(newLine,5)
                        fprintf('\n');
                    end
                end
            end
            if overFlowVals > 0 % values not in complete row
                for n = 1:overFlowVals
                    fprintf(FILE,'%16E',densityMatrix(end-n+1));
                end
                fprintf(FILE,'\n');
            end
            fprintf('\n');
            
            % output file tail
            fprintf( FILE, ['attribute "dep" string "positions"\n' ...
                'object "regular positions regular connections" class field\n' ...
                'component "positions" value 1\n' ...
                'component "connections" value 2\n' ...
                'component "data" value 3\n']);
            
            %fprintf( FILE, 'object "Untitled" call field\n'); % alternate tail
            
            fclose(FILE);
            
        end
        
        function ZCfunGrid = computeZCfunctionsGrid(gridRes, order, scaleOption, chi_coeff_cell, chi_nlm_rst_cell, parOp, poolSize)
            % To speed up reconstructions we can precompute and save the Zernike
            % functions for each voxel in the grid. This scripts will compute all the
            %  N_ZCmoments 3D double precision complex arrays which can be saved for
            %  later use.
            
            % Default scaling is to consider object-voxels to be scaled to fit within 70% of unit ball
            % radius.
            
            
            % INPUT
            % -------------------------------------------------------------
            % <gridRes> : side length of the cubic grid
            
            % <order> : integer
            % The maximum expansion order
            
            % <chi_coeff_cell> : cell array
            % The chi coefficients from ZC.computeChiCoeffs (precomputed and
            % loaded from mat files typically (can be expensive to compute))
            
            % <chi_nlm_rst_cell> : cell array
            % The indices to label the chi coefficients from ZC.computeChiCoeffs
            
            % <scale_option> : 1,2,3
            % 1  =  scaling so that the maximum distance between a filled voxel
            %       and the grid centroid is 70 % of the unit ball radius.
            
            % 2  =  (Used by Novotni & Klein) =  scaling by 2 x Radius of gyration
            %       (average distance between filled voxels and the grid
            %       centroid).
            
            % 3  =  scaling so that the maximum distance between a
            %       filled voxel and the grid centroid is 100 % of the unit ball radius.
            
            % parOp : true/false
            % Use parallell computation
            
            % poolSize : number of workers to use for parallell computation
            
            
            % OUTPUT
            % -------------------------------------------------------------
            % <ZCfunGrid> : NxNxNxM complex array
            % contains all ZC-function values up to order M for each voxel in
            % the NxNxN cubic grid
            
            
            % --------------------------------------------------------------------------
            
            nmoms = round((1/6)*(1+order)*(2+order)*(3+order));
            index_list = zeros(nmoms,3);
            ZCfunGrid = complex(zeros(gridRes, gridRes, gridRes, nmoms) );
            
            [dimX, dimY, dimZ, ~] = size(ZCfunGrid);
            
            xdimvec = (1:dimX);
            xCOG = mean(xdimvec);
            
            ydimvec = (1:dimY);
            yCOG = mean(ydimvec);
            
            zdimvec = (1:dimZ);
            zCOG = mean(zdimvec);
            
            switch scaleOption
                case 1
                    scale = 1/(dimX/2) * 0.7;
                case 2
                    gv = (-dimX/2:dimX/2)/(dimX/2);
                    [xdata, ydata, zdata] = meshgrid(gv, gv, gv);
                    
                    d = sqrt(xdata.^2 + ydata.^2 + zdata.^2);
                    
                    scale = 1/ (2*mean(d(d<=1)));
                    
                case 3
                    scale = 1/(dimX/2);
            end
            
            if parOp
                
                p = gcp('nocreate');
                
                if isempty(p)
                    fprintf('\nSetting up a parpool with %d workers\n', poolSize);
                    parpool(poolSize)
                else
                    poolSize = p.NumWorkers;
                    fprintf('\nParpool running with %d workers\n', poolSize);
                end
            end
            
            fprintf('\n\n\t PRE-PROCESSING DATA');
            fprintf('\n---------------------------------------');
            
            % only reconstruct voxels within unit ball
            [x_ptsMesh, y_ptsMesh, z_ptsMesh] = meshgrid(1:dimX, 1:dimY, 1:dimZ);
            
            x_ptsMesh = (x_ptsMesh - xCOG) * scale;
            y_ptsMesh = (y_ptsMesh - yCOG) * scale;
            z_ptsMesh = (z_ptsMesh - zCOG) * scale;
            
            dists = sqrt( x_ptsMesh.^2 + y_ptsMesh.^2 + z_ptsMesh.^2);
            
            withinBall = dists<1;
            withinBall_tot = sum(sum(sum(withinBall)));
            
            fprintf('\nNumber of Zernike polynomials to compute within unit ball: %5.0f', sum( sum( sum(withinBall) ) ) );
            
            % Pre-compute monomials
            fprintf('\nComputing table of all monomials (x^r * y^s * z^t; r+s+t <= order)');
            
            x_ptsVec = (xdimvec - xCOG) * scale;
            y_ptsVec = (ydimvec - yCOG) * scale;
            z_ptsVec = (zdimvec - zCOG) * scale;
            
            order_vec = 0:order;
            
            monomials_x = power( (ones(order + 1, 1) * x_ptsVec), (ones(dimX, 1) * order_vec)' );
            monomials_y = power( (ones(order + 1, 1) * y_ptsVec), (ones(dimY, 1) * order_vec)' );
            monomials_z = power( (ones(order + 1, 1) * z_ptsVec), (ones(dimZ, 1) * order_vec)' );
            
            
            fprintf('\n\n\t COMPUTING ZC function for a %d grid and up to order %d ', gridRes, order);
            fprintf('\n---------------------------------------');
            
            fprintf('\n\n%5.0f percent completed:  %5.2f min remaining', 0, NaN);
            nvoxc=0;
            
            for x = 1:dimX
                
                start = tic;
                fprintf('\n Doing layer %3.0f (%3.0f)', x, dimX);
                
                for y = 1:dimY
                    
                    for z = 1:dimZ
                        
                        % only process voxels within unit ball
                        if withinBall(x,y,z)
                            
                            momcount=0;
                            
                            for n=0:order
                                
                                l0 = mod(n,2);
                                for l = l0:2:n
                                    
                                    for m = -l:l
                                        % Evaluate zernike polynomial Zp_{n,l,m} at point [x y z]
                                        zp = complex(0);
                                        absM = abs(m);
                                        
                                        chi_values = chi_coeff_cell{n+1, l+1, absM+1};
                                        chi_nlm_rst = chi_nlm_rst_cell{n+1, l+1, absM+1};
                                        nCoeffs = size(chi_nlm_rst, 1);
                                        
                                        % conjugate if m negative
                                        if m < 0
                                            chi_values = conj(chi_values);
                                            % take care of sign
                                            if mod(absM,2)
                                                chi_values = -1*chi_values;
                                            end
                                        end
                                        
                                        % loop over invariants
                                        for i = 1:nCoeffs
                                            
                                            zp = zp + chi_values(i) ...
                                                * monomials_x(chi_nlm_rst(i, 4)+1, x) ...
                                                * monomials_y(chi_nlm_rst(i, 5)+1, y) ...
                                                * monomials_z(chi_nlm_rst(i, 6)+1, z);
                                            
                                        end %i
                                        
                                        momcount = momcount+1;
                                        index_list(momcount,:) = [n l m ];
                                        
                                        ZCfunGrid(x,y,z,momcount) = zp;
                                        
                                    end %m
                                    
                                end % end l
                                
                            end % n
                            
                        end % if within ball
                        
                    end % z
                end % y
                
                stop = toc(start);
                
                nvoxs_in_layer = sum(sum(withinBall(x,:,:)));
                
                nvoxc = nvoxc + nvoxs_in_layer;
                
                t_pervox = stop/nvoxs_in_layer;
                
                tleft = (withinBall_tot-nvoxc)*t_pervox;
                
                fprintf(repmat('\b',1,69))
                %     fprintf('\n\n%5.0f percent completed:  %5.2f min remaining', x/dimX*100, stop/60*(dimX-x));
                fprintf('\n\n%5.0f percent completed:  %5.2f min remaining', nvoxc/withinBall_tot*100, tleft/60);
                
            end % x
            
            stopend = toc(startbeg);
            
            fprintf('\n Total execution time: %2.2f min', stopend/60);
            
            saveFilename = sprintf('ZP_%dgrid_order%d.mat', grid_res, order);
            
            fprintf('\n\n Saving data to file %s\n', saveFilename);
            
            save(saveFilename, 'ZCfunGrid', 'scaleOption', 'scale', '-v7.3');
            
        end
        
        function Descriptors = moments2descriptors(moments)
            % accepts as input either
            % ZC.Moments.CellValues (which is a 3D array and not a cell despite its name)
            % or
            % ZC.Moments.IndicesList (which is a Nx5 array containing the n,l,m
            % indices (col 1-3) and the real (col 4) and imaginary (col 5) ZC moments
            % or 
            % ZC.Moments.Values (which is vector with the complex-valued
            % moments)
            
            inputSize = size(moments);
            
            dim = numel(inputSize);
            
            switch dim
                                
                case 3
                    
                    order = (inputSize(3)-1)/2;
                    
                    nInvariants = ZC.numberOfInvariants(order);
                    
                    Descriptors = zeros(nInvariants,1);
                    
                    inv_count = 0;
                    
                    for n = 0:order
                        
                        % sum_tmp = 0; % NB This is the bug in the original N&K
                        % code that caused the invariants to be cumulatative. Keep
                        % this information for legacy reasons.
                        
                        for l = mod(n,2):2:n
                            
                            sum_tmp = 0;
                            
                            for m=-l:l
                                
                                mom = moments(n+1, l+1, m+l+1);
                                
                                % Legacy code
%                                 absM = abs(m);
%                                 
%                                 %  The ZC_nlm moment
%                                 mom = moments(n+1, l+1, absM+1);
%                                 
%                                 %conjugate if m negative
%                                 if m<0
%                                     mom = conj(mom);
%                                     % take care of sign for odd m
%                                     if mod(absM,2)
%                                         mom = -1*mom;
%                                     end
%                                 end
                                
                                sum_tmp = sum_tmp + norm(mom)^2;
                                % the C++ std:norm function gives the square of the L2
                                %(euclidian norm), which is the so called field norm
                                
                            end
                            
                            inv_count = inv_count + 1;
                            Descriptors(inv_count,1) =  sqrt(sum_tmp);
                            
                        end
                    end
                    
                    
                    
                case 2
                    
                    if any(inputSize) %assume vector
                        % compute the n,l,m indices for
                        order = ZC.expansionOrder(length(moments));                        
                        moments = [ZC.NLMlabels(order) real(moments) imag(moments)];
                    end
                    
                    ZCmoments = complex(moments(:,4), moments(:,5));
                    
                    d = [true; diff(moments(:,2)) ~= 0; true];
                    
                    startRowId = find(d);
                    
                    nInvariants = numel(startRowId)-1;
                    
                    Descriptors = zeros(nInvariants,1);
                    
                    for n = 1:numel(startRowId)-1
                        
                        sumInterval = startRowId(n):(startRowId(n+1)-1);
                        
                        Descriptors(n) = norm(ZCmoments(sumInterval),2);
                        
                    end
              
                    
            end
            
        end
        
    end % method (Static)
    
end % (classDef)

