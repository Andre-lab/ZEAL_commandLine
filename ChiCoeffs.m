classdef ChiCoeffs < handle
    % Class to load the object-independent chi coefficients used for ZC
    % moment computations.
    %             % used for the ZC moment computation. The folder path is taken from
    %             % obj.Settings.ChiCoeffPath and the files themself are assumed to
    %             % be mat-files with name 'chiCoeffs_order_X.mat', where X is the
    %             % ZC order.
    
    properties ( SetAccess = private )
        Values
        Indices
        Order
    end
    
    methods
        
        function obj = ChiCoeffs()
            obj.loadData(); % constructor
        end
        
        
        function loadData(obj, varargin)
            % actual loading of data, decoupled from the construcyor so that we can reload if needed

                        
            persistent nQueries % keep track of how many times we call the load function
            
            if isempty(nQueries)
                nQueries = 1;
            else
                nQueries = nQueries + 1;
            end
            
            default_Order = 20;
            default_ChiCoeffFolder = fullfile(pwd,'chi_coefficients');
            
            p = inputParser;
            
            addOptional(p, 'Order', default_Order);
            addOptional(p, 'ChiCoeffFolder', default_ChiCoeffFolder);
            
            parse(p, varargin{:});
            
            obj.Order = p.Results.Order;
            
            ChiCoeffpath = fullfile(p.Results.ChiCoeffFolder, sprintf('chiCoeffs_order_%d.mat', obj.Order));
            
            fprintf('\n ChiCoeffs class queried %d time(s)\n\t Loading Chi coefficients for order %1.0f', nQueries, obj.Order);
            
            data = load(ChiCoeffpath,'chi_coeff_cell','chi_nlm_rst_cell');
            
            obj.Values = data.chi_coeff_cell;
            obj.Indices = data.chi_nlm_rst_cell;

        end
        
    end
    
end