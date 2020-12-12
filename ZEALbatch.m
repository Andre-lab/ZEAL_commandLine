classdef ZEALbatch
    
    % Class for performing multiple ZEAL computations in serial or parallell
    % from structures specified in a Matlab cell array or in a CSV file
    
    % INPUT
    
    % CELL ARRAY
    % batchList = {{'S11','S12'},{'S21','S22'},{'S31','S32'}}
    
    % If using Python, specify batch list as a list of lists
    % batchList = [["S11","S12"],["S21","S22"],["S31","S32"]]
    
    % Python: get second structure-pair
    % batchList[1]
    %
    % Python: get second structure in second structure pair
    % batchList[1][1]
    
    % CSV FILE
    % CSV format: ZEAL ALign Mode
    %1_pdbID_or_file,2_pdbID_or_file
    
    % By default, the second structure is defined as the "rotating"
    % structure
    
    % Ex. See ./sample_data/fileList_alignMode.csv
    
    % where each line contains the PDB ID (+optional chain ID) or the filepath to a PDB or CIF file for a pair of structures
    
    % CSV format: ZEAL Single Mode
    
    % subsequent lines: 1_pdbID_or_file,2_pdbID_or_file
    
    % where each line contains the PDB ID (+optional chain ID) or the
    % filepath to a PDB or CIF file for a structure
    
    % Ex. See /sample_data/fileList_singleMode.csv in sample_data in this folder.
    
    
    
    %in batch mode, where a ZEAL session is started
    % for
    % C
    % batchFile = 'sample_data/fileList_download_align.csv';
    %
    % parop = false;
    % nworkers = 2;
    %
    % if parop
    %
    %     if isempty(gcp)
    %         parpool(nworkers)
    %     end
    %
    %     pool = gcp;
    %
    % end
    
    properties
        Batch
        AlignMode
        Results
    end
    
    methods
        
        function obj = ZEALbatch(batchList, varargin)
            
            if iscell(batchList)
                
                [obj.Batch.Fixed, obj.Batch.Rotating] = ZEALbatch.importCell(batchList);
                
            elseif exist(batchList, 'file')
                
                [fixed, rotating] = ZEALbatch.importBatchFile(batchList);
                
                obj.Batch.Fixed = fixed;
                obj.Batch.Rotating = rotating;
                    
            else
                
                error('Batch list has to be specified as a matlab cell array or as a csv-file.');
                
            end
            
            obj.Batch.N = length(obj.Batch.Fixed);
            obj.AlignMode = ZEALbatch.checkAlignMode(obj.Batch);
            
            
            % Get options for ZEAL, both default and any supplied in this constructor
            
            ZEAL_options = ZEAL('dummy', 'InputParserOnly', true, varargin{:});
            
            % Preallocate arrays for holding results
            nInvariants = ZC.numberOfInvariants(ZEAL_options.Settings.Order);
            
            fix_descriptors = zeros(obj.Batch.N, nInvariants);
            rot_descriptors = zeros(obj.Batch.N, nInvariants);
            score = zeros(obj.Batch.N, 1);
            
            if isempty(varargin)
                zealOptions = {'Order', 20};
            else                
                zealOptions = varargin;
            end
            
            try
                
                if obj.AlignMode % ALIGN MODE
                    
                    fix = obj.Batch.Fixed;
                    rot = obj.Batch.Rotating;
                    
                    % parfor i = 1:obj.Batch.N
                    for i = 1:obj.Batch.N
                        
                        if ~isempty(fix{i}) && ~isempty(rot{i})
                            shape_i = ZEAL(fix{i}, 'rot', rot{i}, zealOptions{:});
                        else
                            warning('No structure pair specified (input number %d)', i);
                        end
                        
                        fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                        rot_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                        score(i) = shape_i.Score;
                    end
                    
                    % Gather (adapted for parfor)
                    obj.Results.Fixed.Descriptors = fix_descriptors;
                    obj.Results.Rotating.Descriptors = rot_descriptors;
                    obj.Results.Score = score;
                    
                else % SINGLE MODE
                    
                    fix = obj.Batch.Fixed;
                    
                    % parfor i = 1:obj.Batch.N
                    for i = 1:obj.Batch.N
                        
                        if ~isempty(fix{i})
                            shape_i = ZEAL(fix{i}, zealOptions{:});
                            fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                        else
                            warning('No structure specified (input number %d)', i);
                        end
                    end
                    
                    % Gather (adapted for parfor)
                    obj.Results.Fixed.Descriptors = fix_descriptors;
                    
                end
                
                
            catch ME
                
                error(ME.message);
                
            end
            
        end
        
    end
    
    
    methods (Static)
        
        function [fixed, rotating] = importBatchFile(batchFile)
            
            try
                           
                fid = fopen(batchFile);
                tline = fgetl(fid);
                
                fixed = {};
                rotating={};
                count = 0;
                
                while ischar(tline)
                    
                    count = count + 1;
                    
                    if ischar(tline)
                        row_data = split(tline,',');
                        
                        if ~isempty(row_data{1})
                            fixed(count) = row_data(1);
                        end
                        
                        if length(row_data)>1 % rot structure specified
                            
                            if ~isempty(row_data{2})
                                rotating(count) = row_data(2);
                            end
                            
                        else
                            rotating(count) = {''};
                        end
                        
                    end
                    
                    tline = fgetl(fid);
                    
                end
                
                fclose(fid);
                           
            catch ME
                
                fprintf(ME.message);
                
            end
                        
        end
        
        function [fixed, rotating] = importCell(batchList)
            
            fixed = cell(numel(batchList),1);
            rotating = fixed;
            
            try
                for i = 1:length(batchList)
                    
                    i_pair = batchList{i};
                    
                    if ~isempty(i_pair)
                        
                        if numel(batchList(i)) == 1
                            
                            fixed{i} = i_pair;
                            
                            rotating{i} = '';
                            
                        else
                            
                            if ~isempty(i_pair{1})
                                fixed(i) = i_pair(1);
                            else
                                fixed(i) = {''};
                            end
                            
                            if ~isempty(i_pair{2})
                                rotating(i) = i_pair(2);
                            else
                                rotating(i) = {''};
                            end
                            
                        end
                        
                    else
                        
                        fixed{i} = '';
                        rotating{i} = '';
                        
                    end
                    
                end
                
            catch ME
                
                fprintf(ME.message);
                
            end
            
        end
        
        function AlignMode = checkAlignMode(batch)
            
            %fixed_empty_tf = cellfun(@isempty,batch.Rixed);
            rotating_empty_tf = cellfun(@isempty,batch.Rotating);
                
                if sum(rotating_empty_tf) > 0
                    AlignMode = false;
                else
                    AlignMode = true;
                end     
                
        end
        
        
        
        
    end % methods (Static)
    
end