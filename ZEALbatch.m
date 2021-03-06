classdef ZEALbatch
    
    % Class for performing multiple ZEAL computations in serial or parallell
    % from structures specified in a Matlab cell array or in a CSV file
    
    % INPUT
   
    % TEXT FILE : ZEAL ALign Mode
    % ---------------------------------
    % CSV format
    % 1_pdbID_or_file,2_pdbID_or_file
    
    % By default, the second structure is defined as the "rotating"
    % structure
    
    % Ex. See ./sample_data/fileList_alignMode.csv
    
    % where each line contains the PDB ID (+optional chain ID) or the filepath to a PDB or CIF file for a pair of structures
    
    % TEXT FILE : ZEAL Single Mode
    % ---------------------------------
    % subsequent lines: 1_pdbID_or_file,2_pdbID_or_file
    
    % where each line contains the PDB ID (+optional chain ID) or the
    % filepath to a PDB or CIF file for a structure
    
    % Ex. See /sample_data/fileList_singleMode.txt in sample_data in this folder.
    
    % CELL ARRAY
    % ---------------------------------
    % batchList_alignMode = {{'S11','S12'},{'S21','S22'},{'S31','S32'}}
    % batchList_singleMode = {'S1','S2,'S3','S4','S5','S6'}
    
    % If using Python, specify batch list as a list of lists
    
    % batchList_alignMode = [["S11","S12"],["S21","S22"],["S31","S32"]]
    % batchList_singleMode = ["S1","S2","S3","S4","S5",S6"]
    
    % batchList_alignMode = [["sample_data/highTMset/1_2JERA_3H7CX/2JERA.pdb", "sample_data/highTMset/1_2JERA_3H7CX/3H7CX.pdb"],["sample_data/highTMset/2_2O90A_5F3MA/2O90A.pdb", "sample_data/highTMset/2_2O90A_5F3MA/5F3MA.pdb"],["3KTIA","5C90A"]]
   
    % ---------------------------------
    % TODO: "DRY" the code 
    % TODO: Unit test class
    % ---------------------------------
    
    % Filip (Persson) Ljung 
    %
    % filip.persson@gmail.com
    % 2020-12
    
    properties
        Batch
        Results
    end
    
    properties (Hidden)
        pool
        AlignMode
    end
    
    methods
        
        function obj = ZEALbatch(batchList, varargin)
            
            % IMPORT JOB LIST
            if iscell(batchList)
                
                [obj.Batch.Fixed, obj.Batch.Rotating] = ZEALbatch.importCell(batchList);
                
            elseif exist(batchList, 'file')
                
                [fixed, rotating] = ZEALbatch.importBatchFile(batchList);
                
                obj.Batch.Fixed = fixed;
                obj.Batch.Rotating = rotating;
                
            else
                
                error('Batch list has to be specified as a matlab cell array or as a txt or csv file.');
                
            end
            
            obj.Batch.N = length(obj.Batch.Fixed);
            obj.AlignMode = ZEALbatch.checkAlignMode(obj.Batch);
            
            % OPTIONAL ARGUMENTS
            
            % defaults
            default_Parallel = false;
            default_NumCores = 4;
            default_PDBoutput = false;
            default_OutputPath = pwd;
            default_KeepMoments = false;
            
            p = inputParser;
            p.KeepUnmatched = true;
            
            addOptional(p, 'Parallel', default_Parallel);
            addOptional(p, 'NumCores', default_NumCores);
            addOptional(p, 'PDBoutput', default_PDBoutput);
            addOptional(p, 'OutputPath', default_OutputPath);
            addOptional(p, 'KeepMoments', default_KeepMoments);
            
            parse(p, varargin{:});
            
            parOp = p.Results.Parallel;
            numCores = p.Results.NumCores;
            PDBoutput = p.Results.PDBoutput;
            outputPath = p.Results.OutputPath;
            keepMoments = p.Results.KeepMoments;
            
            % set up parpool
            if parOp
                
                if isempty(gcp)
                    parpool(numCores)
                end
                
                obj.pool = gcp;
                
                if obj.pool.NumWorkers ~= numCores
                    fprintf('\n Restarting parallel pool to change number of cores to allocate\n');
                    delete(obj.pool)
                    parpool(numCores)
                end
                
            end
            
            % Get options for ZEAL, both default and any supplied in this constructor
            ZEAL_options = ZEAL('dummy', 'InputParserOnly', true, p.Unmatched);
            
            % Preallocate arrays for holding results
            nInvariants = ZC.numberOfInvariants(ZEAL_options.Settings.Order);
            nMoments = ZC.numberOfMoments(ZEAL_options.Settings.Order);
            
            fix_descriptors = zeros(obj.Batch.N, nInvariants);
            rot_descriptors = zeros(obj.Batch.N, nInvariants);
            
           
            
            if keepMoments
                fix_moments = zeros(obj.Batch.N, nMoments);
                rot_moments = zeros(obj.Batch.N, nMoments);              
            end
            
            fix_Rg = zeros(obj.Batch.N, 1);
            rot_Rg = zeros(obj.Batch.N, 1);
            
            score = zeros(obj.Batch.N, 1);
            
            % set non_default options for ZEAL
            if isempty(varargin)
                % zealOptions.Order = 20;
                zealOptions = struct;
            else
                zealOptions = p.Unmatched;
            end
            
            startTime = tic;
            
            try
                
                N = obj.Batch.N;
                errors = 0;
                           
                if obj.AlignMode % ALIGN MODE
                    
                    fix = obj.Batch.Fixed;
                    rot = obj.Batch.Rotating;
                    
                    if parOp
                        
                        % D = parallel.pool.DataQueue;
                        % afterEach(D, @updateParStatus);
                        % parCount = 1;
                                                                        
                        parfor i = 1:obj.Batch.N
                            
                            try
                                fprintf('Doing id %d (%d):\n\t %s\n\t %s', i, N, fix{i}, rot{i})
                                
                                
                                if ~isempty(fix{i}) && ~isempty(rot{i})
                                    shape_i = ZEAL(fix{i}, 'rot', rot{i}, zealOptions);
                                else
                                    warning('No structure pair specified (input number %d)', i);
                                end
                                
                                fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                                rot_descriptors(i,:) = shape_i.rotating.ZC.Descriptors;
                                
                                if keepMoments
                                    fix_moments(i,:) = shape_i.fixed.ZC.Moments.Values;
                                    rot_moments(i,:) = shape_i.fixed.ZC.Moments.Values;                                  
                                end
                                
                                score(i) = shape_i.Score;
                                
                                fix_Rg(i) = shape_i.fixed.Rg;
                                rot_Rg(i) = shape_i.rotating.Rg;
                                
                                if PDBoutput
                                    [~, name, ~] = fileparts(fix{i});
                                    fix_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    [~, name, ~] = fileparts(rot{i});
                                    rot_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    save2pdb(shape_i, 'fixName', fix_savename, 'rotName', rot_savename, 'folderPath', outputPath);
                                    
                                end
                                
                                % D.send(i);
                                
                            catch ME
                                warning(ME.message)
                                errors = errors + 1;
                            end
                            
                        end
                                                
                    else % parOp false
                        
                        for i = 1:N
                            
                            try
                                fprintf('Doing id %d (%d):\n\t %s\n\t %s', i, N, fix{i}, rot{i})
                                
                                if ~isempty(fix{i}) && ~isempty(rot{i})
                                    shape_i = ZEAL(fix{i}, 'rot', rot{i}, zealOptions);
                                else
                                    warning('No structure pair specified (input number %d)', i);
                                end
                                
                                fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                                rot_descriptors(i,:) = shape_i.rotating.ZC.Descriptors;
                                score(i) = shape_i.Score;
                                
                                if keepMoments
                                    fix_moments(i,:) = shape_i.fixed.ZC.Moments.Values;
                                    rot_moments(i,:) = shape_i.fixed.ZC.Moments.Values;                                   
                                end
                                
                                fix_Rg(i) = shape_i.fixed.Rg;
                                rot_Rg(i) = shape_i.rotating.Rg;
                                
                                if PDBoutput
                                    [~, name, ~] = fileparts(fix{i});
                                    fix_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    [~, name, ~] = fileparts(rot{i});
                                    rot_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    save2pdb(shape_i, 'fixName', fix_savename, 'rotName', rot_savename, 'folderPath', outputPath);
                                    
                                end
                                
                            catch ME
                                warning(ME.message)
                                errors = errors + 1;
                            end
                            
                        end
                        
                    end
                    
                    % Gather (adapted for parfor)
                    obj.Results.Fixed.Descriptors = fix_descriptors;
                    obj.Results.Rotating.Descriptors = rot_descriptors;
                    
                    if keepMoments
                        obj.Results.Fixed.Moments.Values = fix_moments;
                        obj.Results.Rotating.Moments.Values = rot_moments;                        
                    end
                                                       
                    obj.Results.Score = score;
                    
                    obj.Results.Fixed.Rg = fix_Rg;
                    obj.Results.Rotating.Rg = rot_Rg;
                    
                else % SINGLE MODE
                    
                    fix = obj.Batch.Fixed;
                    
                    if parOp
                        
                        % D = parallel.pool.DataQueue;
                        % afterEach(D, @updateParStatus);
                        % parCount = 1;
                        
                        parfor i = 1:N
                            
                            try
                                
                                fprintf('\n Doing id %d (%d): \n\t%s', i, N, fix{i})
                                
                                if ~isempty(fix{i})
                                    shape_i = ZEAL(fix{i}, zealOptions);
                                    fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                                    
                                    if keepMoments
                                        fix_moments(i,:) = shape_i.fixed.ZC.Moments.Values;                                                                            
                                    end
                                
                                               
                                    fix_Rg(i) = shape_i.fixed.Rg;

                                else
                                    warning('No structure specified (input number %d)', i);
                                end
                                
                                if PDBoutput
                                    [~, name, ~] = fileparts(fix{i});
                                    fix_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    save2pdb(shape_i, 'fixName', fix_savename, 'folderPath', outputPath);
                                end
                                
                                % D.send(i);
                                
                            catch ME
                                warning(ME.message)
                                errors = errors + 1;
                            end
                                                                                    
                        end
                                                                        
                    else % parOp false
                        
                        for i = 1:N
                            
                            try
                                
                                fprintf('\n Doing id %d (%d): \n\t%s', i, N, fix{i})
                                
                                if ~isempty(fix{i})
                                    shape_i = ZEAL(fix{i}, zealOptions);
                                    fix_descriptors(i,:) = shape_i.fixed.ZC.Descriptors;
                                    
                                    if keepMoments
                                        fix_moments(i,:) = shape_i.fixed.ZC.Moments.Values;                                                                           
                                    end            
                                    
                                    fix_Rg(i) = shape_i.fixed.Rg;
                                    
                                else
                                    warning('No structure specified (input number %d)', i);
                                end
                                
                                if PDBoutput
                                    [~, name, ~] = fileparts(fix{i});
                                    fix_savename = sprintf('%d_%s_ZEAL.pdb', i, name);
                                    
                                    save2pdb(shape_i, 'fixName', fix_savename, 'folderPath', outputPath);
                                end
                                
                            catch ME
                                warning(ME.message)
                                errors = errors + 1;
                            end
                            
                        end
                        
                    end
                    
                    % Gather (adapted for parfor)
                    obj.Results.Fixed.Descriptors = fix_descriptors;
                    obj.Results.Fixed.Rg = fix_Rg;
                    
                    if keepMoments
                        obj.Results.Fixed.Moments.Values = fix_moments;                       
                    end
                    
                end
                
            catch ME
                
                warning(ME.message);
                
            end

            obj.Results.Errors = errors;
            obj.Results.ComputationTime = toc(startTime);
            
            
            fprintf('\n\n Batch job finished.');
            fprintf('\n Total computation time: %s (HH:MM:SS)', datestr(seconds(obj.Results.ComputationTime),'HH:MM:SS'));
            if obj.AlignMode
                fprintf('\n                         %3.2f s / alignment', obj.Results.ComputationTime/N);
            else
                fprintf('\n                         %3.2f s / structure', obj.Results.ComputationTime/N);
            end
            
            fprintf('\n')
            
            %             function [] = updateParStatus(~)
            %
            %                 fprintf('\n Progress: %2.2f', parCount/obj.Batch.N);
            %                 parCount = parCount +1;
            %
            %             end
                        
            
        end
        
        function status = closeParPool(obj)
           
            if ~isempty(obj.pool)
                
               delete(obj.pool) 
               
            end
                 
            status = obj.pool;
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
            rotating = cell(numel(batchList),1);
            
            try
                for i = 1:length(batchList)
                    
                    i_pair = batchList{i};
                    
                    if ~isempty(i_pair)
                        
                        if ~iscell(i_pair)
                            
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