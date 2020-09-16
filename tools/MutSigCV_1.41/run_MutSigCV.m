function run_MutSigCV(path_mut, path_cov, path_gen, path_out, path_dict, path_ref)
    %   Created on Wed Sep 16 2020.
    %
    %   author: Yoann Pradat
    %
    %   Runs the package MutSigCV v1.41
    %
    %   Parameters
    %   ---------
    %   path_mut:
    %       path to mutation file
    %   path_cov:
    %       path to coveage file
    %   path_gen:
    %       path to gene covariates file
    %	path_out: string
    %	    Path to out folder
    %	path_dict: string
    %	    Path to mutation dictionary for defining effects from MAF's Variant_Classification column.
    % 
    %   Example
    %   -------
    %   # fusion parameters
    %   matlab -nodisplay -nodesktop -r "run_MutSigCV('$path_mut', '$path_cov', '$path_gen', '$path_dict', '$path_ref')"
    %

    %% Display input parameters
    disp(['Input parameters'])
    disp(['path_mut: ', path_mut])
    disp(['path_cov: ', path_cov])
    disp(['path_gen: ', path_gen])
    disp(['path_out: ', path_out])
    disp(['path_dict: ', path_dict])
    disp(['path_ref: ', path_ref])
    
    %% Add paths to source code
    fprintf('Working directory: %s\n', pwd);
    addpath('./src');


    %% Run MutSigCV
    ts = tic;
    
	MutSigCV(...
	    path_mut,...           % mutation counts
	    path_cov,...           % coverage counts
	    path_gen,...           % gene covariates
        path_out,...           % out folder
	    path_dict,...          % mutation effects dict
	    path_ref);             % folder ref genome
    
    te = toc(ts);
    disp(['total running time MutSigCV: ', num2str(te), ' seconds']);

    %% Exit program
    exit;
end
