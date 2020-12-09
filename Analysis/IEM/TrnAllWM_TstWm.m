% MMH 3/11/20- train/test inverted encoding model for oriSpin data. 

% training/testing on SPATIAL POSITION localizer, sanity check to make sure we
% can get nicely tuned responses.
clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));


for ss=1:length(sublist)
    
    allchanresp = [];

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','IEM','IEM_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnAllWM_TstWM_%s.mat',substr));
    
    %% loop over ROIs and run the model for each.

    for vv = 1:length(ROI_names) 
         %% check if there's data here

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s %s\n',substr, ROI_names{vv});
        
        %% pull out the data for localizer (training)
        mainPosIEM = mainSig(vv).targPos;
        mainDatIEM = mainSig(vv).dat_avg;
        
        mainDatIEM = mainDatIEM - repmat(mean(mainDatIEM,2), 1, size(mainDatIEM, 2));
       
        runLabsIEM = mainSig(vv).runLabs;
        condLabs = mainSig(vv).condLabs;
        
        randBoundPos = mainSig(vv).randBoundPos;
        boundPos = mainSig(vv).boundPos;
        %% Do the IEM - train and test within the localizer, cross-validated.

        n_spat_chans = 9; 
        n_deg = 360;
        
        if mod(n_deg,n_spat_chans)~=0
            fprintf('\nn_ori_chans needs to go into 180 evenly\nReturning...\n')
            return;
        end

        basis_pwr = n_spat_chans-1; % good default is - 1, but as noted can play with this to see what happens... critical to understand.

        % define function to make basis funcs
        % HAVE TO DIVIDE BY 2 HERE SO THAT THE SPACE WORKS OUT
        make_basis_function = @(xx,mu) (cosd((xx-mu)/2)).^basis_pwr;
        % in circular space xx would go from 0-pi, mu would be between 0 and
        % pi, and the function would use cos instead of cosd

        % my x-axis (offset by half-degrees so that the wedge masks can be
        % perfectly centered).
        xx = linspace(1,n_deg,n_deg) - 0.5;


        %%  generate stim mask

        stim_mask = zeros(numel(mainPosIEM),length(xx));
    
        for tt = 1:numel(mainPosIEM)   
            % find the nearest index.
            [~,my_ind] = min(abs(mainPosIEM(tt)-xx));
           
            stim_mask(tt,my_ind)=1; % put "1" in the right orientation column
        end


        %% train and test the model

        chan_resp = NaN(size(mainDatIEM,1), length(xx));   

        unvals = unique(runLabsIEM);
        
        for rr = 1:numel(unique(runLabsIEM))

            trninds = runLabsIEM~=unvals(rr);
            tstinds = runLabsIEM==unvals(rr);

            % start the loop from 1:20 (e.g. if num chans == 9). Each time shift the
            % basis set centers (mu's) to completely cover the space after all
            % iterations.

            for b = 1:(n_deg/n_spat_chans)

                basis_set = NaN(n_deg,n_spat_chans); % basis-set can go in here
                chan_center = b:n_deg/n_spat_chans:n_deg; % my channel centers on this itteration of b

                for cc = 1:n_spat_chans
                    basis_set(:,cc) = make_basis_function(xx,chan_center(cc));
                end

                % now generate the design matrix
                trnX = stim_mask(trninds,:)*basis_set;

                if rank(trnX)~=size(trnX,2)
                    fprintf('\nrank deficient training set Design Matrix\nReturning...\n')
                    return;
                end

                % then solve for the weights!
                w = trnX\mainDatIEM(trninds,:); % compute weights
                
                chan_resp(tstinds,chan_center) = ((w*w')\w*mainDatIEM(tstinds,:)').';

            end % loop over basis set centers.

        end

        %% re-center and average

        % where do you want them to be centered? this works well b/c center of
        % xx axis.
        shift_to = 180;

        chan_resp_shift = NaN(size(chan_resp));
        for t = 1:size(chan_resp,1)
            chan_resp_shift(t,:) =  wshift('1D', chan_resp(t,:),ceil(mainPosIEM(t)-shift_to));
        end

        allchanresp(vv).chan_resp = chan_resp;
        allchanresp(vv).chan_resp_shift = chan_resp_shift;
        
        allchanresp(vv).condLabs = condLabs;
        allchanresp(vv).targPos = mainSig(1).targPos;
        allchanresp(vv).boundPos = mainSig(1).boundPos;
        allchanresp(vv).randBoundPos = mainSig(1).randBoundPos; 
        allchanresp(vv).RTLabs = mainSig(1).RTLabs;
        allchanresp(vv).dist_to_real_bound = mainSig(1).dist_to_real_bound;
        allchanresp(vv).dist_to_rand_bound = mainSig(1).dist_to_rand_bound;
        allchanresp(vv).dir_to_real_bound = mainSig(1).dir_to_real_bound;
        allchanresp(vv).dir_to_rand_bound = mainSig(1).dir_to_rand_bound;
        
    end


    fprintf('saving to %s\n',fn2save);
    save(fn2save,'allchanresp','ROI_names','xx','shift_to');

end