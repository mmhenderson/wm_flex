%% Train/test inverted encoding model.
% training on one condition of the main task, testing on the other
% condition. See if cross-condition encoding model is worse than within
% condition. 

%%
clear
close all;

sublist = [2:7];
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

nConds=2;

% setting up for IEM
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

nBins=8;
bin_centers=[0:45:359];
bin_size=diff(bin_centers(1:2));

%% loop over subjects
for ss=1:length(sublist)
    
    allchanresp = [];

    substr = sprintf('S%02d',sublist(ss));
    fn2load = fullfile(exp_path,'Samples',sprintf('MainTaskSignalByTrial_%s.mat',substr));
    load(fn2load);
    save_dir = fullfile(exp_path,'Analysis','IEM','IEM_results');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    fn2save = fullfile(save_dir,sprintf('TrnTestAcrossConds_%s.mat',substr));
    
    %% loop over ROIs and run the model for each.

    for vv = 1:length(ROI_names) 
         %% check if there's data here

        if length(mainSig)<vv || isempty(mainSig(vv).dat_avg) || size(mainSig(vv).dat_avg,2)<9
            fprintf('skipping area %s because not enough voxels\n',ROI_names{vv})
            continue
        end
        fprintf('running model for %s %s\n',substr, ROI_names{vv});
        
       %% pull out the data for main task
        mainPosIEM_all = mainSig(vv).targPos;
        mainDatIEM_all = mainSig(vv).dat_avg;
      
        % subtract mean over voxels (third dim)
        mainDatIEM_all = mainDatIEM_all - repmat(mean(mainDatIEM_all,2), 1, size(mainDatIEM_all, 2));

%         runLabsIEM_all = mainSig(vv).runLabs;
        condLabs_all = mainSig(vv).condLabs;

        chan_resp_all = NaN(size(mainDatIEM_all,1), length(xx));   
        
        % bin these for classifier - want 8 bins that are roughly centered at
        % 0, 45, 90, 135,...
        posLabs=mainPosIEM_all;
        binLabs = zeros(size(posLabs));        
        for bb=1:nBins
            inds_this_bin = abs(posLabs-(bin_centers(bb)-0.0001))<bin_size/2 | abs((posLabs-360)-(bin_centers(bb)-0.0001))<bin_size/2;
            binLabs(inds_this_bin) = bb;
        end
        neach_tst = sum(repmat(binLabs,1,nBins)==repmat((1:nBins),size(binLabs,1),1));
        assert(~any(binLabs==0))

        % create cross-validation labels
        % on each fold, want to leave out one pair of trials (one from
        % each label bin) so that training set stays balanced.
        % basically each bin contributes once to each cross-validation
        % fold.           
        cvLabs=zeros(size(condLabs_all));
        for bb=1:nBins
            for cc=1:nConds
                inds = binLabs==bb & condLabs_all==cc;
                cvLabs(inds) = 1:sum(inds);
            end
        end
        assert(~any(cvLabs==0))
        nCV = numel(unique(cvLabs));
        cvLabsIEM_all  = cvLabs;
        
        % looping over conditions to serve as testing set (train on other
        % condition)
        for cond=1:nConds
           
            inds2use_tst=condLabs_all==cond;
            inds2use_trn=condLabs_all~=cond;
            
            mainPosIEM_tst = mainPosIEM_all(inds2use_tst);
            mainDatIEM_tst = mainDatIEM_all(inds2use_tst,:,:);
            cvLabsIEM_tst = cvLabsIEM_all(inds2use_tst);
       
            mainPosIEM_trn = mainPosIEM_all(inds2use_trn);
            mainDatIEM_trn = mainDatIEM_all(inds2use_trn,:,:);
            cvLabsIEM_trn = cvLabsIEM_all(inds2use_trn);
       
        %%  generate stim mask

            stim_mask = zeros(numel(mainPosIEM_trn),length(xx));

            for tt = 1:numel(mainPosIEM_trn)   
                % find the nearest index.
                [~,my_ind] = min(abs(mainPosIEM_trn(tt)-xx));

                stim_mask(tt,my_ind)=1; % put "1" in the right orientation column
            end


            %% train and test the model

            chan_resp = NaN(size(mainDatIEM_tst,1), length(xx));   

            unvals = unique(cvLabsIEM_tst);
            assert(all(unvals==unique(cvLabsIEM_trn)));
            
            for rr = 1:numel(unique(cvLabsIEM_tst))

                trninds = cvLabsIEM_trn~=unvals(rr);
                tstinds = cvLabsIEM_tst==unvals(rr);

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
                    w = trnX\mainDatIEM_trn(trninds,:); % compute weights

                    chan_resp(tstinds,chan_center) = ((w*w')\w*mainDatIEM_tst(tstinds,:)').';

                end % loop over basis set centers.

            end
            
            chan_resp_all(inds2use_tst,:) = chan_resp;
            
        end
        
        chan_resp = chan_resp_all;
        mainPosIEM_tst = mainPosIEM_all;
        %% re-center and average

        % where do you want them to be centered? this works well b/c center of
        % xx axis.
        shift_to = 180;

        chan_resp_shift = NaN(size(chan_resp));
        for t = 1:size(chan_resp,1)
            chan_resp_shift(t,:) =  wshift('1D', chan_resp(t,:),ceil(mainPosIEM_tst(t)-shift_to));
        end

        allchanresp(vv).chan_resp = chan_resp;
        allchanresp(vv).chan_resp_shift = chan_resp_shift;
        
        allchanresp(vv).condLabs = condLabs_all;
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