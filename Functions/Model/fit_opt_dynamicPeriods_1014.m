%
function misfit = fit_opt_dynamicPeriods_1014(X,sample_data,expo_mids,model_interval,misfit_type,plot_fig,burialfrac_bnds)


% Unpack
model_time = X(1)/1000;
expo_dur = X(2)/1000;

data = sample_data.org_data;
ages14 = sample_data.ages14;


% Generate exposure and burial periods
recent_exp = mean(ages14)/1000;
exposed_or_not = var_periods_logical(model_time,expo_mids,expo_dur,'exposure',[],recent_exp,model_interval);
        

if ~isempty(burialfrac_bnds)
%     % Determine burial fractions for each sample depending on
%     % elevation order
%     burial_frac_int = 0.05; % Interval for burial fractions within given range
%     burial_fracs = burial_frac_bnds(1):burial_frac_int:burial_frac_bnds(2)+burial_frac_int; % Burial ractions
%     
%     [~,pos_idx] = sort(sample_data.position,'descend'); % Get sample position order
%     sample_arr = 1:length(sample_data.position); sample_arr_rev = fliplr(sample_arr);
%     for c = sample_arr % Create array for each sample, staggered by position
%         this_pos_idx = pos_idx(c);
%         sample_burial_fracs(:,this_pos_idx) = burial_fracs(sample_arr(c):end-sample_arr_rev(c));
%     end
    
    % Determine burial fractions for each sample depending on
    % relative elevation difference
    n_bur_frac = 10; % Set number of iterations to test burial fraction
    bur_frac_arr = linspace(burialfrac_bnds(1),burialfrac_bnds(2),n_bur_frac);
    pos_scaled = rescale(sample_data.position);
    for ii = 1:length(bur_frac_arr)
        sample_burial_fracs(ii,:) = bur_frac_arr(ii)*pos_scaled;
    end
    
    % Plot
    %figure; plot(sample_burial_fracs,sample_data.position,'o-');
    
end


% Evaluate each sample
for ii = 1:length(data)
    
    data_s = data{ii};
    
    % If no burial fractions are specified, run once per sample
    if isempty(burialfrac_bnds)
        
        % Calculate predicted concentrations and misfit
        [this_s_predN,this_s_misfit] = fit_sample(exposed_or_not,data_s);
        
        this_scenario_predN{ii} = this_s_predN;
        this_scenario_misfit(ii) = this_s_misfit;
        
    % Otherwise, do range of burial fractions for each sample
    else
        
        for iii = 1:length(sample_burial_fracs(:,1))
            
            this_frac_buried = sample_burial_fracs(iii,ii);
            this_recent_exp = data_s.meanage14/1000;
            
            % Generate exposure and burial periods
            this_exposed_or_not = var_periods_logical(model_time,expo_mids,expo_dur,'exposure',[],this_recent_exp,model_interval,this_frac_buried);

            % Calculate predicted concentrations and misfit
            [this_burfrac_predN,this_s_misfit] = fit_sample(this_exposed_or_not,data_s,misfit_type);
            
            burfrac_predN{iii} = this_burfrac_predN;
            this_burfrac_misfit(iii) = this_s_misfit;
            
        end
        
        % Get the bestfit burial fraction
        if numel(unique(this_burfrac_misfit))==numel(this_burfrac_misfit)
            burfrac_min_misfit = min(this_burfrac_misfit);
            this_s_predN = burfrac_predN(burfrac_min_misfit==this_burfrac_misfit);
            this_s_bestfit_burfrac = sample_burial_fracs((burfrac_min_misfit==this_burfrac_misfit),ii);
        else % Use the first result if no unique bestfit
            burfrac_min_misfit = min(this_burfrac_misfit);
            this_s_predN = burfrac_predN(1);
            this_s_bestfit_burfrac = sample_burial_fracs(1,ii);
        end
        this_scenario_predN{ii} = this_s_predN{1};
        this_scenario_misfit(ii) = burfrac_min_misfit;
        this_scenario_burfrac(ii) = this_s_bestfit_burfrac;
        
    end
end

if strcmpi(misfit_type,'all')
    misfit = mean(this_scenario_misfit,'omitnan');
elseif strcmpi(misfit_type,'minmax')
    %misfit_mean = mean(this_scenario_misfit,'omitnan');
    misfit_min = min(this_scenario_misfit,[],'omitnan');
    misfit_max = max(this_scenario_misfit,[],'omitnan');
    misfit = mean([misfit_min,misfit_max]);
    %misfit = mean([misfit_mean,misfit_minmax]);
elseif strcmpi(misfit_type,'minBe')
    [~,minBe_idx] = min(sample_data.CC.Be10(:,9),[],'omitnan');
    misfit = this_scenario_misfit(minBe_idx);
elseif strcmpi(misfit_type,'maxBe')
    [~,maxBe_idx] = max(sample_data.CC.Be10(:,9),[],'omitnan');
    misfit = this_scenario_misfit(maxBe_idx);
elseif strcmpi(misfit_type,'minC')
    [~,minC_idx] = min(sample_data.CC.C14(:,9),[],'omitnan');
    misfit = this_scenario_misfit(minC_idx);
elseif strcmpi(misfit_type,'maxC')
    [~,maxC_idx] = max(sample_data.CC.C14(:,9),[],'omitnan');
    misfit = this_scenario_misfit(maxC_idx);
end


% Plot concentrations for scenario
if ~isempty(plot_fig)
    if plot_fig(1,2).Children(1).Color == [0,1,0]
        delete(plot_fig(1,2).Children(1:length(this_scenario_predN)));
    end
    plot_Scenario(exposed_or_not,sample_data,this_scenario_predN,plot_fig)
    
    if 0 % Change to zero to not plot concs vs elevation
        plot_14Celev(sample_data,this_scenario_predN)
    end
end
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate each sample %%%%%%%%%%%%%%%%%%%%%%%%%
  
    function [predN,misfit] = fit_sample(exposed_or_not,data_s,misfit_type)
        
        ee = [];
        cover = [];
        
        % Calculate predicted concentrations
        predN = forward_model_1014(exposed_or_not,data_s,ee,cover);
        
        
        % Calculate misfit
        if sample_data.logical_1014
            if strcmpi(misfit_type,'minBe') || strcmpi(misfit_type,'maxBe')
                misfit = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            elseif strcmpi(misfit_type,'minC') || strcmpi(misfit_type,'maxC')
                misfit = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            else
                misfit_10 = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
                misfit_14 = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
                misfit = mean([misfit_10,misfit_14]); % Average misfits of each nuclide
            end
        else
            if data_s.nuclide10 == 1
                misfit = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            end
            if data_s.nuclide14 == 1
                misfit = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            end
        end
        
    end
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    function plot_Scenario(exposed_or_not,data,scenario_predN,plot_fig)
        
        figure(plot_fig(1));
        
        % Plot time-series data
        ts_plot = plot_expbur_time(exposed_or_not,plot_fig);
        
        % Plot predicted concentrations for samples
        for iv = 1:length(scenario_predN)
            this_predN = scenario_predN{iv};
            plot_concs_1014_pred(data,plot_fig,this_predN,[]);
        end
        
        % Plot predicted concentration pathways for samples
        %plot_pred_pathway(sample_data,plot_fig,predN,[]); % FIX
        
        % Plot exposure/burial periods
        plot_ExpBur(exposed_or_not,ts_plot,plot_fig);
        
        drawnow;
        
    end


  %%%%%%%%%%%%%%%%%%%%%%%%%% Plot 14C vs elevation %%%%%%%%%%%%%%%%%%%%%%%%
  
    function plot_14Celev(sample_data,scenario_predN)
        
        fig = figure(20);
        
        pred_col = 'g';
        
        %delete(fig.Children(1:length(this_scenario_predN)));
        
        h = errbar(sample_data.CC.C14(:,9)/1000,sample_data.position,sample_data.CC_uncert.C14(:,9)/1000,'-','horiz'); hold on;
        set(h,'Color','k','LineWidth',1.5);
        plot(sample_data.CC.C14(:,9)/1000,sample_data.position,'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1);
        for iv = 1:length(scenario_predN)
            this_predN = scenario_predN{iv};
            plot(this_predN.sum.N14/1000,sample_data.position(iv),'o','Color',pred_col,'MarkerFaceColor',pred_col,'MarkerSize',8);
        end
        hold off;
        xlabel('^{14}C (k atoms/g)')
        ylabel('Relative elevation (m)')

    end

end
