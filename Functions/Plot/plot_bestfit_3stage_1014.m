%
% plot_bestfit_3stage_1014(sample_data,best_fits,model_interval,expo_intervals,bur_intervals,x_lim,y_lim)
%
%%

function iso_plot = plot_bestfit_3stage_1014(sample_data,best_fits,model_interval,expo_intervals,bur_intervals,x_lim,y_lim)

bur_dur = best_fits.bestfit_bur_dur;
model_time = best_fits.bestfit_model_time;
if isfield(best_fits,'bestfit_burfracs')
    burfracs = best_fits.bestfit_burfracs;
end

if numel(bur_dur)>1
    bur_dur = mean(bur_dur);
    model_time = mean(model_time);
end


% Determine whether samples have multiple nuclide measurements
logical_1014 = any(sample_data.logical_10 & sample_data.logical_14);


% Re-organise data for each sample
for a = 1:length(sample_data.s)
    data{a}.name = sample_data.s{a}.name{1};
    data{a}.cover = sample_data.cover;
    %data{a}.ee = ee;
    data{a}.pp = sample_data.pp;
    data{a}.sf10 = sample_data.sf1026(a);
    data{a}.cp10 = sample_data.cp1026(a);
    data{a}.sf14 = sample_data.sf14(a);
    data{a}.cp14 = sample_data.cp14(a);
    data{a}.top_z_gcm2 = sample_data.s{a}.top_z_gcm2;
    data{a}.bottom_z_gcm2 = sample_data.s{a}.bottom_z_gcm2;
    data{a}.weight10 = sample_data.s{a}.weight10;
    data{a}.weight14 = sample_data.s{a}.weight14;
    data{a}.nuclide10 = sample_data.s{a}.nuclide10;
    data{a}.nuclide14 = sample_data.s{a}.nuclide14;
    if sample_data.s{a}.nuclide10 == 1
        data{a}.N10 = sample_data.s{a}.N10;
        data{a}.dN10 = sample_data.s{a}.dN10;
    end
    if sample_data.s{a}.nuclide14 == 1
        data{a}.N14 = sample_data.s{a}.N14;
        data{a}.dN14 = sample_data.s{a}.dN14;
        data{a}.meanage14 = sample_data.ages.C14(a,1);
    end
    ages14(a)=data{a}.meanage14;
end

if numel(model_time)>1 || numel(bur_dur)>1
    model_time = mean(model_time);
    bur_dur = mean(bur_dur);
end

% Generate exposure and burial periods
recent_exp = mean(ages14);
exposed_or_not = three_stage_logical(recent_exp,bur_dur,model_time,model_interval);


% Evaluate each sample
for c = 1:length(data)
    
    data_s = data{c};
    
    
    % If no burial fractions are specified, run once per sample
    if ~isfield(best_fits,'bestfit_burfracs')
        
        % Calculate predicted concentrations
        [this_s_predN,~] = fit_sample(exposed_or_not,data_s);
        
        scenario_predN{c} = this_s_predN;
        
        
    % Otherwise, do range of burial fractions for each sample
    else
        
        this_frac_buried = burfracs(c);
        this_recent_exp = data_s.meanage14;
        
        % Generate exposure and burial periods
        this_exposed_or_not = three_stage_logical(this_recent_exp,bur_dur,model_time,model_interval,this_frac_buried);
        
        % Calculate predicted concentrations
        [this_s_predN,~] = fit_sample(this_exposed_or_not,data_s);
        
        scenario_predN{c} = this_s_predN;
        
    end
    
end


% Plot concentrations

plot_fig = plot_concs_1014(sample_data,2,1,expo_intervals,bur_intervals,x_lim,y_lim);

% Plot time-series data
ts_plot = plot_expbur_time(exposed_or_not,plot_fig);

% Plot predicted concentrations for samples
for d = 1:length(scenario_predN)
    this_predN = scenario_predN{d};
    plot_concs_1014_pred(sample_data,plot_fig,this_predN,[]);
end

% Plot predicted concentration pathways for samples
%plot_pred_pathway(sample_data,plot_fig,predN,[]); % FIX

% Plot exposure/burial periods
plot_ExpBur(exposed_or_not,ts_plot,plot_fig);



%%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate each sample %%%%%%%%%%%%%%%%%%%%%%%%%

    function [predN,misfit] = fit_sample(exposed_or_not,data_s)
        
        ee = [];
        cover = [];
        
        % Calculate predicted concentrations
        predN = forward_model_1014(exposed_or_not,data_s,ee,cover);
        
        
        % Calculate misfit
        if logical_1014
            
            misfit_10 = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            misfit_14 = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            misfit = mean([misfit_10,misfit_14]); % Average misfits of each nuclide
        else
            if data_s.nuclide10 == 1
                misfit = ((predN.sum.N10 - data_s.N10)./data_s.dN10) .^2;
            end
            if data_s.nuclide14 == 1
                misfit = ((predN.sum.N14 - data_s.N14)./data_s.dN14) .^2;
            end
        end
        
    end

% Export figure handles
iso_plot = gcf;

end
