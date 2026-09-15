function [Flux_t, solution1, D, AUC] = process_data_host_diff(Flux_t, solution1, t, model)
    Flux_t = Flux_t(1:length(model.rxns),:);
    Sol_ref = solution1.v(1:length(model.rxns));
    D = (Flux_t-Sol_ref);
    D(abs(D)<1e-5)=0;
    AUC = cumtrapz(t,D,2);
end