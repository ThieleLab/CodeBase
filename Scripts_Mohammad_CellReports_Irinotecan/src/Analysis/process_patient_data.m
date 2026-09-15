function [Flux_t, solution, t, Delta_vi, D, AUC] = process_patient_data(file)
% Define function to process patient data    
    load(file, 't', 'Flux_t', 'solution', 'solution1');
    solution = solution1;
    Flux_t(abs(Flux_t) <1e-5) =0;
    solution.v(abs(solution.v) <1e-5) =0;
    Delta_vi = Flux_t - solution.v;
    D = Delta_vi;
    D(abs(D)< 1e-5) = 0;
    AUC = cumtrapz(t, D, 2);
end