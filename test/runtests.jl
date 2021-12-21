using Mimi, Test, MimiGIVE

save_list = [(:temperature, :T), (:co2_cycle, :co2)]

# ------------------------------------------------------------------------------
# Run model and mcs with default settings

m = MimiGIVE.get_model()
run(m)

mcs_results = MimiGIVE.run_mcs(m = m, trials = 10, output_dir = nothing, save_trials = false) # results in memory but nothing saved
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10, output_dir = nothing, save_trials = false, save_list = save_list) # results in memory
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10, output_dir = nothing, save_trials = false, save_list = save_list, results_in_memory = false) # results NOT in memory

MimiGIVE.compute_scc(m, year = 2020) # deterministic

MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5); # mcs no save list
MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5, save_list = save_list); # mcs with save list

# ------------------------------------------------------------------------------
# Run model and mcs with specific SSP socioeconomics

m = MimiGIVE.get_model(socioeconomics_source = :SSP,
                SSPmodel = "Benveniste",
                SSP = "SSP5",
                RCPmodel = "Leach",
                RCP = "RCP8.5")
run(m)

mcs_results = MimiGIVE.run_mcs(m = m, trials = 10; output_dir = nothing, save_trials = true) # results in memory but nothing saved
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10; output_dir = nothing, save_trials = true, save_list = save_list) # results in memory
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10, output_dir = nothing, save_trials = false, save_list = save_list, results_in_memory = false) # results NOT in memory

MimiGIVE.compute_scc(m, year = 2020) # deterministic

MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5); # mcs no save list
MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5, save_list = save_list); # mcs with save list


# ------------------------------------------------------------------------------
# Run model and mcs with RFF socioeconomics

m = MimiGIVE.get_model(socioeconomics_source = :RFF)
run(m)

mcs_results = MimiGIVE.run_mcs(m = m, trials = 10; output_dir = nothing, save_trials = true) # results in memory but nothing saved
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10; output_dir = nothing, save_trials = true, save_list = save_list) # results in memory
mcs_results = MimiGIVE.run_mcs(m = m, trials = 10, output_dir = nothing, save_trials = false, save_list = save_list, results_in_memory = false) # results NOT in memory

MimiGIVE.compute_scc(m, year = 2020) # deterministic

MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5); # mcs no save list
MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)], n = 5, save_list = save_list); # mcs with save list

# ------------------------------------------------------------------------------
# Errors
@test_throws ErrorException m = MimiGIVE.get_model(socioeconomics_source= :False) # not a valid socioeconomics source

# ------------------------------------------------------------------------------
# Other

# SCC for a shorter time horizon (end in 2200)
m = MimiGIVE.get_model()

MimiGIVE.compute_scc(m, year = 2020, last_year = 2200)
MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45)], last_year = 2200, n = 5) # mcs
results = MimiGIVE.compute_scc(m, year = 2020, discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45)], last_year = 2200, n = 5, save_md = true) # mcs
