using MimiGIVE
println("NO EW")
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45, ew=nothing, ew_norm_region=nothing)];
result = MimiGIVE.compute_scc(year = 2030, discount_rates = discount_rates, n = 2, compute_sectoral_values = true, save_md = true, save_cpc = true, certainty_equivalent=true);
result[:scc]
# fails for sectors

println("GDP COUNTRY")
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45, ew=:gdp_country, ew_norm_region="USA")];
result = MimiGIVE.compute_scc(year = 2030, discount_rates = discount_rates, n = 2, compute_sectoral_values = true, save_md = true, save_cpc = true, certainty_equivalent=true);
result[:scc]
# works but identical

println("GDP REGION")
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45, ew=:gdp_region, ew_norm_region="USA")];
result = MimiGIVE.compute_scc(year = 2030, discount_rates = discount_rates, n = 2, compute_sectoral_values = true, save_md = true, save_cpc = true, certainty_equivalent=true);
result[:scc]
# works but identical

println("CONSUMPTION COUNTRY")
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45, ew=:consumption_country, ew_norm_region="USA")];
result = MimiGIVE.compute_scc(year = 2030, discount_rates = discount_rates, n = 2, compute_sectoral_values = true, save_md = true, save_cpc = true, certainty_equivalent=true);
result[:scc]
# fails for sectors

println("CONSUMPTION REGION")
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45, ew=:consumption_region, ew_norm_region="USA")];
result = MimiGIVE.compute_scc(year = 2030, discount_rates = discount_rates, n = 2, compute_sectoral_values = true, save_md = true, save_cpc = true, certainty_equivalent=true);
result[:scc]
# fails for sectors

# eta = 0 fixes it