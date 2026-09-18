using PrecompileTools

@compile_workload begin
    # The SSP socioeconomics source reads only data bundled with this package, so
    # the whole headline path can be exercised here. `compute_scc` is the
    # expensive part and the one that matters most: it is what pulls the CIAM
    # model build, and `MimiCIAM.run_timestep_MimiCIAM_slrcost` in particular,
    # into the cache.
    m = get_model(socioeconomics_source = :SSP, SSP_scenario = "SSP245")
    run(m)
    m[:DamageAggregator, :total_damage]
    compute_scc(m, year = 2020)

    # The default source is :RFF, which differs from :SSP in four components
    # (:Socioeconomic, :Agriculture and the two Agriculture aggregators). Its
    # data is a 1.5 GB DataDeps download that must not happen at build time, but
    # building a model never reads component input data -- that happens in `init`
    # at run time -- so we can build the RFF model and compile its components
    # without running it.
    Mimi.precompile_model(get_model())
end
