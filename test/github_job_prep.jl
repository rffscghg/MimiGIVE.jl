using Pkg

Pkg.add("DataDeps")

using DataDeps

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

register(DataDep(
    "rffsps_v5",
    "RFF SPs version v5",
    "https://zenodo.org/record/6016583/files/rffsps_v5.7z",
    "a39b51d7552d198123b1863ea5a131533715a7a7c9ff6fad4b493594ece49497",
    post_fetch_method=unpack
))

datadep"rffsps_v5"
