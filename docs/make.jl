#using Pkg
#Pkg.activate("../")
#push!(LOAD_PATH,"../")

#using G3RT
using Documenter

makedocs(
    sitename="GSFC EM27/SUN Retrievals",
    remotes=nothing,
    pages = [
        "Main" => "index.md",
        "Basic Usage" => "basic_usage.md",
        "Algorithm Details" => "details.md",
    ]
)
