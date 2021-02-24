using Documenter
using SpawnerRecruitModels

makedocs(
    sitename = "Spawner-Recruit Models",
    format = Documenter.HTML(),
    modules = [SpawnerRecruitModels]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
