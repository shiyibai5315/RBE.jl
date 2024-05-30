using RBE
using Documenter

DocMeta.setdocmeta!(RBE, :DocTestSetup, :(using RBE); recursive=true)

makedocs(;
    modules=[RBE],
    authors="shiyibai5315",
    sitename="RBE.jl",
    format=Documenter.HTML(;
        canonical="https://shiyibai5315.github.io/RBE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/shiyibai5315/RBE.jl",
    devbranch="main",
)
