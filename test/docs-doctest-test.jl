@testitem "Documentation doctests" tags=[:docs, :doctest] begin
    using Documenter: DocMeta, doctest
    using ArchimedLight

    DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true)

    doctest(ArchimedLight; manual=false)
end
