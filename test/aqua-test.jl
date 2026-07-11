@testitem "Aqua package quality" tags = [:aqua, :quality, :fast] begin
    using Aqua
    using ArchimedLight

    Aqua.test_all(ArchimedLight)
end
