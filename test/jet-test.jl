@testitem "JET public API smoke" tags = [:jet, :quality, :fast] begin
    using ArchimedLight
    using Dates
    using JET

    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.LightOptions()
    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.translucent(par=0.15, nir=0.90)
    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.virtual_sensor()
    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.emitter(radiance=10.0)
    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.SkyState(135.0, 35.0, 200.0, 180.0, 0.5, 0.5)

    meteo = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            latitude=15.0,
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]

    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.check_meteo(meteo)
    JET.@test_call target_modules=(ArchimedLight,) ArchimedLight.summarize_meteo(meteo)
end
