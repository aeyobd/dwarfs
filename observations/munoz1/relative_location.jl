using LilGuys
using TOML


props_munoz1 = TOML.parsefile("observed_properties.toml")
props_umi = TOML.parsefile("../ursa_minor/observed_properties.toml")

xi, eta = 60 .* LilGuys.to_tangent(props_munoz1["ra"], props_munoz1["dec"], props_umi["ra"], props_umi["dec"])

@info "units = arcmin"
@info "xi, eta = $xi, $eta"

R_ell = LilGuys.calc_R_ell(xi, eta, props_umi["ellipticity"], props_umi["position_angle"])

@info "R_ell = $R_ell"
