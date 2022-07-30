use pmsim::base::{ParamFluids, SampleMeshes};

fn main() {
    // 0. mesh and liquid parameters
    let mesh = SampleMeshes::column();
    let liq = ParamFluids::new();

    // 1. allocate instance to calculate the liquid pressure along the vertical line
    let geo = Geostatic::new(&mesh, &liq);

    // 2. for each point of the mesh, compute the liquid pressure
    let all_pl = geo.calc_points_pl(&mesh, &liq);

    // 3. for each Gauss point, compute the effective stress tensor
    let all_sigma = geo.calc_stress(&mesh, &liq, &porous, &pads, &ips);

    // 4. define boundary conditions by varying the liquid pressure of points at boundaries
    let pl_ebcs = Vec::new();
    let t_fin = 100.0;
    for point in boundary {
        let pl_0 = geo.pl(point);
        pl_ebcs.push((Pl, |t| {
            let tt = f64::max(0.0, f64::min(t_fin, t));
            pl_0 * (1.0 - tt / t_fin)
        }));
    }
}
