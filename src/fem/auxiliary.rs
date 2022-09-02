use super::Data;
use gemlab::mesh::Cell;
use russell_lab::Vector;

pub fn interpolate(_u_local: &mut Vector, _u_global: &Vector, _data: &Data, _cell: &Cell) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::interpolate;
    use crate::base::{Element, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn interpolate_works() {
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let u_global = Vector::new(data.dof_numbers.n_equation);

        let mut u_local = Vector::new(3);
        interpolate(&mut u_local, &u_global, &data, &data.mesh.cells[0]);
    }
}
