use super::{ElementEquations, State};
use crate::base::{Config, DofNumbers, ElementDofs, ParamSolid};
use crate::model::{allocate_stress_strain_model, StressStrain};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::copy_tensor2;

pub struct ElementSolid<'a> {
    pub info: &'a ElementDofs,
    pub config: &'a Config,
    pub cell: &'a Cell,
    pub param: &'a ParamSolid,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Matrix,
    pub model: Box<dyn StressStrain>,
}

impl<'a> ElementSolid<'a> {
    pub fn new(
        mesh: &'a Mesh,
        dn: &'a DofNumbers,
        config: &'a Config,
        cell: &'a Cell,
        param: &'a ParamSolid,
    ) -> Result<Self, StrError> {
        // extract element info
        let info = dn
            .element_dofs
            .get(&(cell.attribute_id, cell.kind))
            .ok_or("cannot extract CellAttributeId to allocate ElementSolid")?;

        // pad and ips
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(mesh.ndim, kind)?;
        set_pad_coords(&mut pad, &points, &mesh);

        // model
        let two_dim = mesh.ndim == 2;
        let plane_stress = config.plane_stress;
        let model = allocate_stress_strain_model(param, two_dim, plane_stress);

        // done
        Ok({
            ElementSolid {
                info,
                config,
                cell,
                param,
                pad,
                ips: config.integ_point_data(cell)?,
                residual: Vector::new(info.n_equation_local),
                jacobian: Matrix::new(info.n_equation_local, info.n_equation_local),
                model,
            }
        })
    }
}

impl<'a> ElementEquations for ElementSolid<'a> {
    fn residual(&mut self, state: &State) -> Result<(), StrError> {
        let sigma = &state.effective_stress[self.cell.id];
        integ::vec_04_tg(&mut self.residual, &mut self.pad, 0, true, self.ips, |sig, p| {
            copy_tensor2(sig, &sigma[p])
        })
    }
    fn jacobian(&mut self, state: &State) -> Result<(), StrError> {
        let sigma = &state.effective_stress[self.cell.id];
        integ::mat_10_gdg(&mut self.jacobian, &mut self.pad, 0, 0, true, self.ips, |dd, p| {
            self.model.stiffness(dd, &sigma[p])
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementSolid;
    use crate::base::{Config, DofNumbers, Element, ParamSolid, ParamStressStrain, SampleParams};
    use crate::fem::element_equations::ElementEquations;
    use crate::fem::State;
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_chk::assert_vec_approx_eq;
    use std::collections::HashMap;

    #[test]
    fn new_handles_errors() {
        let mut mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, &elements).unwrap();
        let mut config = Config::new();
        mesh.cells[0].attribute_id = 100; // << never do this!
        assert_eq!(
            ElementSolid::new(&mesh, &dn, &config, &mesh.cells[0], &p1).err(),
            Some("cannot extract CellAttributeId to allocate ElementSolid")
        );
        mesh.cells[0].attribute_id = 1;
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementSolid::new(&mesh, &dn, &config, &mesh.cells[0], &p1).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn element_solid_works() {
        let mesh = Samples::one_tri3();
        let young = 10_000.0; // kPa
        let poisson = 0.2; // [-]
        let p1 = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
        };
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, &elements).unwrap();
        let config = Config::new();
        let mut elem = ElementSolid::new(&mesh, &dn, &config, &mesh.cells[0], &p1).unwrap();

        // check residual vector
        let (s00, s11, s01) = (1.0, 2.0, 3.0);
        let mut state = State::new(&mesh, &elements, &dn, &config).unwrap();
        for sigma in &mut state.effective_stress[0] {
            sigma.sym_set(0, 0, s00);
            sigma.sym_set(1, 1, s11);
            sigma.sym_set(0, 1, s01);
        }
        elem.residual(&state).unwrap();
        let ana = integ::AnalyticalTri3::new(&elem.pad);
        let correct = ana.vec_04_tg(s00, s11, s01);
        assert_vec_approx_eq!(elem.residual.as_data(), correct, 1e-15);

        // check Jacobian matrix
        elem.jacobian(&state).unwrap();
        let correct = ana
            .mat_10_gdg(young, poisson, config.plane_stress, config.thickness)
            .unwrap();
        assert_vec_approx_eq!(elem.jacobian.as_data(), correct.as_data(), 1e-12);
    }
}
