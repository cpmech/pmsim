use super::{Data, LocalEquations, State};
use crate::base::{Config, ParamSolid};
use crate::model::{allocate_stress_strain_model, StressStrain};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::copy_tensor2;

/// Implements the local Solid Element equations
pub struct ElementSolid<'a> {
    /// Number of space dimensions
    pub ndim: usize,

    /// Local-to-global mapping
    pub local_to_global: &'a Vec<usize>,

    /// Global configuration
    pub config: &'a Config,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamSolid,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub ips: integ::IntegPointData,

    /// Stress-strain model
    pub model: Box<dyn StressStrain>,
}

impl<'a> ElementSolid<'a> {
    /// Allocates new instance
    pub fn new(data: &'a Data, config: &'a Config, cell: &'a Cell, param: &'a ParamSolid) -> Result<Self, StrError> {
        // scratchpad for numerical integration
        let ndim = data.mesh.ndim;
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, data.mesh);

        Ok({
            ElementSolid {
                ndim,
                local_to_global: &data.dof_numbers.local_to_global[cell.id],
                config,
                cell,
                param,
                pad,
                ips: config.integ_point_data(cell)?,
                model: allocate_stress_strain_model(param, ndim == 2, config.plane_stress),
            }
        })
    }
}

impl<'a> LocalEquations for ElementSolid<'a> {
    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &State) -> Result<(), StrError> {
        let sigma = &state.effective_stress[self.cell.id];
        integ::vec_04_tg(residual, &mut self.pad, 0, true, self.ips, |sig, p, _| {
            copy_tensor2(sig, &sigma[p])
        })
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &State) -> Result<(), StrError> {
        let sigma = &state.effective_stress[self.cell.id];
        integ::mat_10_gdg(jacobian, &mut self.pad, 0, 0, true, self.ips, |dd, p, _| {
            self.model.stiffness(dd, &sigma[p])
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementSolid;
    use crate::base::{Config, Element, ParamSolid, ParamStressStrain, SampleParams};
    use crate::fem::{Data, LocalEquations, State};
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementSolid::new(&data, &config, &mesh.cells[0], &p1).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn element_solid_works_2d() {
        // mesh and parameters
        let mesh = Samples::one_tri3();
        let young = 10_000.0; // kPa
        let poisson = 0.2; // [-]
        let p1 = ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic { young, poisson },
        };
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let mut elem = ElementSolid::new(&data, &config, &mesh.cells[0], &p1).unwrap();

        // set stress state
        let (s00, s11, s01) = (1.0, 2.0, 3.0);
        let mut state = State::new(&data, &config).unwrap();
        for sigma in &mut state.effective_stress[0] {
            sigma.sym_set(0, 0, s00);
            sigma.sym_set(1, 1, s11);
            sigma.sym_set(0, 1, s01);
        }

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
        let neq = 3 * 2;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        let correct = ana.vec_04_tg(s00, s11, s01);
        vec_approx_eq(residual.as_data(), &correct, 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct = ana
            .mat_10_gdg(young, poisson, config.plane_stress, config.thickness)
            .unwrap();
        vec_approx_eq(jacobian.as_data(), correct.as_data(), 1e-12);
    }
}
