use super::{ElementDiffusion, ElementRod, ElementSolid, ElementTrait, FemBase, FemState};
use crate::base::{assemble_matrix, assemble_vector, Config, Elem};
use crate::StrError;
use gemlab::mesh::{Cell, Mesh};
use russell_lab::{deriv1_central5, Matrix, Vector};
use russell_sparse::CooMatrix;

/// Defines a generic finite element, wrapping an "actual" implementation
pub struct GenericElement<'a> {
    /// Connects to the "actual" implementation of local equations
    pub actual: Box<dyn ElementTrait + 'a>,

    /// Holds the local vector of "internal forces" (including the dynamical forces)
    pub f_int: Vector,

    /// Holds the local vector of "external forces"
    pub f_ext: Vector,

    /// Holds the Ke matrix (local Jacobian matrix; derivative of ϕ w.r.t u)
    pub kke: Matrix,
}

/// Holds a collection of (generic) finite elements
pub struct Elements<'a> {
    /// Holds configuration parameters
    pub config: &'a Config<'a>,

    /// All elements
    ///
    /// (ncell)
    pub all: Vec<GenericElement<'a>>,
}

/// Holds auxiliary arguments for the computation of numerical Jacobian matrices
struct ArgsForNumericalJacobian<'a> {
    /// Holds the local vector of "internal forces"
    pub f_int: &'a mut Vector,

    /// Holds the current state
    pub state: &'a mut FemState,
}

impl<'a> GenericElement<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, base: &'a FemBase, config: &'a Config, cell: &Cell) -> Result<Self, StrError> {
        let element = base.amap.get(cell.attribute).unwrap(); // already checked
        let actual: Box<dyn ElementTrait> = match element {
            Elem::Diffusion(p) => Box::new(ElementDiffusion::new(mesh, base, config, p, cell.id)?),
            Elem::Rod(p) => Box::new(ElementRod::new(mesh, base, config, p, cell.id)?),
            Elem::Beam(..) => panic!("TODO: Beam"),
            Elem::Solid(p) => Box::new(ElementSolid::new(mesh, base, config, p, cell.id)?),
            Elem::PorousLiq(..) => panic!("TODO: PorousLiq"),
            Elem::PorousLiqGas(..) => panic!("TODO: PorousLiqGas"),
            Elem::PorousSldLiq(..) => panic!("TODO: PorousSldLiq"),
            Elem::PorousSldLiqGas(..) => panic!("TODO: PorousSldLiqGas"),
        };
        let neq = base.n_local_eq(cell).unwrap();
        Ok(GenericElement {
            actual,
            f_int: Vector::new(neq),
            f_ext: Vector::new(neq),
            kke: Matrix::new(neq, neq),
        })
    }

    /// Calculates the local Jacobian matrix using finite differences
    ///
    /// **Note:** The state may be changed temporarily, but it is restored at the end of the function
    pub fn numerical_jacobian(&mut self, state: &mut FemState) -> Result<(), StrError> {
        let neq = self.f_int.dim();
        let mut args = ArgsForNumericalJacobian {
            f_int: &mut self.f_int,
            state,
        };
        for i in 0..neq {
            for j in 0..neq {
                let at_u = args.state.uu[j];
                let res = deriv1_central5(at_u, &mut args, |u, a| {
                    let original_uu = a.state.uu[j];
                    let original_duu = a.state.duu[j];
                    a.state.uu[j] = u;
                    a.state.duu[j] = u - original_uu;
                    self.actual.backup_secondary_values(a.state);
                    self.actual.update_secondary_values(&mut a.state).unwrap();
                    self.actual.calc_f_int(&mut a.f_int, &a.state).unwrap();
                    self.actual.restore_secondary_values(&mut a.state);
                    a.state.uu[j] = original_uu;
                    a.state.duu[j] = original_duu;
                    Ok(a.f_int[i])
                });
                self.kke.set(i, j, res.unwrap());
            }
        }
        Ok(())
    }
}

impl<'a> Elements<'a> {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, base: &'a FemBase, config: &'a Config) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = mesh
            .cells
            .iter()
            .map(|cell| GenericElement::new(mesh, base, config, cell))
            .collect();
        match res {
            Ok(all) => Ok(Elements { config, all }),
            Err(e) => Err(e),
        }
    }

    /// Returns whether all local Jacobian matrices are symmetric or not
    pub fn all_symmetric_jacobians(&self) -> bool {
        for e in &self.all {
            if !e.actual.symmetric_jacobian() {
                return false;
            }
        }
        return true;
    }

    /// Calculates all local f_int vectors and assembles them into the global F_int vector
    ///
    /// `ignore` (n_equation) holds the equation numbers to be ignored in the assembly process;
    /// i.e., it allows for skipping the essential prescribed values and generating the reduced system.
    pub fn assemble_f_int(&mut self, ff_int: &mut Vector, state: &FemState, ignore: &[bool]) -> Result<(), StrError> {
        for e in &mut self.all {
            e.actual.calc_f_int(&mut e.f_int, state)?;
            assemble_vector(ff_int, &e.f_int, &e.actual.local_to_global(), ignore);
        }
        Ok(())
    }

    /// Calculates all local f_ext vectors and assembles them into the global F_ext vector
    ///
    /// `ignore` (n_equation) holds the equation numbers to be ignored in the assembly process;
    /// i.e., it allows for skipping the essential prescribed values and generating the reduced system.
    pub fn assemble_f_ext(&mut self, ff_ext: &mut Vector, time: f64, ignore: &[bool]) -> Result<(), StrError> {
        for e in &mut self.all {
            e.actual.calc_f_ext(&mut e.f_ext, time)?;
            assemble_vector(ff_ext, &e.f_ext, &e.actual.local_to_global(), ignore);
        }
        Ok(())
    }

    /// Calculates all local Ke matrices and assembles them into K
    ///
    /// `ignore` (n_equation) holds the equation numbers to be ignored in the assembly process;
    /// i.e., it allows for skipping the essential prescribed values and generating the reduced system.
    pub fn assemble_kke(&mut self, kk: &mut CooMatrix, state: &FemState, ignore: &[bool]) -> Result<(), StrError> {
        let tol = self.config.symmetry_check_tolerance;
        for e in &mut self.all {
            e.actual.calc_jacobian(&mut e.kke, state)?;
            assemble_matrix(kk, &e.kke, &e.actual.local_to_global(), ignore, tol)?;
        }
        Ok(())
    }

    /// Initializes all internal variables
    pub fn initialize_internal_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        self.all
            .iter_mut()
            .map(|e| e.actual.initialize_internal_values(state))
            .collect()
    }

    /// Updates secondary values such as stresses and internal variables
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    pub fn update_secondary_values(&mut self, state: &mut FemState) -> Result<(), StrError> {
        self.all
            .iter_mut()
            .map(|e| e.actual.update_secondary_values(state))
            .collect()
    }

    /// Creates a copy of the secondary values (e.g., stress, int_vars)
    pub fn backup_secondary_values(&mut self, state: &FemState) {
        self.all
            .iter_mut()
            .map(|e| e.actual.backup_secondary_values(state))
            .collect()
    }

    /// Restores the secondary values (e.g., stress, int_vars) from the backup
    pub fn restore_secondary_values(&self, state: &mut FemState) {
        self.all
            .iter()
            .map(|e| e.actual.restore_secondary_values(state))
            .collect()
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    pub fn reset_algorithmic_variables(&self, state: &mut FemState) {
        self.all
            .iter()
            .map(|e| e.actual.reset_algorithmic_variables(state))
            .collect()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Elements, GenericElement};
    use crate::base::{Conductivity, Config, Elem, Essential, ParamBeam, ParamPorousLiqGas, StressStrain};
    use crate::base::{ParamDiffusion, ParamPorousLiq, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
    use crate::fem::{FemBase, FemState};
    use gemlab::integ;
    use gemlab::mesh::{Mesh, Samples};
    use russell_lab::{mat_approx_eq, vec_add, vec_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let config = Config::new(&mesh);

        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        assert_eq!(
            GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).err(),
            Some("requested number of integration points is not available for Tri class")
        );
        assert_eq!(
            Elements::new(&mesh, &base, &config).err(),
            Some("requested number of integration points is not available for Tri class")
        );

        let mut p1 = ParamDiffusion::sample();
        p1.ngauss = Some(123); // wrong
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        assert_eq!(
            GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).err(),
            Some("requested number of integration points is not available for Tri class")
        );
        assert_eq!(
            Elements::new(&mesh, &base, &config).err(),
            Some("requested number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        let elements = Elements::new(&mesh, &base, &config).unwrap();
        assert_eq!(elements.all.len(), mesh.cells.len());
    }

    #[test]
    fn num_jacobian_diffusion() {
        // mesh
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut ele = GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x, y| 100.0 + 7.0 * x + 3.0 * y;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[1]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[1]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[1]);

        // check
        ele.actual.calc_jacobian(&mut ele.kke, &state).unwrap();
        let jj_ana = ele.kke.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        mat_approx_eq(&jj_ana, &ele.kke, 1e-11);
        // println!("Ke = \n{}", ele.kke);
    }

    #[test]
    fn num_jacobian_diffusion_transient() {
        // mesh
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let mut config = Config::new(&mesh);
        config.set_transient(true);
        let mut ele = GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x, y| 100.0 + 7.0 * x + 3.0 * y;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[1]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[1]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[1]);

        // check
        ele.actual.calc_jacobian(&mut ele.kke, &state).unwrap();
        let jj_ana = ele.kke.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        mat_approx_eq(&jj_ana, &ele.kke, 1e-11);
    }

    #[test]
    fn num_jacobian_diffusion_variable_conductivity() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion {
            rho: 1.0,
            conductivity: Conductivity::IsotropicLinear { kr: 2.0, beta: 10.0 },
            source: None,
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut ele = GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let tt_field = |x, y| 100.0 + 7.0 * x + 3.0 * y;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[1]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[1]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[1]);

        // check
        ele.actual.calc_jacobian(&mut ele.kke, &state).unwrap();
        let jj_ana = ele.kke.clone();
        ele.numerical_jacobian(&mut state).unwrap();
        // println!("ana: J = \n{}", jj_ana);
        // println!("num: J = \n{}", ele.kke);
        mat_approx_eq(&jj_ana, &ele.kke, 1e-7);
        // note that the "stiffness" is now unsymmetric
        // println!("difference = {:?}", jj_ana.get(0, 2) - jj_ana.get(2, 0));
    }

    #[test]
    fn num_jacobian_solid() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut ele = GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();

        // linear displacement field
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        state.duu[0] = 1.0 + mesh.points[0].coords[0];
        state.duu[1] = 2.0 + mesh.points[0].coords[1];
        state.duu[2] = 1.0 + mesh.points[1].coords[0];
        state.duu[3] = 2.0 + mesh.points[1].coords[1];
        state.duu[4] = 1.0 + mesh.points[2].coords[0];
        state.duu[5] = 2.0 + mesh.points[2].coords[1];
        for i in 0..6 {
            state.uu[i] = state.duu[i];
        }
        ele.actual.update_secondary_values(&mut state).unwrap();
        println!("uu =\n{}", state.uu);

        ele.actual.calc_jacobian(&mut ele.kke, &state).unwrap();
        let jj_ana = ele.kke.clone();
        ele.numerical_jacobian(&mut state).unwrap();

        println!("J(ana)=\n{:.2}", jj_ana);
        println!("J(num)=\n{:.2}", ele.kke);

        mat_approx_eq(&jj_ana, &ele.kke, 1e-8);
    }

    // ----------------- temporary ----------------------------------------

    #[test]
    #[should_panic(expected = "TODO: Beam")]
    fn new_panics_beam() {
        let mesh = Samples::one_lin2();
        let p1 = ParamBeam::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Beam(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiq")]
    fn new_panics_porous_liq() {
        let mesh = Samples::one_tri3();
        let p1 = ParamPorousLiq::sample_brooks_corey_constant();
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiq(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousLiqGas")]
    fn new_panics_porous_liq_gas() {
        let mesh = Samples::one_tri3();
        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiqGas(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiq")]
    fn new_panics_porous_sld_liq() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::PorousSldLiq(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    #[should_panic(expected = "TODO: PorousSldLiqGas")]
    fn new_panics_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new(&mesh);
        GenericElement::new(&mesh, &base, &config, &mesh.cells[0]).unwrap();
    }

    #[test]
    fn calc_and_add_methods_work() {
        // mesh and parameters
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let young = 1500.0;
        let poisson = 0.25;
        let p1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic { young, poisson },
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // calculate solution (c vectors = contributions to R) and set state
        let neq = mesh.points.len() * 2; // 2 DOF per node
        let mut rr_correct = Vector::new(neq);
        let mut kk_correct = Matrix::new(neq, neq);
        let l2g = [[0, 1, 2, 3, 8, 9], [2, 3, 6, 7, 8, 9], [2, 3, 4, 5, 6, 7]];
        for e in 0..3 {
            let delta = (e + 1) as f64;
            let sigma = Tensor2::from_matrix(
                &[
                    [1.0 + delta, 3.0 + delta, 0.0],
                    [3.0 + delta, 2.0 + delta, 0.0],
                    [0.0, 0.0, delta],
                ],
                Mandel::Symmetric2D,
            )
            .unwrap();
            for state in &mut state.gauss[e].solid {
                state.stress.set_tensor(1.0, &sigma);
            }
            let mut pad = Mesh::get_pad(&mesh, e);
            let ana = integ::AnalyticalTri3::new(&mut pad);
            let c = ana.vec_04_tb(&sigma, false);
            for l in 0..6 {
                rr_correct[l2g[e][l]] += c[l];
            }
            let kk = ana.mat_10_bdb(young, poisson, false, 1.0).unwrap();
            for l in 0..6 {
                for m in 0..6 {
                    kk_correct.add(l2g[e][l], l2g[e][m], kk.get(l, m));
                }
            }
        }

        // elements
        let mut elements = Elements::new(&mesh, &base, &config).unwrap();
        let neq = base.dofs.size();
        let nnz_sup = 3 * neq * neq;
        let mut f_int = Vector::new(neq);
        let mut f_ext = Vector::new(neq);
        let mut rr = Vector::new(neq);
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::No).unwrap();
        let ignore = vec![false; neq];
        elements.assemble_f_int(&mut f_int, &state, &ignore).unwrap();
        elements.assemble_f_ext(&mut f_ext, state.t, &ignore).unwrap();
        elements.assemble_kke(&mut kk, &state, &ignore).unwrap();
        vec_add(&mut rr, 1.0, &f_int, -1.0, &f_ext).unwrap();
        let kk_mat = kk.as_dense();
        vec_approx_eq(&rr, &rr_correct, 1e-14);
        mat_approx_eq(&kk_mat, &kk_correct, 1e-12);
    }
}
