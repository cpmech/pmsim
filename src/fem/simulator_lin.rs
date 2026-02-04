use super::FemData;
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::fem::callbacks::calc_jac_sps;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_sparse::{CooMatrix, LinSolver, Sym};
use uuid::Uuid;

/// Performs linear finite element simulations
pub struct SimulatorLin<'a> {
    data_uuid: Uuid,
    mm: CooMatrix,
    kk_bar: CooMatrix,
    ls: LinSolver<'a>,
    u: Vector,
    ddu: Vector,
    rhs: Vector,
}

impl<'a> SimulatorLin<'a> {
    /// Allocates a new instance
    ///
    /// Typical usage:
    ///
    /// ```text
    /// let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    /// ```
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
    ) -> Result<(Self, FemData<'a>), StrError> {
        let data = FemData::new(&mesh, &schema, &config, &essential, &natural)?;
        let ls = LinSolver::new(data.config.lin_sol_genie)?;
        let nsys = data.nsys; // system dimension
        let u = Vector::new(nsys);
        let ddu = Vector::new(nsys);
        let rhs = Vector::new(nsys);
        let (mm, kk_bar) = if data.config.lagrange_mult_method {
            (
                CooMatrix::new(nsys, nsys, data.nnz_mm, data.sym).unwrap(),
                CooMatrix::new(1, 1, 1, Sym::No).unwrap(),
            )
        } else {
            (
                CooMatrix::new(1, 1, 1, Sym::No).unwrap(),
                CooMatrix::new(nsys, nsys, data.nnz_kk_bar, data.sym).unwrap(),
            )
        };
        let solver = SimulatorLin {
            data_uuid: data.uuid,
            mm,
            kk_bar,
            ls,
            u,
            ddu,
            rhs,
        };
        Ok((solver, data))
    }

    // Runs a steady-state simulation
    ///
    /// Typical usage:
    ///
    /// ```text
    /// let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    /// sim.steady(&mut data, true)?;
    //  let state = data.get_state();
    /// ```
    pub fn steady(&mut self, data: &mut FemData<'a>, post_compute_second_values: bool) -> Result<(), StrError> {
        // Check UUID
        if data.uuid != self.data_uuid {
            return Err("The solver requires FemData with matching UUID");
        }

        // Print information about the system and the header
        if data.config.verbose {
            data.print_system_info("N/A");
        }

        // Update loading factor
        let l0 = data.state.lambda;
        let l1 = l0 + 1.0;
        data.state.lambda = l1;

        // Calculate Pᵤ(t), lambda-free part of the prescribed values Ǔ = λ Pᵤ(t)
        data.calc_ppu();

        // Calculate F(t), external forces
        data.calc_ff()?;

        // Calculate Y: internal forces
        data.calc_yy()?;

        // Check the symmetry of the local stiffness matrices
        if let Some(tol) = data.config.enable_symmetry_check {
            data.elements.check_symmetry_kke(&data.state, tol)?;
            data.boundaries.check_symmetry_kke(&data.state, tol)?;
        }

        // Lagrange multipliers method
        let mut empty = Vector::new(0);
        if data.config.lagrange_mult_method {
            // Assemble the right-hand side vector:
            //       ┌       ┐
            //       │ λ₁ F  │  (neq)
            // RHS = │       │
            //       │ λ₁ Pᵤ │  (np)
            //       └       ┘
            for i in 0..data.ndof {
                self.rhs[i] = l1 * data.ff[i];
            }
            for ip in 0..data.np {
                let j = data.ndof + ip;
                self.rhs[j] = l1 * data.ppu[ip];
            }

            // Calculate the stiffness matrix K and assemble it into M
            self.mm.reset();
            data.elements.assemble_kk_lmm(&mut self.mm, &mut data.state)?;
            data.boundaries.assemble_kk_lmm(&mut self.mm, &mut data.state)?;

            // Add constraint matrix to M
            //     ┌         ┐
            //     │  K   Cᵀ │
            // M = │         │
            //     │  C   0  │
            //     └         ┘
            let sym = self.mm.get_info().3;
            match sym {
                Sym::YesLower => {
                    for ip in 0..data.np {
                        let i = data.eq_handler.prescribed()[ip];
                        let j = data.ndof + ip;
                        self.mm.put(j, i, 1.0).unwrap(); // C
                    }
                }
                Sym::YesUpper => {
                    for ip in 0..data.np {
                        let i = data.eq_handler.prescribed()[ip];
                        let j = data.ndof + ip;
                        self.mm.put(i, j, 1.0).unwrap(); // Cᵀ
                    }
                }
                Sym::YesFull | Sym::No => {
                    for ip in 0..data.np {
                        let i = data.eq_handler.prescribed()[ip];
                        let j = data.ndof + ip;
                        self.mm.put(i, j, 1.0).unwrap(); // Cᵀ
                        self.mm.put(j, i, 1.0).unwrap(); // C
                    }
                }
            }

            // Factorize M
            self.ls
                .actual
                .factorize(&mut self.mm, Some(data.config.lin_sol_params))?;

            // Solve the linear system: ΔU = M⁻¹ RHS
            self.ls
                .actual
                .solve(&mut self.ddu, &self.rhs, data.config.verbose_lin_sys_solve)?;

            // Update U and ΔU in the state
            for eq in 0..data.ndof {
                data.state.uu[eq] += self.ddu[eq];
                data.state.dduu[eq] = self.ddu[eq];
            }
        }
        // System partitioning method
        else {
            // Initialize the right-hand side vector: RHS := λ₁ F̄
            for iu in 0..data.nu {
                let eq = data.eq_handler.unknown()[iu];
                self.rhs[iu] = l1 * data.ff[eq];
            }

            // Calculate the stiffness matrix
            self.kk_bar.reset();
            data.kk_check.reset();
            calc_jac_sps(&mut self.kk_bar, &mut empty, l1, &self.u, data)?;

            // Fix the right-hand side: RHS := RHS - Ǩ Ǔ = RHS - λ Ǩ Cᵤ
            data.kk_check.mat_vec_mul_update(&mut self.rhs, -l1, &data.ppu)?;

            // Factorize the stiffness matrix
            self.ls
                .actual
                .factorize(&mut self.kk_bar, Some(data.config.lin_sol_params))?;

            // Solve the linear system for ΔU
            self.ls
                .actual
                .solve(&mut self.ddu, &self.rhs, data.config.verbose_lin_sys_solve)?;

            // Update the solution: U₁ = U₀ + ΔU
            for eq in 0..data.ndof {
                if data.eq_handler.is_unknown(eq) {
                    let iu = data.eq_handler.iu(eq);
                    data.state.uu[eq] += self.ddu[iu];
                    data.state.dduu[eq] = self.ddu[iu];
                }
            }

            // Update prescribed values Ǔ
            data.eq_handler.prescribed().iter().for_each(|&eq| {
                let ip = data.eq_handler.ip(eq);
                data.state.uu[eq] = l1 * data.ppu[ip];
                data.state.dduu[eq] = (l1 - l0) * data.ppu[ip];
            });
        }

        // Compute secondary state variables such as stresses
        if post_compute_second_values {
            data.elements.update_secondary_values(&mut data.state)?;
        }

        // Last output
        data.files.execute(&data.schema, &data.config, &data.state, &data.yy)?;
        data.files.stop(&data.config)?;
        Ok(())
    }
}
