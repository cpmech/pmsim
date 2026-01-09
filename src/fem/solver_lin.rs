use super::FemData;
use super::{calc_gg_lmm, calc_gg_sps, calc_ggu_lmm, calc_ggu_sps};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_sparse::{CooMatrix, LinSolver, Sym};
use uuid::Uuid;

pub struct SolverLin<'a> {
    data_uuid: Uuid,
    kk: CooMatrix,
    kk_bar: CooMatrix,
    ls: LinSolver<'a>,
    u: Vector,
    mdu: Vector,
    gg: Vector,
}

impl<'a> SolverLin<'a> {
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
    ) -> Result<(Self, FemData<'a>), StrError> {
        let data = FemData::new(&mesh, &schema, &config, &essential, &natural)?;
        let ls = LinSolver::new(data.config.lin_sol_genie)?;
        let u = Vector::new(data.ndim);
        let mdu = Vector::new(data.ndim);
        let gg = Vector::new(data.ndim);
        let (kk, kk_bar) = if data.config.lagrange_mult_method {
            (
                CooMatrix::new(data.ndim, data.ndim, data.nnz_kk, data.sym).unwrap(),
                CooMatrix::new(1, 1, 1, Sym::No).unwrap(),
            )
        } else {
            (
                CooMatrix::new(1, 1, 1, Sym::No).unwrap(),
                CooMatrix::new(data.ndim, data.ndim, data.nnz_kk_bar, data.sym).unwrap(),
            )
        };
        let solver = SolverLin {
            data_uuid: data.uuid,
            kk,
            kk_bar,
            ls,
            u,
            mdu,
            gg,
        };
        Ok((solver, data))
    }

    pub fn steady(&mut self, data: &mut FemData<'a>) -> Result<(), StrError> {
        if data.uuid != self.data_uuid {
            return Err("The solver requires FemData with matching UUID");
        }
        data.state.time = 1.0;
        if data.config.lagrange_mult_method {
            self.kk.reset();
            calc_gg_lmm(&mut self.gg, 1.0, &self.u, data)?;
            calc_ggu_lmm(&mut self.kk, 1.0, &self.u, data)?;
            self.ls
                .actual
                .factorize(&mut self.kk, Some(data.config.lin_sol_params))?;
            self.ls
                .actual
                .solve(&mut self.mdu, &self.gg, data.config.verbose_lin_sys_solve)?;
            for eq in 0..data.neq {
                data.state.u[eq] -= self.mdu[eq];
            }
        } else {
            self.kk_bar.reset();
            calc_gg_sps(&mut self.gg, 1.0, &self.u, data)?;
            calc_ggu_sps(&mut self.kk_bar, 1.0, &self.u, data)?;
            self.ls
                .actual
                .factorize(&mut self.kk_bar, Some(data.config.lin_sol_params))?;
            self.ls
                .actual
                .solve(&mut self.mdu, &self.gg, data.config.verbose_lin_sys_solve)?;
            for eq in 0..data.neq {
                if data.eq_handler.is_unknown(eq) {
                    let iu = data.eq_handler.iu(eq);
                    data.state.u[eq] -= self.mdu[iu];
                }
            }
        }
        Ok(())
    }
}
