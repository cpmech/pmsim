use super::{FemState, FileIo, SolverImplicit};
use crate::StrError;
use russell_lab::vec_add;

impl<'a> SolverImplicit<'a> {
    /// Solves the system of equations
    pub fn solve_steady(&mut self, state: &mut FemState, file_io: &mut FileIo) -> Result<(), StrError> {
        // stages loop
        state.step = 0;
        for stage in 0..self.config.nstage {
            // initialize stage
            self.time.initialize_stage(stage, state);
            self.print.stage(state);

            // time loop
            state.lambda = 1.0;
            while state.step < self.config.n_max_timesteps {
                // done if last timestep
                if self.time.last() {
                    self.print.footer();
                    break;
                }

                // update time-related variables
                self.time.update(state)?;

                // update external forces vector F_ext (also updates the load reversal flag)
                state.reverse = self.com.assemble_ff_ext(state.stage, state.lambda, state.t)?;

                // transient/dynamics: old state variables
                if self.config.transient {
                    vec_add(&mut state.u_star, state.beta1, &state.u, state.beta2, &state.v).unwrap();
                };

                // trial displacement u, displacement increment Δu, and trial loading factor ℓ
                if self.config.arc_length_method {
                    self.arc.trial(state)?;
                } else {
                    // the trial displacement is the displacement at the old time (unchanged)
                    state.ddu.fill(0.0);
                    state.lambda = 1.0;
                }

                // reset algorithmic variables
                if !self.config.linear_problem {
                    self.com.elements.reset_algorithmic_variables(state);
                }

                // print time information
                self.print.step(state);

                // iteration loop
                for iteration in 0..self.config.n_max_iterations {
                    self.iterate(iteration, state)?;
                    if self.res.converged() {
                        self.res.add_converged();
                        break;
                    }
                    if !self.config.arc_length_method {
                        if iteration == self.config.n_max_iterations - 1 {
                            return Err("Newton-Raphson did not converge");
                        }
                    }
                }

                // arc-length step adaptation
                if self.config.arc_length_method {
                    self.arc.adapt(state, self.res.converged(), &self.com.ls.ff_ext)?;
                }

                // check if many iterations failed to converge in a single time step
                self.res.add_failed();
                if self.res.too_many_failures() {
                    return Err("too many iterations failed to converge");
                }

                // perform output
                if self.time.out(state) && self.res.converged() {
                    file_io.write_state(state)?;
                }
                state.step += 1;
            }
        }
        Ok(())
    }
}
