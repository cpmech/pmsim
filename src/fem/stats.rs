pub(crate) struct Stats {
    /// Number of accepted steps
    n_step_accepted: usize,

    /// Number of rejected steps
    n_step_rejected: usize,

    /// Number of reductions on Δλ
    n_ddl_reduction: usize,

    /// Current (step) number of iterations
    n_iteration_current: usize,

    /// Minimum number of iterations per step
    n_iteration_min: usize,

    /// Maximum number of iterations per step
    n_iteration_max: usize,

    /// Total number of iterations
    n_iteration_total: usize,

    /// Current (step) number of iterations with large norm(δu)
    n_large_du_current: usize,

    /// Minimum number of iterations with large norm(δu) per step
    n_large_du_min: usize,

    /// Maximum number of iterations with large norm(δu) per step
    n_large_du_max: usize,

    /// Total number of iterations with large norm(δu)
    n_large_du_total: usize,

    /// Indicates whether to record iterations
    recording: bool,
}

impl Stats {
    pub fn new() -> Self {
        Self {
            n_step_accepted: 0,
            n_step_rejected: 0,
            n_ddl_reduction: 0,
            n_iteration_current: 0,
            n_iteration_min: usize::MAX,
            n_iteration_max: 0,
            n_iteration_total: 0,
            n_large_du_current: 0,
            n_large_du_min: usize::MAX,
            n_large_du_max: 0,
            n_large_du_total: 0,
            recording: false,
        }
    }

    pub fn add_step_accepted(&mut self) {
        self.n_step_accepted += 1;
    }

    pub fn add_step_rejected(&mut self) {
        self.n_step_rejected += 1;
    }

    pub fn add_ddl_reduction(&mut self) {
        self.n_ddl_reduction += 1;
    }

    pub fn start_recording(&mut self) {
        self.n_iteration_current = 0;
        self.n_large_du_current = 0;
        self.recording = true;
    }

    pub fn stop_recording(&mut self) {
        if self.n_iteration_current < self.n_iteration_min {
            self.n_iteration_min = self.n_iteration_current;
        }
        if self.n_iteration_current > self.n_iteration_max {
            self.n_iteration_max = self.n_iteration_current;
        }
        if self.n_large_du_current < self.n_large_du_min {
            self.n_large_du_min = self.n_large_du_current;
        }
        if self.n_large_du_current > self.n_large_du_max {
            self.n_large_du_max = self.n_large_du_current;
        }
        self.recording = false;
    }

    pub fn add_iteration(&mut self, iteration: usize) {
        if self.recording && iteration > 0 {
            self.n_iteration_current += 1;
            self.n_iteration_total += 1;
        }
    }

    pub fn add_large_du(&mut self) {
        if self.recording {
            self.n_large_du_current += 1;
            self.n_large_du_total += 1;
        }
    }

    pub fn n_step_accepted(&self) -> usize {
        self.n_step_accepted
    }

    pub fn n_step_rejected(&self) -> usize {
        self.n_step_rejected
    }

    pub fn n_ddl_reduction(&self) -> usize {
        self.n_ddl_reduction
    }

    pub fn n_iteration(&self) -> String {
        format!(
            "(min: {}, max: {}, total: {})",
            self.n_iteration_min, self.n_iteration_max, self.n_iteration_total
        )
    }

    pub fn n_large_du(&self) -> String {
        if self.n_large_du_total > 0 {
            format!(
                "(min: {}, max: {}, total: {})",
                self.n_large_du_min, self.n_large_du_max, self.n_large_du_total
            )
        } else {
            "None".to_string()
        }
    }
}
