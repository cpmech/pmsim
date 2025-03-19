pub(crate) struct Stats {
    /// Number of accepted steps
    n_step_accepted: usize,

    /// Number of rejected steps
    n_step_rejected: usize,

    /// Current (step) number of iterations
    n_iteration_current: usize,

    /// Minimum number of iterations per step
    n_iteration_min: usize,

    /// Maximum number of iterations per step
    n_iteration_max: usize,

    /// Total number of iterations
    n_iteration_total: usize,

    /// Indicates if there are failed iterations
    has_failed_iteration: bool,

    /// Current (step) number of failed iterations
    n_failed_iteration_current: usize,

    /// Minimum number of failed iterations per step
    n_failed_iteration_min: usize,

    /// Maximum number of failed iterations per step
    n_failed_iteration_max: usize,

    /// Total number of failed iterations
    n_failed_iteration_total: usize,

    /// Indicates whether to record iterations
    recording: bool,
}

impl Stats {
    pub fn new() -> Self {
        Self {
            n_step_accepted: 0,
            n_step_rejected: 0,
            n_iteration_current: 0,
            n_iteration_min: usize::MAX,
            n_iteration_max: 0,
            n_iteration_total: 0,
            has_failed_iteration: false,
            n_failed_iteration_current: 0,
            n_failed_iteration_min: usize::MAX,
            n_failed_iteration_max: 0,
            n_failed_iteration_total: 0,
            recording: false,
        }
    }

    pub fn add_step_accepted(&mut self) {
        self.n_step_accepted += 1;
    }

    pub fn add_step_rejected(&mut self) {
        self.n_step_rejected += 1;
    }

    pub fn start_recording(&mut self) {
        self.n_iteration_current = 0;
        self.n_failed_iteration_current = 0;
        self.recording = true;
    }

    pub fn stop_recording(&mut self) {
        if self.n_iteration_current < self.n_iteration_min {
            self.n_iteration_min = self.n_iteration_current;
        }
        if self.n_iteration_current > self.n_iteration_max {
            self.n_iteration_max = self.n_iteration_current;
        }
        if self.has_failed_iteration {
            if self.n_failed_iteration_current < self.n_failed_iteration_min {
                self.n_failed_iteration_min = self.n_failed_iteration_current;
            }
            if self.n_failed_iteration_current > self.n_failed_iteration_max {
                self.n_failed_iteration_max = self.n_failed_iteration_current;
            }
        }
        self.recording = false;
    }

    pub fn add_iteration(&mut self, iteration: usize) {
        if self.recording && iteration > 0 {
            self.n_iteration_current += 1;
            self.n_iteration_total += 1;
        }
    }

    pub fn add_iteration_fail(&mut self) {
        if self.recording {
            self.has_failed_iteration = true;
            self.n_failed_iteration_current += 1;
            self.n_failed_iteration_total += 1;
        }
    }

    pub fn n_accepted_steps(&self) -> usize {
        self.n_step_accepted
    }

    pub fn n_rejected_steps(&self) -> usize {
        self.n_step_rejected
    }

    pub fn n_iteration(&self) -> String {
        format!(
            "(min: {}, max: {}, total: {})",
            self.n_iteration_min, self.n_iteration_max, self.n_iteration_total
        )
    }

    pub fn n_iteration_failed(&self) -> String {
        if self.has_failed_iteration {
            format!(
                "(min: {}, max: {}, total: {})",
                self.n_failed_iteration_min, self.n_failed_iteration_max, self.n_failed_iteration_total
            )
        } else {
            "None".to_string()
        }
    }
}
