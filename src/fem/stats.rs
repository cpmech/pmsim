pub(crate) struct Stats {
    /// Number of accepted steps
    n_step_accepted: usize,

    /// Number of rejected steps
    n_step_rejected: usize,

    /// Number of iterations per step
    ///
    /// (n_step * n_substep)
    n_iteration_success: Vec<usize>,

    /// Number of failed iterations per step
    ///
    /// (n_step * n_substep)
    n_iteration_fail: Vec<usize>,

    /// Indicates whether to record iterations
    recording: bool,
}

impl Stats {
    pub fn new() -> Self {
        Self {
            n_step_accepted: 0,
            n_step_rejected: 0,
            n_iteration_success: Vec::new(),
            n_iteration_fail: Vec::new(),
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
        self.n_iteration_success.push(0);
        self.n_iteration_fail.push(0);
        self.recording = true;
    }

    pub fn stop_recording(&mut self) {
        self.recording = false;
    }

    pub fn add_iteration_success(&mut self) {
        if self.recording {
            let n = self.n_iteration_success.len();
            self.n_iteration_success[n - 1] += 1;
        }
    }

    pub fn add_iteration_fail(&mut self) {
        if self.recording {
            let n = self.n_iteration_fail.len();
            self.n_iteration_fail[n - 1] += 1;
        }
    }

    pub fn n_accepted_steps(&self) -> usize {
        self.n_step_accepted
    }

    pub fn n_rejected_steps(&self) -> usize {
        self.n_step_rejected
    }

    pub fn n_iteration_succeeded(&self) -> usize {
        self.n_iteration_success.iter().sum()
    }

    pub fn n_iteration_failed(&self) -> usize {
        self.n_iteration_fail.iter().sum()
    }
}
