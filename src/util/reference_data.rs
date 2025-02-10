pub(crate) trait ReferenceDataTrait {
    /// Returns the number of load increments or timesteps
    fn nstep(&self) -> usize;

    /// Returns the number of points
    fn npoint(&self) -> usize;

    /// Returns the number of cells/elements
    fn ncell(&self) -> usize;

    /// Returns the displacement component of point p, dimension i
    ///
    /// index is the index of the load increment or timestep
    fn displacement(&self, step: usize, p: usize, i: usize) -> f64;

    /// Returns the number of Gauss points of element/cell e
    ///
    /// step is the index of the load increment or timestep
    /// Returns the number of Gauss points of element/cell e
    fn ngauss(&self, step: usize, e: usize) -> usize;

    /// Returns the stress component of element/cell e, gauss point ip, component i
    ///
    /// step is the index of the load increment or timestep
    /// Returns the stress component of element/cell e, gauss point ip, component i
    fn stresses(&self, step: usize, e: usize, ip: usize, i: usize) -> f64;
}
