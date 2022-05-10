/// Holds an option to initialize stresses
pub enum IniOption {
    /// Geostatic initial state with data = (overburden,total_stress)
    ///
    /// # Note
    ///
    /// * The argument is the overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial isotropic stress state with σ_xx = σ_yy = σ_zz = value
    Isotropic(f64),

    /// Zero initial state
    Zero,
}
