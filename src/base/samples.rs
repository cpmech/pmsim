use super::{
    ParamBeam, ParamConductivity, ParamFluids, ParamLiquidRetention, ParamPorous, ParamRealDensity, ParamRod,
    ParamSeepage, ParamSolid, ParamStressStrain,
};
use gemlab::mesh::Mesh;

/// Holds some sample material parameters and meshes
pub struct Samples;

impl Samples {
    /// Returns sample parameters for the density of water (SI units)
    pub fn param_density_water(incompressible: bool) -> ParamRealDensity {
        let cc = if incompressible { 1e-12 } else { 4.53e-7 }; // Mg/(m³ kPa)
        ParamRealDensity {
            cc,           // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        }
    }

    /// Returns sample parameters for the density of dry air (SI units)
    pub fn param_density_dry_air() -> ParamRealDensity {
        ParamRealDensity {
            cc: 1.17e-5,     // Mg/(m³ kPa)
            p_ref: 0.0,      // kPa
            rho_ref: 0.0012, // Mg/m³
            tt_ref: 25.0,    // ℃
        }
    }

    /// Returns sample parameters for water (SI units)
    pub fn param_water(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: Samples::param_density_water(incompressible),
            density_gas: None,
        }
    }

    /// Returns sample parameters for water and dry air (SI units)
    pub fn param_water_and_dry_air(incompressible: bool) -> ParamFluids {
        ParamFluids {
            density_liquid: Samples::param_density_water(incompressible),
            density_gas: Some(Samples::param_density_dry_air()),
        }
    }

    /// Returns sample parameters for a linear-elastic rod
    pub fn param_rod() -> ParamRod {
        ParamRod {
            density: 2.0,
            young: 1000.0,
            area: 1.0,
        }
    }

    /// Returns sample parameters for an Euler-Bernoulli beam
    pub fn param_beam() -> ParamBeam {
        ParamBeam {
            density: 2.0,
            young: 1000.0,
            shear: 2000.0,
            area: 1.0,
            ii_11: 1.0,
            ii_22: 1.0,
            jj_tt: 1.0,
        }
    }

    /// Returns sample parameters for a solid medium
    pub fn param_solid() -> ParamSolid {
        ParamSolid {
            density: 2.7, // Mg/m²
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: 0.2,    // [-]
            },
            n_integ_point: None,
        }
    }

    /// Returns sample parameters for a porous medium with solid, liquid and gas
    pub fn param_porous_sol_liq_gas(porosity_initial: f64, k_iso: f64) -> ParamPorous {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorous {
            earth_pres_coef_ini: kk0,
            porosity_initial,
            density_solid: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: nu,     // [-]
            },
            retention_liquid: ParamLiquidRetention::PedrosoWilliams {
                with_hysteresis: true,
                lambda_d: 3.0,
                lambda_w: 3.0,
                beta_d: 6.0,
                beta_w: 6.0,
                beta_1: 6.0,
                beta_2: 6.0,
                x_rd: 2.0,
                x_rw: 2.0,
                y_0: 0.95,
                y_r: 0.005,
            },
            conductivity_liquid: ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            conductivity_gas: Some(ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 2.0,
                lambda_1: 0.001,
                alpha: 0.01,
                beta: 10.0,
            }),
            n_integ_point: None,
        }
    }

    /// Returns sample parameters for a porous medium with solid and liquid
    pub fn param_porous_sol_liq(porosity_initial: f64, k_iso: f64) -> ParamPorous {
        let nu = 0.2;
        let kk0 = nu / (1.0 - nu);
        ParamPorous {
            earth_pres_coef_ini: kk0,
            porosity_initial,
            density_solid: 2.7, // Mg/m³
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0, // kPa
                poisson: nu,     // [-]
            },
            retention_liquid: ParamLiquidRetention::PedrosoWilliams {
                with_hysteresis: true,
                lambda_d: 3.0,
                lambda_w: 3.0,
                beta_d: 6.0,
                beta_w: 6.0,
                beta_1: 6.0,
                beta_2: 6.0,
                x_rd: 2.0,
                x_rw: 2.0,
                y_0: 1.0,
                y_r: 0.005,
            },
            conductivity_liquid: ParamConductivity::PedrosoZhangEhlers {
                kx: k_iso, // m/s
                ky: k_iso, // m/s
                kz: k_iso, // m/s
                lambda_0: 0.001,
                lambda_1: 1.2,
                alpha: 0.01,
                beta: 10.0,
            },
            conductivity_gas: None,
            n_integ_point: None,
        }
    }

    /// Returns sample parameters for seepage models with liquid only
    pub fn param_seepage_liq() -> ParamSeepage {
        ParamSeepage {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::BrooksCorey {
                lambda: 0.1,
                pc_ae: 0.1,
                sl_min: 0.1,
                sl_max: 1.0,
            },
            conductivity_liquid: ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            },
            conductivity_gas: None,
            n_integ_point: None,
        }
    }

    /// Returns sample parameters for seepage models with liquid and gas
    pub fn param_seepage_liq_gas() -> ParamSeepage {
        ParamSeepage {
            porosity_initial: 0.4,
            retention_liquid: ParamLiquidRetention::BrooksCorey {
                lambda: 0.1,
                pc_ae: 0.1,
                sl_min: 0.1,
                sl_max: 1.0,
            },
            conductivity_liquid: ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            },
            conductivity_gas: Some(ParamConductivity::Constant {
                kx: 0.1,
                ky: 0.1,
                kz: 0.1,
            }),
            n_integ_point: None,
        }
    }

    /// Returns a sample mesh with a single segment
    ///
    /// ```text
    ///    1
    ///   /
    ///  /
    /// 0
    /// ```
    pub fn mesh_segment() -> Mesh {
        Mesh::from_text(
            r"
            # header
            # space_ndim npoint ncell
                       2      2     1

            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 1.0

            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        1     2  0 1",
        )
        .unwrap()
    }

    /// Returns a sample mesh with a single square
    ///
    /// ```text
    /// 3------2
    /// |      |
    /// |      |
    /// 0------1
    /// ```
    pub fn mesh_square() -> Mesh {
        Mesh::from_text(
            r"
            # header
            # space_ndim npoint ncell
                       2      4     1

            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0

            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3",
        )
        .unwrap()
    }

    /// Returns a sample mesh with two squares
    ///
    /// ```text
    /// 3--------2--------5
    /// |        |        |
    /// |        |        |
    /// |        |        |
    /// 0--------1--------4
    /// ```
    pub fn mesh_two_quads() -> Mesh {
        Mesh::from_text(
            r"
            # space_ndim npoint ncell
                       2      6     2

            # points
            # id    x   y
               0  0.0 0.0
               1  1.0 0.0
               2  1.0 1.0
               3  0.0 1.0
               4  2.0 0.0
               5  2.0 1.0

            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3
               1   2        2     4  1 4 5 2",
        )
        .unwrap()
    }

    /// Returns a sample mesh with a single cube
    ///
    /// ```text
    ///     4-----------7
    ///    /.          /|
    ///   / .         / |
    ///  5-----------6  |
    ///  |  .        |  |
    ///  |  0--------|--3
    ///  | /         | /
    ///  |/          |/
    ///  1-----------2
    /// ```
    pub fn mesh_cube() -> Mesh {
        Mesh::from_text(
            r"
            # header
            # space_ndim npoint ncell
                       3      8     1

            # points
            # id    x   y   z
               0  0.0 0.0 0.0
               1  1.0 0.0 0.0
               2  1.0 1.0 0.0
               3  0.0 1.0 0.0
               4  0.0 0.0 1.0
               5  1.0 0.0 1.0
               6  1.0 1.0 1.0
               7  0.0 1.0 1.0

            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        3     8  0 1 2 3 4 5 6 7",
        )
        .unwrap()
    }

    /// Returns the mesh used in Figure 5.2 of Smith-Griffiths-Margetts book
    ///
    /// ```text
    /// Example shown in Figure 5.2 from [@sgm] page 173
    ///
    ///          0.25       0.5      0.25 kN/m
    ///            ↓         ↓         ↓
    ///    ---    ▷0---------1---------2   Plane-Strain
    ///     |      |       ,'|       ,'|   Young = 1e6 kN/m²
    ///     |      |  0  ,'  |  2  ,'  |   Poisson = 0.3
    ///     |      |   ,'    |   ,'    |
    ///            | ,'   1  | ,'  3   |   connectivity:
    ///    1 m    ▷3'--------4'--------5     0 : 1 0 3
    ///            |       ,'|       ,'|     1 : 3 4 1
    ///     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
    ///     |      |   ,'    |   ,'    |     3 : 4 5 2
    ///     |      | ,'   5  | ,'   7  |     4 : 4 3 6
    ///    ---    ▷6'--------7'--------8     5 : 6 7 4
    ///            △         △         △     6 : 5 4 7
    ///                                      7 : 7 8 5
    ///            |------- 1 m -------|
    ///
    /// Note: the x-y origin is at the top-left (Point #0)
    ///
    /// References
    ///
    /// [@sgm] Smith, Griffiths and Margetts (5th ed) Figure 5.2 p173
    /// ```
    pub fn mesh_sgm_5_2() -> Mesh {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2      9     8

            # points
            # id    x    y
               0  0.0  0.0
               1  0.5  0.0
               2  1.0  0.0
               3  0.0 -0.5
               4  0.5 -0.5
               5  1.0 -0.5
               6  0.0 -1.0
               7  0.5 -1.0
               8  1.0 -1.0

            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     3  1 0 3
               1   1        2     3  3 4 1
               2   1        2     3  2 1 4
               3   1        2     3  4 5 2
               4   1        2     3  4 3 6
               5   1        2     3  6 7 4
               6   1        2     3  5 4 7
               7   1        2     3  7 8 5
        "#,
        )
        .unwrap()
    }

    /// Returns the mesh used in Example 1.6 of Bhatti's book
    ///
    /// ```text
    /// Example 1.6 from [@bhatti] page 32
    ///
    /// Solid bracket with thickness = 0.25 (plane-stress)
    ///
    /// The load a the top is normal to the slanted edge
    /// and has a value of 20 kN/m; thus Qn = -20 kN/m
    ///
    ///              1    load                connectivity:
    /// y=2.0  fixed *'-,__                    eid : vertices
    ///              |     '-,_  3   load        0 :  0, 2, 3
    /// y=1.5 - - -  |        ,'*-,__            1 :  3, 1, 0
    ///              |  1   ,'  |    '-,_  5     2 :  2, 4, 5
    /// y=1.0 - - -  |    ,'    |  3   ,-'*      3 :  5, 3, 2
    ///              |  ,'  0   |   ,-'   |
    ///              |,'        |,-'   2  |   constraints:
    /// y=0.0  fixed *----------*---------*     fixed on x and y
    ///              0          2         4
    ///             x=0.0     x=2.0     x=4.0
    ///
    /// References
    ///
    /// [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
    ///           and Applications, Wiley, 700p.
    /// ```
    pub fn mesh_bhatti_1_6() -> Mesh {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2      6     4
            
            # points
            # id    x   y
               0  0.0 0.0
               1  0.0 2.0
               2  2.0 0.0
               3  2.0 1.5
               4  4.0 0.0
               5  4.0 1.0
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     3  0 2 3
               1   1        2     3  3 1 0
               2   1        2     3  2 4 5
               3   1        2     3  5 3 2
            "#,
        )
        .unwrap()
    }

    /// Returns a sample 2D mesh with triangles and quadrilaterals
    ///
    /// ```text
    /// 3.1  14---------15
    ///       |  [11]333 |    
    /// 3.0  10---------11---------12--------------------13
    ///       |        .' '.        |                     |
    ///       | [8]  .'     '.  [9] |                     |
    ///       | 222.'         '.222 |        [10]         |  L
    ///       |  .'             '.  |         222         |  A
    ///       |.'                 '.|                     |  Y
    /// 2.0   7         [5]         8---------------------9  E
    ///       |'.       222       .'|                     |  R
    ///       |  '.             .'  |                     | 222
    ///       | [4]'.         .'[6] |         [7]         |
    ///       | 222  '.     .'  222 |         222         |
    ///       |        '. .'        |                     |
    /// 1.0   3----------4----------5---------------------6  <-- layer separation
    ///       |        .' '.        |                     |  L
    ///       | [0]  .'     '.  [2] |                     |  A
    ///       | 111.'   [1]   '.111 |         [3]         |  Y
    ///       |  .'     111     '.  |         111         |  E
    ///       |.'                 '.|                     |  R
    /// 0.0   0---------------------1---------------------2 111
    ///
    ///      0.0        1.0        2.0                   4.0
    /// ```
    pub fn mesh_rectangle_tris_quads() -> Mesh {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2     16    12
            
            # points
            # id    x   y
               0  0.0 0.0
               1  2.0 0.0
               2  4.0 0.0
            
               3  0.0 1.0
               4  1.0 1.0
               5  2.0 1.0
               6  4.0 1.0
            
               7  0.0 2.0
               8  2.0 2.0
               9  4.0 2.0
            
              10  0.0 3.0
              11  1.0 3.0
              12  2.0 3.0
              13  4.0 3.0
            
              14  0.0 3.1
              15  1.0 3.1
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0 111        2     3  0 4 3
               1 111        2     3  0 1 4
               2 111        2     3  1 5 4
               3 111        2     4  1 2 6 5
            
               4 222        2     3  3 4 7
               5 222        2     4  4 8 11 7
               6 222        2     3  4 5 8
               7 222        2     4  5 6 9 8
            
               8 222        2     3  7 11 10
               9 222        2     3  8 12 11
              10 222        2     4  8 9 13 12
            
              11 333        2     4  10 11 15 14
            "#,
        )
        .unwrap()
    }

    /// Returns a mesh with quadrilaterals representing a column
    ///
    /// ```text
    /// 3.0   6---------13
    ///       |    [5]   |
    ///       |          |   L
    /// 2.5   5---------12   A
    ///       |    [4]   |   Y
    ///       |          |   E
    /// 2.0   4---------11   R
    ///       |    [3]   |    
    ///       |          |   2
    /// 1.5   3---------10    
    ///       |    [2]   |
    ///       |          |
    /// 1.0   2----------9   <-- layer separation
    ///       |    [1]   |   L
    ///       |          |   A
    /// 0.5   1----------8   Y
    ///       |    [0]   |   E
    ///       |          |   R
    /// 0.0   0----------7   1
    ///
    ///      0.0        1.0
    /// ```
    pub fn mesh_column_two_layers_quads() -> Mesh {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2     14     6
            
            # points
            # id    x   y
               0  0.0 0.0
               1  0.0 0.5
               2  0.0 1.0
               3  0.0 1.5
               4  0.0 2.0
               5  0.0 2.5
               6  0.0 3.0
            
               7  1.0 0.0
               8  1.0 0.5
               9  1.0 1.0
              10  1.0 1.5
              11  1.0 2.0
              12  1.0 2.5
              13  1.0 3.0
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0  7  8 1
               1   1        2     4  1  8  9 2
            
               2   2        2     4  2  9 10 3
               3   2        2     4  3 10 11 4
               4   2        2     4  4 11 12 5
               5   2        2     4  5 12 13 6
            "#,
        )
        .unwrap()
    }

    /// Returns a mesh with triangles and distorted quadrilaterals representing a column
    ///
    /// ```text
    /// 3.1   6---------12
    ///       |    [6]   |       SOLID
    /// 3.0   5---------11       <-- layer separation
    ///       | [5]  ,-' |       POROUS
    ///       |   ,-'    |
    /// 2.5   4.-'       |
    ///       | '.   [4] |       L
    ///       | [3].     |       A
    /// 2.0   3.    '.   |       Y
    ///       | '--__ '. |       E
    ///       |      '--10  1.8  R
    ///       |          |       2
    ///       |   [2]    |
    ///       |          |
    /// 1.0   2----------9       <-- layer separation
    ///       |          |       L
    ///       |    [1]   |       A
    /// 0.5   1.__       |       Y
    ///       |   '--..  |       E
    ///       |  [0]   '-8  0.2  R
    /// 0.0   0----------7       1
    ///
    ///      0.0        1.0
    /// ```
    pub fn mesh_column_distorted_tris_quads() -> Mesh {
        Mesh::from_text(
            r#"
            # header
            # space_ndim npoint ncell
                       2     13     7
            
            # points
            # id    x   y
               0  0.0 0.0
               1  0.0 0.5
               2  0.0 1.0
               3  0.0 2.0
               4  0.0 2.5
               5  0.0 3.0
               6  0.0 3.1
            
               7  1.0 0.0
               8  1.0 0.2
               9  1.0 1.0
              10  1.0 1.8
              11  1.0 3.0
              12  1.0 3.1
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 7 8 1
               1   1        2     4  1 8 9 2
            
               2   2        2     4  2 9 10 3
               3   2        2     3  3 10 4
               4   2        2     3  4 10 11
               5   2        2     3  4 11 5
            
               6   3        2     4  5 11 12 6
            "#,
        )
        .unwrap()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Samples;

    #[test]
    fn sample_params_work() {
        let p = Samples::param_density_water(true);
        assert_eq!(p.cc, 1e-12);

        let p = Samples::param_density_water(false);
        assert_eq!(p.cc, 4.53e-7);

        let p = Samples::param_density_dry_air();
        assert_eq!(p.cc, 1.17e-5);

        let p = Samples::param_water(true);
        assert_eq!(p.density_liquid.cc, 1e-12);

        let p = Samples::param_water(false);
        assert_eq!(p.density_liquid.cc, 4.53e-7);

        let p = Samples::param_water_and_dry_air(true);
        assert_eq!(p.density_liquid.cc, 1e-12);

        let p = Samples::param_water_and_dry_air(false);
        assert_eq!(p.density_liquid.cc, 4.53e-7);

        let p = Samples::param_rod();
        assert_eq!(p.density, 2.0);

        let p = Samples::param_beam();
        assert_eq!(p.density, 2.0);

        let p = Samples::param_solid();
        assert_eq!(p.density, 2.7);

        let p = Samples::param_porous_sol_liq(0.4, 1.0);
        assert_eq!(p.porosity_initial, 0.4);

        let p = Samples::param_porous_sol_liq_gas(0.4, 1.0);
        assert_eq!(p.porosity_initial, 0.4);

        let p = Samples::param_seepage_liq();
        assert_eq!(p.porosity_initial, 0.4);

        let p = Samples::param_seepage_liq_gas();
        assert_eq!(p.porosity_initial, 0.4);
    }

    #[test]
    fn sample_meshes_work() {
        let m = Samples::mesh_segment();
        assert_eq!(m.cells.len(), 1);

        let m = Samples::mesh_square();
        assert_eq!(m.cells.len(), 1);

        let m = Samples::mesh_two_quads();
        assert_eq!(m.cells.len(), 2);

        let m = Samples::mesh_cube();
        assert_eq!(m.cells.len(), 1);
        assert_eq!(m.space_ndim, 3);

        let m = Samples::mesh_sgm_5_2();
        assert_eq!(m.cells.len(), 8);

        let m = Samples::mesh_bhatti_1_6();
        assert_eq!(m.cells.len(), 4);

        let m = Samples::mesh_rectangle_tris_quads();
        assert_eq!(m.cells.len(), 12);

        let m = Samples::mesh_column_two_layers_quads();
        assert_eq!(m.cells.len(), 6);

        let m = Samples::mesh_column_distorted_tris_quads();
        assert_eq!(m.cells.len(), 7);
    }
}
