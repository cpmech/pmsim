use super::{Dof, Element};
use gemlab::mesh::CellAttributeId;
use gemlab::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

/// Holds all attributes/elements; maps CellAttributeId to Element type
///
/// # Examples
///
/// ```
/// use pmsim::base::{display_attr_element, AttrElement, Element};
///
/// let attr_element = AttrElement::from([
///     (1, Element::Solid),
///     (2, Element::PorousSldLiq),
/// ]);
/// assert_eq!(format!("{}", display_attr_element(&attr_element)),
///     "1 → Solid\n\
///      2 → PorousSldLiq\n"
/// );
/// ```
pub type AttrElement = HashMap<CellAttributeId, Element>;

/// Holds all attributes/DOFs; maps (CellAttributeId,GeoKind) to a (nnode,ndof) table
///
/// # Examples
///
/// ```
/// use gemlab::shapes::GeoKind;
/// use pmsim::base::{display_attr_dofs, Dof, Element};
/// use std::collections::HashMap;
///
/// let attr_dofs = HashMap::from([
///     (
///         (1, GeoKind::Tri3),
///         vec![
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///         ],
///     ),
///     (
///         (2, GeoKind::Qua4),
///         vec![
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///             vec![Dof::Ux, Dof::Uy],
///         ],
///     ),
/// ]);
/// assert_eq!(
///     format!("{}", display_attr_dofs(&attr_dofs)),
///     "(1, Tri3) → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
///      (2, Qua4) → [[Ux, Uy], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
/// );
/// ```
pub type AttrDofs = HashMap<(CellAttributeId, GeoKind), Vec<Vec<Dof>>>;

/// Holds all point DOFs (npoint); maps PointId (index of point) to a set of DOFs
///
/// # Examples
///
/// ```
/// use pmsim::base::{display_point_dofs, Dof};
/// use std::collections::HashSet;
///
/// //           {Ux}
/// //           {Uy}
/// //           {Pl}
/// //            2
/// //           / \
/// //    {Ux}  /   \  {Ux}
/// //    {Uy} 5     4 {Uy}
/// //        /       \
/// // {Ux}  /         \  {Ux}
/// // {Uy} 0-----3-----1 {Uy}
/// // {Pl}      {Ux}     {Pl}
/// //           {Uy}
///
/// let point_dofs = vec![
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy, Dof::Pl]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
/// ];
/// assert_eq!(
///     format!("{}", display_point_dofs(&point_dofs)),
///     "0 → [Ux, Uy, Pl]\n\
///      1 → [Ux, Uy, Pl]\n\
///      2 → [Ux, Uy, Pl]\n\
///      3 → [Ux, Uy]\n\
///      4 → [Ux, Uy]\n\
///      5 → [Ux, Uy]\n"
/// );
/// ```
pub type PointDofs = Vec<HashSet<Dof>>;

/// Holds all point equation numbers (npoint); maps PointId (index of point) to a set of equation numbers
///
/// # Examples
///
/// ```
/// use pmsim::base::display_point_equations;
///
/// //            {Ux → 6}
/// //            {Uy → 7}
/// //            {Pl → 8}
/// //                2
/// //               / \
/// //   {Ux → 13}  /   \  {Ux → 11}
/// //   {Uy → 14} 5     4 {Uy → 12}
/// //            /       \
/// // {Ux → 0}  /         \  {Ux → 3}
/// // {Uy → 1} 0-----3-----1 {Uy → 4}
/// // {Pl → 2}   {Ux → 9}    {Pl → 5}
/// //            {Uy → 10}
///
/// let point_equations = vec![
///     vec![ 0, 1, 2],
///     vec![ 3, 4, 5],
///     vec![ 6, 7, 8],
///     vec![ 9, 10],
///     vec![11, 12],
///     vec![13, 14],
/// ];
/// assert_eq!(
///     format!("{}", display_point_equations(&point_equations)),
///     "0 → [0, 1, 2]\n\
///      1 → [3, 4, 5]\n\
///      2 → [6, 7, 8]\n\
///      3 → [9, 10]\n\
///      4 → [11, 12]\n\
///      5 → [13, 14]\n"
/// );
/// ```
pub type PointEquations = Vec<Vec<usize>>;

/// Holds all local-to-global mappings (ncell,n_local_equation)
///
/// # Examples
///
/// ```
/// use pmsim::base::display_local_to_global;
///
/// //       {8} 4---.__
/// //       {9}/ \     `--.___3 {6}   [#] indicates id
/// //         /   \          / \{7}   (#) indicates attribute_id
/// //        /     \  [1]   /   \     {#} indicates equation number
/// //       /  [0]  \ (1)  / [2] \
/// // {0}  /   (1)   \    /  (1)  \
/// // {1} 0---.__     \  /      ___2 {4}
/// //            `--.__\/__.---'     {5}
/// //                   1 {2}
/// //                     {3}
/// let local_to_global = vec![
///     vec![0, 1, 2, 3, 8, 9],
///     vec![2, 3, 6, 7, 8, 9],
///     vec![2, 3, 4, 5, 6, 7],
/// ];
/// assert_eq!(
///     format!("{}", display_local_to_global(&local_to_global)),
///     "0 → [0, 1, 2, 3, 8, 9]\n\
///      1 → [2, 3, 6, 7, 8, 9]\n\
///      2 → [2, 3, 4, 5, 6, 7]\n"
/// );
/// ```
pub type LocalToGlobal = Vec<Vec<usize>>;

/// Returns a string representing an AttrElement data structure
pub fn display_attr_element(attr_element: &AttrElement) -> String {
    let mut buffer = String::new();
    let mut keys: Vec<_> = attr_element.keys().collect();
    keys.sort();
    for key in keys {
        write!(&mut buffer, "{:?} → {:?}\n", key, attr_element.get(key).unwrap()).unwrap();
    }
    buffer
}

/// Returns a string representing an AttrDofs data structure
pub fn display_attr_dofs(attr_dofs: &AttrDofs) -> String {
    let mut buffer = String::new();
    let mut keys: Vec<_> = attr_dofs.keys().collect();
    keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for key in keys {
        write!(&mut buffer, "{:?} → {:?}\n", key, attr_dofs.get(key).unwrap()).unwrap();
    }
    buffer
}

/// Returns a string representing a PointDofs data structure
pub fn display_point_dofs(point_dofs: &PointDofs) -> String {
    let mut buffer = String::new();
    for point_id in 0..point_dofs.len() {
        let mut dofs: Vec<_> = point_dofs[point_id].iter().collect();
        dofs.sort();
        write!(&mut buffer, "{:?} → {:?}\n", point_id, dofs).unwrap();
    }
    buffer
}

/// Returns a string representing a PointEquations data structure
pub fn display_point_equations(point_equations: &PointEquations) -> String {
    let mut buffer = String::new();
    for point_id in 0..point_equations.len() {
        let eqs = &point_equations[point_id];
        write!(&mut buffer, "{:?} → {:?}\n", point_id, eqs).unwrap();
    }
    buffer
}

/// Returns a string representing a LocalToGlobal data structure
pub fn display_local_to_global(local_to_global: &LocalToGlobal) -> String {
    let mut buffer = String::new();
    for cell_id in 0..local_to_global.len() {
        let eqs = &local_to_global[cell_id];
        write!(&mut buffer, "{:?} → {:?}\n", cell_id, eqs).unwrap();
    }
    buffer
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{
        display_attr_dofs, display_attr_element, display_local_to_global, display_point_dofs, display_point_equations,
        AttrElement, Dof, Element,
    };
    use gemlab::shapes::GeoKind;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn display_attr_element_works() {
        let attr_element = AttrElement::from([(1, Element::Solid)]);
        assert_eq!(format!("{}", display_attr_element(&attr_element)), "1 → Solid\n");
    }

    #[test]
    fn display_attr_dofs_works() {
        let attr_dofs = HashMap::from([
            (
                (1, GeoKind::Tri3),
                vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]],
            ),
            (
                (2, GeoKind::Qua4),
                vec![
                    vec![Dof::Ux, Dof::Uy],
                    vec![Dof::Ux, Dof::Uy],
                    vec![Dof::Ux, Dof::Uy],
                    vec![Dof::Ux, Dof::Uy],
                ],
            ),
        ]);
        assert_eq!(
            format!("{}", display_attr_dofs(&attr_dofs)),
            "(1, Tri3) → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             (2, Qua4) → [[Ux, Uy], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
        );
    }

    #[test]
    fn display_point_dofs_works() {
        let point_dofs = vec![
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
            HashSet::from([Dof::Ux, Dof::Uy]),
        ];
        assert_eq!(
            format!("{}", display_point_dofs(&point_dofs)),
            "0 → [Ux, Uy]\n\
             1 → [Ux, Uy]\n\
             2 → [Ux, Uy]\n\
             3 → [Ux, Uy]\n"
        );
    }

    #[test]
    fn display_point_equations_works() {
        let point_equations = vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7]];
        assert_eq!(
            format!("{}", display_point_equations(&point_equations)),
            "0 → [0, 1]\n\
             1 → [2, 3]\n\
             2 → [4, 5]\n\
             3 → [6, 7]\n"
        );
    }

    #[test]
    fn display_local_to_global_works() {
        let local_to_global = vec![vec![0, 1, 2, 3, 8, 9], vec![2, 3, 6, 7, 8, 9], vec![2, 3, 4, 5, 6, 7]];
        assert_eq!(
            format!("{}", display_local_to_global(&local_to_global)),
            "0 → [0, 1, 2, 3, 8, 9]\n\
             1 → [2, 3, 6, 7, 8, 9]\n\
             2 → [2, 3, 4, 5, 6, 7]\n"
        );
    }
}
