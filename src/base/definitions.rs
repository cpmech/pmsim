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
/// use pmsim::base::{AttrElement, Element};
///
/// let attr_element = AttrElement::from([(1, Element::Solid)]);
/// ```
pub type AttrElement = HashMap<CellAttributeId, Element>;

/// Holds all attributes/DOFs; maps (CellAttributeId,GeoKind) to a (nnode,ndof) table
///
/// # Examples
///
/// ```
/// use gemlab::shapes::GeoKind;
/// use pmsim::base::{display_attr_dofs, AttrElement, Dof, Element};
/// use std::collections::HashMap;
///
/// let attr_element = AttrElement::from([(1, Element::Solid), (2, Element::Solid)]);
/// let attr_dofs = HashMap::from([
///     (
///         (1, GeoKind::Tri3),
///         vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]],
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
///     format!("{}", display_attr_dofs(&attr_element, &attr_dofs)),
///     "(1, Tri3) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
///      (2, Qua4) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
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
/// let point_dofs = vec![
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
///     HashSet::from([Dof::Ux, Dof::Uy]),
/// ];
/// assert_eq!(
///     format!("{}", display_point_dofs(&point_dofs)),
///     "0 → [Ux, Uy]\n\
///      1 → [Ux, Uy]\n\
///      2 → [Ux, Uy]\n\
///      3 → [Ux, Uy]\n"
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
/// let point_equations = vec![
///     vec![0, 1],
///     vec![2, 3],
///     vec![4, 5],
///     vec![6, 7],
/// ];
/// assert_eq!(
///     format!("{}", display_point_equations(&point_equations)),
///     "0 → [0, 1]\n\
///      1 → [2, 3]\n\
///      2 → [4, 5]\n\
///      3 → [6, 7]\n"
/// );
/// ```
pub type PointEquations = Vec<Vec<usize>>;

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

/// Returns a string representing a pair of (AttrElement,AttrDofs)
///
/// # Panics
///
/// The `attr_element` map must have the keys used in the `attr_dofs` map,
/// otherwise a panic will occur.
pub fn display_attr_dofs(attr_element: &AttrElement, attr_dofs: &AttrDofs) -> String {
    let mut buffer = String::new();
    let mut keys: Vec<_> = attr_dofs.keys().collect();
    keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for key in keys {
        let elem = attr_element.get(&key.0).unwrap();
        let dofs = attr_dofs.get(key).unwrap();
        write!(&mut buffer, "{:?} → {:?} → {:?}\n", key, elem, dofs).unwrap();
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{
        display_attr_dofs, display_attr_element, display_point_dofs, display_point_equations, AttrElement, Dof, Element,
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
        let attr_element = AttrElement::from([(1, Element::Solid), (2, Element::Solid)]);
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
            format!("{}", display_attr_dofs(&attr_element, &attr_dofs)),
            "(1, Tri3) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             (2, Qua4) → Solid → [[Ux, Uy], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
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
}
