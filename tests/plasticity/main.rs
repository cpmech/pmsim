mod general_vm_single_elem_2d;
mod general_vm_soft_single_elem_2d;
mod general_vm_soft_single_elem_3d;
mod von_mises_2x2_elements_2d;
// The kinematic suite lives in its own folder (kinematic_vm), which also holds the Python
// reference material of the model; the path attribute keeps the test names unchanged.
#[path = "kinematic_vm/von_mises_kinematic.rs"]
mod von_mises_kinematic;
mod von_mises_single_element_2d;
