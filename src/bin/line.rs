use lattice_qft::lattice::Lattice3d;

const CUBE: usize = 3;

fn main() {
    let lattice: Lattice3d<CUBE, CUBE, CUBE> = Lattice3d::new();
}
