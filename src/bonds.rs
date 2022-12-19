use std::{collections::VecDeque, ops::Deref};

use crate::{field::Field, lattice::Lattice, pause};

const VERBOSE: bool = false;

/// Models a field that has values on the bonds between sites.
pub struct BondsField<'a, T, const D: usize, const SIZE: usize>(
    Field<'a, [T; D * 2_usize], D, SIZE>,
)
where
    [(); D * 2_usize]:;

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    T: Default,
    [(); D * 2_usize]:,
{
    /// Constructor for a new BondsField with the default value.
    fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        let values: Vec<()> = vec![(); SIZE];
        let values: Vec<[T; D * 2_usize]> = values
            .into_iter()
            .map(|_| {
                let ary: [(); D * 2_usize] = [(); D * 2_usize];
                let ary: [T; D * 2_usize] = ary.map(|_| T::default());
                ary
            })
            .collect();

        BondsField(Field::<'a, [T; D * 2_usize], D, SIZE> { values, lattice })
    }
}

impl<'a, T, const D: usize, const SIZE: usize> Deref for BondsField<'a, T, D, SIZE>
where
    [(); D * 2_usize]:,
{
    type Target = Field<'a, [T; D * 2_usize], D, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    T: Clone,
    [(); D * 2_usize]:,
{
    /// Activate a link of the LinksField
    pub fn set_value(&mut self, index: usize, direction: usize, value: T) {
        let neighbour = self.lattice.get_neighbours_array(index)[direction];
        self.0.values[index][direction] = value.clone();
        self.0.values[neighbour][(direction + D) % (D * 2_usize)] = value;
    }
}

impl<'a, T, const D: usize, const SIZE: usize> BondsField<'a, T, D, SIZE>
where
    T: PartialEq + std::fmt::Debug,
    [(); D * 2_usize]:,
{
    /// A method to check if the links, which are saved on both of the
    /// sites, are correct at both points.
    #[allow(dead_code)]
    pub fn check_coherence(&self) {
        for (index, links) in self.values.iter().enumerate() {
            for (direction, link) in links.iter().enumerate() {
                assert_eq!(
                    *link,
                    self.values[self.lattice.get_neighbours_array(index)[direction]]
                        [(direction + D) % (D * 2)]
                );
            }
        }
    }
}

/// Models the activation of outgoing links from a lattice site
pub struct LinksField<'a, const D: usize, const SIZE: usize>(BondsField<'a, bool, D, SIZE>)
where
    [(); D * 2_usize]:;

impl<'a, const D: usize, const SIZE: usize> LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Constructor for a new LinksField on the lattice initialized to be false everywhere.
    pub fn new(lattice: &'a Lattice<D, SIZE>) -> Self {
        assert_eq!(bool::default(), false);
        LinksField(BondsField::<bool, D, SIZE>::new(lattice))
    }
}

impl<'a, const D: usize, const SIZE: usize> Deref for LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    type Target = BondsField<'a, bool, D, SIZE>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a, const D: usize, const SIZE: usize> LinksField<'a, D, SIZE>
where
    [(); D * 2_usize]:,
{
    /// Activate a link of the LinksField
    pub fn activate(&mut self, index: usize, direction: usize) {
        self.0.set_value(index, direction, true);
    }

    pub fn is_active(&self, index: usize, direction: usize) -> bool {
        self.values[index][direction]
    }

    /// Return collection of all clusters
    pub fn collect_clusters(&self) -> Vec<Vec<usize>> {
        // Remember if a given site has already been considered for a cluster
        let mut unvisited: Vec<bool> = vec![true; SIZE];
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for index in 0..SIZE {
            if unvisited[index] {
                // Build a new cluster if you come across a site not yet in one.
                clusters.push(self.build_cluster(index, &mut unvisited));
            }
        }
        clusters
    }

    /// Build a cluster from activated links
    fn build_cluster(&self, index: usize, unvisited: &mut Vec<bool>) -> Vec<usize> {
        let mut cluster: Vec<usize> = Vec::with_capacity(SIZE / 2);
        let mut checklist: VecDeque<usize> = VecDeque::with_capacity(SIZE / 2);
        // Add the new site to a checklist for sites to check for
        // activated, unvisited neighbours to be added to the cluster
        if VERBOSE {
            println!("Starting to build the cluster around {index}");
        }
        checklist.push_back(index);
        unvisited[index] = false;

        // Continuously working through the checklist, for each entry...
        while let Some(check) = checklist.pop_front() {
            // ...add them to the cluster list and...
            cluster.push(check);
            if VERBOSE {
                pause();
                println!(
                    "Looking at {check}, with links {:?} and neighbours {:?}",
                    self.values[check],
                    self.lattice.get_neighbours_array(check)
                );
            }

            // ...check each neighbour...
            for (direction, neighbour) in self
                .lattice
                .get_neighbours_array(check)
                .into_iter()
                .enumerate()
            {
                // ...if it hasn't been visited and the link is active,
                // add it to the checklist and mark it as visited.
                if unvisited[neighbour] && self.values[check][direction] {
                    if VERBOSE {
                        println!("{neighbour} gets added to the checklist :)");
                    }
                    checklist.push_back(neighbour);
                    unvisited[neighbour] = false;
                } else if VERBOSE {
                    println!("{neighbour} is not added to the checklist :(")
                }
            }
        }

        cluster.shrink_to_fit();
        cluster
    }
}

#[test]
fn test_links_field() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let mut links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    // Setting the link at [0, 2, 1] (index = 28) in the third direction
    // (direction = 2).
    assert_eq!(links_field.values[28][2], false);
    assert_eq!(links_field.values[48][5], false);
    links_field.activate(28, 2);
    //links_field.print_values_formated(SIZE_ARY);
    assert_eq!(links_field.values[28][2], true);
    assert_eq!(links_field.values[48][5], true);

    // Setting the link at [2, 1, 0] (index = 6) in the negatice first direction
    // (direction = 3).
    assert_eq!(links_field.values[6][3], false);
    assert_eq!(links_field.values[5][0], false);
    links_field.activate(6, 3);
    assert_eq!(links_field.values[6][3], true);
    assert_eq!(links_field.values[5][0], true);
}

#[test]
fn test_all_links_deactivated() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    let clusters: Vec<Vec<usize>> = links_field.collect_clusters();

    let mut test: Vec<Vec<usize>> = Vec::new();
    for index in 0..SIZE {
        let mut vect: Vec<usize> = Vec::new();
        vect.push(index);
        test.push(vect);
    }

    assert_eq!(test, clusters);
}

#[test]
fn test_all_links_active() {
    const D: usize = 3;
    const SIZE_ARY: [usize; D] = [4, 5, 3];
    const SIZE: usize = SIZE_ARY[0] * SIZE_ARY[1] * SIZE_ARY[2];

    let lattice: Lattice<D, SIZE> = Lattice::new(SIZE_ARY);
    let mut links_field: LinksField<D, SIZE> = LinksField::new(&lattice);

    for links in links_field.0 .0.values.iter_mut() {
        for link in links.iter_mut() {
            *link = true;
        }
    }

    let clusters: Vec<Vec<usize>> = links_field.collect_clusters();

    let mut test: Vec<Vec<usize>> = Vec::new();
    let mut cluster: Vec<usize> = Vec::with_capacity(SIZE);
    for index in 0..SIZE {
        cluster.push(index);
    }
    test.push(cluster);

    assert_eq!(test.len(), clusters.len());
}
