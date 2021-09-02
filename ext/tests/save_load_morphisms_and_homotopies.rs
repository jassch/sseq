use algebra::module::homomorphism::ModuleHomomorphism;
use ext::chain_complex::{ChainComplex, ChainHomotopy};
use ext::utils::construct;
use ext::resolution_homomorphism::ResolutionHomomorphism;
use saveload::{Load, Save};
use std::io::{Cursor, Read, Seek, SeekFrom};
use std::sync::Arc;
use fp::vector::FpVector;

#[test]
fn test_save_load_resolution_hom() {
    let resolution = Arc::new(construct("S_2", None).unwrap());
    resolution.compute_through_bidegree(20, 20);

    let class_0_1_0 = vec![0];
    let hom = ResolutionHomomorphism::from_class("x_(0,1,0)".into(), resolution.clone(), resolution.clone(), 1, 1, &class_0_1_0);
    hom.extend(10,10);

    let mut cursor: Cursor<Vec<u8>> = Cursor::new(Vec::new());
    hom.save(&mut cursor).unwrap();

    cursor.seek(SeekFrom::Start(0)).unwrap();

    let hom2 = ResolutionHomomorphism::load(&mut cursor, &(resolution.clone(), resolution.clone())).unwrap();
    assert_eq!(0, cursor.bytes().count());

    // TODO: figure out which quasi_inverses should be nontrivial, so that this is a better test
    assert_eq!(
        hom.get_map(2).quasi_inverse(7),
        hom2.get_map(2).quasi_inverse(7)
    );

    hom.extend(20, 20);
    hom2.extend(20, 20);

    assert_eq!(
        hom.get_map(11).quasi_inverse(15),
        hom2.get_map(11).quasi_inverse(15)
    );

    assert_eq!(
        hom.get_map(13).quasi_inverse(12),
        hom2.get_map(13).quasi_inverse(12)
    );
}


#[test]
fn test_save_load_chain_homotopy() {
    let resolution = Arc::new(construct("S_2", None).unwrap());
    resolution.compute_through_bidegree(20, 20);

    let class_0_1_0 = vec![0];
    let hom_0_1 = ResolutionHomomorphism::from_class("x_(0,1,0)".into(), resolution.clone(), resolution.clone(), 1, 1, &class_0_1_0);
    hom_0_1.extend(20,20);

    let class_1_1_0 = vec![0];
    let hom_1_1 = ResolutionHomomorphism::from_class("x_(1,1,0)".into(), resolution.clone(), resolution.clone(), 1, 1, &class_1_1_0);
    hom_1_1.extend(20,20);

    // htpy for hom_1_1 in the middle, hom_0_1 on the right

    let f = |source_s, source_t, idx, row: &mut FpVector| {
            let mid_s = source_s - 1;

            hom_0_1.get_map(source_s)
                .compose(hom_1_1.get_map(mid_s))
                .apply_to_basis_element(row.as_slice_mut(), 1, source_t, idx)
        };

    let htpy = ChainHomotopy::new(
        &*resolution,
        &*resolution,
        2, // total s
        3, // total t
        f.clone(),
    );

    htpy.extend(10,10);

    let mut cursor: Cursor<Vec<u8>> = Cursor::new(Vec::new());
    htpy.save(&mut cursor).unwrap();

    cursor.seek(SeekFrom::Start(0)).unwrap();

    let htpy2 = ChainHomotopy::load(&mut cursor, &(&*resolution, &*resolution, 2, 3, f.clone())).unwrap();
    assert_eq!(0, cursor.bytes().count());


    for s in 2..=10 {
        let htpy_map = htpy.homotopy(s);
        let htpy2_map = htpy2.homotopy(s);
        assert_eq!(
            htpy_map.min_degree(),
            htpy2_map.min_degree()
        );
        assert_eq!(
            htpy_map.next_degree(),
            htpy2_map.next_degree()
        );
        for t in htpy_map.min_degree()..htpy_map.next_degree() {
            for idx in 0..resolution.number_of_gens_in_bidegree(s,t) {
                assert_eq!(
                    htpy_map.output(t, idx),
                    htpy2_map.output(t, idx)
                );
            }
        }
    }

    /*
    assert_eq!(
        hom.get_map(2).quasi_inverse(7),
        hom2.get_map(2).quasi_inverse(7)
    );

    hom.extend(20, 20);
    hom2.extend(20, 20);

    assert_eq!(
        hom.get_map(11).quasi_inverse(15),
        hom2.get_map(11).quasi_inverse(15)
    );

    assert_eq!(
        hom.get_map(13).quasi_inverse(12),
        hom2.get_map(13).quasi_inverse(12)
    );
    */
}