#![cfg_attr(rustfmt, rustfmt_skip)]
use crate::module::{Module, QuotientModule};
use fp::matrix::{QuasiInverse, Subspace};
use fp::vector::{FpVector, SliceMut};
use std::sync::Arc;

use crate::module::homomorphism::ModuleHomomorphism;

pub struct QuotientHomomorphism<F: ModuleHomomorphism> {
    f: Arc<F>,
    s: Arc<QuotientModule<F::Source>>,
    t: Arc<QuotientModule<F::Target>>,
}

impl<F: ModuleHomomorphism> QuotientHomomorphism<F> {
    pub fn new(
        f: Arc<F>,
        s: Arc<QuotientModule<F::Source>>,
        t: Arc<QuotientModule<F::Target>>,
    ) -> Self {
        QuotientHomomorphism { f, s, t }
    }
}

impl<F: ModuleHomomorphism> ModuleHomomorphism for QuotientHomomorphism<F> {
    type Source = QuotientModule<F::Source>;
    type Target = QuotientModule<F::Target>;

    fn source(&self) -> Arc<Self::Source> {
        Arc::clone(&self.s)
    }
    fn target(&self) -> Arc<Self::Target> {
        Arc::clone(&self.t)
    }
    fn degree_shift(&self) -> i32 {
        self.f.degree_shift()
    }

    fn apply_to_basis_element(
        &self,
        result: SliceMut,
        coeff: u32,
        input_degree: i32,
        input_idx: usize,
    ) {
        let output_degree = input_degree - self.degree_shift();
        let mut result_ = FpVector::new(self.prime(), self.t.module.dimension(output_degree));
        self.f.apply_to_basis_element(
            result_.as_slice_mut(),
            coeff,
            input_degree,
            self.s.basis_list[input_degree][input_idx],
        );

        self.t.reduce(output_degree, result_.as_slice_mut());
        self.t.old_basis_to_new(output_degree, result, result_.as_slice());
    }

    fn kernel(&self, _degree: i32) -> &Subspace {
        unimplemented!();
    }

    fn quasi_inverse(&self, _degree: i32) -> &QuasiInverse {
        unimplemented!();
    }

    fn compute_kernels_and_quasi_inverses_through_degree(&self, _degree: i32) {
        unimplemented!();
    }
}

pub struct QuotientHomomorphismSource<F: ModuleHomomorphism> {
    f: Arc<F>,
    s: Arc<QuotientModule<F::Source>>,
    t: Arc<F::Target>,
}

impl<F: ModuleHomomorphism> QuotientHomomorphismSource<F> {
    pub fn new(f: Arc<F>, s: Arc<QuotientModule<F::Source>>, t: Arc<F::Target>) -> Self {
        QuotientHomomorphismSource { f, s, t }
    }
}
impl<F: ModuleHomomorphism> ModuleHomomorphism for QuotientHomomorphismSource<F> {
    type Source = QuotientModule<F::Source>;
    type Target = F::Target;

    fn source(&self) -> Arc<Self::Source> {
        Arc::clone(&self.s)
    }
    fn target(&self) -> Arc<Self::Target> {
        Arc::clone(&self.t)
    }
    fn degree_shift(&self) -> i32 {
        self.f.degree_shift()
    }

    fn apply_to_basis_element(
        &self,
        result: SliceMut,
        coeff: u32,
        input_degree: i32,
        input_idx: usize,
    ) {
        self.f.apply_to_basis_element(
            result,
            coeff,
            input_degree,
            self.s.basis_list[input_degree][input_idx],
        );
    }

    fn kernel(&self, _degree: i32) -> &Subspace {
        unimplemented!();
    }

    fn quasi_inverse(&self, _degree: i32) -> &QuasiInverse {
        unimplemented!();
    }

    fn compute_kernels_and_quasi_inverses_through_degree(&self, _degree: i32) {
        unimplemented!();
    }
}