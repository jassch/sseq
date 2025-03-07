use crate::algebra::{
    adem_algebra::AdemBasisElement,
    milnor_algebra::{MilnorBasisElement, PPartEntry},
    AdemAlgebra, MilnorAlgebra, SteenrodAlgebraBorrow, SteenrodAlgebraT,
};
use crate::module::{BoundedModule, Module, ZeroModule};
use fp::prime::{Binomial, ValidPrime};
use fp::vector::SliceMut;

use std::sync::Arc;

#[cfg(feature = "json")]
use {serde::Deserialize, serde_json::Value};

/// This is $\mathbb{RP}_{\mathrm{min}}^{\mathrm{max}}$. The cohomology is the subquotient of
/// $\mathbb{F}_2[x^\pm]$ given by elements of degree between min and max (inclusive)
///
/// The `clear_bottom` option, if selected, mods out by the elements in the *$A(2)$ submodule*
/// generated by degrees less than `min`. This is useful for approximating $\mathrm{tmf} \wedge
/// \mathbb{RP}_{-\infty}^n$, c.f. Proposition 2.2 of Bailey and Ricka. Note that this quotient
/// always has minimum degree -1 mod 8.
pub struct RealProjectiveSpace<A: SteenrodAlgebraT> {
    algebra: Arc<A>,
    pub min: i32,
    pub max: Option<i32>, // If None,  then RP^oo
    pub clear_bottom: bool,
}

impl<A: SteenrodAlgebraT> std::fmt::Display for RealProjectiveSpace<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let clear = if self.clear_bottom {
            " (clear_bottom)"
        } else {
            ""
        };

        if let Some(max) = self.max {
            write!(f, "RP^{}_{}{}", max, self.min, clear)
        } else {
            write!(f, "RP_{}{}", self.min, clear)
        }
    }
}

impl<A: SteenrodAlgebraT> PartialEq for RealProjectiveSpace<A> {
    fn eq(&self, other: &Self) -> bool {
        self.min == other.min && self.max == other.max
    }
}

impl<A: SteenrodAlgebraT> Eq for RealProjectiveSpace<A> {}

impl<A: SteenrodAlgebraT> Module for RealProjectiveSpace<A> {
    type Algebra = A;

    fn algebra(&self) -> Arc<A> {
        Arc::clone(&self.algebra)
    }

    fn min_degree(&self) -> i32 {
        self.min
    }

    fn max_computed_degree(&self) -> i32 {
        i32::max_value()
    }

    fn dimension(&self, degree: i32) -> usize {
        if degree < self.min {
            return 0;
        }
        if let Some(m) = self.max {
            if degree > m {
                return 0;
            }
        }

        if self.clear_bottom
            && (degree == self.min + 1
                || degree == self.min + 1 + 1
                || degree == self.min + 1 + 2
                || degree == self.min + 1 + 4
                || degree == self.min + 1 + 8)
        {
            return 0;
        }
        1
    }

    fn basis_element_to_string(&self, degree: i32, _idx: usize) -> String {
        // It is an error to call the function if self.dimension(degree) == 0
        format!("x^{{{}}}", degree)
    }

    fn act_on_basis(
        &self,
        mut result: SliceMut,
        coeff: u32,
        op_degree: i32,
        op_index: usize,
        mod_degree: i32,
        mod_index: usize,
    ) {
        assert!(op_index < self.algebra().dimension(op_degree, mod_degree));
        assert!(mod_index < self.dimension(mod_degree));

        let output_degree = mod_degree + op_degree;

        if op_degree == 0 || coeff == 0 || self.dimension(output_degree) == 0 {
            return;
        }

        if match self.algebra.steenrod_algebra() {
            SteenrodAlgebraBorrow::BorrowAdem(a) => coef_adem(a, op_degree, op_index, mod_degree),
            SteenrodAlgebraBorrow::BorrowMilnor(a) => {
                coef_milnor(a, op_degree, op_index, mod_degree)
            }
        } {
            result.add_basis_element(0, 1);
        }
    }
}

// Compute the coefficient of the operation on x^j.
fn coef_adem(algebra: &AdemAlgebra, op_deg: i32, op_idx: usize, mut j: i32) -> bool {
    let elt: &AdemBasisElement = algebra.basis_element_from_index(op_deg, op_idx);
    // Apply Sq^i to x^j and see if it is zero
    for i in elt.ps.iter().rev() {
        let c = if j >= 0 {
            i32::binomial(ValidPrime::new(2), j, *i as i32)
        } else {
            i32::binomial(ValidPrime::new(2), -j + (*i as i32) - 1, *i as i32)
        };
        if c == 0 {
            return false;
        }
        // Somehow j += 1 produces the same answer...
        j += *i as i32;
    }
    true
}

fn coef_milnor(algebra: &MilnorAlgebra, op_deg: i32, op_idx: usize, mut mod_degree: i32) -> bool {
    if mod_degree == 0 {
        return false;
    }

    let elt: &MilnorBasisElement = algebra.basis_element_from_index(op_deg, op_idx);

    let sum: PPartEntry = elt.p_part.iter().sum();
    if mod_degree < 0 {
        mod_degree = sum as i32 - mod_degree - 1;
    } else if mod_degree < sum as i32 {
        return false;
    }

    let mod_degree = mod_degree as PPartEntry;

    let mut list = Vec::with_capacity(elt.p_part.len() + 1);
    list.push(mod_degree - sum);
    list.extend_from_slice(&elt.p_part);

    PPartEntry::multinomial2(&list) == 1
}

impl<A: SteenrodAlgebraT> ZeroModule for RealProjectiveSpace<A> {
    fn zero_module(algebra: Arc<A>, min_degree: i32) -> Self {
        Self::new(algebra, min_degree, Some(min_degree - 1), false)
    }
}

impl<A: SteenrodAlgebraT> RealProjectiveSpace<A> {
    pub fn new(algebra: Arc<A>, min: i32, max: Option<i32>, clear_bottom: bool) -> Self {
        assert_eq!(*algebra.prime(), 2);
        if let Some(max) = max {
            assert!(max >= min);
        }
        Self {
            algebra,
            min,
            max,
            clear_bottom,
        }
    }
}

#[cfg(feature = "json")]
#[derive(Deserialize, Debug)]
struct RPSpec {
    min: i32,
    clear_bottom: Option<bool>,
    max: Option<i32>,
}

#[cfg(feature = "json")]
impl<A: SteenrodAlgebraT> RealProjectiveSpace<A> {
    pub fn from_json(algebra: Arc<A>, json: &Value) -> error::Result<Self> {
        let spec: RPSpec = RPSpec::deserialize(json)?;
        let clear_bottom = spec.clear_bottom.unwrap_or(false);
        let mut min = spec.min;
        if clear_bottom {
            let x = (spec.min + 1).rem_euclid(8);
            if x != 0 {
                min += 8 - x;
            }
        }

        Ok(Self {
            algebra,
            min,
            clear_bottom,
            max: spec.max,
        })
    }

    pub fn to_json(&self, json: &mut Value) {
        json["name"] = Value::String(self.to_string());
        json["type"] = Value::from("real projective space");
        json["min"] = Value::from(self.min);
        if let Some(max) = self.max {
            json["max"] = Value::from(max);
        }
        if self.clear_bottom {
            json["clear_bottom"] = Value::Bool(true);
        }
    }
}

impl<A: SteenrodAlgebraT> BoundedModule for RealProjectiveSpace<A> {
    /// `max_degree` is the a degree such that if t > `max_degree`, then `self.dimension(t) = 0`.
    fn max_degree(&self) -> i32 {
        self.max.unwrap_or(i32::MAX)
    }
}
